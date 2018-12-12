#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pthread.h>

#include "graph.h"
#include "stack.h"
#include "tour.h"

#define DEBUG 0

// ==================================================
//
//  constants & important globals
//
// ==================================================

// constants
const int MASTER = 0;
const int START = 0;
const int DEFAULT_STACK_SIZE = 64;

// MPI
int np;   // MPI number of processes
int rank; // MPI rank

// pthreads
pthread_mutex_t best_tour_mutex;
pthread_mutex_t global_stack_mutex;
pthread_cond_t global_stack_full;
pthread_cond_t global_stack_empty;
Stack* global_stack;

// important globals
Tour* best = NULL;       // best tour
int ncities;             // number of cities
int waiting_threads = 0; // number of threads waiting on global_stack_empty
int nthreads;            // number of threads

// ==================================================
//
//  auxiliary macros
//
// ==================================================

#define MASTER 0

#define MPI_TAG_TOUR_N      11
#define MPI_TAG_TOUR_COST   12
#define MPI_TAG_TOUR_COUNT  13
#define MPI_TAG_TOUR_CITIES 14

// ==================================================
//
//  auxiliary functions
//
// ==================================================

static void mpiassert(int result);

static void updatebest(Tour*);

static void globalpush(Stack*, Tour*);
static Tour* globalpop(Stack*);
static void pushwork(Stack*, Tour*);
static Tour* popwork(Stack*);

#define pushcopywork(s, t) (pushwork(s, tour_copy(t)))

static pthread_t newthread(int id, Stack* stacks);

static Stack** dividework(Tour* tour, int n);

// ==================================================
//
//  functions
//
// ==================================================

Stack* BFS(Tour* beginning, int nt) {
    Stack* stack = stack_new(ncities * ncities * ncities);
    stack_push_copy(stack, beginning);

    while (stack->size < nt) {
        Tour* tour = stack_pop(stack);
        if (!tour) {
            goto END;
        }

        for (int neighbor = 0; neighbor < ncities; neighbor++) {
            if (!tour_visits_city(tour, neighbor)) {
                tour_add_city(tour, neighbor);
                stack_push_copy(stack, tour);
                tour_remove_last_city(tour);
            }
        }

        tour_free(tour);
    }

    END: {
        return stack;
    }
}

static bool feasible(Tour* tour, int city) {
    // if the city has already been tour_visits_city
    if (tour_visits_city(tour, city)) {
        return false;
    }

    // if it can lead to a least cost tour
    tour_add_city(tour, city);
    int newcost = tour->cost;
    tour_remove_last_city(tour);
    if (newcost > best->cost) {
        return false;
    }

    return true;
}

void findbest(int id, Stack* stack) {
    #if DEBUG
        printf("STARTING THREAD %d\n", id);
    #endif

    Tour* tour;
    while ((tour = popwork(stack))) {
        if (tour->count == ncities && tour->cost < best->cost) {
            updatebest(tour);
        } else {
            for (int neighbor = 0; neighbor < ncities; neighbor++) {
                if (feasible(tour, neighbor)) {
                    tour_add_city(tour, neighbor);
                    pushcopywork(stack, tour);
                    tour_remove_last_city(tour);
                }
            }
        }

        tour_free(tour);
    }

    stack_free(stack);
    #if DEBUG
        printf("ENDING THREAD %d\n", id);
    #endif
}

void* threadfindbest(void* arguments) {
    int id = *((int*)((void**)arguments)[0]);
    Stack* stack = (Stack*)((void**)arguments)[1];
    findbest(id, stack);
    return NULL;
}

static void sendtours(Tour** tours, int n, int to) {
    MPI_Send(&n, 1, MPI_INT, to, MPI_TAG_TOUR_N, MPI_COMM_WORLD);
    for (int i = 0; i < n; i++) {
        MPI_Send(&tours[i]->cost, 1, MPI_INT, to, MPI_TAG_TOUR_COST,
            MPI_COMM_WORLD);
        int count = tours[i]->count;
        MPI_Send(&count, 1, MPI_INT, to, MPI_TAG_TOUR_COUNT, MPI_COMM_WORLD);
        for (int j = 0; j < count; j++) {
            MPI_Send(&tours[i]->cities[j], 1, MPI_INT, to, MPI_TAG_TOUR_CITIES,
                MPI_COMM_WORLD);
        }
    }
}

static Tour** receivetours(int* n) {
    MPI_Status s;

    MPI_Recv(n, 1, MPI_INT, MASTER, MPI_TAG_TOUR_N, MPI_COMM_WORLD, &s);

    Tour** tours = malloc(*n * sizeof(Tour*));
    for (int i = 0; i < *n; i++) {
        tours[i] = tour_new(ncities);
        MPI_Recv(&tours[i]->cost, 1, MPI_INT, MASTER, MPI_TAG_TOUR_COST,
            MPI_COMM_WORLD, &s);
        MPI_Recv(&tours[i]->count, 1, MPI_INT, MASTER, MPI_TAG_TOUR_COUNT,
            MPI_COMM_WORLD, &s);
        for (int j = 0; j < tours[i]->count; j++) {
            MPI_Recv(&tours[i]->cities[j], 1, MPI_INT, MASTER,
                MPI_TAG_TOUR_CITIES, MPI_COMM_WORLD, &s);
        }
    }

    return tours;
}

void master(void) {
    // starting tour
    Tour* beginning = tour_new(ncities);
    tour_add_city(beginning, START);

    // dividing tasks
    Stack** nodestacks = dividework(beginning, np);

    // sending tasks
    for (int i = 0; i < np - 1; i++) {
        Stack* stack = nodestacks[i];
        int ntours = stack->size;
        Tour** tours = malloc(ntours * sizeof(Tour*));

        int j = 0;
        while (!stack_empty(stack)) {
            tours[j++] = stack_pop(stack);
        }

        sendtours(tours, ntours, i + 1);

        // freeing memory
        stack_free(nodestacks[i]);
        for (int k = 0; k < ntours; k++) {
            tour_free(tours[k]);
        }
    }

    free(nodestacks);
}

void worker(void) {
    int ntours;
    Tour** tours = receivetours(&ntours);

    // global stack
    global_stack = stack_new(ncities * ncities * ncities);

    // pthreads
    pthread_mutex_init(&best_tour_mutex, NULL);
    pthread_mutex_init(&global_stack_mutex, NULL);
    pthread_cond_init(&global_stack_full, NULL);
    pthread_cond_init(&global_stack_empty, NULL);

    // best tour
    best = tour_new(ncities);
    best->cost = INFINITY;

    // dividing tasks
    Stack** stacks = dividework(tours[0], nthreads);
    int j = 0;
    for (int i = 1; i < ntours; i++) {
        if (j == nthreads) {
            j = 0;
        }
        stack_push(stacks[nthreads - 1 - j++], tours[i]);
    }

    // running the threads
    pthread_t tojoin[nthreads];
    for (int i = 0; i < nthreads; i++) {
        tojoin[i] = newthread(i, stacks[i]);
    }

    // waiting for threads to finish
    for (int i = 0; i < nthreads; i++) {
        pthread_join(tojoin[i], NULL);
    }

    tour_print(best);

    tour_free(best);
    graph_free();
    pthread_mutex_destroy(&best_tour_mutex);
    pthread_mutex_destroy(&global_stack_mutex);
    pthread_cond_destroy(&global_stack_full);
    pthread_cond_destroy(&global_stack_empty);
    stack_free(global_stack);
}

int main(int argc, char** argv) {
    // arguments (graph path & number of threads)
    graph_load(argv[1]);
    nthreads = atoi(argv[2]);
    ncities = graph_size;

    // MPI
    mpiassert(MPI_Init(&argc, &argv));
    mpiassert(MPI_Comm_size(MPI_COMM_WORLD, &np));
    mpiassert(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

    #if DEBUG
        printf("INITIALIZED RANK %d\n", rank);
    #endif

    if (rank == MASTER) {
        master();
    } else {
        worker();
    }
    
    #if DEBUG
        printf("FINILIZED RANK %d\n", rank);
    #endif

    mpiassert(MPI_Finalize());
    return 0;
}

// ==================================================
//
//  auxiliary
//
// ==================================================

static void mpiassert(int result) {
    assert(result == MPI_SUCCESS);
}

static void updatebest(Tour* t) {
    pthread_mutex_lock(&best_tour_mutex);
    tour_free(best);
    best = tour_copy(t);
    #if DEBUG
        printf("UPDATED BEST ");
        tour_print(t);
    #endif
    pthread_mutex_unlock(&best_tour_mutex);
}

static void globalpush(Stack* stack, Tour* tour) {
    pthread_mutex_lock(&global_stack_mutex);
    while (stack_full(global_stack)) {
        pthread_cond_wait(&global_stack_full, &global_stack_mutex);
    }
    int half = (stack->size / 2) - 1;
    for (int i = 0; i < half; i++) {
        stack_push(global_stack, stack_pop(stack));
    }
    stack_push(global_stack, tour);
    pthread_cond_broadcast(&global_stack_empty);
    pthread_mutex_unlock(&global_stack_mutex);
}

static Tour* globalpop(Stack* stack) {
    pthread_mutex_lock(&global_stack_mutex);
    while (stack_empty(global_stack)) {
        if (++waiting_threads == nthreads) {
            pthread_cond_broadcast(&global_stack_empty);
            pthread_mutex_unlock(&global_stack_mutex);
            return NULL;
        }
        pthread_cond_wait(&global_stack_empty, &global_stack_mutex);
        waiting_threads--;
    }
    Tour* tour = stack_pop(global_stack);
    int half = (stack->size / 2) - 1;
    for (int i = 0; i < half; i++) {
        stack_push(stack, stack_pop(global_stack));
    }
    pthread_cond_broadcast(&global_stack_full);
    pthread_mutex_unlock(&global_stack_mutex);
    return tour;
}

static void pushwork(Stack* stack, Tour* tour) {
    if (stack_full(stack)) {
        globalpush(stack, tour);
    } else {
        stack_push(stack, tour);
    }
}

static Tour* popwork(Stack* stack) {
    return (stack_empty(stack)) ? globalpop(stack) : stack_pop(stack);
}


static pthread_t newthread(int id, Stack* stack) {
    int* idpointer = malloc(sizeof(int));
    assert(idpointer);
    *idpointer = id;

    void** arguments = malloc(2 * sizeof(void*));
    assert(arguments);

    arguments[0] = idpointer;
    arguments[1] = stack;

    pthread_t thread;
    pthread_create(&thread, NULL, threadfindbest, arguments);
    return thread;
}

// TODO: move
static Stack** dividework(Tour* tour, int n) {
    // dividing tasks
    Stack* stack = BFS(tour, n);
    Stack** divided = malloc(n * sizeof(Stack*));
    for (int i = 0; i < n; i++) {
        divided[i] = stack_new(DEFAULT_STACK_SIZE);
    }
    int i = 0;
    while (!stack_empty(stack)) {
        stack_push(divided[i], stack_pop(stack));
        i = (i < n - 1) ? i + 1 : 0;
    }
    stack_free(stack);
    return divided;
}
