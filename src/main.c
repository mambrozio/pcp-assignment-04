#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pthread.h>

#include "graph.h"
#include "stack.h"
#include "tour.h"

#define DEBUG 1

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

#define MPI_TAG_TOUR_N       11
#define MPI_TAG_TOUR_COST    12
#define MPI_TAG_TOUR_COUNT   13
#define MPI_TAG_TOUR_CITIES  14
#define MPI_TAG_SENDING_TOUR 15
#define MPI_TAG_DONE         16

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

// send & receive

#define send_int(data, destination, tag) \
    (MPI_Send(&data, 1, MPI_INT, destination, tag, MPI_COMM_WORLD)) \

#define receive_int(address, source, tag, status) \
    (MPI_Recv(address, 1, MPI_INT, source, tag, MPI_COMM_WORLD, status)) \

static void send_tour(Tour*, int destination);
static Tour* receive_tour(int source);
static void send_tours(Tour**, int n, int destination);
static Tour** receive_tours(int source, int* n);

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

        send_tours(tours, ntours, i + 1);

        // freeing memory
        stack_free(nodestacks[i]);
        for (int k = 0; k < ntours; k++) {
            tour_free(tours[k]);
        }
    }
    free(nodestacks);

    // best tour
    best = tour_new(ncities);
    best->cost = INFINITY;

    MPI_Status status;
    for (int all = 1; all < np;) {
        int source;
        receive_int(&source, MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
        switch (status.MPI_TAG) {
            case MPI_TAG_SENDING_TOUR: {
                Tour* received = receive_tour(source);
                assert(received);
                assert(received->count == ncities);
                if (received->cost < best->cost) {
                    tour_free(best);
                    best = received;
                } else {
                    tour_free(received);
                }
                send_tour(best, source);
                break;
            }
            case MPI_TAG_DONE: {
                all++;
                break;
            }
        default:
            assert(false);
        }
    }
}

void worker(void) {
    int ntours;
    Tour** tours = receive_tours(MASTER, &ntours);

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

    int data = 1;
    send_int(data, MASTER, MPI_TAG_DONE);

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
    double t1 = 0.0;
    double t2 = 0.0;

    // MPI
    mpiassert(MPI_Init(&argc, &argv));
    mpiassert(MPI_Comm_size(MPI_COMM_WORLD, &np));
    mpiassert(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

    #if DEBUG
        printf("INITIALIZED RANK %d\n", rank);
    #endif

    if (rank == MASTER) {
        t1 = MPI_Wtime();
        master();
        t2 = MPI_Wtime();
    } else {
        worker();
    }

    #if DEBUG
        printf("FINILIZED RANK %d\n", rank);
    #endif

    if (rank == MASTER) {
        printf("Elapsed time: %.25f\n", t2 - t1);
        printf("Best route");
        tour_print(best);
    }

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

static void lock(pthread_mutex_t* mutex, const char* id) {
    #if DEBUG
        printf("TRYING TO ACQUIRE LOCK <%s>\n", id);
    #endif
    pthread_mutex_lock(mutex);
    #if DEBUG
        printf("ACQUIRED LOCK <%s>\n", id);
    #endif
}

static void unlock(pthread_mutex_t* mutex, const char* id) {
    #if DEBUG
        printf("TRYING TO REALEASE LOCK <%s>\n", id);
    #endif
    pthread_mutex_unlock(mutex);
    #if DEBUG
        printf("RELEASED LOCK <%s>\n", id);
    #endif
}

// important: does not free tour
static void updatebest(Tour* tour) {
    lock(&best_tour_mutex, "updatebest");

    // sending a local best to master and potentially receiving a global best
    send_int(rank, MASTER, MPI_TAG_SENDING_TOUR);
    send_tour(tour, MASTER);
    tour = receive_tour(MASTER);

    // replacing the best tour
    tour_free(best);
    best = tour_copy(tour);

    #if DEBUG
        printf("UPDATED BEST ");
        tour_print(tour);
    #endif
    unlock(&best_tour_mutex, "updatebest");
}

static void globalpush(Stack* stack, Tour* tour) {
    lock(&global_stack_mutex, "globalpush");
    while (stack_full(global_stack)) {
        #if DEBUG
            printf("-- RANK %d WAITING (global_stack_full)\n", rank);
        #endif
        pthread_cond_wait(&global_stack_full, &global_stack_mutex);
        #if DEBUG
            printf("-- RANK %d DONE (global_stack_full)\n", rank);
        #endif
    }
    // int half = (stack->size / 2) - 1;
    // for (int i = 0; i < half; i++) {
        stack_push(global_stack, stack_pop(stack));
    // }
    stack_push(global_stack, tour);
    pthread_cond_broadcast(&global_stack_empty);
    unlock(&global_stack_mutex, "globalpush");
}

static Tour* globalpop(Stack* stack) {
    lock(&global_stack_mutex, "globalpop");
    while (stack_empty(global_stack)) {
        if (++waiting_threads == nthreads) {
            pthread_cond_broadcast(&global_stack_empty);
            unlock(&global_stack_mutex, "globalpop-done");
            return NULL;
        }
        #if DEBUG
            printf("-- RANK %d WAITING (global_stack_empty)\n", rank);
        #endif
        pthread_cond_wait(&global_stack_empty, &global_stack_mutex);
        #if DEBUG
            printf("-- RANK %d DONE (global_stack_empty)\n", rank);
        #endif
        waiting_threads--;
    }
    Tour* tour = stack_pop(global_stack);
    // int half = (stack->size / 2) - 1;
    // for (int i = 0; i < half; i++) {
        stack_push(stack, stack_pop(global_stack));
    // }
    pthread_cond_broadcast(&global_stack_full);
    unlock(&global_stack_mutex, "globalpop");
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

static void send_tour(Tour* tour, int to) {
    send_int(tour->cost, to, MPI_TAG_TOUR_COST);
    send_int(tour->count, to, MPI_TAG_TOUR_COUNT);
    for (int i = 0; i < tour->count; i++) {
        send_int(tour->cities[i], to, MPI_TAG_TOUR_CITIES);
    }
}

static Tour* receive_tour(int source) {
    MPI_Status status;
    Tour* tour = tour_new(ncities);
    receive_int(&tour->cost, source, MPI_TAG_TOUR_COST, &status);
    receive_int(&tour->count, source, MPI_TAG_TOUR_COUNT, &status);
    for (int i = 0; i < tour->count; i++) {
        receive_int(&tour->cities[i], source, MPI_TAG_TOUR_CITIES, &status);
    }
    return tour;
}

static void send_tours(Tour** tours, int n, int to) {
    send_int(n, to, MPI_TAG_TOUR_N);
    for (int i = 0; i < n; i++) {
        send_tour(tours[i], to);
    }
}

static Tour** receive_tours(int from, int* n) {
    MPI_Status status;
    receive_int(n, from, MPI_TAG_TOUR_N, &status);
    Tour** tours = malloc(*n * sizeof(Tour*));
    for (int i = 0; i < *n; i++) {
        tours[i] = receive_tour(from);
    }
    return tours;
}
