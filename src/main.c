#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pthread.h>

#define DEBUG 0

// ==================================================
//
//  data structures & typedefs
//
// ==================================================

typedef int* Graph;

typedef struct Tour {
   int* cities; // cities in partial tour
   int count; // number of cities in partial tour
   int cost; // cost of partial tour
} Tour;

typedef struct Stack {
   int size;        // size of the stack
   int capacity;    // capacity of the stack (size^3)
   Tour** array;    // array of tours (the tasks)
} Stack;

// ==================================================
//
//  constants & important globals
//
// ==================================================

// constants
const int MASTER = 0;
const int START = 0;
const int NO_CITY = -1;
const int INFINITY = INT_MAX;

// MPI
int np;   // MPI number of processes
int rank; // MPI rank
MPI_Datatype mpi_stack;

// pthreads
pthread_mutex_t best_tour_mutex;
pthread_mutex_t global_stack_mutex;
pthread_cond_t global_stack_full;
pthread_cond_t global_stack_empty;
Stack* global_stack;

// important globals
Tour* best = NULL;       // best tour
int ncities;             // number of cities
Graph graph;             // the cities' graph
int waiting_threads = 0; // number of threads waiting on global_stack_empty
int nthreads;            // number of threads

// ==================================================
//
//  auxiliary macros
//
// ==================================================

#define ERROR(msg) do { \
    fprintf(stderr, "\terror: " msg "\n"); \
    exit(1); \
} while (0); \

#define GRAPH(i, j) (graph[i*ncities+j])

#define FIRST_CITY(t) (t->cities[0])
#define LAST_CITY(t) (t->cities[t->count - 1])

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

static void loadgraph(const char* path);
static void printgraph(void);

static Tour* newtour(void);
static Tour* copytour(Tour*);
static void freetour(Tour*);
static void addcity(Tour* t, int city);
static void removelastcity(Tour*);
static bool visited(Tour* tour, int city);
static void printtour(Tour*);

static void updatebest(Tour*);

static Stack* newstack(void);
static void freestack(Stack*);
static void push(Stack*, Tour*);
static Tour* pop(Stack*);
#define pushcopy(s, t) (push(s, copytour(t)))
static bool empty(Stack*);
static bool full(Stack*);
static void printstack(Stack*);

static void pushwork(Stack*, Tour*);
static Tour* popwork(Stack*);
#define pushcopywork(s, t) (pushwork(s, copytour(t)))

static pthread_t newthread(int id, Stack* tasks);

static Stack** dividework(Tour* tour, int n);

// ==================================================
//
//  functions
//
// ==================================================

Stack* BFS(Tour* beginning, int nt) {
    Stack* stack = newstack();
    pushcopy(stack, beginning);

    while (stack->size < nt) {
        Tour* tour = pop(stack);
        if (!tour) {
            goto END;
        }

        for (int neighbor = 0; neighbor < ncities; neighbor++) {
            if (!visited(tour, neighbor)) {
                addcity(tour, neighbor);
                pushcopy(stack, tour);
                removelastcity(tour);
            }
        }

        freetour(tour);
    }

    END: {
        return stack;
    }
}

static bool feasible(Tour* tour, int city) {
    // if the city has already been visited
    if (visited(tour, city)) {
        return false;
    }

    // if it can lead to a least cost tour
    addcity(tour, city);
    int newcost = tour->cost;
    removelastcity(tour);
    if (newcost > best->cost) {
        return false;
    }

    return true;
}

void findbest(int id, Stack* stack) {
    #if DEBUG
        printf("STARTING THREAD %d\n", id);
    #endif

    while (!empty(stack)) {
        Tour* tour = popwork(stack);
        if (!tour) {
            goto END;
        }

        if (tour->count == ncities && tour->cost < best->cost) {
            updatebest(tour);
        } else {
            for (int neighbor = 0; neighbor < ncities; neighbor++) {
                if (feasible(tour, neighbor)) {
                    addcity(tour, neighbor);
                    pushcopywork(stack, tour);
                    removelastcity(tour);
                }
            }
        }

        freetour(tour);
    }

    END: {
        freestack(stack);
    }
}

void* threadfindbest(void* arguments) {
    int id = *((int*)((void**)arguments)[0]);
    Stack* tasks = (Stack*)((void**)arguments)[1];
    findbest(id, tasks);
    return NULL;
}

static void sendtours(Tour** tours, int n, int to) {
    MPI_Send(&n, 1, MPI_INT, to, MPI_TAG_TOUR_N, MPI_COMM_WORLD);
    for (int i = 0; i < n; i++) {
        MPI_Send(&tours[i]->cost, 1, MPI_INT, to, MPI_TAG_TOUR_COST, MPI_COMM_WORLD);
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
        tours[i] = newtour();
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
    Tour* beginning = newtour();
    addcity(beginning, START);

    // dividing tasks
    Stack** nodetasks = dividework(beginning, np);

    // sending tasks
    for (int i = 0; i < np - 1; i++) {
        Stack* stack = nodetasks[i];
        int ntours = stack->size;
        Tour** tours = malloc(ntours * sizeof(Tour*));

        int j = 0;
        while (!empty(stack)) {
            tours[j++] = pop(stack);
        }

        sendtours(tours, ntours, i + 1);

        // freeing memory
        freestack(nodetasks[i]);
        for (int k = 0; k < ntours; k++) {
            freetour(tours[k]);
        }
    }

    free(nodetasks);
}

void worker(void) {
    int ntours;
    Tour** tours = receivetours(&ntours);

    // FIXME
    // Receive one or more Tours from MPI

    // pthreads
    pthread_mutex_init(&best_tour_mutex, NULL);
    pthread_mutex_init(&global_stack_mutex, NULL);
    pthread_cond_init(&global_stack_full, NULL);
    pthread_cond_init(&global_stack_empty, NULL);

    // best tour
    best = newtour();
    best->cost = INFINITY;

    // dividing tasks
    Stack** stacks = dividework(tours[0], nthreads);
    int j = 0;
    for (int i = 1; i < ntours; i++) {
        if (j == nthreads) {
            j = 0;
        }
        push(stacks[nthreads - 1 - j++], tours[i]);
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

    printtour(best);

    freetour(best);
    free(graph);
    pthread_mutex_destroy(&best_tour_mutex);
    pthread_mutex_destroy(&global_stack_mutex);
    pthread_cond_destroy(&global_stack_full);
    pthread_cond_destroy(&global_stack_empty);
}

int main(int argc, char** argv) {
    // arguments (graph path & number of threads)
    loadgraph(argv[1]);
    nthreads = atoi(argv[2]);

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

// loads the global variables "n" and "graph"
static void loadgraph(const char* path) {
    // file
    FILE* file = fopen(path, "r");
    if (!file) {
        ERROR("can't open graph file");
    }

    // number of vertices
    fscanf(file, "%d", &ncities);
    if (ncities <= 1) {
        ERROR("number of vertices in the graph must greater than 1");
    }

    // weights
    graph = malloc(ncities*ncities*sizeof(int));
    assert(graph);
    for (int i = 0; i < ncities; i++) {
        for (int j = 0; j < ncities; j++) {
            fscanf(file, "%d", &graph[i*ncities + j]);
            if (i == j && GRAPH(i, j) != 0) {
                ERROR("diagonal entries must be zero");
            }
            if (i != j && GRAPH(i, j) <= 0) {
                ERROR("off-diagonal entries must be positive");
            }
        }
    }

    fclose(file);
}

// prints the global variable "graph"
static void printgraph(void) {
    printf("Number of cities = %d\n", ncities);
    printf("Matrix = \n");
    for (int i = 0; i < ncities; i++) {
        for (int j = 0; j < ncities; j++) {
            printf("%2d ", GRAPH(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

static Tour* newtour(void) {
    Tour* t = malloc(sizeof(Tour));
    assert(t);
    t->cities = malloc(ncities * sizeof(int));
    assert(t->cities);
    for (int i = 0; i < ncities; i++) {
        t->cities[i] = NO_CITY;
    }
    t->count = 0;
    t->cost = 0;
    return t;
}

static Tour* copytour(Tour* t1) {
    Tour* t2 = newtour();
    for (int i = 0; i < t1->count; i++) {
        t2->cities[i] = t1->cities[i];
    }
    for (int i = t1->count; i < ncities; i++) {
        t2->cities[i] = NO_CITY;
    }
    t2->count = t1->count;
    t2->cost = t1->cost;
    return t2;
}

static void freetour(Tour* t) {
    free(t->cities);
    free(t);
}

static void addcity(Tour* t, int city) {
    if (t->count == ncities) {
        ERROR("invalid number of cities in tour (add)");
    }

    if (t->count > 0) {
        t->cost += GRAPH(LAST_CITY(t), city);
    }
    t->cities[t->count++] = city;

    if (t->count == ncities) {
        t->cost += GRAPH(LAST_CITY(t), FIRST_CITY(t));
    }
}

static void removelastcity(Tour* t) {
    if (t->count <= 1) {
        ERROR("invalid number of cities in tour (remove)");
    }

    if (t->count == ncities) {
        t->cost -= GRAPH(LAST_CITY(t), FIRST_CITY(t));
    }

    if (t->count == 2) {
        t->cost = 0;
    } else {
        t->cost -= GRAPH(t->cities[t->count - 2], LAST_CITY(t));
    }
    
    t->cities[--t->count] = NO_CITY;
}

static bool visited(Tour* tour, int city) {
    for (int i = 0; i < tour->count; i++) {
        if (tour->cities[i] == city) {
            return true;
        }
    }
    return false;
}

static void printtour(Tour* tour) {
    printf("Tour = {Cost = ");
    if (tour->cost == INFINITY) {
        printf("INFINITY");
    } else {
        printf("%d", tour->cost);
    }
    printf(", Cities (%d) = {", tour->count);
    for (int i = 0; i < tour->count; i++) {
        printf("%d", tour->cities[i]);
        if (i == tour->count - 1) {
            break;
        }
        printf(", ");
    }
    printf("}}\n");
}

static void updatebest(Tour* t) {
    pthread_mutex_lock(&best_tour_mutex);
    freetour(best);
    best = copytour(t);
    #if DEBUG
        printf("UPDATED BEST ");
        printtour(t);
    #endif
    pthread_mutex_unlock(&best_tour_mutex);
}

static Stack* newstack(void) {
    Stack* s = malloc(sizeof(Stack));
    assert(s);
    s->size = 0;
    s->capacity = ncities*ncities;
    s->array = malloc(s->capacity * sizeof(Tour));
    assert(s->array);
    return s;
}

static void freestack(Stack* s) {
    assert(empty(s));
    free(s->array);
    free(s);
}

static void push(Stack* stack, Tour* tour) {
    assert(!full(stack));
    stack->array[stack->size++] = tour;
    #if DEBUG
        printf("PUSH ");
        printtour(tour);
    #endif
}

static Tour* pop(Stack* stack) {
    if (empty(stack)) {
        #if DEBUG
            printf("POP NULL\n");
        #endif
        return NULL;
    }

    Tour* tour = stack->array[--stack->size];
    #if DEBUG
        printf("POP ");
        printtour(tour);
    #endif
    return tour;
}

static bool empty(Stack* s) {
    return s->size == 0;
}

static bool full(Stack* s) {
    return s->size == s->capacity;
}

static void printstack(Stack* stack) {
    printf("Stack = {\n");
    printf("\tSize = %d\n", stack->size);
    printf("\tCapacity = %d\n", stack->capacity);
    printf("\tTours = [[\n");

    for (int i = stack->size - 1; i >= 0; i--) {
        printf("\t\t");
        printtour(stack->array[i]);
    }

    printf("\t]]\n}\n");
}

static void pushwork(Stack* tasks, Tour* tour) {
    #if DEBUG
        printf("GLOBAL\n");
    #endif
    if (full(tasks)) {
        pthread_mutex_lock(&global_stack_mutex);
        while (full(global_stack)) {
            pthread_cond_wait(&global_stack_full, &global_stack_mutex);
        }
        for (int i = 0; i < (tasks->size / 2) - 1; i++) {
            push(global_stack, pop(tasks));
        }
        push(global_stack, tour);
        pthread_mutex_unlock(&global_stack_mutex);
    } else {
        push(tasks, tour);
    }
}

static Tour* popwork(Stack* tasks) {
    #if DEBUG
        printf("GLOBAL\n");
    #endif
    Tour* tour;
    if (empty(tasks)) {
        pthread_mutex_lock(&global_stack_mutex);
        while (empty(global_stack)) {
            if (++waiting_threads == nthreads) {
                return NULL;
            }
            pthread_cond_wait(&global_stack_empty, &global_stack_mutex);
            waiting_threads--;
        }
        tour = pop(global_stack);
        for (int i = 0; i < (global_stack->size / 2) - 1; i++) {
            push(tasks, pop(global_stack));
        }
        pthread_mutex_unlock(&global_stack_mutex);
    } else {
        tour = pop(tasks);
    }
    return tour;
}

static pthread_t newthread(int id, Stack* tasks) {
    int* idpointer = malloc(sizeof(int));
    assert(idpointer);
    *idpointer = id;

    void** arguments = malloc(2 * sizeof(void*));
    assert(arguments);

    arguments[0] = idpointer;
    arguments[1] = tasks;

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
        divided[i] = newstack();
    }
    int i = 0;
    while (!empty(stack)) {
        push(divided[i], pop(stack));
        i = (i < n - 1) ? i + 1 : 0;
    }
    freestack(stack);
    return divided;
}
