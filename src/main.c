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

// TODO
// typedef void*(ThreadFunction)(void*);

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
const int START = 0;
const int NO_CITY = -1;
const int INFINITY = INT_MAX;

// MPI
int np;   // MPI number of processes
int rank; // MPI rank

// pthreads
pthread_mutex_t best_tour_mutex;

// important globals
Tour* best = NULL;  // best tour
int n;              // number of cities
Graph graph;        // the cities' graph

// ==================================================
//
//  auxiliary macros
//
// ==================================================

#define ERROR(msg) do { \
    fprintf(stderr, "\terror: " msg "\n"); \
    exit(1); \
} while (0); \

#define IS_MAIN_PROCESS (rank == 0)

#define GRAPH(i, j) (graph[i*n+j])

#define FIRST_CITY(t) (t->cities[0])
#define LAST_CITY(t) (t->cities[t->count - 1])

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
static void printstack(Stack*);

static pthread_t newthread(int id, Stack* tasks);

// ==================================================
//
//  functions
//
// ==================================================

// // TODO
// int main(int argc, char** argv) {
//     mpiassert(MPI_Init(&argc, &argv));
//     mpiassert(MPI_Comm_size(MPI_COMM_WORLD, &np));
//     mpiassert(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

//     printf("--- rank: %d\n", rank);

//     if (IS_MAIN_THREAD) {
//         mpiassert(MPI_Finalize());
//     }

//     return 0;
// }

Stack* BFS(Tour* beginning, int nt) {
    // stack
    Stack* stack = newstack();
    pushcopy(stack, beginning); // the tour that visits only the home town

    while (stack->size < nt) {
        Tour* tour = pop(stack);

        // for each neighboring city
        for (int neighbor = 0; neighbor < n; neighbor++) {
            if (!visited(tour, neighbor)) {
                addcity(tour, neighbor);
                pushcopy(stack, tour);
                removelastcity(tour);
            }
        }

        freetour(tour);
    }

    return stack;
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

    // TODO: access to best->cost without locking will probably lead to errors

    while (!empty(stack)) {
        Tour* tour = pop(stack);
        if (tour->count == n && tour->cost < best->cost) {
            updatebest(tour);
        } else {
            // for each neighboring city
            for (int neighbor = n - 1; neighbor >= 0; neighbor--) {
                if (feasible(tour, neighbor)) {
                    addcity(tour, neighbor);
                    pushcopy(stack, tour);
                    removelastcity(tour);
                }
            }
        }

        freetour(tour);
    }

    // TODO: steal work and goto top

    freestack(stack);
}

void* threadfindbest(void* arguments) {
    int id = *((int*)((void**)arguments)[0]);
    Stack* tasks = (Stack*)((void**)arguments)[1];
    findbest(id, tasks);
    return NULL;
}

int main(int argc, char** argv) {
    // arguments (graph path & number of threads)
    const char* path = argv[1];
    const int nthreads = atoi(argv[2]);
    loadgraph(path);

    // pthreads
    pthread_mutex_init(&best_tour_mutex, NULL);

    // starting tour
    Tour* beginning = newtour();
    addcity(beginning, START);

    // best tour
    best = copytour(beginning);
    best->cost = INFINITY;

    // dividing tasks
    Stack* tasks = BFS(beginning, nthreads);
    Stack* stacks[nthreads];
    for (int i = 0; i < nthreads; i++) {
        stacks[i] = newstack();
    }
    {
        int i = 0;
        while (!empty(tasks)) {
            push(stacks[i], pop(tasks));
            i = (i < nthreads - 1) ? i + 1 : 0;
        }
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

    printf("--- BEST TOUR ---\n");
    printtour(best);


    freetour(beginning);
    freetour(best);
    pthread_mutex_destroy(&best_tour_mutex);
    free(graph);
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
    fscanf(file, "%d", &n);
    if (n <= 1) {
        ERROR("number of vertices in the graph must greater than 1");
    }

    // weights
    graph = malloc(n*n*sizeof(int));
    assert(graph);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%d", &graph[i*n + j]);
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
    printf("Number of cities = %d\n", n);
    printf("Matrix = \n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%2d ", GRAPH(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

static Tour* newtour(void) {
    Tour* t = malloc(sizeof(Tour));
    assert(t);
    t->cities = malloc(n * sizeof(int));
    assert(t->cities);
    for (int i = 0; i < n; i++) {
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
    for (int i = t1->count; i < n; i++) {
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
    if (t->count == n) {
        ERROR("invalid number of cities in tour");
    }

    if (t->count > 0) {
        t->cost += GRAPH(LAST_CITY(t), city);
    }
    t->cities[t->count++] = city;

    if (t->count == n) {
        t->cost += GRAPH(LAST_CITY(t), FIRST_CITY(t));
    }
}

static void removelastcity(Tour* t) {
    if (t->count <= 1) {
        ERROR("invalid number of cities in tour");
    }

    if (t->count == n) {
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
    s->capacity = n*n*n;
    s->array = malloc(s->capacity * sizeof(Tour));
    assert(s->array);
    return s;
}

static void freestack(Stack* s) {
    assert(empty(s));
    free(s->array);
    free(s);
}

static void push(Stack* s, Tour* t) {
    s->array[s->size++] = t;
    assert(s->size != s->capacity);
    #if DEBUG
        printf("PUSH ");
        printtour(t);
    #endif
}

static Tour* pop(Stack* s) {
    assert(!empty(s));
    Tour* t = s->array[--s->size];
    #if DEBUG
        printf("POP ");
        printtour(t);
    #endif
    return t;
}

static bool empty(Stack* s) {
    return s->size == 0;
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
