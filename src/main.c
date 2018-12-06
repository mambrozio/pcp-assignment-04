#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#define ERROR(msg) do { \
    fprintf(stderr, "\terror: " msg "\n"); \
    exit(1); \
} while (0); \

#define IS_MAIN_THREAD (rank == 0)

// graph
typedef int* Graph;
#define GRAPH(i, j) (graph[i*n+j])

// constants
const int NO_CITY = -1;

// important globals
int n;          // number of cities
Graph graph;    // the cities' graph
int np;         // MPI number of processes
int rank;       // MPI rank

// auxiliary functions
// static void mpiassert(int result);
static void loadgraph(const char* path);
static void printgraph(void);

// ==================================================
//
//  data structures
//
// ==================================================

typedef struct Tour {
   int* cities; // cities in partial tour
   int count; // number of cities in partial tour
   int cost; // cost of partial tour
} Tour;

typedef struct Stack {
   Tour** list; // list of tours (the tasks)
   int list_sz; // TODO
   int list_alloc; // TODO
} Stack;

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

Tour* best_tour = NULL;

static bool feasible(Tour* tour, int city) {
    // if it can lead to a least cost tour
    int lastcity = tour->cities[tour->count - 1];
    int newcost = tour->cost + GRAPH(lastcity, city);
    if (newcost > best_tour->cost) {
        return false;
    }

    // if the vertex has already been visited
    for (int i = 0; i < tour->count; i++) {
        if (tour->cities[i] == city) {
            return false;
        }
    }

    return true;
}

void DFS(Tour* tour) {
    if (tour->count == n && tour->cost < best_tour->cost) {
        best_tour = tour;
    } else {
        int lastcity = tour->cities[tour->count - 1];

        // for each neighboring city
        for (int newcity = 0; newcity < n; newcity++) {
            if (GRAPH(lastcity, newcity)) {
                if (feasible(tour, newcity)) {
                    tour->cities[tour->count++] = newcity; // add new city
                }
                DFS(tour);
                tour->cities[--tour->count] = NO_CITY; // remove last city
            }
        }
    }
}

int main(int argc, char** argv) {
    // graph
    const char* path = argv[1];
    loadgraph(path);
    printgraph();

    // TODO
    int start = 0;
    // int finish = 2;

    Tour* tour = malloc(sizeof(Tour));
    tour->cities = malloc(n * sizeof(int));
    tour->cities[0] = start;
    for (int i = 1; i < n; i++) {
        tour->cities[i] = NO_CITY;
    }
    tour->count = 1;
    tour->cost = 0;

    best_tour = tour;
    DFS(tour);

    printf("best_tour =\n");
    printf("cities =\n");
    for (int i = 0; i < n; i++) {
        printf("%d ", tour->cities[i]);
    }
    printf("count = %d\n", best_tour->count);
    printf("cost = %d\n", best_tour->cost);

    free(tour->cities);
    free(tour);
    free(graph);
    return 0;
}

// ==================================================
//
//  auxiliary
//
// ==================================================

// static void mpiassert(int result) {
//     assert(result == MPI_SUCCESS);
// }

// loads the global variables "n" and "graph"
static void loadgraph(const char* path) {
    // file
    FILE* file = fopen(path, "r");
    if (!file) {
        ERROR("can't open graph file");
    }

    // number of vertices
    fscanf(file, "%d", &n);
    if (n <= 0) {
        ERROR("number of vertices in the graph must be positive");
    }

    // weights
    graph = malloc(n*n*sizeof(int));
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
