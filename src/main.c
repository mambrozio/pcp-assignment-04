#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define ERROR(msg) do { \
    fprintf(stderr, "\terror: " msg "\n"); \
    exit(1); \
} while (0); \

#define IS_MAIN_THREAD (rank == 0)

// graph
typedef int* Graph;
#define GRAPH(i, j) (graph[i*n+j])
#define GRAPH_SET(i, j, v) do { GRAPH(i, j) = v } while (0);

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

int main(int argc, char** argv) {
    // graph
    const char* path = argv[1];
    loadgraph(path);
    printgraph();

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
    graph = malloc(n*n*sizeof(int));

    // weights
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
