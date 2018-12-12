#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "graph.h"

Graph graph;
int graph_size = -1;

void graph_load(const char* path) {
    FILE* file = fopen(path, "r");
    assert(file);

    fscanf(file, "%d", &graph_size);
    assert(graph_size > 1);

    graph = malloc(graph_size * graph_size * sizeof(int));
    assert(graph);

    for (int i = 0; i < graph_size; i++) {
        for (int j = 0; j < graph_size; j++) {
            fscanf(file, "%d", &graph[i * graph_size + j]);
            if (i == j) {
                assert(path(i, j) == 0);
            } else { 
                assert(path(i, j) > 0);
            }
        }
    }

    fclose(file);
}

void graph_print(void) {
    printf("Number of cities = %d\n", graph_size);
    printf("Matrix = \n");
    for (int i = 0; i < graph_size; i++) {
        for (int j = 0; j < graph_size; j++) {
            printf("%2d ", path(i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void graph_free(void) {
    free(graph);
}
