#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "tour.h"

#define firstcity(t) (t->cities[0])
#define lastcity(t)  (t->cities[t->count - 1])

const int NO_CITY = -1;

const int INFINITY = INT_MAX;

Tour* tour_new(int capacity) {
    Tour* tour = malloc(sizeof(Tour));
    assert(tour);
    tour->cities = malloc(capacity * sizeof(int));
    assert(tour->cities);
    for (int i = 0; i < capacity; i++) {
        tour->cities[i] = NO_CITY;
    }
    tour->capacity = capacity;
    tour->count = 0;
    tour->cost = 0;
    return tour;
}

Tour* tour_copy(Tour* tour) {
    Tour* copy = tour_new(tour->capacity);
    for (int i = 0; i < tour->count; i++) {
        copy->cities[i] = tour->cities[i];
    }
    for (int i = tour->count; i < tour->capacity; i++) {
        copy->cities[i] = NO_CITY;
    }
    copy->count = tour->count;
    copy->cost = tour->cost;
    return copy;
}

void tour_free(Tour* tour) {
    free(tour->cities);
    free(tour);
}

void tour_print(Tour* tour) {
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

void tour_add_city(Tour* tour, City city) {
    assert(tour->count != tour->capacity);

    if (tour->count > 0) {
        tour->cost += path(lastcity(tour), city);
    }
    tour->cities[tour->count++] = city;

    if (tour->count == tour->capacity) {
        tour->cost += path(lastcity(tour), firstcity(tour));
    }
}

void tour_remove_last_city(Tour* tour) {
    assert(tour->count > 1);

    if (tour->count == tour->capacity) {
        tour->cost -= path(lastcity(tour), firstcity(tour));
    }

    if (tour->count == 2) {
        tour->cost = 0;
    } else {
        tour->cost -= path(tour->cities[tour->count - 2], lastcity(tour));
    }
    
    tour->cities[--tour->count] = NO_CITY;
}

bool tour_visits_city(Tour* tour, City city) {
    for (int i = 0; i < tour->count; i++) {
        if (tour->cities[i] == city) {
            return true;
        }
    }
    return false;
}
