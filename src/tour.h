#if !defined(tour_h)
#define tour_h

#include <stdbool.h>

#include "graph.h"

typedef int City;

typedef struct Tour {
   int* cities;  // cities in the tour
   int capacity; // capacity of the array of cities
   int count;    // number of cities in the tour
   int cost;     // cost of the tour (partial cost or full cost)
} Tour;

extern const int NO_CITY;
extern const int INFINITY;

extern Tour* tour_new(int capacity);
extern Tour* tour_copy(Tour*);
extern void  tour_free(Tour*);
extern void  tour_print(Tour*);

extern void tour_add_city(Tour*, City);
extern void tour_remove_last_city(Tour*);
extern bool tour_visits_city(Tour*, City);

#endif
