#if !defined(stack_h)
#define stack_h

#include "tour.h"

typedef struct Stack {
   int size;        // size of the stack
   int capacity;    // capacity of the stack
   Tour** array;    // array of tours (the tasks)
} Stack;

extern Stack* stack_new(int capacity);
extern void   stack_free(Stack*);
extern void   stack_push(Stack*, Tour*);
extern Tour*  stack_pop(Stack*);
extern void   stack_transfer(Stack* from, Stack* to);
extern void   stack_print(Stack*);

#define stack_push_copy(s, t) (stack_push(s, tour_copy(t)))

#define stack_empty(s) (s->size == 0)
#define stack_full(s)  (s->size == s->capacity)

#endif
