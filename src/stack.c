#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "stack.h"
#include "tour.h"

Stack* stack_new(int capacity) {
    Stack* stack = malloc(sizeof(Stack));
    assert(stack);
    stack->size = 0;
    stack->capacity = capacity;
    stack->array = malloc(capacity * sizeof(Tour));
    assert(stack->array);
    return stack;
}

void stack_free(Stack* stack) {
    assert(stack_empty(stack));
    free(stack->array);
    free(stack);
}

void stack_push(Stack* stack, Tour* tour) {
    assert(!stack_full(stack));
    stack->array[stack->size++] = tour;
}

Tour* stack_pop(Stack* stack) {
    return (stack_empty(stack)) ? NULL : stack->array[--stack->size];
}

void stack_transfer(Stack* from, Stack* to) {
    for (int i = 0; i < from->size; i++) {
        stack_push(to, from->array[i]);
    }
    from->size = 0;
}

void stack_print(Stack* stack) {
    printf("Stack = {\n");
    printf("\tSize = %d\n", stack->size);
    printf("\tCapacity = %d\n", stack->capacity);
    printf("\tTours = [[\n");

    for (int i = stack->size - 1; i >= 0; i--) {
        printf("\t\t");
        tour_print(stack->array[i]);
    }

    printf("\t]]\n}\n");
}
