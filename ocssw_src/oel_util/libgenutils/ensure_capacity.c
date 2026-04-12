#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <limits.h>  

void* ensure_capacity(void *ptr, size_t *current_capacity, size_t required_capacity, size_t element_size) {
    if (required_capacity > *current_capacity) {
        if(required_capacity > (SIZE_MAX / element_size /2)) {
            fprintf(stderr, "-E-: %s:%d required capacity too large. Exiting.", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        size_t new_capacity = required_capacity * 2;
        void *new_ptr  = NULL;
        if (ptr != NULL) {
            new_ptr = realloc(ptr, new_capacity * element_size);
        } else {
            new_ptr = malloc(new_capacity * element_size);
        }
        if (new_ptr == NULL) {
            fprintf(stderr, "-E-: %s:%d allocation failed. Exiting.", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        *current_capacity = new_capacity;
        // printf("-I-: allocated %ld elements for required %ld elements\n", *current_capacity,
        //        required_capacity);
        return new_ptr;
    }
    return ptr;
}