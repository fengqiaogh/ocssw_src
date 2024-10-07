#include <genutils.h>
#include <stdlib.h>

/* -----------------------------------------------------------------------
 * allocate and zero out memory.  Print decent error message.
 * ----------------------------------------------------------------------- */
void* allocateMemory(size_t numBytes, const char* name) {
    void* ptr = calloc(numBytes, 1);
    if (ptr == NULL) {
        printf("-E- %s line %d : error allocating memory for %s.\n",
                __FILE__, __LINE__, name);
        exit(1);
    }
    return ptr;
}
