#include <genutils.h>
#include <stdlib.h>

size_t fread_swap(int little_endian, void *ptr, size_t size,
        size_t nmemb, FILE *stream) {
    size_t result;

    result = fread(ptr, size, nmemb, stream);
    if (little_endian != endianess()) {
        swapc_bytes(ptr, size, nmemb);
    }
    return result;
}

size_t fwrite_swap(int little_endian, const void *ptr, size_t size,
        size_t nmemb, FILE *stream) {
    size_t result;
    void *ptr2;

    if (little_endian != endianess()) {
        ptr2 = malloc(size * nmemb);
        swapc_bytes2(ptr, ptr2, size, nmemb);
        result = fwrite(ptr2, size, nmemb, stream);
        free(ptr2);
    } else {
        result = fwrite(ptr, size, nmemb, stream);
    }
    return result;
}
