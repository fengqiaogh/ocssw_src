#include <string.h>
#include <stdlib.h>

int swapc_bytes(char *in, int nbyte, int ntime) {
    char *tmpbuf, *ptr;
    int i, j, k;

    tmpbuf = (char *) malloc(nbyte);
    ptr = in;

    for (j = 0; j < ntime; j++) {
        for (i = 0, k = nbyte - 1; i < nbyte; i++, k--)
            tmpbuf[i] = ptr[k];

        memcpy(ptr, tmpbuf, nbyte);

        ptr += nbyte;
    }
    free(tmpbuf);

    return 0;
}

int swapc_bytes2(const char *in, char *out, int nbyte, int ntime) {
    const char *ptr;
    char *ptr2;
    int i, j, k;

    ptr = in;
    ptr2 = out;

    for (j = 0; j < ntime; j++) {
        for (i = 0, k = nbyte - 1; i < nbyte; i++, k--)
            ptr2[i] = ptr[k];

        ptr += nbyte;
        ptr2 += nbyte;
    }

    return 0;
}
