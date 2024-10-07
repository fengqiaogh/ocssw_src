#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

void make_psurf(signed short *psurf, char *file, int Nlines, int Npixels, float psurf0) {
    int fd, idx;

    fd = open(file, O_RDONLY);
    read(fd, psurf, Nlines * Npixels * sizeof (signed short));
    close(fd);

    for (idx = 0; idx < Nlines * Npixels; idx++) {
        if (psurf[idx] >= 0) psurf[idx] = (signed short) (psurf0 - psurf[idx] / 10. + 0.5);
        else psurf[idx] = (signed short) (psurf0 + 0.5);
    }

}
