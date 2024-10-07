/*
      get red, green and blue look up table from an ascii file
 */
#include <genutils.h>

#include <stdio.h>

int getlut_file(char *lut_file, short *rlut, short *glut, short * blut) {
    int i;
    FILE *file_ptr;

    if ((file_ptr = fopen(lut_file, "r")) == NULL) {
        printf("[getlut_file] error opening the lut file\n");
        return (1);
    }

    /* read the luts from file */

    for (i = 0; i < 256; i++) {
        fscanf(file_ptr, "%hd%hd%hd\n", &rlut[i], &glut[i], &blut[i]);
    }

    fclose(file_ptr);

    return (0);
}

