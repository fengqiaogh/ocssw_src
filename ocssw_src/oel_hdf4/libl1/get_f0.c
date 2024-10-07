#include <stdio.h>
#include <stdlib.h>
#include "l1.h"

int get_f0_neckel(int32_t wl, int32_t width, float *f0) {
    static int firstCall = 1;
    static int nf0table = 600;
    static float f0_table[600];
    static int min_wl = 330;
    static int del_wl = 1;

    FILE *fp;
    char *filedir;
    char filename[FILENAME_MAX];
    char line [80];
    int32_t i, imin, imax;

    if (firstCall) {
        float tmp_wl;
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            return (0);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/f0_table_neckel.dat");

        if ((fp = fopen(filename, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",
                    __FILE__, __LINE__, filename);
            return (0);
        }

        for (i = 0; i < nf0table; i++) {
            if (fgets(line, 80, fp) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: error reading %s at line %d\n",
                        __FILE__, __LINE__, filename, i);
                return (0);
            }
            sscanf(line, "%f %f", &tmp_wl, &f0_table[i]);
        }

        firstCall = 0;
    }

    i = (wl - min_wl) / del_wl;
    imin = MAX(i - width / 2 / del_wl, 0);
    imax = MIN(i + width / 2 / del_wl, nf0table);

    if (i < 0 || i > nf0table) {
        fprintf(stderr,
                "-E- %s line %d: wavelength of %d outside F0 table range.\n",
                __FILE__, __LINE__, wl);
        return (0);
    }

    *f0 = 0.0;
    for (i = imin; i <= imax; i++)
        *f0 += f0_table[i];

    *f0 /= (imax - imin + 1);

    return (1);
}

int get_f0_neckel_(int32_t *wl, int32_t *width, float *f0) {
    return (get_f0_neckel(*wl, *width, f0));
}

int get_f0_thuillier(int32_t wl, int32_t width, float *f0) {
    static int firstCall = 1;
    static int nf0table = 771;
    static float f0_table[771];
    static int min_wl = 380;
    static int del_wl = 1;

    FILE *fp;
    char *filedir;
    char filename[FILENAME_MAX];
    char line [80];
    int32_t i, imin, imax;

    if (firstCall) {
        float tmp_wl;
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            return (0);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/f0_table_thuillier.dat");

        if ((fp = fopen(filename, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",
                    __FILE__, __LINE__, filename);
            return (0);
        }

        for (i = 0; i < nf0table; i++) {
            if (fgets(line, 80, fp) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: error reading %s at line %d\n",
                        __FILE__, __LINE__, filename, i);
                return (0);
            }
            sscanf(line, "%f %f", &tmp_wl, &f0_table[i]);
        }

        firstCall = 0;
    }

    i = (wl - min_wl) / del_wl;
    imin = MAX(i - width / 2 / del_wl, 0);
    imax = MIN(i + width / 2 / del_wl, nf0table);

    if (i < 0 || i > nf0table) {
        fprintf(stderr,
                "-E- %s line %d: wavelength of %d outside F0 table range.\n",
                __FILE__, __LINE__, wl);
        return (0);
    }

    *f0 = 0.0;
    for (i = imin; i <= imax; i++)
        *f0 += f0_table[i];

    *f0 /= (imax - imin + 1);

    return (1);
}

int get_f0_thuillier_(int32_t *wl, int32_t *width, float *f0) {
    return (get_f0_thuillier(*wl, *width, f0));
}

void get_f0_thuillier_ext(int32_t wl, int32_t width, float *f0) {
    static int firstCall = 1;
    static int nheader = 15;
    static int nf0table = 2198;
    static float f0_table[2198];
    static int min_wl = 200;
    static int del_wl = 1;

    FILE *fp;
    char *filedir;
    char filename[FILENAME_MAX];
    char line [80];
    int32_t i, imin, imax;

    if (firstCall) {
        float tmp_wl;
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(EXIT_FAILURE);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/Thuillier_F0.dat");

        if ((fp = fopen(filename, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",
                    __FILE__, __LINE__, filename);
            exit(EXIT_FAILURE);
        }

        if (want_verbose)
            printf("Reading Thuillier_F0.dat\n");
        for (i = 0; i < nheader; i++) {
            fgets(line, 80, fp);
            /*
            printf("%s",line);
             */
        }

        for (i = 0; i < nf0table; i++) {
            if (fgets(line, 80, fp) == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: error reading %s at line %d\n",
                        __FILE__, __LINE__, filename, i);
                exit(EXIT_FAILURE);
            }
            sscanf(line, "%f %f", &tmp_wl, &f0_table[i]);
        }

        firstCall = 0;
    }

    i = (wl - min_wl) / del_wl;
    imin = MAX(i - width / 2 / del_wl, 0);
    imax = MIN(i + width / 2 / del_wl, nf0table);

    if (i < 0 || i > nf0table) {
        fprintf(stderr,
                "-E- %s line %d: wavelength of %d outside F0 table range.\n",
                __FILE__, __LINE__, wl);
        exit(EXIT_FAILURE);
    }

    *f0 = 0.0;
    for (i = imin; i <= imax; i++)
        *f0 += f0_table[i];

    *f0 /= (imax - imin + 1);

}

int get_f0_thuillier_ext_(int32_t *wl, int32_t *width, float *f0) {
    return (get_f0_thuillier(*wl, *width, f0));
}

