
/*
  cc -64  -fullwarn -O2 hybrid.c -o hybrid -I$HDFINC -L$HDFLIB -lmfhdf -ldf -lz
  cc -64  -fullwarn -O2 hybrid.c -o /land2/jack/TEST/hybrid -I$HDFINC 
-L$HDFLIB -lmfhdf -ldf -lz
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mfhdf.h"

#define MAXNAMELENGTH 200
#define RADIUS 2
#define SDSNAME1 { "refl_443", "refl_555", "refl_670", "refl_865", "maximum_NDVI", "EVI", "input_file" }
#define SDSNAME2 { "minimum_refl_443", "refl_555", "refl_670", "refl_865", "NDVI", "EVI", "input_file" }
#define SDSNAME3 { "refl_443", "refl_555", "refl_670", "refl_865", "NDVI", "EVI", "input_file" }

enum {
    maxNDVI, minBlue, Hybrid, Nfiles
};

enum {
    BLUE, GREEN, RED, NIR, NDVI, EVI, INPUT, Nbands
};

typedef struct {
    char name[MAXNAMELENGTH];
    int32 file_id, id, index, num_type, rank, n_attr, Nl, Np, *plane, Nplanes;
    int32 start[H4_MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dim_sizes[MAX_VAR_DIMS];
    void *data, *fillvalue;
} SDS;

double R3(double x);

void main(int argc, char *argv[]) {
    char filename[Nfiles][MAXNAMELENGTH];
    SDS sds[Nfiles][Nbands];
    char *sdsname1[Nbands] = SDSNAME1;
    char *sdsname2[Nbands] = SDSNAME2;
    char *sdsname3[Nbands] = SDSNAME3;
    char openmode[Nfiles];
    int k, j, ib, irow, jcol, idx, source;
    float NIRreldiff, NIRness;
    int npix, nexam, nrepl;
    int Nrows, Ncols;
    double weigh, sum, val1, val2;

    if (argc < 4) {
        fprintf(stderr, "Usage: hybrid <max NDVI file> <min Blue file> <hybrid file>\n");
        exit(-1);
    }

    for (j = 0; j < Nfiles; j++)
        strcpy(filename[j], argv[j + 1]);

    for (ib = 0; ib < Nbands; ib++) {
        strcpy(sds[maxNDVI][ib].name, sdsname1[ib]);
        strcpy(sds[minBlue][ib].name, sdsname2[ib]);
        strcpy(sds[Hybrid][ib].name, sdsname3[ib]);
    }

    openmode[maxNDVI] = openmode[minBlue] = DFACC_RDONLY;
    openmode[Hybrid] = DFACC_CREATE;

    for (j = 0; j < Nfiles; j++) {
        if ((sds[j][0].file_id = SDstart(filename[j], openmode[j])) == -1) {
            fprintf(stderr, "Cannot open file %s\n", filename[j]);
            exit(-1);
        }
        for (ib = 0; ib < Nbands; ib++)
            sds[j][ib].file_id = sds[j][0].file_id;
        printf("File %s:\n", filename[j]);
        if (j == maxNDVI || j == minBlue)
            for (ib = 0; ib < Nbands; ib++) {
                if ((sds[j][ib].index = SDnametoindex(sds[j][ib].file_id,
                        sds[j][ib].name)) == -1) {
                    fprintf(stderr, "Cannot locate SDS %s in file %s\n", sds[j][ib].name,
                            filename[j]);
                    exit(-1);
                }
                if ((sds[j][ib].id = SDselect(sds[j][ib].file_id, sds[j][ib].index)) ==
                        -1) {
                    fprintf(stderr, "Cannot select SDS %s in file %s\n", sds[j][ib].name,
                            filename[j]);
                    exit(-1);
                }
                if (SDgetinfo(sds[j][ib].id, sds[j][ib].name, &sds[j][ib].rank,
                        sds[j][ib].dim_sizes, &sds[j][ib].num_type, &sds[j][ib].n_attr) == -1) {
                    fprintf(stderr, "Cannot get info from SDS %s in file %s\n",
                            sds[j][ib].name, filename[j]);
                    exit(-1);
                }
            } else
            for (ib = 0; ib < Nbands; ib++) {
                sds[j][ib].num_type = sds[maxNDVI][ib].num_type;
                sds[j][ib].rank = 2;
                sds[j][ib].dim_sizes[0] = sds[maxNDVI][ib].dim_sizes[0];
                sds[j][ib].dim_sizes[1] = sds[maxNDVI][ib].dim_sizes[1];
                if ((sds[j][ib].id = SDcreate(sds[j][ib].file_id, sds[j][ib].name,
                        sds[j][ib].num_type, sds[j][ib].rank, sds[j][ib].dim_sizes)) == -1) {
                    fprintf(stderr, "Cannot create SDS %s in file %s\n", sds[j][ib].name,
                            filename[j]);
                    exit(-1);
                }
            }
        for (ib = 0; ib < Nbands; ib++) {
            sds[j][ib].Nl = sds[j][ib].dim_sizes[0];
            sds[j][ib].Np = sds[j][ib].dim_sizes[1];
            printf("  SDS %d: \"%s\" %dx%d\n", ib + 1, sds[j][ib].name, sds[j][ib].Np,
                    sds[j][ib].Nl);
            sds[j][ib].data = malloc((2 * RADIUS + 1) * sds[j][ib].Np *
                    DFKNTsize(sds[j][ib].num_type));
        }
    }

    for (j = 0; j < Nfiles; j++)
        for (ib = 0; ib < Nbands; ib++) {
            sds[j][ib].start[1] = 0;
            sds[j][ib].edges[0] = 1;
            sds[j][ib].edges[1] = sds[j][ib].Np;
        }

    npix = nexam = nrepl = 0;
    Nrows = sds[maxNDVI][0].Nl;
    Ncols = sds[maxNDVI][0].Np;

    for (irow = 0; irow < Nrows; irow++) {

        if (irow % 100 == 0)
            printf("row %d\n", irow);

        for (j = 0; j < Nfiles; j++)
            for (ib = 0; ib < Nbands; ib++) {
                sds[j][ib].start[0] = MAX(0, irow - RADIUS);
                sds[j][ib].edges[0] = MIN(Nrows - sds[j][ib].start[0], irow + RADIUS -
                        sds[j][ib].start[0] + 1);
                if (SDreaddata(sds[j][ib].id, sds[j][ib].start, NULL, sds[j][ib].edges,
                        (char *) sds[j][ib].data + (2 * RADIUS + 1 - sds[j][ib].edges[0]) * Ncols *
                        DFKNTsize(sds[j][ib].num_type)) == -1) {
                    printf("  Can't read SDS \"%s\"\n", sds[j][ib].name);
                    exit(-1);
                }
            }

        for (jcol = 0; jcol < Ncols; jcol++) {

            idx = RADIUS * Ncols + jcol;

            if (((int16 *) sds[minBlue][NIR].data)[idx] != 0)
                NIRreldiff = (((int16 *) sds[maxNDVI][NIR].data)[idx] - ((int16
                    *) sds[minBlue][NIR].data)[idx]) / fabs(((int16 *) sds[minBlue][NIR].data)[idx])
                ;
            else
                NIRreldiff = 1;

            source = minBlue;
            if (((int16 *) sds[minBlue][BLUE].data)[idx] != 0 && /* Not water */
                    ((int16 *) sds[maxNDVI][BLUE].data)[idx] != 0) { /* Not fill value */
                npix++;
                if (NIRreldiff > 0.10) { /* Signed difference */
                    nexam++;
                    val1 = 0;
                    val2 = 0;
                    sum = 0;
                    for (k = -RADIUS; k <= RADIUS; k++)
                        if (irow + k >= 0 &&
                                irow + k <= Nrows - 1)
                            for (j = -RADIUS; j <= RADIUS; j++)
                                if (jcol + j >= 0 &&
                                        jcol + j <= Ncols - 1) {
                                    weigh = R3(fabs(0.5 * k)) * R3(fabs(0.5 * j));
                                    sum += weigh;
                                    val1 += weigh * ((int16 *) sds[maxNDVI][NIR].data)[idx + k *
                                            Ncols + j];
                                    val2 += weigh * ((int16 *) sds[minBlue][NIR].data)[idx + k *
                                            Ncols + j];
                                }
                    if (sum > 0) {
                        val1 /= sum;
                        val2 /= sum;
                    }
                    NIRness = (((int16 *) sds[maxNDVI][NIR].data)[idx] - ((int16
                            *) sds[maxNDVI][BLUE].data)[idx]) / fabs(((int16 *) sds[maxNDVI][BLUE].data)[idx
                            ]);
                    if (fabs(NIRness) > 0.2 && /* Too much spectral variation for a 
cloud */
                            ((int16 *) sds[maxNDVI][BLUE].data)[idx] < 2500 && /* Not bright 
enough for most clouds, but enough for most deserts (would be a better test at 
412nm) */
                            abs(val1 - ((int16 *) sds[maxNDVI][NIR].data)[idx]) <= abs(val2 -
                            ((int16 *) sds[minBlue][NIR].data)[idx]) &&
                            (((int16 *) sds[maxNDVI][NDVI].data)[idx] < 3000 || /* Either 
not vegetation */
                            ((int16 *) sds[maxNDVI][BLUE].data)[idx] < 500 || /* Or blue is 
dark enough (filters out smoke over dense vegetation) */
                            ((int16 *) sds[minBlue][NIR].data)[idx] == 0)) {
                        /*
                        printf("irow %d  jcol %d  NIRreldiff %g  NIRness %g\n", irow, jcol, 
                        NIRreldiff, NIRness);
                         */
                        source = maxNDVI;
                        nrepl++;
                    }
                }
            }

            for (ib = 0; ib < Nbands; ib++)
                switch (sds[Hybrid][ib].num_type) {
                case DFNT_UINT8: ((uint8 *) sds[Hybrid][ib].data)[idx] = ((uint8
                            *) sds[source][ib].data)[idx];
                    break;
                case DFNT_INT16: ((int16 *) sds[Hybrid][ib].data)[idx] = ((int16
                            *) sds[source][ib].data)[idx];
                    break;
                }

        }

        for (ib = 0; ib < Nbands; ib++)
            if (SDwritedata(sds[Hybrid][ib].id, sds[Hybrid][ib].start, NULL,
                    sds[Hybrid][ib].edges, sds[Hybrid][ib].data) == -1) {
                fprintf(stderr, "Cannot write row %d of SDS %s\n", irow,
                        sds[Hybrid][ib].name);
                exit(-1);
            }

    }

    for (j = 0; j < Nfiles; j++) {
        for (ib = 0; ib < Nbands; ib++)
            SDendaccess(sds[j][ib].id);
        SDend(sds[j][ib].file_id);
    }

    printf("%d land pixels\n", npix);
    printf("%d pixels challenged\n", nexam);
    printf("%d pixels replaced\n", nrepl);

}

double R3(double x) {
    if (x < 0 || x >= 2)
        return 0;
    else if (x <= 1) /* 0 <= x <= 1 */
        return 2 / 3. + x * x * x / 2. - x * x;
    else /* 1 <= x <= 2 */
        return (2 - x) * (2 - x) * (2 - x) / 6;
}


