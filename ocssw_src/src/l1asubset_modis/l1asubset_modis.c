/*
    Modification history:
    Programmer       Organization      Date      Description of change
    --------------   ------------    --------    ---------------------
    Joel Gales       Futuretech      03/28/03    Original Development
    Joel Gales       Futuretech      10/18/07    Rescaled pixels with
                                                 value less than 4095
                                                 are set to 4095.
    Joel Gales       Futuretech      10/26/07    Back out prev change
    Joel Gales       Futuretech      10/26/07    Put back prev change
                                                 put check if < 0 before
                                                 setting to 4095
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "hdf.h"
#include "mfhdf.h"

#include <timeutils.h>

#define basename(s) (strrchr((s), '/') == NULL ? (s) : strrchr((s), '/') + 1)


#define VERSION "1.45"

void usage(char *progname) {
    printf("This is version %s of %s (compiled on %s %s)\n",
            VERSION, progname, __DATE__, __TIME__);

    printf("\nUsage: %s ifile ofile\n", progname);
    printf("   ifile      = input filename\n");
    printf("   ofile      = output filename\n");
    exit(0);
}

int main(int argc, char *argv[]) {

    int32 i;
    int32 j;
    int32 k;
    int32 l;
    int32 n;
    int32 sds_id;
    int32 i32;
    int32 dims[8];
    int32 rank;
    int32 nelem;
    int32 dtype;
    int32 nattrs;
    int32 status;

    int32 vdid;
    int32 vdid_w;
    int32 nflds;
    int32 nrec;
    int32 nread;
    int32 nwrite;
    int32 rescale;
    int32 last_bad = 0;
    int32 daynotneeded[] = {11, 12, 13, -1}; /* 0-based */
    int32 ngtnotneeded[] = {1, 4, 10, 13, 14, 15, 16, -1}; /* 0-based */

    int32 nlines;
    int32 npix;

    int16 idata_HR[16];
    int16 sz;
    int16 *data_1km;

    float32 *median_250_500;
    float32 sumx, sumy, sumx2, sumxy, a, b, scale_parm[2];
    float32 f32;
    float32 rescale_parm[8];

    static char buffer[2 * 2048];
    static char buf2[2048];

    char sat[8];

    char *data;
    char *cptr;

    int32 HDFfid_r;
    int32 sd_id_r;
    int32 HDFfid_w;
    int32 sd_id_w;
    int32 sds_id_w;
    int32 ndatasets;
    int32 nglobal_attr;
    int32 count;
    int32 dim_id_r;
    int32 dim_id_w;
    int32 start[3] = {0, 0, 0};
    int32 start_w[3] = {0, 0, 0};
    int32 maxsize;
    int32 *lonebuf;
    int32 coremeta = 0;

    char *tmp_str;

    int compare_int16(const void *, const void *);
    int compare_float32(const void *, const void *);

    FILE *fp;

    struct tm *tmnow;
    time_t tnow;

    char prodtime[80];


    if (argc == 1) {
        usage(argv[0]);
        exit(0);
    }

    time(&tnow);
    tmnow = gmtime(&tnow);
    strftime(prodtime, 80, "%Y-%m-%dT%XZ", tmnow);

    HDFfid_r = Hopen(argv[1], DFACC_READ, 0);
    status = Vstart(HDFfid_r);
    sd_id_r = SDstart(argv[1], DFACC_RDONLY);


    /* Create output HDF file */
    HDFfid_w = Hopen(argv[2], DFACC_CREATE, 0);
    status = Vstart(HDFfid_w);
    sd_id_w = SDstart(argv[2], DFACC_RDWR);

    printf("argc: %d\n", argc);

    /*
      rescale =  0 (Perform rescaling)
      rescale =  1 (Perform rescaling with granule-based parameters)
      rescale = -1 (No rescaling)
     */
    rescale = 0;
    if (argc >= 4) {
        rescale = atoi(argv[3]);

    }

    /* Determine number of datasets in input files */
    SDfileinfo(sd_id_r, &ndatasets, &nglobal_attr);

    /* For each dataset ... */
    for (j = 0; j < ndatasets; j++) {
        sds_id = SDselect(sd_id_r, j);
        status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

        if (strcmp(buffer, "SD_250m") == 0) continue;
        if (strcmp(buffer, "SD_500m") == 0) continue;

        if (strcmp(buffer, "SRCA_250m") == 0) continue;
        if (strcmp(buffer, "SRCA_500m") == 0) continue;

        if (strcmp(buffer, "BB_250m") == 0) continue;
        if (strcmp(buffer, "BB_500m") == 0) continue;

        if (strcmp(buffer, "SV_250m") == 0) continue;
        if (strcmp(buffer, "SV_500m") == 0) continue;

        if (strcmp(buffer, "EV_250m") == 0) continue;
        if (strcmp(buffer, "EV_500m") == 0) continue;


        /*
        printf("Name :%s  rank: %d  type: %d\n", buffer, rank, dtype);
        printf("dims: %d %d %d\n", dims[0], dims[1], dims[2]);
         */

        sds_id_w = SDcreate(sd_id_w, buffer, dtype, rank, dims);
        if (sds_id_w == -1) {
            printf("Field: %s cannot be created\n", buffer);
            exit(-1);
        }

        for (i = 0; i < nattrs; i++) {
            status = SDreadattr(sds_id, i, (VOIDP) buf2);
            status = SDattrinfo(sds_id, i, buffer, &dtype, &count);
            status = SDsetattr(sds_id_w, buffer, dtype, count, (VOIDP) buf2);
        }

        for (i = 0; i < rank; i++) {
            dim_id_r = SDgetdimid(sds_id, i);
            dim_id_w = SDgetdimid(sds_id_w, i);
            SDdiminfo(dim_id_r, buffer, &count, &dtype, &nattrs);
            SDsetdimname(dim_id_w, buffer);
        }

        SDendaccess(sds_id_w);
        SDendaccess(sds_id);
    } /* End Create SD */


    /* Copy Datasets */
    /* ------------- */
    n = 0;
    for (j = 0; j < ndatasets; j++) {
        sds_id = SDselect(sd_id_r, j);
        status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

        if (strcmp(buffer, "SD_250m") == 0) continue;
        if (strcmp(buffer, "SD_500m") == 0) continue;

        if (strcmp(buffer, "SRCA_250m") == 0) continue;
        if (strcmp(buffer, "SRCA_500m") == 0) continue;

        if (strcmp(buffer, "BB_250m") == 0) continue;
        if (strcmp(buffer, "BB_500m") == 0) continue;

        if (strcmp(buffer, "SV_250m") == 0) continue;
        if (strcmp(buffer, "SV_500m") == 0) continue;

        if (strcmp(buffer, "EV_250m") == 0) continue;
        if (strcmp(buffer, "EV_500m") == 0) continue;

        sds_id_w = SDselect(sd_id_w, n++);

        /*
        printf("% 3d Name :%s  rank: %d  type: %d\n", j, buffer, rank, dtype);
        for (k=0; k<rank; k++) printf("dim %d: %d\n", k, dims[k]);
         */

        nelem = dims[0];
        for (k = 1; k < rank; k++) nelem *= dims[k];

        data = (char *) calloc(nelem, DFKNTsize(dtype));

        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) data);


#if 0
        /* Check whether last scan is bad */
        /* ------------------------------ */
        if (strcmp(buffer, "Frame count array") == 0 && last_bad == 0) {
            memcpy(&i16, &data[(6 * 202 + 1)*2], 2);
            if (i16 == 0) last_bad = 1;
            printf("last_bad: %d\n", last_bad);
        }
#endif


        /* Flag unneeded day bands */
        /* ----------------------- */
        if (strcmp(buffer, "EV_1km_day") == 0) {
            printf("Masking Unneeded Day Bands\n");
            for (k = 0; k < dims[0]; k++) {
                l = 0;
                while (daynotneeded[l] != -1) {
                    if (k == 0) printf("%d\n", daynotneeded[l]);
                    memset(&data[k * 2 * dims[1] * dims[2] + daynotneeded[l]*2 * dims[2]],
                            255, 2 * dims[2]);
                    l++;
                }
            }
        }

        /* Flag unneeded night bands */
        /* ------------------------- */
        if (strcmp(buffer, "EV_1km_night") == 0) {
            printf("Masking Unneeded Night Bands\n");
            for (k = 0; k < dims[0]; k++) {
                l = 0;
                while (ngtnotneeded[l] != -1) {
                    if (k == 0) printf("%d\n", ngtnotneeded[l]);
                    memset(&data[k * 2 * dims[1] * dims[2] + ngtnotneeded[l]*2 * dims[2]],
                            255, 2 * dims[2]);
                    l++;
                }
            }
        }


        status = SDwritedata(sds_id_w, start_w, NULL, dims, (VOIDP) data);
        if (status == -1) {
            printf("write status: %d\n\n", status);
            exit(-1);
        }

        free(data);

        SDendaccess(sds_id);
        SDendaccess(sds_id_w);
    }


    if (SDselect(sd_id_r, SDnametoindex(sd_id_r, "EV_250m")) == -1 ||
            SDselect(sd_id_r, SDnametoindex(sd_id_r, "EV_500m")) == -1) {
        printf("No 250/500 EV field found.\n");
        rescale = -1;
    }

    if (rescale == -1) {
        printf("No rescaling performed\n");
        goto skip_rescale;
    }


    /* Read rescaling parameters */
    /* ------------------------- */
    if (rescale == 0) {
        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            return (1);
        }

        strcpy(buf2, tmp_str);
        strcat(buf2, "/modisa/l1a_scaling.dat");

        if ((fp = fopen(buf2, "r")) == NULL) {
            printf("Error: Unable to open rescaling parameter file: %s.\n", buf2);
            return (1);
        }

        for (i = 0; i < 4; i++) {
            if (fscanf(fp, "%f %f\n",
                    &rescale_parm[2 * i], &rescale_parm[2 * i + 1]) == EOF) {
                printf("Error: Insufficent number of scaling parameters: %s.\n",
                        buf2);
                return (1);
            }
        }
    }


    /* Rescale "4095" 1km B pixels (B10) using 500m B pixels (B3) */
    /* ---------------------------------------------------------- */
    printf("Rescale saturated 1km B pixels (B10) using 500m B pixels (B3)\n");
    sz = 2;
    sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "EV_1km_day"));
    status = SDgetinfo(sds_id_w, buffer, &rank, dims, &dtype, &nattrs);
    start[1] = 2;
    dims[1] = 1;

    nelem = dims[0];
    for (k = 1; k < rank; k++) nelem *= dims[k];

    nlines = dims[0];
    npix = dims[2];

    data_1km = (int16 *) calloc(nelem, sz);
    median_250_500 = (float32 *) calloc(nelem, sizeof (float32));

    status = SDreaddata(sds_id_w, start, NULL, dims, (VOIDP) data_1km);

    sds_id = SDselect(sd_id_r, SDnametoindex(sd_id_r, "EV_500m"));
    status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

    dims[0] = 2;
    start[1] = 0;
    dims[1] = 1;

    nelem = dims[0];
    for (k = 1; k < rank; k++) nelem *= dims[k];

    data = (char *) calloc(nelem, sz);

    sumx = 0;
    sumy = 0;
    sumx2 = 0;
    sumxy = 0;
    n = 0;
    for (i = 0; i < nlines; i++) {
        start[0] = dims[0] * i;
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) data);

        if (i % 400 == 0) printf("i: %d\n", i);

        for (j = 0; j < npix; j++) {

            memcpy(&idata_HR[0], &data[sz * 2 * (j + 0 * npix)], 2 * sz);
            memcpy(&idata_HR[2], &data[sz * 2 * (j + 1 * npix)], 2 * sz);

            qsort((void *) idata_HR, dims[0] * dims[0], sizeof (int16),
                    compare_int16);

            for (k = 0; k < 4; k++) {
                if (idata_HR[k] > 0) {

                    if ((k % 2) == 0)
                        median_250_500[i * npix + j] =
                            0.5 * (idata_HR[(2 + k) / 2] + idata_HR[(4 + k) / 2]);
                    else
                        median_250_500[i * npix + j] = (float32) idata_HR[(3 + k) / 2];

                    if (data_1km[i * npix + j] > 0 && data_1km[i * npix + j] < 4095) {
                        n++;
                        sumx += median_250_500[i * npix + j];
                        sumy += data_1km[i * npix + j];
                        sumx2 += median_250_500[i * npix + j] * median_250_500[i * npix + j];
                        sumxy += median_250_500[i * npix + j] * data_1km[i * npix + j];
                    }
                    break;
                }
            }
        }
    }
    a = (sumxy - sumx * sumy / n) / (sumx2 - sumx * sumx / n);
    b = (sumy / n) - a * (sumx / n);
    scale_parm[0] = b;
    scale_parm[1] = a;
    status = SDsetattr(sd_id_w, "B_10_scale_parm_computed", DFNT_FLOAT32, 2,
            (VOIDP) scale_parm);

    if (rescale == 0) {
        b = rescale_parm[0];
        a = rescale_parm[1];
        status = SDsetattr(sd_id_w, "B_10_scale_parm", DFNT_FLOAT32, 2,
                (VOIDP) & rescale_parm[0]);
    } else {
        status = SDsetattr(sd_id_w, "B_10_scale_parm", DFNT_FLOAT32, 2,
                (VOIDP) scale_parm);
    }
    printf("a: %f  b: %f\n", a, b);

    for (i = 0; i < nlines; i++)
        for (j = 0; j < npix; j++) {
            if (data_1km[i * npix + j] == 4095) {
                data_1km[i * npix + j] = a * median_250_500[i * npix + j] + b;
                if (data_1km[i * npix + j] >= 0 && data_1km[i * npix + j] < 4095)
                    data_1km[i * npix + j] = 4095;
            }
            if (data_1km[i * npix + j] < 0) data_1km[i * npix + j] = 32767;
        }

    status = SDgetinfo(sds_id_w, buffer, &rank, dims, &dtype, &nattrs);
    start[0] = 0;
    start[1] = 2;
    dims[1] = 1;
    status = SDwritedata(sds_id_w, start, NULL, dims, (VOIDP) data_1km);

    free(data);
    free(data_1km);
    free(median_250_500);

    SDendaccess(sds_id);
    SDendaccess(sds_id_w);



    /* Rescale "4095" 1km G pixels (B12) using 500m G pixels (B4) */
    /* ---------------------------------------------------------- */
    printf("Rescale saturated 1km G pixels (B12) using 500m G pixels (B4)\n");
    sz = 2;
    sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "EV_1km_day"));
    status = SDgetinfo(sds_id_w, buffer, &rank, dims, &dtype, &nattrs);
    start[1] = 4;
    dims[1] = 1;

    nelem = dims[0];
    for (k = 1; k < rank; k++) nelem *= dims[k];

    nlines = dims[0];
    npix = dims[2];

    data_1km = (int16 *) calloc(nelem, sz);
    median_250_500 = (float32 *) calloc(nelem, sizeof (float32));

    status = SDreaddata(sds_id_w, start, NULL, dims, (VOIDP) data_1km);

    sds_id = SDselect(sd_id_r, SDnametoindex(sd_id_r, "EV_500m"));
    status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

    dims[0] = 2;
    start[1] = 1;
    dims[1] = 1;

    nelem = dims[0];
    for (k = 1; k < rank; k++) nelem *= dims[k];

    data = (char *) calloc(nelem, sz);

    sumx = 0;
    sumy = 0;
    sumx2 = 0;
    sumxy = 0;
    n = 0;
    for (i = 0; i < nlines; i++) {
        start[0] = dims[0] * i;
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) data);

        if (i % 400 == 0) printf("i: %d\n", i);

        for (j = 0; j < npix; j++) {

            memcpy(&idata_HR[0], &data[sz * 2 * (j + 0 * npix)], 2 * sz);
            memcpy(&idata_HR[2], &data[sz * 2 * (j + 1 * npix)], 2 * sz);

            qsort((void *) idata_HR, dims[0] * dims[0], sizeof (int16),
                    compare_int16);

            for (k = 0; k < 4; k++) {
                if (idata_HR[k] > 0) {

                    if ((k % 2) == 0)
                        median_250_500[i * npix + j] =
                            0.5 * (idata_HR[(2 + k) / 2] + idata_HR[(4 + k) / 2]);
                    else
                        median_250_500[i * npix + j] = (float32) idata_HR[(3 + k) / 2];

                    if (data_1km[i * npix + j] > 0 && data_1km[i * npix + j] < 4095) {
                        n++;
                        sumx += median_250_500[i * npix + j];
                        sumy += data_1km[i * npix + j];
                        sumx2 += median_250_500[i * npix + j] * median_250_500[i * npix + j];
                        sumxy += median_250_500[i * npix + j] * data_1km[i * npix + j];
                    }
                    break;
                }
            }
        }
    }
    a = (sumxy - sumx * sumy / n) / (sumx2 - sumx * sumx / n);
    b = (sumy / n) - a * (sumx / n);
    scale_parm[0] = b;
    scale_parm[1] = a;
    status = SDsetattr(sd_id_w, "B_12_scale_parm_computed", DFNT_FLOAT32, 2,
            (VOIDP) scale_parm);

    if (rescale == 0) {
        b = rescale_parm[2];
        a = rescale_parm[3];
        status = SDsetattr(sd_id_w, "B_12_scale_parm", DFNT_FLOAT32, 2,
                (VOIDP) & rescale_parm[2]);
    } else {
        status = SDsetattr(sd_id_w, "B_12_scale_parm", DFNT_FLOAT32, 2,
                (VOIDP) scale_parm);
    }
    printf("a: %f  b: %f\n", a, b);

    for (i = 0; i < nlines; i++)
        for (j = 0; j < npix; j++) {
            if (data_1km[i * npix + j] == 4095) {
                data_1km[i * npix + j] = a * median_250_500[i * npix + j] + b;
                if (data_1km[i * npix + j] >= 0 && data_1km[i * npix + j] < 4095)
                    data_1km[i * npix + j] = 4095;
            }
            if (data_1km[i * npix + j] < 0) data_1km[i * npix + j] = 32767;
        }

    status = SDgetinfo(sds_id_w, buffer, &rank, dims, &dtype, &nattrs);
    start[0] = 0;
    start[1] = 4;
    dims[1] = 1;
    status = SDwritedata(sds_id_w, start, NULL, dims, (VOIDP) data_1km);

    free(data);
    free(data_1km);
    free(median_250_500);

    SDendaccess(sds_id);
    SDendaccess(sds_id_w);


    /* Rescale "4095" 1km R pixels (B13) using 250m R pixels (B1) */
    /* ---------------------------------------------------------- */
    printf("Rescale saturated 1km R pixels (B13) using 250m R pixels (B1)\n");
    sz = 2;
    sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "EV_1km_day"));
    status = SDgetinfo(sds_id_w, buffer, &rank, dims, &dtype, &nattrs);
    start[1] = 5;
    dims[1] = 1;

    nelem = dims[0];
    for (k = 1; k < rank; k++) nelem *= dims[k];

    nlines = dims[0];
    npix = dims[2];

    data_1km = (int16 *) calloc(nelem, sz);
    median_250_500 = (float32 *) calloc(nelem, sizeof (float32));

    status = SDreaddata(sds_id_w, start, NULL, dims, (VOIDP) data_1km);

    sds_id = SDselect(sd_id_r, SDnametoindex(sd_id_r, "EV_250m"));
    status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

    dims[0] = 4;
    start[1] = 0;
    dims[1] = 1;

    nelem = dims[0];
    for (k = 1; k < rank; k++) nelem *= dims[k];

    data = (char *) calloc(nelem, sz);

    sumx = 0;
    sumy = 0;
    sumx2 = 0;
    sumxy = 0;
    n = 0;
    for (i = 0; i < nlines; i++) {
        start[0] = dims[0] * i;
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) data);

        if (i % 400 == 0) printf("i: %d\n", i);

        for (j = 0; j < npix; j++) {

            memcpy(&idata_HR[0], &data[sz * 4 * (j + 0 * npix)], 4 * sz);
            memcpy(&idata_HR[4], &data[sz * 4 * (j + 1 * npix)], 4 * sz);
            memcpy(&idata_HR[8], &data[sz * 4 * (j + 2 * npix)], 4 * sz);
            memcpy(&idata_HR[12], &data[sz * 4 * (j + 3 * npix)], 4 * sz);

            qsort((void *) idata_HR, dims[0] * dims[0], sizeof (int16),
                    compare_int16);

            for (k = 0; k < dims[0] * dims[0]; k++) {
                if (idata_HR[k] > 0) {

                    if ((k % 2) == 0)
                        median_250_500[i * npix + j] =
                            0.5 * (idata_HR[(14 + k) / 2] + idata_HR[(16 + k) / 2]);
                    else
                        median_250_500[i * npix + j] = (float32) idata_HR[(15 + k) / 2];

                    if (data_1km[i * npix + j] > 0 && data_1km[i * npix + j] < 4095) {
                        n++;
                        sumx += median_250_500[i * npix + j];
                        sumy += data_1km[i * npix + j];
                        sumx2 += median_250_500[i * npix + j] * median_250_500[i * npix + j];
                        sumxy += median_250_500[i * npix + j] * data_1km[i * npix + j];
                    }
                    break;
                }
            }
        }
    }
    a = (sumxy - sumx * sumy / n) / (sumx2 - sumx * sumx / n);
    b = (sumy / n) - a * (sumx / n);
    scale_parm[0] = b;
    scale_parm[1] = a;
    status = SDsetattr(sd_id_w, "B_13_scale_parm_computed", DFNT_FLOAT32, 2,
            (VOIDP) scale_parm);

    if (rescale == 0) {
        b = rescale_parm[4];
        a = rescale_parm[5];
        status = SDsetattr(sd_id_w, "B_13_scale_parm", DFNT_FLOAT32, 2,
                (VOIDP) & rescale_parm[4]);
    } else {
        status = SDsetattr(sd_id_w, "B_13_scale_parm", DFNT_FLOAT32, 2,
                (VOIDP) scale_parm);
    }
    printf("a: %f  b: %f\n", a, b);

    for (i = 0; i < nlines; i++)
        for (j = 0; j < npix; j++) {
            if (data_1km[i * npix + j] == 4095) {
                data_1km[i * npix + j] = a * median_250_500[i * npix + j] + b;
                if (data_1km[i * npix + j] >= 0 && data_1km[i * npix + j] < 4095)
                    data_1km[i * npix + j] = 4095;
            }
            /*      if (data_1km[i*npix+j] == -32768) data_1km[i*npix+j] = 32767;*/
            if (data_1km[i * npix + j] < 0) data_1km[i * npix + j] = 32767;
        }

    status = SDgetinfo(sds_id_w, buffer, &rank, dims, &dtype, &nattrs);
    start[0] = 0;
    start[1] = 5;
    dims[1] = 1;
    status = SDwritedata(sds_id_w, start, NULL, dims, (VOIDP) data_1km);

    free(data);
    free(data_1km);
    free(median_250_500);

    SDendaccess(sds_id);
    SDendaccess(sds_id_w);



    /* Rescale "4095" 1km IR pixels (B16) using 250m IR pixels (B2) */
    /* ------------------------------------------------------------ */
    printf("Rescale saturated 1km IR pixels (B16) using 250m IR pixels (B2)\n");
    sz = 2;
    sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "EV_1km_day"));
    status = SDgetinfo(sds_id_w, buffer, &rank, dims, &dtype, &nattrs);
    start[1] = 10;
    dims[1] = 1;

    nelem = dims[0];
    for (k = 1; k < rank; k++) nelem *= dims[k];

    nlines = dims[0];
    npix = dims[2];

    data_1km = (int16 *) calloc(nelem, sz);
    median_250_500 = (float32 *) calloc(nelem, sizeof (float32));

    status = SDreaddata(sds_id_w, start, NULL, dims, (VOIDP) data_1km);

    sds_id = SDselect(sd_id_r, SDnametoindex(sd_id_r, "EV_250m"));
    status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

    dims[0] = 4;
    start[1] = 1;
    dims[1] = 1;

    nelem = dims[0];
    for (k = 1; k < rank; k++) nelem *= dims[k];

    data = (char *) calloc(nelem, sz);

    sumx = 0;
    sumy = 0;
    sumx2 = 0;
    sumxy = 0;
    n = 0;
    for (i = 0; i < nlines; i++) {
        start[0] = dims[0] * i;
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) data);

        if (i % 400 == 0) printf("i: %d\n", i);

        for (j = 0; j < npix; j++) {

            memcpy(&idata_HR[0], &data[sz * 4 * (j + 0 * npix)], 4 * sz);
            memcpy(&idata_HR[4], &data[sz * 4 * (j + 1 * npix)], 4 * sz);
            memcpy(&idata_HR[8], &data[sz * 4 * (j + 2 * npix)], 4 * sz);
            memcpy(&idata_HR[12], &data[sz * 4 * (j + 3 * npix)], 4 * sz);

            qsort((void *) idata_HR, dims[0] * dims[0], sizeof (int16),
                    compare_int16);

            for (k = 0; k < dims[0] * dims[0]; k++) {
                if (idata_HR[k] > 0) {

                    if ((k % 2) == 0)
                        median_250_500[i * npix + j] =
                            0.5 * (idata_HR[(14 + k) / 2] + idata_HR[(16 + k) / 2]);
                    else
                        median_250_500[i * npix + j] = (float32) idata_HR[(15 + k) / 2];

                    if (data_1km[i * npix + j] > 0 && data_1km[i * npix + j] < 4095) {
                        n++;
                        sumx += median_250_500[i * npix + j];
                        sumy += data_1km[i * npix + j];
                        sumx2 += median_250_500[i * npix + j] * median_250_500[i * npix + j];
                        sumxy += median_250_500[i * npix + j] * data_1km[i * npix + j];
                    }
                    break;
                }
            }
        }
    }
    a = (sumxy - sumx * sumy / n) / (sumx2 - sumx * sumx / n);
    b = (sumy / n) - a * (sumx / n);
    scale_parm[0] = b;
    scale_parm[1] = a;
    status = SDsetattr(sd_id_w, "B_16_scale_parm_computed", DFNT_FLOAT32, 2,
            (VOIDP) scale_parm);

    if (rescale == 0) {
        b = rescale_parm[6];
        a = rescale_parm[7];
        status = SDsetattr(sd_id_w, "B_16_scale_parm", DFNT_FLOAT32, 2,
                (VOIDP) & rescale_parm[6]);
    } else {
        status = SDsetattr(sd_id_w, "B_16_scale_parm", DFNT_FLOAT32, 2,
                (VOIDP) scale_parm);
    }
    printf("a: %f  b: %f\n", a, b);

    for (i = 0; i < nlines; i++)
        for (j = 0; j < npix; j++) {
            if (data_1km[i * npix + j] == 4095) {
                data_1km[i * npix + j] = a * median_250_500[i * npix + j] + b;
                if (data_1km[i * npix + j] >= 0 && data_1km[i * npix + j] < 4095)
                    data_1km[i * npix + j] = 4095;
            }
            if (data_1km[i * npix + j] < 0) data_1km[i * npix + j] = 32767;
        }

    status = SDgetinfo(sds_id_w, buffer, &rank, dims, &dtype, &nattrs);
    start[0] = 0;
    start[1] = 10;
    dims[1] = 1;
    status = SDwritedata(sds_id_w, start, NULL, dims, (VOIDP) data_1km);

    free(data);
    free(data_1km);
    free(median_250_500);

    SDendaccess(sds_id);
    SDendaccess(sds_id_w);


skip_rescale:

    /* Write lone Vdatas */
    /* ----------------- */
    maxsize = VSlone(HDFfid_r, NULL, 0);
    lonebuf = (int32 *) calloc(maxsize, sizeof (int32));
    maxsize = VSlone(HDFfid_r, lonebuf, maxsize);


    /* For each lone vdata ... */
    /* ----------------------- */
    for (i = 0; i < maxsize; i++) {
        vdid = VSattach(HDFfid_r, lonebuf[i], "r");
        vdid_w = VSattach(HDFfid_w, -1, "w");
        VSgetname(vdid, buffer);
        VSsetname(vdid_w, buffer);
        /*    printf("Vlone name: %s\n", buffer);*/

        VSsetinterlace(vdid_w, VSgetinterlace(vdid));


        /* Replace "," in fieldname string by 0 */
        /* ------------------------------------ */
        nflds = VSgetfields(vdid, buf2);
        k = strlen(buf2);
        for (j = 0; j < k; j++) if (buf2[j] == ',') buf2[j] = 0;


        /* Parse fieldname string and define fields */
        /* ---------------------------------------- */
        cptr = &buf2[0];
        for (j = 0; j < nflds; j++) {
            VSfdefine(vdid_w, cptr, VFfieldtype(vdid, j), VFfieldorder(vdid, j));

            while (1) {
                cptr++;
                if (*cptr == 0) {
                    cptr++;
                    break;
                }
            }
        }

        for (j = 0; j < k - 1; j++) if (buf2[j] == 0) buf2[j] = ',';
        status = VSsetfields(vdid, buf2);
        status = VSsetfields(vdid_w, buf2);

        VSQuerycount(vdid, &nrec);
        /*    printf("nrec: %d  vsize: %d\n", nrec, VSsizeof(vdid, buf2));*/
        data = (char *) calloc(nrec * VSsizeof(vdid, buf2), 1);

        nread = VSread(vdid, (uint8*) data, nrec, VSgetinterlace(vdid));
        nwrite = VSwrite(vdid_w, (uint8*) data, nrec, VSgetinterlace(vdid_w));
        /*    printf("nread: %d\n", nread);*/

        if (nread != nwrite) {
            printf("%d %d %d\n", i, nread, nwrite);
        }

        free(data);
        VSdetach(vdid_w);
        VSdetach(vdid);
    }
    free(lonebuf);


    /* Write Global Attributes */
    /* ----------------------- */
    for (i = 0; i < nglobal_attr; i++) {
        status = SDattrinfo(sd_id_r, i, buffer, &dtype, &count);

        if (strcmp(buffer, "CoreMetadata.0") == 0) {
            coremeta = 1;
        }

        if (strcmp(buffer, "Number of scans") == 0) {
            strcpy(buffer, "Number of Scans");
        }

        /* printf("%d %s\n", count, buffer);*/
        data = (char *) calloc(count, DFKNTsize(dtype));

        status = SDreadattr(sd_id_r, i, (VOIDP) data);

        if (strcmp(buffer, "Number of Scans") == 0 && last_bad == 1) {
            i32 = 202;
            memcpy(data, &i32, sizeof (i32));
        }


        if (strcmp(buffer, "Processing Time") == 0) strcpy(data, ydhmsf(now(), 'L'));
        status = SDsetattr(sd_id_w, buffer, dtype, count, (VOIDP) data);
        /*    printf("%d %d\n", i, status);*/

        if (strcmp(buffer, "Satellite") == 0) {
            strcpy(sat, data);
        }

        free(data);
    }


    if (coremeta == 0) {

#if 0
        strcpy(buf2, getenv("SIMBIOS_ROOT"));
        strcat(buf2, "/data/modis/static/");
        strcat(buf2, sat);
        strcat(buf2, "/coremeta.txt");
        printf("%s\n", buf2);
        fp = fopen(buf2, "r");

        if (fp == 0x0) {
            printf("\"coremeta.txt\" cannot be found.\n");
        } else {
            buffer[0] = 10;
            buffer[1] = 0;

            while (fgets(buf2, 80, fp) != NULL) {
                strcat(buffer, buf2);
            }
            if (feof(fp) != 0)
                status = SDsetattr(sd_id_w, "CoreMetadata.0", DFNT_CHAR,
                    strlen(buffer) + 1, (VOIDP) buffer);
        }
#endif

        /* Write CoreMetadata (jmg) */
        strcpy(buffer, "\n");
        strcat(buffer, "GROUP                  = INVENTORYMETADATA\n");
        strcat(buffer, "  GROUPTYPE            = MASTERGROUP\n");

        strcat(buffer, "    GROUP              = ECSDATAGRANULE\n");

        strcat(buffer, "      OBJECT           = LOCALGRANULEID\n");
        strcat(buffer, "        NUM_VAL        = 1\n");


        memset(buf2, 0, sizeof (buf2));
        strcpy(buf2, getenv("L0File"));
        strcpy(buf2, basename(buf2));

        strcat(buffer, "        VALUE          = ");
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "      END_OBJECT       = LOCALGRANULEID\n");

        strcat(buffer, "      OBJECT           = PRODUCTIONDATETIME\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        strcat(buffer, "\"");
        strcat(buffer, prodtime);
        strcat(buffer, "\"\n");
        strcat(buffer, "      END_OBJECT        = PRODUCTIONDATETIME\n");

        strcat(buffer, "      OBJECT           = DAYNIGHTFLAG\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "DAYNIGHTFLAG")) != -1)
            status = SDreadattr(sd_id_r, i, buf2);
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "      END_OBJECT       = DAYNIGHTFLAG\n");

        strcat(buffer, "      OBJECT           = REPROCESSINGACTUAL\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        strcat(buffer, "\"processed once\"");
        strcat(buffer, "\n");
        strcat(buffer, "      END_OBJECT       = REPROCESSINGACTUAL\n");

        strcat(buffer, "    END_GROUP          = ECSDATAGRANULE\n");

        strcat(buffer, "    GROUP              = ORBITCALCULATEDSPATIALDOMAIN\n");

        strcat(buffer, "      OBJECT           = ORBITCALCULATEDSPATIALDOMAINCONTAINER\n");
        strcat(buffer, "        CLASS          = \"1\"\n");
        strcat(buffer, "        OBJECT         = EQUATORCROSSINGDATE\n");
        strcat(buffer, "          CLASS        = \"1\"\n");
        strcat(buffer, "            NUM_VAL    = 1\n");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "EQUATORCROSSINGDATE")) != -1)
            status = SDreadattr(sd_id_r, i, buf2);
        strcat(buffer, "            VALUE      = ");
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "        END_OBJECT      = EQUATORCROSSINGDATE\n");

        strcat(buffer, "        OBJECT         = EQUATORCROSSINGTIME\n");
        strcat(buffer, "          CLASS        = \"1\"\n");
        strcat(buffer, "            NUM_VAL    = 1\n");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "EQUATORCROSSINGTIME")) != -1)
            status = SDreadattr(sd_id_r, i, buf2);
        strcat(buffer, "            VALUE      = ");
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "        END_OBJECT      = EQUATORCROSSINGTIME\n");

        strcat(buffer, "        OBJECT         = ORBITNUMBER\n");
        strcat(buffer, "          CLASS        = \"1\"\n");
        strcat(buffer, "            NUM_VAL    = 1\n");
        strcat(buffer, "            VALUE      = ");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "ORBITNUMBER")) != -1) {
            status = SDreadattr(sd_id_r, i, &i32);
            sprintf(buf2, "%6d\n", i32);
        } else sprintf(buf2, "%6d\n", -1);
        strcat(buffer, buf2);
        strcat(buffer, "        END_OBJECT     = ORBITNUMBER\n");

        strcat(buffer, "        OBJECT         = EQUATORCROSSINGLONGITUDE\n");
        strcat(buffer, "          CLASS        = \"1\"\n");
        strcat(buffer, "            NUM_VAL    = 1\n");
        strcat(buffer, "            VALUE      = ");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "EQUATORCROSSINGLONGITUDE")) != -1) {
            status = SDreadattr(sd_id_r, i, &f32);
            sprintf(buf2, "%8.5f\n", f32);
        } else sprintf(buf2, "%8.5f\n", 999.0);
        strcat(buffer, buf2);
        strcat(buffer, "        END_OBJECT     = EQUATORCROSSINGLONGITUDE\n");

        strcat(buffer, "      END_OBJECT       = ORBITCALCULATEDSPATIALDOMAINCONTAINER\n");
        strcat(buffer, "    END_GROUP          = ORBITCALCULATEDSPATIALDOMAIN\n");

        strcat(buffer, "    GROUP              = PGEVERSIONCLASS\n");
        strcat(buffer, "      OBJECT           = PGEVERSION\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        strcat(buffer, "\"4.0.3\"\n");
        strcat(buffer, "      END_OBJECT       = PGEVERSION\n");
        strcat(buffer, "    END_GROUP          = PGEVERSIONCLASS\n");

        strcat(buffer, "    GROUP              = SPATIALDOMAINCONTAINER\n");
        strcat(buffer, "      GROUP            = HORIZONTALSPATIALDOMAINCONTAINER\n");
        strcat(buffer, "        GROUP          = GPOLYGON\n");
        strcat(buffer, "          OBJECT       = GPOLYGONCONTAINER\n");
        strcat(buffer, "            CLASS      = \"1\"\n");

        strcat(buffer, "            GROUP      = GRINGPOINT\n");
        strcat(buffer, "              CLASS    = \"1\"\n");
        strcat(buffer, "              OBJECT   = GRINGPOINTLONGITUDE\n");
        strcat(buffer, "                NUM_VAL= 4\n");
        strcat(buffer, "                CLASS  = \"1\"\n");
        strcat(buffer, "                VALUE  = ");
        memset(buf2, 0, sizeof (buf2));
        sprintf(buf2, "(%8.5f, %8.5f, %8.5f, %8.5f)\n", -1., -1., -1., -1.);
        strcat(buffer, buf2);
        strcat(buffer, "              END_OBJECT= GRINGPOINTLONGITUDE\n");

        strcat(buffer, "              OBJECT   = GRINGPOINTLATITUDE\n");
        strcat(buffer, "                NUM_VAL= 4\n");
        strcat(buffer, "                CLASS  = \"1\"\n");
        strcat(buffer, "                VALUE  = ");
        memset(buf2, 0, sizeof (buf2));
        sprintf(buf2, "(%8.5f, %8.5f, %8.5f, %8.5f)\n", -1., -1., -1., -1.);
        strcat(buffer, buf2);
        strcat(buffer, "              END_OBJECT= GRINGPOINTLATITUDE\n");

        strcat(buffer, "              OBJECT   = GRINGPOINTSEQUENCENO\n");
        strcat(buffer, "                NUM_VAL= 4\n");
        strcat(buffer, "                CLASS  = \"1\"\n");
        strcat(buffer, "                VALUE  = ");
        memset(buf2, 0, sizeof (buf2));
        sprintf(buf2, "(%4d, %4d, %4d, %4d)\n", 1, 2, 3, 4);
        strcat(buffer, buf2);
        strcat(buffer, "              END_OBJECT= GRINGPOINTSEQUENCENO\n");
        strcat(buffer, "            END_GROUP  = GRINGPOINT\n");

        strcat(buffer, "            GROUP      = GRING\n");
        strcat(buffer, "              CLASS    = \"1\"\n");
        strcat(buffer, "              OBJECT   = EXCLUSIONGRINGFLAG\n");
        strcat(buffer, "                NUM_VAL= 1\n");
        strcat(buffer, "                CLASS  = \"1\"\n");
        strcat(buffer, "                VALUE  = ");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "EXCLUSIONGRINGFLAG")) != -1)
            status = SDreadattr(sd_id_r, i, buf2);
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "              END_OBJECT= EXCLUSIONGRINGFLAG\n");
        strcat(buffer, "            END_GROUP   = GRING\n");
        strcat(buffer, "          END_OBJECT    = GPOLYGONCONTAINER\n");
        strcat(buffer, "        END_GROUP       = GPOLYGON\n");
        strcat(buffer, "      END_GROUP         = HORIZONTALSPATIALDOMAINCONTAINER\n");
        strcat(buffer, "    END_GROUP           = SPATIALDOMAINCONTAINER\n");

        strcat(buffer, "    GROUP              = RANGEDATETIME\n");

        strcat(buffer, "      OBJECT           = RANGEENDINGDATE\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "RANGEENDINGDATE")) != -1)
            status = SDreadattr(sd_id_r, i, buf2);
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "      END_OBJECT        = RANGEENDINGDATE\n");

        strcat(buffer, "      OBJECT           = RANGEENDINGTIME\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "RANGEENDINGTIME")) != -1)
            status = SDreadattr(sd_id_r, i, buf2);
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "      END_OBJECT        = RANGEENDINGTIME\n");

        strcat(buffer, "      OBJECT           = RANGEBEGINNINGDATE\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "RANGEBEGINNINGDATE")) != -1)
            status = SDreadattr(sd_id_r, i, buf2);
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "      END_OBJECT        = RANGEBEGINNINGDATE\n");

        strcat(buffer, "      OBJECT           = RANGEBEGINNINGTIME\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        memset(buf2, 0, sizeof (buf2));
        if ((i = SDfindattr(sd_id_r, "RANGEBEGINNINGTIME")) != -1)
            status = SDreadattr(sd_id_r, i, buf2);
        strcat(buffer, "\"");
        strcat(buffer, buf2);
        strcat(buffer, "\"\n");
        strcat(buffer, "      END_OBJECT        = RANGEBEGINNINGTIME\n");

        strcat(buffer, "    END_GROUP           = RANGEDATETIME\n");

        strcat(buffer, "    GROUP              = ADDITIONALATTRIBUTES\n");

        strcat(buffer, "      OBJECT           = ADDITIONALATTRIBUTESCONTAINER\n");
        strcat(buffer, "        CLASS          = \"1\"\n");

        strcat(buffer, "        OBJECT         = ADDITIONALATTRIBUTENAME\n");
        strcat(buffer, "          CLASS        = \"1\"\n");
        strcat(buffer, "            NUM_VAL    = 1\n");
        strcat(buffer, "            VALUE      = GRANULENUMBER\n");
        strcat(buffer, "        END_OBJECT     = ADDITIONALATTRIBUTENAME\n");

        strcat(buffer, "      END_OBJECT       = ADDITIONALATTRIBUTESCONTAINER\n");
        strcat(buffer, "    END_GROUP          = ADDITIONALATTRIBUTES\n");


        strcat(buffer, "    GROUP              = CollectionDescriptionClass\n");
        strcat(buffer, "      OBJECT           = SHORTNAME\n");
        strcat(buffer, "        NUM_VAL        = 1\n");
        strcat(buffer, "        VALUE          = ");
        strcat(buffer, "\"MYD01SS\"\n");
        strcat(buffer, "      END_OBJECT        = SHORTNAME\n");
        strcat(buffer, "    END_GROUP           = CollectionDescriptionClass\n");

        strcat(buffer, "END_GROUP               = INVENTORYMETADATA\n");
        strcat(buffer, "END\n");
        strcat(buffer, "");

        status = SDsetattr(sd_id_w, "CoreMetadata.0", DFNT_CHAR, strlen(buffer) + 1,
                (VOIDP) buffer);

    }


    SDend(sd_id_r);
    Hclose(HDFfid_r);
    Vend(HDFfid_r);

    SDend(sd_id_w);
    Hclose(HDFfid_w);
    Vend(HDFfid_w);

    return 0;
}

int compare_int16(const void * p1, const void * p2) {
    int16 a = *((int16 *) p1);
    int16 b = *((int16 *) p2);

    if (a > b)
        return 1;
    else if (a < b)
        return -1;
    else
        return 0;
}

int compare_float32(const void * p1, const void * p2) {
    float32 a = *((float32 *) p1);
    float32 b = *((float32 *) p2);

    if (a > b)
        return 1;
    else if (a < b)
        return -1;
    else
        return 0;
}




/*
qsort ((void *) data_1km, nelem, sizeof(int16), compare_int16);
for (k=0; k<nelem; k++) if (data_1km[k] > 0) break;
for (l=0; l<nelem; l++) if (data_1km[nelem-l-1] < 32767) break;
printf("%f\n", (float32) data_1km[(k+(nelem-l-1))/2]);

qsort ((void *) median_250_500, nelem, sizeof(float32), compare_float32);
for (k=0; k<nelem; k++) if (median_250_500[k] > 0) break;
for (l=0; l<nelem; l++) if (median_250_500[nelem-l-1] < 32767) break;
printf("%f\n", median_250_500[(k+(nelem-l-1))/2]);

f32 = 1. / (res_b[0]*data_1km[(k+(nelem-l-1))/2]/median_250_500[(k+(nelem-l-1))/2]+ 
            1000*res_b[1]/median_250_500[(k+(nelem-l-1))/2] + const_b);
status = SDsetattr(sd_id_w, "B_TC_RS", DFNT_FLOAT32, 1, (VOIDP) &f32);
printf("%f\n",f32);
 */



/*
if (masking == 1) {
  if (argc == 6) {
    scan_incr = atoi(argv[4]);
    pixl_incr = atoi(argv[5]);
  } else {
    printf("Scan and Pixel Increments must be defined for masking.\n");
    exit(-1);
  }
}
 */


#if 0
/* Flag skipped pixels and scans (if masking on) */
/* --------------------------------------------- */
if (strncmp(buffer, "EV", 2) == 0 && masking == 1) {
    printf("Masking %s\n", buffer);
    for (k = 0; k < dims[0]; k++) {
        if (k % scan_incr != 0) {
            memset(&data[k * 2 * dims[1] * dims[2]], 255, 2 * dims[1] * dims[2]);
        } else {

            for (l = 0; l < dims[1]; l++)
                for (i = 0; i < dims[2]; i++)
                    if ((i % pixl_incr) != 0)
                        memset(&data[k * 2 * dims[1] * dims[2] + l * 2 * dims[2] + i * 2], 255, 2);
        }
    }
}
#endif
