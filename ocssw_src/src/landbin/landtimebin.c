/*
    Modification history:
    Programmer       Organization      Date      Description of change
    --------------   ------------    --------    ---------------------
    Joel Gales       Futuretech      05/07/01    Original Development

 */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "hdf.h"
#include "mfhdf.h"

int main(int argc, char *argv[]) {

    int32 i;
    int32 j;
    int32 k;
    int32 sds_id;
    int32 dims[8];
    int32 rank;
    int32 dtype;
    int32 nattrs;

    static char buffer[2048];
    static char attr_buf[2048];
    char *data;
    char *data_w;
    char *str_ptr;
    char *comma;
    char *tmp_str;
    int16 *eval_flds;
    int16 * fileindex[3];
    unsigned char *n_water_pix;

    int32 nfiles;
    int32 scan;
    int32 nscan;
    int32 ncol;
    int32 offset;
    int32 *HDFfid_r;
    int32 *sd_id_r;
    int32 HDFfid_w;
    int32 sd_id_w;
    int32 sds_id_w;
    int32 ndatasets;
    int32 nglobal_attr;
    int32 count;
    int32 dim_id_r;
    int32 dim_id_w;
    int32 start[2] = {0, 0};
    int32 edges[2] = {0, 0};
    int32 mask_row_off = 0;
    int32 mask_col_off = 0;
    int32 sub_scan[2] = {0, 0};
    int32 sub_pixl[2] = {0, 0};
    int32 water_min = 1000;

    int16 fillvalue = -32767;
    int16 water_val = -16384;
    int16 eval_flds_flag=0;

    int16 f_index;
    int16 f_ndvi;
    int16 f_blue;

    int16 icol, irow;
    float32 NIRreldiff, NIRness;
    float64 val1, val2, sum, weigh;
    float64 v1, v2;


    int16 maxNDVI = 0;
    int16 minBLUE = 1;

    int16 RADIUS = 2;
    int16 NEVALFLD = 4;
    int16 FILE_BSIZ = NEVALFLD * (2 * RADIUS + 1);
    int16 FIELD_BSIZ = 2 * RADIUS + 1;

    int16 NDVI = 0 * FIELD_BSIZ;
    int16 BLUE = 1 * FIELD_BSIZ;
    int16 RED = 2 * FIELD_BSIZ;
    int16 NIR = 3 * FIELD_BSIZ;

    int32 eval_flds_index[NEVALFLD];

    time_t tnow;
    struct tm *tmnow;

    FILE *fp, *mask_4km;

    double R3(double x);

    setlinebuf(stdout);

    printf("landtimebin V2.31 (02/12/02)\n");

    if (argc == 1) {
        printf("landtimebin <File containing list of input landbin files> <Output landbin file> <NDVI/BLUE/HYBRID>");
        return 0;
    }

    /* Determine number of input files */
    /* ------------------------------- */
    nfiles = 0;
    fp = fopen(argv[1], "r");
    if (fp == NULL) {
        printf("Input listing file: \"%s\" not found.\n", argv[1]);
        return -1;
    }
    while (fgets(buffer, 256, fp) != NULL) nfiles++;
    fclose(fp);
    printf("%d input files\n", nfiles);


    /* Open input files and initialize HDF interface */
    /* --------------------------------------------- */
    HDFfid_r = (int32 *) calloc(nfiles, sizeof (int32));
    sd_id_r = (int32 *) calloc(nfiles, sizeof (int32));

    fp = fopen(argv[1], "r");
    for (i = 0; i < nfiles; i++) {
        fgets(buffer, 256, fp);
        buffer[strlen(buffer) - 1] = 0;
        HDFfid_r[i] = Hopen(buffer, DFACC_READ, 0);
        if (HDFfid_r[i] == -1) {
            printf("Input file: \"%s\" not found.\n", buffer);
            return -1;
        }
        Vstart(HDFfid_r[i]);
        sd_id_r[i] = SDstart(buffer, DFACC_RDONLY);

        j = SDfindattr(sd_id_r[i], "Subset Scan Range");
        if (j != -1) SDreadattr(sd_id_r[i], j, (VOIDP) & sub_scan);

        j = SDfindattr(sd_id_r[i], "Subset Pixel Range");
        if (j != -1) SDreadattr(sd_id_r[i], j, (VOIDP) & sub_pixl);
    }
    fclose(fp);


    /* Create output HDF file */
    /* ---------------------- */
    HDFfid_w = Hopen(argv[2], DFACC_CREATE, 0);
    Vstart(HDFfid_w);
    sd_id_w = SDstart(argv[2], DFACC_RDWR);


    /* Set compositing flag */
    /* -------------------- */
    if (argc == 3 || (strcmp(argv[3], "NDVI") == 0))
        eval_flds_flag = 0;
    else if (strcmp(argv[3], "BLUE") == 0)
        eval_flds_flag = 1;
    else if (strcmp(argv[3], "HYBRID") == 0)
        eval_flds_flag = 2;

    printf("eval_flds_flag: %d\n", eval_flds_flag);


    /* Get # of scans and allocate dynamic arrays */
    /* ------------------------------------------ */
    sds_id = SDselect(sd_id_r[0], 0);
    SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);
    SDendaccess(sds_id);
    nscan = dims[0];
    ncol = dims[1];
    data = (char *) calloc(nfiles * ncol * 2, sizeof (char));
    data_w = (char *) calloc(ncol * 2, sizeof (char));
    eval_flds = (int16 *) calloc(nfiles * ncol * FILE_BSIZ, sizeof (int16));

    n_water_pix = (unsigned char *) calloc(8640 * FIELD_BSIZ, sizeof (unsigned char));

    printf("%d\n", argc);
    if (argc > 4) {
        water_min = atol(argv[4]);
    }
    printf("water_min: %4d\n", water_min);

    fileindex[0] = (int16 *) calloc(ncol*FIELD_BSIZ, sizeof (int16));
    fileindex[1] = (int16 *) calloc(ncol*FIELD_BSIZ, sizeof (int16));
    fileindex[2] = (int16 *) calloc(ncol*FIELD_BSIZ, sizeof (int16));


    /* Determine number of datasets in input files */
    /* ------------------------------------------- */
    SDfileinfo(sd_id_r[0], &ndatasets, &nglobal_attr);


    /* For each dataset in input file, create corresponding dataset in output file */
    /* --------------------------------------------------------------------------- */
    for (j = 0; j < ndatasets; j++) {
        sds_id = SDselect(sd_id_r[0], j);
        SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);
        sds_id_w = SDcreate(sd_id_w, buffer, dtype, rank, dims);


        /* Find index number of compositing dataset */
        /* ---------------------------------------- */
        if (strcmp(buffer, "NDVI") == 0) eval_flds_index[0] = j;
        if (strcmp(buffer, "refl_443") == 0) eval_flds_index[1] = j;
        if (strcmp(buffer, "refl_670") == 0) eval_flds_index[2] = j;
        if (strcmp(buffer, "refl_865") == 0) eval_flds_index[3] = j;
        //if (strcmp(buffer, "input_file") == 0) input_file_index = j;


        /* Copy dataset attributes */
        /* ----------------------- */
        for (i = 0; i < nattrs; i++) {
            SDreadattr(sds_id, i, (VOIDP) attr_buf);
            SDattrinfo(sds_id, i, buffer, &dtype, &count);
            SDsetattr(sds_id_w, buffer, dtype, count, (VOIDP) attr_buf);
        }


        /* Set output file dim names */
        /* ------------------------- */
        for (i = 0; i < rank; i++) {
            dim_id_r = SDgetdimid(sds_id, i);
            dim_id_w = SDgetdimid(sds_id_w, i);
            SDdiminfo(dim_id_r, buffer, &count, &dtype, &nattrs);
            SDsetdimname(dim_id_w, buffer);
        }

        SDendaccess(sds_id_w);
        SDendaccess(sds_id);
    }


    /* Open 4km water pixel mask */
    /* ------------------------- */
    i = SDfindattr(sd_id_r[0], "Projection");
    SDattrinfo(sd_id_r[0], i, buffer, &dtype, &count);
    SDreadattr(sd_id_r[0], i, (VOIDP) attr_buf);

    tmp_str = getenv("OCDATAROOT");
    if (tmp_str == 0x0) {
        printf("Environment variable OCDATAROOT not defined.\n");
        exit(1);
    }
    strcpy(buffer, tmp_str);

    if (strcmp(attr_buf, "Plate Carree") == 0) {
        strcat(buffer, "/common/water_pix_4km_platte_carre.dat");
        mask_4km = fopen(buffer, "rb");
        if (mask_4km == 0x0) {
            printf("Water mask file: %s not found.\n", buffer);
            exit(1);
        }
    }

    if (strcmp(attr_buf, "Sinusoidal") == 0) {
        strcat(buffer, "/common/water_pix_4km_sinusoidal.dat");
        mask_4km = fopen(buffer, "rb");
        if (mask_4km == 0x0) {
            printf("Water mask file: %s not found.\n", buffer);
            exit(1);
        }
    }


    i = SDfindattr(sd_id_r[0], "Command line");
    SDattrinfo(sd_id_r[0], i, buffer, &dtype, &count);
    SDreadattr(sd_id_r[0], i, (VOIDP) attr_buf);
    str_ptr = strstr(attr_buf, "-box");
    if (str_ptr != 0x0) {
        str_ptr += 5;
        comma = strstr(str_ptr, ",");
        *comma = 0;
        mask_col_off = atoi(str_ptr);
        str_ptr = comma + 1;
        comma = strstr(str_ptr, ",");
        *comma = 0;
        mask_row_off = atoi(str_ptr);
    }


    /* For each scan ... */
    /* ----------------- */
    offset = RADIUS;
    edges[1] = ncol;

    for (scan = 0; scan < nscan; scan++) {

        if ((scan % 200) == 0) {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("Scan:%6d %s", scan, asctime(tmnow));
        }

        /* Read compositing dataset from each input file */
        /* --------------------------------------------- */
        if (scan == 0) {
            start[0] = scan;
            edges[0] = RADIUS + 1;
        } else if (scan < nscan - 2) {
            start[0] = scan + RADIUS;
            edges[0] = 1;
        } else {
            goto DATA_WRITE;
        }

        for (i = 0; i < nfiles; i++) {
            for (j = 0; j < NEVALFLD; j++) {
                sds_id = SDselect(sd_id_r[i], eval_flds_index[j]);
                SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

                SDreaddata(sds_id, start, NULL, edges,
                        (VOIDP) & eval_flds[(FILE_BSIZ * i + FIELD_BSIZ * j + offset) * ncol]);
            }
        } /* file loop */



        /* Read 4km land/water/mixed mask */
        /* ------------------------------ */
        fseek(mask_4km, (mask_row_off + start[0]) * 8640, SEEK_SET);
        fread(&n_water_pix[8640 * offset], edges[0] * 8640, 1, mask_4km);


        /* Determine file index for compositing */
        /* ------------------------------------ */
        for (k = offset; k < FIELD_BSIZ; k++) {
            for (j = 0; j < ncol; j++) {
                fileindex[maxNDVI][ncol * k + j] = 0;
                fileindex[minBLUE][ncol * k + j] = 0;
            }
        }


        /* Max NDVI fileindex */
        /* ------------------ */
        for (k = offset; k < FIELD_BSIZ; k++) {
            for (i = 1; i < nfiles; i++) {
                for (j = 0; j < ncol; j++) {
                    f_index = fileindex[maxNDVI][ncol * k + j];
                    if (eval_flds[ncol * (FILE_BSIZ * i + NDVI + k) + j] >
                            eval_flds[ncol * (FILE_BSIZ * f_index + NDVI + k) + j])
                        fileindex[maxNDVI][ncol * k + j] = i;
                }
            }
        }

        /* Min Blue fileindex */
        /* ------------------ */
        for (k = offset; k < FIELD_BSIZ; k++) {
            for (i = 1; i < nfiles; i++) {
                for (j = 0; j < ncol; j++) {
                    f_index = fileindex[minBLUE][ncol * k + j];

                    if (eval_flds[ncol * (FILE_BSIZ * f_index + BLUE + k) + j] < 0 &&
                            eval_flds[ncol * (FILE_BSIZ * i + BLUE + k) + j] >= 0)
                        fileindex[minBLUE][ncol * k + j] = i;

                    if (eval_flds[ncol * (FILE_BSIZ * f_index + BLUE + k) + j] >= 0 &&
                            eval_flds[ncol * (FILE_BSIZ * i + BLUE + k) + j] >= 0 &&
                            eval_flds[ncol * (FILE_BSIZ * i + RED + k) + j] >= 0 &&
                            (eval_flds[ncol * (FILE_BSIZ * i + BLUE + k) + j] <
                            eval_flds[ncol * (FILE_BSIZ * f_index + BLUE + k) + j]))
                        fileindex[minBLUE][ncol * k + j] = i;

                    if (eval_flds[ncol * (FILE_BSIZ * f_index + RED + k) + j] < 0 &&
                            eval_flds[ncol * (FILE_BSIZ * i + RED + k) + j] < 0 &&
                            (eval_flds[ncol * (FILE_BSIZ * i + RED + k) + j] >
                            eval_flds[ncol * (FILE_BSIZ * f_index + RED + k) + j]))
                        fileindex[minBLUE][ncol * k + j] = i;

                }
            }
        }


        /* Fill Value Handling */
        /* ------------------- */
        for (k = offset; k < FIELD_BSIZ; k++) {
            for (j = 0; j < ncol; j++) {

                f_index = fileindex[maxNDVI][ncol * k + j];
                if (eval_flds[ncol * (FILE_BSIZ * f_index + NDVI + k) + j] == -32767)
                    fileindex[maxNDVI][ncol * k + j] = fillvalue;

                f_index = fileindex[minBLUE][ncol * k + j];
                if (eval_flds[ncol * (FILE_BSIZ * f_index + BLUE + k) + j] < 0)
                    fileindex[minBLUE][ncol * k + j] = fillvalue;
            }
        }



        if (eval_flds_flag != 2) goto DATA_WRITE;


        /* HYBRID */
        /* ------ */
        for (j = 0; j < ncol; j++) {

            f_ndvi = fileindex[maxNDVI][ncol * RADIUS + j];
            f_blue = fileindex[minBLUE][ncol * RADIUS + j];

            fileindex[2][ncol * RADIUS + j] = f_blue;

            if (f_blue != fillvalue &&
                    f_ndvi != fillvalue &&
                    n_water_pix[(8640 * offset)+(mask_col_off + j)] == 0 &&
                    eval_flds[ncol * (FILE_BSIZ * f_ndvi + NIR + RADIUS) + j] != fillvalue &&
                    eval_flds[ncol * (FILE_BSIZ * f_blue + NIR + RADIUS) + j] != fillvalue &&
                    eval_flds[ncol * (FILE_BSIZ * f_ndvi + BLUE + RADIUS) + j] != fillvalue) {

                if (eval_flds[ncol * (FILE_BSIZ * f_blue + NIR + RADIUS) + j] != 0) {
                    NIRreldiff = ((eval_flds[ncol * (FILE_BSIZ * f_ndvi + NIR + RADIUS) + j] -
                            eval_flds[ncol * (FILE_BSIZ * f_blue + NIR + RADIUS) + j]) /
                            abs(eval_flds[ncol * (FILE_BSIZ * f_blue + NIR + RADIUS) + j]));
                } else {
                    NIRreldiff = 1.0;
                }

                if (eval_flds[ncol * (FILE_BSIZ * f_blue + BLUE + RADIUS) + j] != 0 &&
                        eval_flds[ncol * (FILE_BSIZ * f_ndvi + BLUE + RADIUS) + j] != 0) {

                    if (NIRreldiff > 0.10) {
                        val1 = 0;
                        val2 = 0;
                        sum = 0;

                        for (irow = 0; irow < 2 * RADIUS + 1; irow++) {

                            for (icol = j - RADIUS; icol <= j + RADIUS; icol++) {
                                if (icol < 0 || icol > ncol) break;

                                weigh = R3(fabs(0.5 * (irow - RADIUS))) * R3(fabs(0.5 * (icol - j)));
                                sum += weigh;

                                f_ndvi = fileindex[maxNDVI][ncol * irow + icol];
                                f_blue = fileindex[minBLUE][ncol * irow + icol];

                                if (f_ndvi != fillvalue)
                                    v1 = weigh * eval_flds[ncol * (FILE_BSIZ * f_ndvi + NIR + irow) + icol];
                                else
                                    v1 = 0.0;

                                if (f_blue != fillvalue)
                                    v2 = weigh * eval_flds[ncol * (FILE_BSIZ * f_blue + NIR + irow) + icol];
                                else
                                    v2 = 0.0;

                                if (v1 < 0.0) v1 = 0.0;
                                if (v2 < 0.0) v2 = 0.0;

                                val1 += v1;
                                val2 += v2;

                            } /* icol loop */
                        } /* irow loop */

                        if (sum > 0) {
                            val1 /= sum;
                            val2 /= sum;
                        }


                        f_ndvi = fileindex[maxNDVI][ncol * RADIUS + j];
                        f_blue = fileindex[minBLUE][ncol * RADIUS + j];

                        NIRness = (eval_flds[ncol * (FILE_BSIZ * f_ndvi + NIR + RADIUS) + j] -
                                eval_flds[ncol * (FILE_BSIZ * f_ndvi + BLUE + RADIUS) + j]) /
                                abs(eval_flds[ncol * (FILE_BSIZ * f_ndvi + BLUE + RADIUS) + j]);

                        if (fabs(NIRness) > 0.2 &&
                                eval_flds[ncol * (FILE_BSIZ * f_ndvi + BLUE + RADIUS) + j] < 2500 &&
                                fabs(val1 - eval_flds[ncol * (FILE_BSIZ * f_ndvi + NIR + RADIUS) + j]) <=
                                fabs(val2 - eval_flds[ncol * (FILE_BSIZ * f_blue + NIR + RADIUS) + j]) &&

                                (eval_flds[ncol * (FILE_BSIZ * f_ndvi + NDVI + RADIUS) + j] < 3000 ||
                                eval_flds[ncol * (FILE_BSIZ * f_ndvi + BLUE + RADIUS) + j] < 500 ||
                                eval_flds[ncol * (FILE_BSIZ * f_blue + NIR + RADIUS) + j] == 0)) {
                            fileindex[2][ncol * RADIUS + j] = f_ndvi;
                        }

                    } /* (NIRreldiff > 0.10) */
                }
            }

        } /* HYBRID col loop */



    DATA_WRITE:

        start[0] = scan;
        edges[0] = 1;

        /* For each dataset ... */
        /* -------------------- */
        for (j = 0; j < ndatasets; j++) {
            sds_id = SDselect(sd_id_r[0], j);
            SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);
            SDendaccess(sds_id);

            if (strcmp(buffer, "input_file") == 0) {
                memcpy(&data_w[0], &fileindex[eval_flds_flag][ncol * RADIUS], 2 * ncol);
            } else {
                /* For each file read dataset */
                /* -------------------------- */
                for (i = 0; i < nfiles; i++) {
                    sds_id = SDselect(sd_id_r[i], j);
                    edges[1] = ncol;
                    SDreaddata(sds_id, start, NULL, edges,
                            (VOIDP) & data[2 * i * ncol]);
                    SDendaccess(sds_id);
                }

                /* Store compositing data in write buffer */
                /* -------------------------------------- */
                for (i = 0; i < ncol; i++) {
                    if (fileindex[eval_flds_flag][ncol * RADIUS + i] == fillvalue) {
                        memcpy(&data_w[2 * i], &fillvalue, 2);
                    }
                        /* mask water and mixed with water val = -16384 (JMG) */

                    else if (n_water_pix[(8640 * offset)+(mask_col_off + i)] >= water_min) {
                        memcpy(&data_w[2 * i], &water_val, 2);
                    } else {
                        memcpy(&data_w[2 * i], &data[2 * (fileindex[eval_flds_flag][ncol * RADIUS + i] * ncol + i)], 2);
                    }
                }
            } /* if "input_file" */

            /* Write to output file */
            sds_id_w = SDselect(sd_id_w, j);
            SDwritedata(sds_id_w, start, NULL, edges, (VOIDP) data_w);
            SDendaccess(sds_id_w);
        } /* dataset loop */


        /* Rotate Arrays */
        /* ------------- */
        for (k = 1; k < FIELD_BSIZ; k++) {
            for (i = 0; i < nfiles; i++) {
                memcpy(&eval_flds[ncol * (FILE_BSIZ * i + NDVI + (k - 1))],
                        &eval_flds[ncol * (FILE_BSIZ * i + NDVI + k)],
                        ncol * sizeof (int16));

                memcpy(&eval_flds[ncol * (FILE_BSIZ * i + BLUE + (k - 1))],
                        &eval_flds[ncol * (FILE_BSIZ * i + BLUE + k)],
                        ncol * sizeof (int16));

                memcpy(&eval_flds[ncol * (FILE_BSIZ * i + NIR + (k - 1))],
                        &eval_flds[ncol * (FILE_BSIZ * i + NIR + k)],
                        ncol * sizeof (int16));
            }

            memcpy(&fileindex[maxNDVI][ncol * (k - 1)],
                    &fileindex[maxNDVI][ncol * k],
                    ncol * sizeof (int16));

            memcpy(&fileindex[minBLUE][ncol * (k - 1)],
                    &fileindex[minBLUE][ncol * k],
                    ncol * sizeof (int16));

            memcpy(&n_water_pix[8640 * (k - 1)],
                    &n_water_pix[8640 * k],
                    8640 * sizeof (int16));
        }

        offset = 2 * RADIUS;


    } /* scan loop */


    /* Write Global Attributes */
    /* ----------------------- */
    i = SDfindattr(sd_id_r[0], "Projection");
    SDattrinfo(sd_id_r[0], i, buffer, &dtype, &count);
    SDreadattr(sd_id_r[0], i, (VOIDP) attr_buf);
    SDsetattr(sd_id_w, buffer, dtype, count, (VOIDP) attr_buf);

    i = SDfindattr(sd_id_r[0], "Pixel size (meters)");
    SDattrinfo(sd_id_r[0], i, buffer, &dtype, &count);
    SDreadattr(sd_id_r[0], i, (VOIDP) attr_buf);
    SDsetattr(sd_id_w, buffer, dtype, count, (VOIDP) attr_buf);

    i = SDfindattr(sd_id_r[0], "Compositing criterion");
    SDattrinfo(sd_id_r[0], i, buffer, &dtype, &count);
    SDreadattr(sd_id_r[0], i, (VOIDP) attr_buf);
    SDsetattr(sd_id_w, buffer, dtype, count, (VOIDP) attr_buf);

    if ((sub_scan[0] != 0) || (sub_scan[1] != 0) ||
            (sub_pixl[0] != 0) || (sub_pixl[1] != 0)) {
        SDsetattr(sd_id_w, "Subset Scan Range", DFNT_INT32, 2, sub_scan);
        SDsetattr(sd_id_w, "Subset Pixel Range", DFNT_INT32, 2, &sub_pixl);
    }

    fp = fopen(argv[1], "r");
    attr_buf[0] = 0;
    while (fgets(buffer, 256, fp) != NULL) {
        buffer[strlen(buffer) - 1] = 0;
        strcat(buffer, "\n,");
        strcat(attr_buf, buffer);
    }
    fclose(fp);
    attr_buf[strlen(attr_buf) - 1] = 0;
    SDsetattr(sd_id_w, "File list", DFNT_CHAR8, (int) strlen(attr_buf) + 1, attr_buf);



    /* Close input files */
    /* ----------------- */
    for (i = 0; i < nfiles; i++) {
        SDend(sd_id_r[i]);
        Hclose(HDFfid_r[i]);
        Vend(HDFfid_r[i]);
    }
    free(HDFfid_r);
    free(sd_id_r);
    free(data);
    free(data_w);
    free(eval_flds);
    free(fileindex[0]);
    free(fileindex[1]);
    free(fileindex[2]);
    free(n_water_pix);


    /* Close output file */
    /* ----------------- */
    SDend(sd_id_w);
    Hclose(HDFfid_w);
    Vend(HDFfid_w);

    fclose(mask_4km);

    return 0;
}

double R3(double x) {
    if (x < 0 || x >= 2)
        return 0;
    else if (x <= 1) /* 0 <= x <= 1 */
        return 2 / 3. + x * x * x / 2. - x * x;
    else /* 1 <= x <= 2 */
        return (2 - x) * (2 - x) * (2 - x) / 6;
}


