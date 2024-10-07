/*
    Modification history:
    Programmer       Organization      Date      Description of change
    --------------   ------------    --------    ---------------------
    Joel Gales       Futuretech      08/10/03    Original Development

 */

/*
  Revision 2.0.2 09/13/22
  Add support for CLO
  B. Yang

  Revision 2.0.2 12/06/17
  Add support for prodlist for netcdf4 files
  J. Gales

  Revision 1.000 01/08/13
  Don't read control points if lonlatinterp == 0
  J. Gales

  Revision 0.993 05/21/10
  Fix east/westmost lon for extracts crossing dateline
  J. Gales

  Revision 0.992 09/26/08
  Only print prodlist (arg 9) if it exists
  J. Gales

  Revision 0.991 06/18/08
  Don't extract on "pxls" dimention for OCTS
  J. Gales

  Revision 0.990 12/05/06
  Fix slat,clat,...,slon,clon,...,& lon/lat boundary metadata
  J. Gales

  Revision 0.983 08/08/06
  Make cntl_pt_cols 1-based
  J. Gales

  Revision 0.982 06/15/06
  Set octs only if Level 1
  J. Gales

  Revision 0.981 06/01/06
  Make sure sscan is odd and escan is even if OCTS
  J. Gales

  Revision 0.98 05/30/06
  Fix problem with OCTS L1A extraction
  J. Gales
 */


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "mfhdf.h"
#include <genutils.h>

#include "l2extract.h"
#include "hdf.h"
#include <timeutils.h>

clo_optionList_t *optionList = NULL;

int main(int argc, char *argv[]) {
    int32 i;
    int32 j;
    int32 k;
    int32 n;
    int32 jds;
    int32 sds_id;
    int32 dims[8];
    int32 dims_w[256][8];
    int32 dims_lonlat[2];
    int32 rank;
    int32 nelem;
    int32 dtype;
    int32 nattrs;
    int32 status;
    int32 LAC_start;
    int32 i32;
    int32 lonlatinterp = 1;

    int32 zero = 0;
    // int32 one = 1;

    int32 vg_ref;
    int32 vgid_r;
    int32 vgid_w;
    int32 tag;
    int32 ref;
    // int32 listlen = 0;
    int32 npixl;
    int32 nscans;
    // int32 nread;
    // int32 nwrite;
    int32 sscan, escan, nscan;
    int32 spix, epix;
    // int32 maskout;
    // int32 fldsize;
    int32 octs = 0;
    int32 seawifs_l1a = 0;

    int32 startyear = 0, cenyear = 0, endyear = 0;
    int32 startday = 0, cenday = 0, endday = 0;
    int32 startmsec = 0, cenmsec = 0, endmsec = 0;
    float64 utime;

    float32 northern_lat = -90.0;
    float32 southern_lat = +90.0;
    float32 western_lon = +180.0;
    float32 eastern_lon = -180.0;

    static char buffer[2048];
    static char buf2[2048];
    static char prodlist[2048];
    static char geophyslist[2048];
    static char wavelist[2048];

    static char fieldname[FIELDNAMELENMAX];
    static char attrname[FIELDNAMELENMAX];
    char *infile, *outfile;  // , parm_list[1024];
    char *prodlist_str = "";
    char *waves_str = "";

    char *data;
    char *data2;

    int32 HDFfid_r;
    int32 sd_id_r;
    int32 HDFfid_w;
    int32 sd_id_w;
    int32 sds_id_w;
    int32 ndatasets;
    int32 nglobal_attr;
    int32 numPositionOptions;
    int32 count;
    int32 dim_id_r;
    int32 dim_id_w;
    int32 start[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int32 edge[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int32 start_w[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // int32 maxsize;
    int32 n_cntl_pnts;
    int32 pix_sub, sc_sub;
    float32 f32, delta;

    int32 ntilts;
    int16 old_tilt_flags[20];
    int16 old_tilt_ranges[20][2];
    int16 new_tilt_flags[20];
    int16 new_tilt_ranges[20][2];
    int16 stilt, etilt;

    int32 nvgrp = 3;
    char *vgrpname[] = {"Scan-Line Attributes", "Geophysical Data", "Navigation Data"};

    float32 lllon, lllat, lrlon, lrlat, ullon, ullat, urlon, urlat;
    float32 sclon, sclat, cclon, cclat, eclon, eclat;
    float32 csolz;

    char cntl_pt_cols_buf[4 * 4096];
    char *ver = "3.1";
    // char proc_con[2048];
    char title[128];

    float32 f32_buf[4096];
    float32 spline_arr[4096];

    uint8 datelinecross = 0;

    printf("This is version %s of %s (compiled on %s %s)\n", ver, "l2extract", __DATE__, __TIME__);

    /*** load input parameters into local variables */

    clo_setEnablePositionOptions(1);

    optionList = clo_createList();

    l2extract_init_options(optionList, ver);
    if (argc == 1) {
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }

    l2extract_read_options(optionList, argc, argv);
    numPositionOptions = clo_getPositionNumOptions(optionList);

    if (numPositionOptions == 0) {
        infile = clo_getString(optionList, "ifile");
        outfile = clo_getString(optionList, "ofile");
        spix = clo_getInt(optionList, "spix");
        epix = clo_getInt(optionList, "epix");
        sscan = clo_getInt(optionList, "sline");
        escan = clo_getInt(optionList, "eline");
        // pix_sub = clo_getInt(optionList, "pix_sub");
        // sc_sub = clo_getInt(optionList, "sc_sub");
        pix_sub = 1;
        sc_sub = 1;

        if (clo_isSet(optionList, "product")) {
            prodlist_str = clo_getRawString(optionList, "product");
        }
        if (clo_isSet(optionList, "wavelist")) {
            waves_str = clo_getRawString(optionList, "wavelist");
        }
    } else if (numPositionOptions == 8) {
        infile = clo_getPositionString(optionList, 0);
        spix = atoi(clo_getPositionString(optionList, 1));
        epix = atoi(clo_getPositionString(optionList, 2));
        sscan = atoi(clo_getPositionString(optionList, 3));
        escan = atoi(clo_getPositionString(optionList, 4));
        pix_sub = atoi(clo_getPositionString(optionList, 5));
        sc_sub = atoi(clo_getPositionString(optionList, 6));
        outfile = clo_getPositionString(optionList, 7);
    } else if (numPositionOptions == 9) {
        infile = clo_getPositionString(optionList, 0);
        spix = atoi(clo_getPositionString(optionList, 1));
        epix = atoi(clo_getPositionString(optionList, 2));
        sscan = atoi(clo_getPositionString(optionList, 3));
        escan = atoi(clo_getPositionString(optionList, 4));
        pix_sub = atoi(clo_getPositionString(optionList, 5));
        sc_sub = atoi(clo_getPositionString(optionList, 6));
        outfile = clo_getPositionString(optionList, 7);
        prodlist_str = clo_getPositionString(optionList, 8);
    } else {
        printf("\nNote: Number of Positional Arguments should be 0, 8 or 9!\n");
        exit(EXIT_FAILURE);
    }

    if (pix_sub != 1 || sc_sub != 1) {
        printf("Subsampling not yet implemented.\n");
        exit(EXIT_FAILURE);
    }

    prodlist[0] = 0;
    if (prodlist_str[0] != 0) {
        strcpy(prodlist, ",l2_flags,");
        strcat(prodlist, prodlist_str);
        strcat(prodlist, ",");
        printf("product: %s\n", prodlist_str);
    }

    if (clo_getBool(optionList, "verbose")) {
        want_verbose = 1;
    }
    if (waves_str[0] != 0) {
        strcat(wavelist, waves_str);
        if (want_verbose)
            printf("Wavelength requested by the user: %s\n", waves_str);
    }
    // if it is not HDF4 then try netCDF
    if (!Hishdf(infile)) {
        return extractNetCDF(infile, outfile, spix, epix, sscan, escan, prodlist, wavelist);
    }

    HDFfid_r = Hopen(infile, DFACC_READ, 0);
    status = Vstart(HDFfid_r);
    sd_id_r = SDstart(infile, DFACC_RDONLY);

    /* Get number of scan lines */
    /* ------------------------ */
    status = SDreadattr(sd_id_r, SDfindattr(sd_id_r, "Number of Scan Lines"), (VOIDP)&i32);

    /* Read Title, Check if OCTS L1 */
    /* ---------------------------- */
    status = SDreadattr(sd_id_r, SDfindattr(sd_id_r, "Title"), (VOIDP)title);
    if (strncmp(title, "OCTS", 4) == 0)
        octs = 1;
    if (strstr(title, "Level-2") != 0)
        octs = 0;
    if (strcmp(title, "SeaWiFS Level-1A Data") == 0)
        seawifs_l1a = 1;

    /* If OCTS make sure sscan is odd and escan is even */
    /* ------------------------------------------------ */
    if (octs == 1) {
        if ((sscan % 2) == 0)
            sscan -= 1;
        if ((escan % 2) == 1)
            escan -= 1;
        printf("Actual OCTS sscan: %d\n", sscan);
        printf("Actual OCTS escan: %d\n", escan);
    }

    /* Determine number of datasets in input files */
    SDfileinfo(sd_id_r, &ndatasets, &nglobal_attr);

    if (sscan == -1 && spix == -1)
        lonlatinterp = 0;

    if (sscan == -1) {
        sscan = 1;
        escan = i32;
    }

    if (escan > i32 * (octs + 1)) {
        printf("escan: %d greater than # of scan lines: %d\n", escan, i32);
        exit(BOUNDS_ERROR);
    }

    status = SDreadattr(sd_id_r, SDfindattr(sd_id_r, "Pixels per Scan Line"), (VOIDP)&i32);

    if (spix == -1) {
        spix = 1;
        epix = i32;
    }

    if (epix > i32) {
        printf("epix: %d greater than # of pixels per scan: %d\n", epix, i32);
        exit(BOUNDS_ERROR);
    }

    printf("sscan: %d  escan: %d\n", sscan, escan);
    nscan = escan - sscan + 1;

    /* Create output HDF file */
    /* ---------------------- */
    HDFfid_w = Hopen(outfile, DFACC_CREATE, 0);
    status = Vstart(HDFfid_w);
    sd_id_w = SDstart(outfile, DFACC_RDWR);

    /* For each dataset (CREATE) ... */
    /* ----------------------------- */
    vg_ref = Vfind(HDFfid_r, "Geophysical Data");
    vgid_r = Vattach(HDFfid_r, vg_ref, "r");
    strcpy(geophyslist, ",");
    for (i = 0; i < Vntagrefs(vgid_r); i++) {
        status = Vgettagref(vgid_r, i, &tag, &ref);
        sds_id = SDselect(sd_id_r, SDreftoindex(sd_id_r, ref));
        status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);
        SDendaccess(sds_id);
        strcat(geophyslist, buffer);
        strcat(geophyslist, ",");
    }
    Vdetach(vgid_r);

    for (jds = 0; jds < ndatasets; jds++) {
        sds_id = SDselect(sd_id_r, jds);
        status = SDgetinfo(sds_id, fieldname, &rank, dims, &dtype, &nattrs);

        printf("Name: %s  rank: %d  type: %d\n", fieldname, rank, dtype);

        /* Check if L2 product field & in prodlist */
        /* --------------------------------------- */
        if (prodlist[0] != 0) {
            strcpy(buffer, ",");
            strcat(buffer, fieldname);
            strcat(buffer, ",");

            if (strstr(geophyslist, buffer) == NULL)
                goto LBL1;
            if (strstr(prodlist, buffer) == NULL)
                continue;
        }

    LBL1:

        /* Determine extract dimensions */
        /* ---------------------------- */
        for (i = 0; i < rank; i++) {
            dim_id_r = SDgetdimid(sds_id, i);
            dim_id_w = SDgetdimid(sds_id_w, i);
            SDdiminfo(dim_id_r, buffer, &count, &i32, &i32);

            dims_w[jds][i] = dims[i];

            if (octs == 1) {
                /* OCTS */

                if (strcmp(fieldname, "samp_table") == 0)
                    continue;

                if (strcmp(buffer, "rec") == 0) {
                    dims_w[jds][i] = (escan - sscan + 1) / 2;
                }

                if (strcmp(buffer, "lines") == 0) {
                    dims_w[jds][i] = escan - sscan + 1;
                }

                if (strcmp(buffer, "blines") == 0) {
                    dims_w[jds][i] = 5 * (escan - sscan + 1);
                }

                if (strcmp(buffer, "nsamp") == 0) {
                    dims_w[jds][i] = epix - spix + 1;
                }

                //      if (strcmp(buffer, "pxls") == 0) {
                //  dims_w[jds][i] = (epix - spix + 1) / 20 + 1;
                //}
            } else {
                if (strcmp(buffer, "Number of Scan Lines") == 0) {
                    dims_w[jds][i] = escan - sscan + 1;
                }

                if (strcmp(buffer, "Pixels per Scan Line") == 0) {
                    dims_w[jds][i] = epix - spix + 1;
                }

                if (strcmp(buffer, "Number of Pixel Control Points") == 0) {
                    if (lonlatinterp) {
                        n_cntl_pnts = count;
                        dims_w[jds][i] = epix - spix + 1;
                    }
                }
            }
        }

        for (k = 0; k < rank; k++) printf("%d\n", dims_w[jds][k]);
        printf("---------------\n");

        /* Create extract SDS */
        /* ------------------ */
        sds_id_w = SDcreate(sd_id_w, fieldname, dtype, rank, dims_w[jds]);
        if (sds_id_w == -1) {
            printf("Field: %s cannot be created\n", fieldname);
            exit(EXIT_FAILURE);
        }

        /* Set dimension names */
        /* ------------------- */
        for (i = 0; i < rank; i++) {
            dim_id_r = SDgetdimid(sds_id, i);
            dim_id_w = SDgetdimid(sds_id_w, i);
            SDdiminfo(dim_id_r, buffer, &count, &i32, &i32);
            SDsetdimname(dim_id_w, buffer);
        }

        /* Read & Write SDS attributes */
        /* --------------------------- */
        for (i = 0; i < nattrs; i++) {
            status = SDreadattr(sds_id, i, (VOIDP)buf2);
            status = SDattrinfo(sds_id, i, attrname, &dtype, &count);
            status = SDsetattr(sds_id_w, attrname, dtype, count, (VOIDP)buf2);
        }

        /* Get lon/lat control points */
        /* -------------------------- */
        if (strcmp(fieldname, "cntl_pt_cols") == 0 && lonlatinterp == 1) {
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP)cntl_pt_cols_buf);

            /* Convert cntl pnts from I32 to F32 */
            /* --------------------------------- */
            for (i = 0; i < n_cntl_pnts; i++) {
                memcpy(&i32, &cntl_pt_cols_buf[4 * i], sizeof(int32));
                f32 = (float32)i32;
                memcpy(&cntl_pt_cols_buf[4 * i], &f32, sizeof(float32));
            }
        }

        SDendaccess(sds_id_w);
        SDendaccess(sds_id);
    } /* End Create SD */
    /*  printf("-----------------------------\n");*/

    /* Get #scan lines, #pixels per scan */
    /* --------------------------------- */
    for (i = 0; i < nglobal_attr; i++) {
        status = SDattrinfo(sd_id_r, i, buffer, &dtype, &count);
        data = (char *)calloc(count, DFKNTsize(dtype));
        status = SDreadattr(sd_id_r, i, (VOIDP)data);

        if (strcmp(buffer, "Number of Scan Lines") == 0) {
            memcpy(&nscans, data, 4);
        }

        if (strcmp(buffer, "Pixels per Scan Line") == 0) {
            memcpy(&npixl, data, 4);
        }
        free(data);
    }

    /* Copy Datasets */
    /* ------------- */
    n = 0;

    for (jds = 0; jds < ndatasets; jds++) {
        sds_id = SDselect(sd_id_r, jds);
        status = SDgetinfo(sds_id, fieldname, &rank, dims, &dtype, &nattrs);

        /* Check if L2 product field & in prodlist */
        /* --------------------------------------- */
        if (prodlist[0] != 0) {
            strcpy(buffer, ",");
            strcat(buffer, fieldname);
            strcat(buffer, ",");

            if (strstr(geophyslist, buffer) == NULL)
                goto LBL2;
            if (strstr(prodlist, buffer) == NULL)
                continue;
        }

    LBL2:

        printf("[Copy] %s\n", fieldname);

        sds_id_w = SDselect(sd_id_w, n++);

        /* Determine start & nelem */
        /* ----------------------- */
        for (i = 0; i < rank; i++) {
            start[i] = 0;

            dim_id_w = SDgetdimid(sds_id_w, i);
            SDdiminfo(dim_id_w, buffer, &count, &i32, &nattrs);
            SDsetdimname(dim_id_w, buffer);

            if (octs == 1) {
                /* OCTS */
                if (strcmp(fieldname, "samp_table_") != 0) {
                    if (strcmp(buffer, "rec") == 0) {
                        start[i] = (sscan - 1) / 2;
                    }

                    if (strcmp(buffer, "lines") == 0) {
                        start[i] = sscan - 1;
                    }

                    if (strcmp(buffer, "blines") == 0) {
                        start[i] = 5 * (sscan - 1);
                    }

                    if (strcmp(buffer, "nsamp") == 0) {
                        start[i] = spix - 1;
                    }

                    //    if (strcmp(buffer, "pxls") == 0) {
                    //  start[i] = (spix - 1) / 20;
                    //}
                }
            } else {
                if (strcmp(buffer, "Number of Scan Lines") == 0) {
                    start[i] = sscan - 1;
                }

                if (strcmp(buffer, "Pixels per Scan Line") == 0) {
                    start[i] = spix - 1;
                }
            }

            if (i == 0)
                nelem = dims_w[jds][i];
            else
                nelem *= dims_w[jds][i];
        }

        if (seawifs_l1a == 0) {
            if (strcmp(fieldname, "slat") == 0 || strcmp(fieldname, "slon") == 0 ||
                strcmp(fieldname, "elon") == 0 || strcmp(fieldname, "elat") == 0 ||
                strcmp(fieldname, "clon") == 0 || strcmp(fieldname, "clat") == 0) {
                continue;
            }
        }

        /* Allocate data buffer */
        /* -------------------- */
        if (strcmp(fieldname, "latitude") == 0 || strcmp(fieldname, "longitude") == 0)
            nelem = nscans * npixl;

        data = (char *)calloc(nelem, DFKNTsize(dtype));

        if (strcmp(fieldname, "cntl_pt_cols") == 0 && lonlatinterp == 1) {
            for (i = 0; i < dims_w[jds][0]; i++) {
                i32 = i + 1;
                memcpy(&data[i * 4], &i32, 4);
            }
            /* Latitude */
            /* -------- */
        } else if (strcmp(fieldname, "latitude") == 0 && lonlatinterp == 1) {
            dims_lonlat[0] = dims_w[jds][0];
            dims_lonlat[1] = n_cntl_pnts;

            for (i = 0; i < dims_w[jds][0]; i++) {
                start[0] = i + sscan - 1;
                edge[0] = 1;
                edge[1] = dims_lonlat[1];
                status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)&data[4 * i * npixl]);
            }

            /* Interpolate Latitude */
            /* -------------------- */
            for (j = 0; j < dims_w[jds][0]; j++) {
                spline((float *)cntl_pt_cols_buf, (float *)&data[j * 4 * npixl], n_cntl_pnts, 1e30, 1e30,
                       spline_arr);
                for (i = spix - 1; i < epix; i++) {
                    splint((float *)cntl_pt_cols_buf, (float *)&data[j * 4 * npixl], spline_arr, n_cntl_pnts,
                           i + 1.0, &f32_buf[i - spix + 1]);

                    /* Determine northmost & southmost lat */
                    /* ----------------------------------- */
                    if (f32_buf[i - spix + 1] > northern_lat)
                        northern_lat = f32_buf[i - spix + 1];
                    if (f32_buf[i - spix + 1] < southern_lat)
                        southern_lat = f32_buf[i - spix + 1];
                }

                memcpy(&data[j * 4 * dims_w[jds][1]], f32_buf, 4 * dims_w[jds][1]);
            }

            /* Longitude */
            /* --------- */
        } else if (strcmp(fieldname, "longitude") == 0 && lonlatinterp == 1) {
            dims_lonlat[0] = dims_w[jds][0];
            dims_lonlat[1] = n_cntl_pnts;

            for (i = 0; i < dims_w[jds][0]; i++) {
                start[0] = i + sscan - 1;
                edge[0] = 1;
                edge[1] = dims_lonlat[1];
                status = SDreaddata(sds_id, start, NULL, edge, (VOIDP)&data[4 * i * npixl]);
            }

            /* Interpolate Longitude */
            /* --------------------- */
            for (j = 0; j < dims_w[jds][0]; j++) {
                /* Remove any dateline discontinuity in the longitudes */
                /* --------------------------------------------------- */
                for (i = 1; i < dims_lonlat[1]; i++) {
                    memcpy(&f32_buf[1], &data[4 * (j * npixl + i)], 4);
                    memcpy(&f32_buf[0], &data[4 * (j * npixl + i - 1)], 4);
                    delta = f32_buf[1] - f32_buf[0];
                    if (delta < -180) {
                        f32_buf[1] += 360;
                        datelinecross = 1;
                    } else if (delta > 180) {
                        f32_buf[1] -= 360;
                        datelinecross = 1;
                    }
                    memcpy(&data[4 * (j * npixl + i)], &f32_buf[1], 4);
                }

                spline((float *)cntl_pt_cols_buf, (float *)&data[j * 4 * npixl], n_cntl_pnts, 1e30, 1e30,
                       spline_arr);
                for (i = spix - 1; i < epix; i++) {
                    splint((float *)cntl_pt_cols_buf, (float *)&data[j * 4 * npixl], spline_arr, n_cntl_pnts,
                           i + 1.0, &f32_buf[i - spix + 1]);
                    if (f32_buf[i - spix + 1] > +180)
                        f32_buf[i - spix + 1] -= 360;
                    if (f32_buf[i - spix + 1] < -180)
                        f32_buf[i - spix + 1] += 360;

                    /* Determine eastmost & westmost lon */
                    /* --------------------------------- */
                    if (f32_buf[i - spix + 1] > eastern_lon)
                        eastern_lon = f32_buf[i - spix + 1];
                    if (f32_buf[i - spix + 1] < western_lon)
                        western_lon = f32_buf[i - spix + 1];
                }
                memcpy(&data[j * 4 * dims_w[jds][1]], f32_buf, 4 * dims_w[jds][1]);
            }
        } else {
            status = SDreaddata(sds_id, start, NULL, dims_w[jds], (VOIDP)data);
            if (status == -1) {
                printf("read status: %d for %s\n\n", status, fieldname);
                exit(EXIT_FAILURE);
            }
        }

        /* Modify time/date fields */
        /* ----------------------- */
        if (strcmp(fieldname, "year") == 0) {
            memcpy(&startyear, &data[0], 4);
            memcpy(&cenyear, &data[4 * (dims_w[jds][0] / 2)], 4);
            memcpy(&endyear, &data[4 * (dims_w[jds][0] - 1)], 4);
        }

        if (strcmp(fieldname, "day") == 0) {
            memcpy(&startday, &data[0], 4);
            memcpy(&cenday, &data[4 * (dims_w[jds][0] / 2)], 4);
            memcpy(&endday, &data[4 * (dims_w[jds][0] - 1)], 4);
        }

        if (strcmp(fieldname, "msec") == 0) {
            memcpy(&startmsec, &data[0], 4);
            memcpy(&cenmsec, &data[4 * (dims_w[jds][0] / 2)], 4);
            memcpy(&endmsec, &data[4 * (dims_w[jds][0] - 1)], 4);
        }

        if (strcmp(fieldname, "longitude") == 0) {
            memcpy(&ullon, &data[0], 4);
            memcpy(&lllon, &data[4 * (dims_w[jds][1] * (dims_w[jds][0] - 1))], 4);

            memcpy(&urlon, &data[4 * (dims_w[jds][1] - 1)], 4);
            memcpy(&lrlon, &data[4 * (dims_w[jds][0] * dims_w[jds][1] - 1)], 4);

            memcpy(&cclon, &data[4 * (dims_w[jds][1] * dims_w[jds][0] / 2 - dims_w[jds][1] / 2)], 4);
            memcpy(&sclon, &data[4 * (dims_w[jds][1] / 2)], 4);
            memcpy(&eclon, &data[4 * (dims_w[jds][1] * dims_w[jds][0] - dims_w[jds][1] / 2)], 4);
        }

        if (strcmp(fieldname, "latitude") == 0) {
            memcpy(&ullat, &data[0], 4);
            memcpy(&lllat, &data[4 * (dims_w[jds][1] * (dims_w[jds][0] - 1))], 4);

            memcpy(&urlat, &data[4 * (dims_w[jds][1] - 1)], 4);
            memcpy(&lrlat, &data[4 * (dims_w[jds][0] * dims_w[jds][1] - 1)], 4);

            memcpy(&cclat, &data[4 * (dims_w[jds][1] * dims_w[jds][0] / 2 - dims_w[jds][1] / 2)], 4);
            memcpy(&sclat, &data[4 * (dims_w[jds][1] / 2)], 4);
            memcpy(&eclat, &data[4 * (dims_w[jds][1] * dims_w[jds][0] - dims_w[jds][1] / 2)], 4);
        }

        if (octs == 1) {
            if (strcmp(fieldname, "lon") == 0) {
                memcpy(&ullon, &data[0], 4);
                memcpy(&lllon, &data[4 * (dims_w[jds][0] - 1) * dims_w[jds][1]], 4);

                memcpy(&urlon, &data[4 * (dims_w[jds][1] - 1)], 4);
                memcpy(&lrlon, &data[4 * ((dims_w[jds][0] - 1) * dims_w[jds][1] + dims_w[jds][1] - 1)], 4);

                memcpy(&cclon, &data[4 * ((dims_w[jds][0] / 2) * dims_w[jds][1] + dims_w[jds][1] / 2)], 4);
                memcpy(&sclon, &data[4 * (dims_w[jds][1] / 2)], 4);
                memcpy(&eclon, &data[4 * ((dims_w[jds][0] - 1) * dims_w[jds][1] + dims_w[jds][1] / 2)], 4);
            }

            if (strcmp(fieldname, "lat") == 0) {
                memcpy(&ullat, &data[0], 4);
                memcpy(&lllat, &data[4 * (dims_w[jds][0] - 1) * dims_w[jds][1]], 4);

                memcpy(&urlat, &data[4 * (dims_w[jds][1] - 1)], 4);
                memcpy(&lrlat, &data[4 * ((dims_w[jds][0] - 1) * dims_w[jds][1] + dims_w[jds][1] - 1)], 4);

                memcpy(&cclat, &data[4 * ((dims_w[jds][0] / 2) * dims_w[jds][1] + dims_w[jds][1] / 2)], 4);
                memcpy(&sclat, &data[4 * (dims_w[jds][1] / 2)], 4);
                memcpy(&eclat, &data[4 * ((dims_w[jds][0] - 1) * dims_w[jds][1] + dims_w[jds][1] / 2)], 4);
            }
        }

        if (strcmp(fieldname, "csol_z") == 0) {
            memcpy(&csolz, &data[4 * (dims_w[jds][0] / 2)], 4);
        }

        if (strcmp(fieldname, "ntilts") == 0) {
            memcpy(&ntilts, data, 4);
        }

        if (strcmp(fieldname, "tilt_flags") == 0) {
            memcpy(old_tilt_flags, data, 2 * 20);
        }

        if (strcmp(fieldname, "tilt_ranges") == 0) {
            memset(new_tilt_flags, 0, 2 * 20);
            memset(new_tilt_ranges, 0, 2 * 2 * 20);
            memcpy(old_tilt_ranges, data, 2 * 2 * 20);

            for (i = 0; i < ntilts; i++) {
                if (sscan >= old_tilt_ranges[i][0] && sscan <= old_tilt_ranges[i][1]) {
                    stilt = i;
                }

                if (escan >= old_tilt_ranges[i][0] && escan <= old_tilt_ranges[i][1]) {
                    etilt = i;
                }
            }
            ntilts = etilt - stilt + 1;

            for (i = 0; i < ntilts; i++) {
                new_tilt_ranges[i][0] = old_tilt_ranges[stilt + i][0] - sscan + 1;
                new_tilt_ranges[i][1] = old_tilt_ranges[stilt + i][1] - sscan + 1;

                if (new_tilt_ranges[i][0] < 1)
                    new_tilt_ranges[i][0] = 1;
                if (new_tilt_ranges[i][1] < 1)
                    new_tilt_ranges[i][0] = 1;

                if (new_tilt_ranges[i][0] > nscan)
                    new_tilt_ranges[i][0] = nscan;
                if (new_tilt_ranges[i][1] > nscan)
                    new_tilt_ranges[i][1] = nscan;

                new_tilt_flags[i] = old_tilt_flags[stilt + i];
                new_tilt_flags[i] = old_tilt_flags[stilt + i];
            }

            memcpy(data, new_tilt_ranges, 2 * 2 * 20);
        }

        /* Write extract field */
        /* ------------------- */
        status = SDwritedata(sds_id_w, start_w, NULL, dims_w[jds], (VOIDP)data);
        if (status == -1) {
            printf("write status: %d\n\n", status);
            exit(EXIT_FAILURE);
        }

        /* Write slon,elon,clon */
        /* -------------------- */
        if (strcmp(fieldname, "longitude") == 0) {
            data2 = (char *)calloc(dims_w[jds][0], DFKNTsize(dtype));

            /* slon */
            for (j = 0; j < dims_w[jds][0]; j++) {
                k = DFKNTsize(dtype);
                memcpy(&data2[k * j], &data[k * j * dims_w[jds][1]], DFKNTsize(dtype));
            }
            sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "slon"));
            status = SDwritedata(sds_id_w, start_w, NULL, dims_w[jds], (VOIDP)data2);
            if (status == -1) {
                printf("write status: %d\n\n", status);
                exit(EXIT_FAILURE);
            }

            /* elon */
            for (j = 0; j < dims_w[jds][0]; j++) {
                k = DFKNTsize(dtype);
                memcpy(&data2[k * j], &data[k * (j * dims_w[jds][1] + dims_w[jds][1] - 1)], DFKNTsize(dtype));
            }
            sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "elon"));
            status = SDwritedata(sds_id_w, start_w, NULL, dims_w[jds], (VOIDP)data2);
            if (status == -1) {
                printf("write status: %d\n\n", status);
                exit(EXIT_FAILURE);
            }

            /* clon */
            for (j = 0; j < dims_w[jds][0]; j++) {
                k = DFKNTsize(dtype);
                memcpy(&data2[k * j], &data[k * (j * dims_w[jds][1] + dims_w[jds][1] / 2)], DFKNTsize(dtype));
            }
            sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "clon"));
            status = SDwritedata(sds_id_w, start_w, NULL, dims_w[jds], (VOIDP)data2);
            if (status == -1) {
                printf("write status: %d\n\n", status);
                exit(EXIT_FAILURE);
            }
            free(data2);
        }

        /* Write slat,elat,clat */
        /* -------------------- */
        if (strcmp(fieldname, "latitude") == 0) {
            data2 = (char *)calloc(dims_w[jds][0], DFKNTsize(dtype));

            /* slat */
            for (j = 0; j < dims_w[jds][0]; j++) {
                k = DFKNTsize(dtype);
                memcpy(&data2[k * j], &data[k * j * dims_w[jds][1]], DFKNTsize(dtype));
            }
            sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "slat"));
            status = SDwritedata(sds_id_w, start_w, NULL, dims_w[jds], (VOIDP)data2);
            if (status == -1) {
                printf("write status: %d\n\n", status);
                exit(EXIT_FAILURE);
            }

            /* elat */
            for (j = 0; j < dims_w[jds][0]; j++) {
                k = DFKNTsize(dtype);
                memcpy(&data2[k * j], &data[k * (j * dims_w[jds][1] + dims_w[jds][1] - 1)], DFKNTsize(dtype));
            }
            sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "elat"));
            status = SDwritedata(sds_id_w, start_w, NULL, dims_w[jds], (VOIDP)data2);
            if (status == -1) {
                printf("write status: %d\n\n", status);
                exit(EXIT_FAILURE);
            }

            /* clat */
            for (j = 0; j < dims_w[jds][0]; j++) {
                k = DFKNTsize(dtype);
                memcpy(&data2[k * j], &data[k * (j * dims_w[jds][1] + dims_w[jds][1] / 2)], DFKNTsize(dtype));
            }
            sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "clat"));
            status = SDwritedata(sds_id_w, start_w, NULL, dims_w[jds], (VOIDP)data2);
            if (status == -1) {
                printf("write status: %d\n\n", status);
                exit(EXIT_FAILURE);
            }
            free(data2);
        }

        free(data);

        SDendaccess(sds_id);
        SDendaccess(sds_id_w);

    } /* Dataset Loop */

    /* Correct "tilt" fields */
    /* --------------------- */
    sds_id = SDselect(sd_id_w, SDnametoindex(sd_id_w, "ntilts"));
    if (sds_id != -1) {
        status = SDgetinfo(sds_id, fieldname, &rank, dims, &dtype, &nattrs);
        status = SDwritedata(sds_id, &zero, NULL, dims, (VOIDP)&ntilts);
        if (status == -1) {
            printf("write ntilts status: %d\n\n", status);
            exit(EXIT_FAILURE);
        }
    }

    sds_id = SDselect(sd_id_w, SDnametoindex(sd_id_w, "tilt_flags"));
    if (sds_id != -1) {
        status = SDgetinfo(sds_id, fieldname, &rank, dims, &dtype, &nattrs);
        status = SDwritedata(sds_id, &zero, NULL, dims, (VOIDP)new_tilt_flags);
        if (status == -1) {
            printf("write tilt_flags status: %d for %d\n\n", status, jds);
            exit(EXIT_FAILURE);
        }
    }

    /* If Seawifs L1A determine geo limits */
    /* ----------------------------------- */
    if (seawifs_l1a) {
        start_w[0] = 0;
        dims_w[0][0] = epix - spix + 1;
        data2 = (char *)calloc(dims_w[0][0], sizeof(float));

        sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "slon"));
        status = SDreaddata(sds_id_w, start_w, NULL, dims_w[0], (VOIDP)data2);

        memcpy(&ullon, &data2[0], sizeof(float));
        memcpy(&lllon, &data2[sizeof(float) * (epix - spix)], sizeof(float));

        sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "clon"));
        status = SDreaddata(sds_id_w, start_w, NULL, dims_w[0], (VOIDP)data2);

        memcpy(&sclon, &data2[0], sizeof(float));
        memcpy(&eclon, &data2[sizeof(float) * (epix - spix)], sizeof(float));
        memcpy(&cclon, &data2[sizeof(float) * (epix - spix) / 2], sizeof(float));

        sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "elon"));
        status = SDreaddata(sds_id_w, start_w, NULL, dims_w[0], (VOIDP)data2);

        memcpy(&urlon, &data2[0], sizeof(float));
        memcpy(&lrlon, &data2[sizeof(float) * (epix - spix)], sizeof(float));

        eastern_lon = -180.0;
        western_lon = +180.0;
        if (ullon > eastern_lon)
            eastern_lon = ullon;
        if (lllon > eastern_lon)
            eastern_lon = lllon;
        if (sclon > eastern_lon)
            eastern_lon = sclon;
        if (eclon > eastern_lon)
            eastern_lon = eclon;
        if (urlon > eastern_lon)
            eastern_lon = urlon;
        if (lrlon > eastern_lon)
            eastern_lon = lrlon;

        if (ullon < western_lon)
            western_lon = ullon;
        if (lllon < western_lon)
            western_lon = lllon;
        if (sclon < western_lon)
            western_lon = sclon;
        if (eclon < western_lon)
            western_lon = eclon;
        if (urlon < western_lon)
            western_lon = urlon;
        if (lrlon < western_lon)
            western_lon = lrlon;

        sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "slat"));
        status = SDreaddata(sds_id_w, start_w, NULL, dims_w[0], (VOIDP)data2);

        memcpy(&ullat, &data2[0], sizeof(float));
        memcpy(&lllat, &data2[sizeof(float) * (epix - spix)], sizeof(float));

        sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "clat"));
        status = SDreaddata(sds_id_w, start_w, NULL, dims_w[0], (VOIDP)data2);

        memcpy(&sclat, &data2[0], sizeof(float));
        memcpy(&eclat, &data2[sizeof(float) * (epix - spix)], sizeof(float));
        memcpy(&cclat, &data2[sizeof(float) * (epix - spix) / 2], sizeof(float));

        sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, "elat"));
        status = SDreaddata(sds_id_w, start_w, NULL, dims_w[0], (VOIDP)data2);

        memcpy(&urlat, &data2[0], sizeof(float));
        memcpy(&lrlat, &data2[sizeof(float) * (epix - spix)], sizeof(float));

        northern_lat = -90.0;
        southern_lat = +90.0;
        if (ullat > northern_lat)
            northern_lat = ullat;
        if (lllat > northern_lat)
            northern_lat = lllat;
        if (sclat > northern_lat)
            northern_lat = sclat;
        if (eclat > northern_lat)
            northern_lat = eclat;
        if (urlat > northern_lat)
            northern_lat = urlat;
        if (lrlat > northern_lat)
            northern_lat = lrlat;

        if (ullat < southern_lat)
            southern_lat = ullat;
        if (lllat < southern_lat)
            southern_lat = lllat;
        if (sclat < southern_lat)
            southern_lat = sclat;
        if (eclat < southern_lat)
            southern_lat = eclat;
        if (urlat < southern_lat)
            southern_lat = urlat;
        if (lrlat < southern_lat)
            southern_lat = lrlat;

        /*
        if (strcmp(buffer, "Scene Center Solar Zenith") == 0) {
          memcpy(data, &csolz, 4);
        }

        if (strcmp(buffer, "Easternmost Longitude") == 0 && lonlatinterp == 1) {
          memcpy(data, &eastern_lon, 4);
        }
        if (strcmp(buffer, "Westernmost Longitude") == 0 && lonlatinterp == 1) {
          memcpy(data, &western_lon, 4);
        }
         */
    } else {
        if (datelinecross == 1) {
            eastern_lon = lllon;
            western_lon = lrlon;
        }
    }

    /* Write Global Attributes */
    /* ----------------------- */
    for (i = 0; i < nglobal_attr; i++) {
        status = SDattrinfo(sd_id_r, i, buffer, &dtype, &count);

        data = (char *)calloc(count, DFKNTsize(dtype));

        status = SDreadattr(sd_id_r, i, (VOIDP)data);

        if (strcmp(buffer, "Number of Scan Lines") == 0) {
            if (octs == 1)
                i32 = (escan - sscan + 1) / 2;
            else
                i32 = escan - sscan + 1;
            memcpy(data, &i32, 4);
        }

        if (strcmp(buffer, "Pixels per Scan Line") == 0) {
            i32 = epix - spix + 1;
            memcpy(data, &i32, 4);
        }

        if (strcmp(buffer, "Number of Scan Control Points") == 0 && lonlatinterp == 1) {
            i32 = escan - sscan + 1;
            memcpy(data, &i32, 4);
        }

        if (strcmp(buffer, "Number of Pixel Control Points") == 0 && lonlatinterp == 1) {
            i32 = epix - spix + 1;
            memcpy(data, &i32, 4);
        }

        if (strcmp(buffer, "Scene Center Scan Line") == 0) {
            i32 = (escan - sscan + 1) / 2 + 1; /*((escan - sscan + 1) % 2);*/
            memcpy(data, &i32, 4);
        }

        if (strcmp(buffer, "Start Time") == 0) {
            if (startyear != 0) {
                utime = yds2unix(startyear, startday, startmsec / ((double)1000.0));
                memcpy(data, ydhmsf(utime, 'G'), 17);
            }
        }
        if (strcmp(buffer, "Start Year") == 0) {
            if (startyear != 0)
                memcpy(data, &startyear, 4);
        }
        if (strcmp(buffer, "Start Day") == 0) {
            if (startday != 0)
                memcpy(data, &startday, 4);
        }
        if (strcmp(buffer, "Start Millisec") == 0) {
            if (startmsec != 0)
                memcpy(data, &startmsec, 4);
        }

        if (strcmp(buffer, "End Time") == 0) {
            if (endyear != 0) {
                utime = yds2unix(endyear, endday, endmsec / ((double)1000.0));
                memcpy(data, ydhmsf(utime, 'G'), 17);
            }
        }
        if (strcmp(buffer, "End Year") == 0) {
            if (endyear != 0)
                memcpy(data, &endyear, 4);
        }
        if (strcmp(buffer, "End Day") == 0) {
            if (endday != 0)
                memcpy(data, &endday, 4);
        }
        if (strcmp(buffer, "End Millisec") == 0) {
            if (endmsec != 0)
                memcpy(data, &endmsec, 4);
        }

        if (strcmp(buffer, "Scene Center Time") == 0) {
            if (cenyear != 0) {
                utime = yds2unix(cenyear, cenday, cenmsec / ((double)1000.0));
                memcpy(data, ydhmsf(utime, 'G'), 17);
            }
        }

        if (strcmp(buffer, "Upper Left Longitude") == 0) {
            memcpy(data, &ullon, 4);
        }
        if (strcmp(buffer, "Upper Right Longitude") == 0) {
            memcpy(data, &urlon, 4);
        }
        if (strcmp(buffer, "Upper Left Latitude") == 0) {
            memcpy(data, &ullat, 4);
        }
        if (strcmp(buffer, "Upper Right Latitude") == 0) {
            memcpy(data, &urlat, 4);
        }
        if (strcmp(buffer, "Start Center Longitude") == 0) {
            memcpy(data, &sclon, 4);
        }
        if (strcmp(buffer, "Start Center Latitude") == 0) {
            memcpy(data, &sclat, 4);
        }
        if (strcmp(buffer, "Scene Center Longitude") == 0) {
            memcpy(data, &cclon, 4);
        }
        if (strcmp(buffer, "Scene Center Latitude") == 0) {
            memcpy(data, &cclat, 4);
        }
        if (strcmp(buffer, "Scene Center Solar Zenith") == 0) {
            memcpy(data, &csolz, 4);
        }
        if (strcmp(buffer, "Lower Left Longitude") == 0) {
            memcpy(data, &lllon, 4);
        }
        if (strcmp(buffer, "Lower Right Longitude") == 0) {
            memcpy(data, &lrlon, 4);
        }
        if (strcmp(buffer, "Lower Left Latitude") == 0) {
            memcpy(data, &lllat, 4);
        }
        if (strcmp(buffer, "Lower Right Latitude") == 0) {
            memcpy(data, &lrlat, 4);
        }
        if (strcmp(buffer, "End Center Longitude") == 0) {
            memcpy(data, &eclon, 4);
        }
        if (strcmp(buffer, "End Center Latitude") == 0) {
            memcpy(data, &eclat, 4);
        }

        if (strcmp(buffer, "Northernmost Latitude") == 0 && lonlatinterp == 1) {
            memcpy(data, &northern_lat, 4);
        }
        if (strcmp(buffer, "Southernmost Latitude") == 0 && lonlatinterp == 1) {
            memcpy(data, &southern_lat, 4);
        }
        if (strcmp(buffer, "Easternmost Longitude") == 0 && lonlatinterp == 1) {
            memcpy(data, &eastern_lon, 4);
        }
        if (strcmp(buffer, "Westernmost Longitude") == 0 && lonlatinterp == 1) {
            memcpy(data, &western_lon, 4);
        }

        if (strcmp(buffer, "LAC Pixel Start Number") == 0) {
            status = SDreadattr(sd_id_r, i + 1, (VOIDP)&i32);
            memcpy(&LAC_start, data, 4);
            LAC_start += i32 * (spix - 1);
            memcpy(data, &LAC_start, 4);
        }

        status = SDsetattr(sd_id_w, buffer, dtype, count, (VOIDP)data);

        /* Set Creation Time to current time
           ("Processing Time" is copied from input file */
        if (strcmp(buffer, "Processing Time") == 0) {
            get_time(data);
            status = SDsetattr(sd_id_w, "Creation Time", dtype, count, (VOIDP)data);
        }

        /*    printf("%d %d\n", i, status);*/
        free(data);
    }

    /* Write extract attributes if Level 1 */
    /* ----------------------------------- */
    for (i = 0; i < nglobal_attr; i++) {
        status = SDattrinfo(sd_id_r, i, buffer, &dtype, &count);

        data = (char *)calloc(count, DFKNTsize(dtype));

        status = SDreadattr(sd_id_r, i, (VOIDP)data);

        if (strcmp(buffer, "Title") == 0) {
            if (strstr(data, "Level-1A") != NULL) {
                /* Set Extract Pixel Offset Attribute */
                /* ---------------------------------- */
                i32 = spix - 1;
                status = SDsetattr(sd_id_w, "Extract Pixel Offset", DFNT_INT32, 1, (VOIDP)&i32);

                /* Set Extract Pixel Count Attribute */
                /* --------------------------------- */
                i32 = epix - spix + 1;
                status = SDsetattr(sd_id_w, "Extract Pixel Count", DFNT_INT32, 1, (VOIDP)&i32);

                /* Set Extract Line Offset Attribute */
                /* --------------------------------- */
                i32 = sscan - 1;
                status = SDsetattr(sd_id_w, "Extract Line Offset", DFNT_INT32, 1, (VOIDP)&i32);

                /* Set Extract Line Count Attribute */
                /* -------------------------------- */
                i32 = escan - sscan + 1;
                status = SDsetattr(sd_id_w, "Extract Line Count", DFNT_INT32, 1, (VOIDP)&i32);
            }
            free(data);
            break;
        }
    }

    /* Build and fill Vgroups */
    /* ---------------------- */
    for (j = 0; j < nvgrp; j++) {
        vg_ref = Vfind(HDFfid_r, vgrpname[j]);
        vgid_r = Vattach(HDFfid_r, vg_ref, "r");
        vgid_w = Vattach(HDFfid_w, -1, "w");
        Vgetname(vgid_r, buffer);
        Vsetname(vgid_w, buffer);
        Vgetclass(vgid_r, buffer);
        Vsetclass(vgid_w, buffer);

        for (i = 0; i < Vntagrefs(vgid_r); i++) {
            status = Vgettagref(vgid_r, i, &tag, &ref);
            sds_id = SDselect(sd_id_r, SDreftoindex(sd_id_r, ref));
            status = SDgetinfo(sds_id, fieldname, &rank, dims, &dtype, &nattrs);

            /* Check if L2 product field & in prodlist */
            /* --------------------------------------- */
            if (prodlist[0] != 0 && strcmp(vgrpname[j], "Geophysical Data") == 0) {
                strcpy(buf2, ",");
                strcat(buf2, fieldname);
                strcat(buf2, ",");

                if (strstr(prodlist, buf2) == NULL)
                    goto LBL3;
            }

            sds_id_w = SDselect(sd_id_w, SDnametoindex(sd_id_w, fieldname));
            Vaddtagref(vgid_w, DFTAG_NDG, SDidtoref(sds_id_w));

            SDendaccess(sds_id_w);
        LBL3:
            SDendaccess(sds_id);
        }
        Vdetach(vgid_r);
        Vdetach(vgid_w);
    }

    SDend(sd_id_r);
    Hclose(HDFfid_r);
    Vend(HDFfid_r);

    SDend(sd_id_w);
    Hclose(HDFfid_w);
    Vend(HDFfid_w);

    return 0;
}
