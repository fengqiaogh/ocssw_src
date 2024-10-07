/* ====================================================== */
/* Module l1_octs_hdf.c                                   */
/*                                                        */
/* Functions to open and read a OCTS HDF l1b file.        */
/*                                                        */
/* Written By:                                            */
/*    JT                                                  */
/*    NASA/SIMBIOS Project                                */
/*    1/99                                                */
/* Modifications By:                                      */
/*    J. Gales                                            */
/*    Futuretech                                          */
/*    NASA/SIMBIOS Project                                */
/*    6/99                                                */
/*                                                        */
/*    10/00                                               */
/*                                                        */
/*    Add support for OCTS L1A                            */
/*                                                        */
/*    Sean Baiey                                          */
/*    Futuretech                                          */
/*    6/08                                                */
/*    Numerous minor mods to get extract processing       */
/*    to work properly                                    */
/* ====================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <hdf4utils.h>
#include "l1.h"
#include "filehdr_struc.h"
#include "l1_octs.h"

#include <hdf.h>
#include <mfhdf.h>

#define RADIUS 6378.137    /* Earth radius in km */
#define FL 1.0/298.257     /* Earth flatening factor */
#define NR 560             /* # rows */      
#define NC 30              /* max # columns */     

#define MAXOCLIN 6700      /* max # lines */
#define MAXOCPIX 2218      /* max # pixels */
#define MAXOCARR 10000
#define NOCBANDS 8

void reform_octs_time(char *time);
int navigation(int32_t fileID);
int CalcViewAngle(float lon1, float lat1, float pos[3], float usun[]);
int LeapCheck(int yr);

static int32_t msec_start;
static int32_t msec[MAXOCLIN];
static float tilt[MAXOCLIN];
static int16_t year, day, nline, npix, sline;
static int32_t spix;
static int16_t maxday;
static float solz1, sola1, senz1, sena1;
static float *lon, *lat, *senz, *sena, *solz, *sola;
static float m[NOCBANDS], b[NOCBANDS];
static int ncol;
static int32_t scansPerScene;
static int32_t linesPerScan;
static int16_t *tilt_flag;
static int16_t *gain;
static float *inst_temp;
static int16_t *miss_qual;
static int16_t num_tab;
static int32_t extr_line_offset = 0;
static int16_t *start_scan;
static int16_t last_table = 0;
static int32_t samp_start[8];
static int32_t samp_edges[8];
static int16_t sample_table[3][8][2][400][2];


/* =========================================================================  */
/*                                                                            */
/*  Module get_octs_cal.c                                                     */
/*                                                                            */
/*  Functions to open the given calibration HDF file and apply apply to       */
/*  L1A counts, returns TOA radiances.                                        */
/*                                                                            */
/*  Modification history:                                                     */
/*                                                                            */
/*  Sean Bailey	Futuretech Corporation	16 Dec 2005	Original development  */
/*									      */
/* =========================================================================  */

#define NTILT   3
#define NDET    10
#define NBND    12
#define NSEG 2
#define NGAIN   4

int32_t get_octs_cal(char *file, int16_t year, int16_t day, int32_t msec[MAXOCLIN],
        int32_t recnum, int16_t npix, int32_t spix, int32_t tilt,
        int16_t gainset[MAXOCLIN], float inst_temp[MAXOCLIN],
        int16_t sample_table[3][8][2][400][2], int32_t scansPerScene,
        uint16_t l1acnts[NOCBANDS][MAXOCPIX], float l1brads[MAXOCPIX][NOCBANDS]) {

    static float basic_offset [NBND][NDET][NGAIN];
    static float calib_coeff_gain [NBND][NDET][NGAIN];
    static float det_norm_slope [NOCBANDS][NDET];
    static float det_norm_offset [NOCBANDS][NDET];
    static float L_Border [NOCBANDS][NDET];
    static float tilt_refl [NBND][NTILT];
    static float time_dep [NBND][NSEG][3];

    static float unit_conv = 1.0e-1;
    static int firstCall = 1;
    static double refjsec;
    static int32_t refyear = 1996;
    static int refday = 306;

    int16_t subpix;
    int i, j, idet, iscn, ign, itilt, mside, igscnmod, rdet;
    float eta, Lt, A, tcorr;
    double jsec = yds2unix(year, day, ((double) (msec[recnum])) / 1000.0);
    int segment = 0;
    if (jsec >= 854769600.)
        segment = 1;

    if (firstCall) {

        char name [H4_MAX_NC_NAME] = "";
        char sdsname[H4_MAX_NC_NAME] = "";
        int32_t sd_id;
        int32_t sds_id;
        int32_t rank;
        int32_t nt;
        int32_t dims[H4_MAX_VAR_DIMS];
        int32_t nattrs;
        int32_t start[3] = {0, 0, 0};
        int32_t status;

        /* get the current calibration */
        if (file == NULL) {
            fprintf(stderr,
                    "-E %s Line %d: No calibration file specified.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        if (want_verbose)
            printf("Loading caltable: %s\n", file);

        /* Open the file */
        sd_id = SDstart(file, DFACC_RDONLY);
        if (sd_id == FAIL) {
            fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                    __FILE__, __LINE__, file, DFACC_RDONLY);
            exit(1);
        }

        strcpy(sdsname, "basic_offset");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) basic_offset);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "calib_coeff_gain");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) calib_coeff_gain);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "det_norm_slope");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) det_norm_slope);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "det_norm_offset");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) det_norm_offset);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "L_Border");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) L_Border);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "tilt_refl");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) tilt_refl);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "time_dep");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);

        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) time_dep);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        /* terminate access to the SD interface and close the file */
        status = SDend(sd_id);
        refjsec = yds2unix(refyear, refday, (double) 0.0);

        firstCall = 0;
    }

    /* Apply the appropriate coefficients */

    /* Tilt index for reflectance correction */
    itilt = 1;
    if (tilt < -10) itilt = 0;
    if (tilt > 10) itilt = 2;

    for (j = 0; j < npix; j++) {
        subpix = spix + j;
        for (i = 0; i < NOCBANDS; i++) {

            iscn = recnum / 2;
            mside = recnum % 2;
            //rdet = sample_table[0][itilt][i][mside][j][1];
            rdet = sample_table[itilt][i][mside][subpix][1];
            idet = (rdet + 200) % 10;
            igscnmod = 1;

            /* account for detectors from previous or following scan (or two) */
            if (rdet < -9) igscnmod -= 2;
            if (rdet >= -9 && rdet < 0) igscnmod -= 1;
            if (rdet >= 10 && rdet < 20) igscnmod += 1;
            if (rdet >= 20) igscnmod += 2;

            if (idet < 0.) idet += 9;
            if (idet < 0 || idet > 9) {
                printf("-W- idet %d outside bounds [0,9]\n", idet);
            }

            if (iscn < 0 || iscn > scansPerScene - 1) {
                printf("-W- iscn %d outside bounds [0,%d]\n", iscn, scansPerScene);
                l1brads[j][i] = 0.0;
            } else if (l1acnts[i][j] == 0) {
                l1brads[j][i] = 0.0;
            } else if (inst_temp[iscn] > 30) {
                printf("inst_temp: %8.5f\n", inst_temp[iscn]);
                l1brads[j][i] = 0.0;
            } else if (iscn + igscnmod < 0 || iscn + igscnmod > scansPerScene - 1) {
                // Some of the detectors are out of range
                l1brads[j][i] = 0.0;
            } else {
                ign = gainset[iscn + igscnmod];
                eta = tilt_refl[i][itilt];

                Lt = unit_conv * (l1acnts[i][j] - basic_offset[i][idet][ign]) /
                        calib_coeff_gain[i][4][ign];

                if (Lt <= L_Border[i][idet]) {
                    A = det_norm_slope[i][idet] * Lt + det_norm_offset[i][idet];
                } else {
                    A = det_norm_slope[i][idet] * L_Border[i][idet] +
                            det_norm_offset[i][idet];
                }

                /* Calculate time dependent correction... */
                tcorr = time_dep[i][segment][1]*(jsec - refjsec) + time_dep[i][segment][0];
                l1brads[j][i] = (A * Lt / eta) * tcorr;
            }

        } /* end i */

    } /* end j */

    return 0;
}


/* ------------------------------------------------------ */
/* openl1_octs() - opens a OCTS HDF level-1 file */
/*                    for reading, if not already opened. */
/*                    Reads the global attributes.        */
/*                    Calculates the subscene info, if    */
/*                    necessary, ane gets the navigation  */
/*                    info.                               */
/*                                                        */
/* CALLS: getHDFattr()                                    */
/*        navigation()                                    */
/*        rdSDS()                                         */
/*                                                        */

/* ------------------------------------------------------ */
int openl1_octs(filehandle *l1file) {
    int i, j;
    int32_t status;
    int32_t pixPerScan;
    int32_t linesPerScene;
    int32_t msec_temp1[MAXOCARR], msec_temp2[MAXOCARR];
    char buf[32];
    int32_t fileID;
    int32_t sds_index, sds_id;
    int32_t dims[H4_MAX_VAR_DIMS];

    /* open the file and get the file ID */
    fileID = SDstart(l1file->name, DFACC_RDONLY);

    /* get attributes pertaining to the full scene */
    status = getHDFattr(fileID, "Start Year", "", (VOIDP) & year);
    status = getHDFattr(fileID, "Start Day", "", (VOIDP) & day);
    status = getHDFattr(fileID, "Number of Scan Lines", "", (VOIDP) & scansPerScene);
    status = getHDFattr(fileID, "Pixels per Scan Line", "", (VOIDP) & pixPerScan);
    status = getHDFattr(fileID, "Lines per Scan", "", (VOIDP) & linesPerScan);
    status = getHDFattr(fileID, "Orbit Node Longitude", "", (VOIDP) & l1file->orbit_node_lon);
    status = getHDFattr(fileID, "Orbit Number", "", (VOIDP) & l1file->orbit_number);
    status = getHDFattr(fileID, "Node Crossing Time", "", (VOIDP) & buf);
    reform_octs_time(buf);
    l1file->node_crossing_time = zulu2unix(buf);


    if (getDims(fileID, "pxl", dims) != 0) {
        printf("-E- %s: Error reading control-point column dimension.\n",
                "__FILE__");
        exit(FATAL_ERROR);
    }
    ncol = dims[0];


    linesPerScene = scansPerScene * linesPerScan;

    /* Reformat OCTS Times */
    //   reform_octs_time(l1file->node_crossing_time);

    /* Adjust spix if this is an extract file from l1extract */
    /* 19 May 2008, BAF                                      */
    status = getHDFattr(fileID, "Extract Pixel Offset", "", (VOIDP) & spix);
    if (status != -1) {
        status = getHDFattr(fileID, "Extract Line Offset", "", (VOIDP) & extr_line_offset);
        printf("File is Level-1A extract starting on line %d.\n", extr_line_offset + 1);
    } else {
        extr_line_offset = 0;
        spix = 0;
    }
    sline = 0;
    nline = linesPerScene;
    npix = pixPerScan;


    /* Check that number of scene lines not greater than allocation */
    if (nline > MAXOCLIN) {
        printf("-E- %s: Number of scene lines: %d, greater than array allocation: %d.\n",
                "openl1_octs", nline, MAXOCLIN);
        status = SDend(fileID);
        exit(FATAL_ERROR);
    }


    /* Check that number of pixels per scan not greater than allocation */
    if (npix > MAXOCPIX) {
        printf("-E- %s: Number of scan pixels: %d, greater than array allocation: %d.\n",
                "openl1_octs", npix, MAXOCPIX);
        status = SDend(fileID);
        exit(FATAL_ERROR);
    }


    lon = (float *) calloc(npix*nline, sizeof (float));
    lat = (float *) calloc(npix*nline, sizeof (float));
    senz = (float *) calloc(npix*nline, sizeof (float));
    sena = (float *) calloc(npix*nline, sizeof (float));
    solz = (float *) calloc(npix*nline, sizeof (float));
    sola = (float *) calloc(npix*nline, sizeof (float));

    l1file->npix = npix;
    l1file->nscan = nline;
    l1file->sd_id = fileID;
    strcpy(l1file->spatialResolution, "3.5 km");

    /* compute and store the view angles for the scene
    or subscene */
    status = navigation(fileID);

    /* get the time for each line of full scene */
    status = rdSDS(fileID, "msec", 0, 0, 0, 0, (VOIDP) & msec_temp1);

    for (i = 0; i < scansPerScene; i++) {
        for (j = 0; j < linesPerScan; j++) {
            msec_temp2[i * linesPerScan + j] = msec_temp1[i];
        }
    }
    msec_start = msec_temp2[0];
    j = 0;
    for (i = sline; i < (sline + nline); i++) {
        msec[j] = msec_temp2[i];
        j++;
    }

    /* get the tilt for each line of full scene */
    status = rdSDS(fileID, "tilt", 0, 0, 0, 0, (VOIDP) & tilt);

    maxday = 365 + abs(LeapCheck(year));

    /* get image data scaling parameters (y=mx+b) */
    if (l1file->format == FT_OCTSL1B) {
        status = getHDFattr(fileID, "slope", "l1b_b1_data", (VOIDP) & m[0]);
        status = getHDFattr(fileID, "slope", "l1b_b2_data", (VOIDP) & m[1]);
        status = getHDFattr(fileID, "slope", "l1b_b3_data", (VOIDP) & m[2]);
        status = getHDFattr(fileID, "slope", "l1b_b4_data", (VOIDP) & m[3]);
        status = getHDFattr(fileID, "slope", "l1b_b5_data", (VOIDP) & m[4]);
        status = getHDFattr(fileID, "slope", "l1b_b6_data", (VOIDP) & m[5]);
        status = getHDFattr(fileID, "slope", "l1b_b7_data", (VOIDP) & m[6]);
        status = getHDFattr(fileID, "slope", "l1b_b8_data", (VOIDP) & m[7]);
        status = getHDFattr(fileID, "intercept", "l1b_b1_data", (VOIDP) & b[0]);
        status = getHDFattr(fileID, "intercept", "l1b_b2_data", (VOIDP) & b[1]);
        status = getHDFattr(fileID, "intercept", "l1b_b3_data", (VOIDP) & b[2]);
        status = getHDFattr(fileID, "intercept", "l1b_b4_data", (VOIDP) & b[3]);
        status = getHDFattr(fileID, "intercept", "l1b_b5_data", (VOIDP) & b[4]);
        status = getHDFattr(fileID, "intercept", "l1b_b6_data", (VOIDP) & b[5]);
        status = getHDFattr(fileID, "intercept", "l1b_b7_data", (VOIDP) & b[6]);
        status = getHDFattr(fileID, "intercept", "l1b_b8_data", (VOIDP) & b[7]);
    }


    /* L1A GAC Processing */
    if (l1file->format == FT_OCTSL1A) {

        /* Get tilt flag, instrument temperature, and gain index */
        tilt_flag = (int16_t *) calloc(scansPerScene, sizeof (int16_t));
        status = rdSDS(fileID, "tilt_flag", 0, 0, 0, 0, (VOIDP) & tilt_flag[0]);
        inst_temp = (float *) calloc(4 * scansPerScene, sizeof (float));
        status = rdSDS(fileID, "inst_temp", 0, 0, 0, 0, (VOIDP) & inst_temp[0]);
        gain = (int16_t *) calloc(8 * scansPerScene, sizeof (int16_t));
        status = rdSDS(fileID, "gain", 0, 0, 0, 0, (VOIDP) & gain[0]);
        miss_qual = (int16_t *) calloc(scansPerScene, sizeof (int16_t));
        status = rdSDS(fileID, "miss_qual", 0, 0, 0, 0, (VOIDP) & miss_qual[0]);

        /* Read GAC subsampling table */
        status = rdSDS(fileID, "num_tables", 0, 0, 1, 1, (VOIDP) & num_tab);
        start_scan = (int16_t *) calloc(num_tab, sizeof (int16_t));
        status = rdSDS(fileID, "start_line", 0, 0, 0, 0, (VOIDP) & start_scan[0]);

        for (i = 0; i < num_tab; i++) {
            if ((extr_line_offset + sline) / linesPerScan < start_scan[i]) {
                last_table = i - 1;
                break;
            }
        }

        samp_start[0] = last_table;
        samp_edges[0] = 1;
        samp_start[1] = 0;
        samp_edges[1] = 3;
        samp_start[2] = 0;
        samp_edges[2] = 8;
        samp_start[3] = 0;
        samp_edges[3] = 2;
        samp_start[4] = 0;
        samp_edges[4] = 400;
        samp_start[5] = 0;
        samp_edges[5] = 2;

        sds_index = SDnametoindex(fileID, "samp_table");
        sds_id = SDselect(fileID, sds_index);

        status = SDreaddata(sds_id, samp_start, NULL, samp_edges, (VOIDP) sample_table);
    }

    return (0);
}

/* ------------------------------------------------------ */
/* readl1_octs() - reads a OCTS HDF level-1 record.   */
/*                                                        */

/* ------------------------------------------------------ */
int readl1_octs(filehandle *l1file, int32_t recnum, l1str *l1rec) {
    int i, j, status;
    uint16_t dataarr[NOCBANDS][MAXOCPIX];
    float scanarr[MAXOCPIX][NOCBANDS];
    int32_t fileID;
    int32_t start[8], edges[8];
    int32_t sds_index, sds_id;
    int32_t tilt_deg;
    int32_t scan;
    int32_t ip, ib, iw;
    int32_t nwave = l1rec->l1file->nbands;
    int32_t *bindx = l1rec->l1file->bindx;
    char *cal_path = l1_input->calfile;
    int32_t rank;
    int32_t nt;
    int32_t dims[H4_MAX_VAR_DIMS];
    int32_t nattrs;
    char name [H4_MAX_NC_NAME] = "";
    static int FirstCall = 1;
    int current_scan = recnum / linesPerScan;
    fileID = l1file->sd_id;

    /* load scan times */
    /* if the day changed between the start of the scene
    and the start of the subscene, adjust the day/year
    accordingly */
    static int dayIncrimented = 0;
    if (msec[recnum] < msec_start) {
        if(!dayIncrimented) {
            dayIncrimented = 1;
            day = day + 1;
            if (day > maxday) {
                year = year + 1;
                day = day - maxday;
            }
        }
    }
    if (msec[recnum] >= 86400000L) {
        msec[recnum] = msec[recnum] - 86400000L;
        if(!dayIncrimented) {
            dayIncrimented = 1;
            day = day + 1;
            if (day > maxday) {
                year = year + 1;
                day = day - maxday;
            }
        }
    }
    l1rec->scantime = yds2unix(year, day, (double) (msec[recnum] / 1.e3));

    /* load standard navigation */
    memcpy(l1rec->lon, &lon [recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->lat, &lat [recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->solz, &solz[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->sola, &sola[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->senz, &senz[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->sena, &sena[recnum * npix], sizeof (float)*npix);


    if (l1file->format == FT_OCTSL1A) {

        /* Read the l1a data */
        sds_index = SDnametoindex(fileID, "l1a_data");
        sds_id = SDselect(fileID, sds_index);

        start[1] = recnum;
        start[2] = 0;
        edges[0] = 1;
        edges[1] = 1;
        edges[2] = npix;
        for (i = 0; i < NOCBANDS; i++) {
            start[0] = i;
            status = SDreaddata(sds_id, start, NULL, edges, &dataarr[i][0]);
        }

        scan = (recnum + extr_line_offset) / linesPerScan;

        l1rec->tilt = tilt[current_scan];

        tilt_deg = 0;
        if (tilt_flag[current_scan] == 1) tilt_deg = +20;
        if (tilt_flag[current_scan] == 2) tilt_deg = 0;
        if (tilt_flag[current_scan] == 4) tilt_deg = -20;

        /* Read the next subsample table if necessary */
        if (FirstCall == 1 && recnum > 0) {
            for (i = 0; i < num_tab; i++) {
                if (scan / linesPerScan < start_scan[i]) {
                    last_table = i;
                    FirstCall = 0;
                    break;
                }
            }
        }


        if (scan >= start_scan[last_table + 1]) {
            last_table++;
            samp_start[0] = last_table;

            sds_index = SDnametoindex(fileID, "samp_table");
            sds_id = SDselect(fileID, sds_index);
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, samp_start, NULL, samp_edges, (VOIDP) sample_table);
        }

        /* Calibrate the l1a data */
        status = get_octs_cal(cal_path, year, day, msec, recnum, npix, spix, tilt_deg,
                gain, inst_temp, sample_table, scansPerScene,
                dataarr, scanarr);

        if (status < 0) {
            fprintf(stderr,
                    "-E- %s line %d: Error applying calibration table \"%s\".\n",
                    __FILE__, __LINE__, cal_path);
            exit(status);
        }

    }


    /* load the L1B image data arrays (Lt)*/
    if (l1file->format == FT_OCTSL1B) {

        l1rec->tilt = tilt[current_scan];

        status = rdSDS(fileID, "l1b_b1_data", recnum, 0, 1, npix, (VOIDP) & dataarr[0][0]);
        status = rdSDS(fileID, "l1b_b2_data", recnum, 0, 1, npix, (VOIDP) & dataarr[1][0]);
        status = rdSDS(fileID, "l1b_b3_data", recnum, 0, 1, npix, (VOIDP) & dataarr[2][0]);
        status = rdSDS(fileID, "l1b_b4_data", recnum, 0, 1, npix, (VOIDP) & dataarr[3][0]);
        status = rdSDS(fileID, "l1b_b5_data", recnum, 0, 1, npix, (VOIDP) & dataarr[4][0]);
        status = rdSDS(fileID, "l1b_b6_data", recnum, 0, 1, npix, (VOIDP) & dataarr[5][0]);
        status = rdSDS(fileID, "l1b_b7_data", recnum, 0, 1, npix, (VOIDP) & dataarr[6][0]);
        status = rdSDS(fileID, "l1b_b8_data", recnum, 0, 1, npix, (VOIDP) & dataarr[7][0]);
        for (i = 0; i < NOCBANDS; i++) {
            for (j = 0; j < npix; j++) {
                dataarr[i][j] = dataarr[i][j] & 8191;
            }
        }

        for (i = 0; i < npix; i++) {
            scanarr[i][0] = (float) dataarr[0][i] * m[0] + b[0];
            scanarr[i][1] = (float) dataarr[1][i] * m[1] + b[1];
            scanarr[i][2] = (float) dataarr[2][i] * m[2] + b[2];
            scanarr[i][3] = (float) dataarr[3][i] * m[3] + b[3];
            scanarr[i][4] = (float) dataarr[4][i] * m[4] + b[4];
            scanarr[i][5] = (float) dataarr[5][i] * m[5] + b[5];
            scanarr[i][6] = (float) dataarr[6][i] * m[6] + b[6];
            scanarr[i][7] = (float) dataarr[7][i] * m[7] + b[7];
        }
    }
    /* copy the scan array over all bands, and npix pixels */
    for (ip = 0; ip < npix; ip++) {

        // if solz NaN set navfail
        if (isnan(l1rec->solz[ip]))
            l1rec->navfail[ip] = 1;

        for (iw = 0; iw < nwave; iw++) {
            ib = bindx[iw];
            l1rec->Lt [ip * nwave + ib] = scanarr[ip][iw];

            /* Also set flags */
            if (gain[iw * scansPerScene + current_scan] > 0)
                l1rec->navwarn[ip] = 1;
            if (miss_qual[current_scan] > 0)
                l1rec->navwarn[ip] = 1;
            if (dataarr[iw][ip] > 1022)
                l1rec->hilt[ip] = 1;
            if (scanarr[ip][iw] <= 0.0)
                l1rec->navfail[ip] = 1;
        }
    }

    l1rec->npix = l1file->npix;

    return (0);
}

/* ------------------------------------------------------ */
/* navigation() - generates navigation info interpolated  */
/*              for each pixel, using info from a NASDA   */
/*              OCTS HDF file                             */
/*                                                        */

/* ------------------------------------------------------ */
int navigation(int32_t fileID) {
    int16_t det;

    int16_t *pxl;
    float *xctl;
    float *inlon, *inlat;
    float *insolz, *insola;
    float *insenz, *insena;
    float *usun, *pos;
    int *row;
    int *indx;
    float *yctl;
    float *in1, *in2;

    float iusun[3], ipos[3];
    float usun_sum;
    int i, j, k, l;
    int eline, nlon, nlat;
    int ilat[MAXOCLIN], ilon[MAXOCLIN];
    float out1[MAXOCLIN], out2[MAXOCLIN];
    int jout;
    float *spl_aux;
    float *x_ctl_ll, *y_ctl_ll, *z_ctl_ll;
    float *x_ctl_sol, *y_ctl_sol, *z_ctl_sol;
    float *x_ctl_sen, *y_ctl_sen, *z_ctl_sen;
    float * in_ptr[3][3], *out_ptr[3][2];

    inlon = (float *) calloc(ncol*scansPerScene, sizeof (float));
    inlat = (float *) calloc(ncol*scansPerScene, sizeof (float));
    insolz = (float *) calloc(ncol*scansPerScene, sizeof (float));
    insola = (float *) calloc(ncol*scansPerScene, sizeof (float));
    insenz = (float *) calloc(ncol*scansPerScene, sizeof (float));
    insena = (float *) calloc(ncol*scansPerScene, sizeof (float));
    usun = (float *) calloc(3 * scansPerScene, sizeof (float));
    pos = (float *) calloc(3 * scansPerScene, sizeof (float));
    row = (int *) calloc(scansPerScene, sizeof (int));
    pxl = (int16_t *) calloc(ncol, sizeof (int16_t));
    xctl = (float *) calloc(ncol, sizeof (float));
    indx = (int *) calloc(scansPerScene, sizeof (int));
    yctl = (float *) calloc(scansPerScene, sizeof (float));
    in1 = (float *) calloc(scansPerScene, sizeof (float));
    in2 = (float *) calloc(scansPerScene, sizeof (float));

    /* read the data sets */
    rdSDS(fileID, "det", 0, 0, 1, 1, (VOIDP) & det);
    rdSDS(fileID, "pxl", 0, 0, 0, 0, (VOIDP) pxl); /* control pt. pixels */
    rdSDS(fileID, "lon", 0, 0, 0, 0, (VOIDP) inlon); /* geodetic longitude */
    rdSDS(fileID, "lat", 0, 0, 0, 0, (VOIDP) inlat); /* geodetic latitude  */
    rdSDS(fileID, "sun_ref", 0, 0, 0, 0, (VOIDP) usun); /* geocentric ECR solar refercence vector  */
    rdSDS(fileID, "orb_vec", 0, 0, 0, 0, (VOIDP) pos); /* geocentric ECR S/C position           */

    det = det - 1;
    for (i = 0; i < ncol; i++) {
        pxl[i] = pxl[i] - 1;
    }

    /* convert sun vector to unit vector */
    for (i = 0; i < scansPerScene; i++) {
        usun_sum = 0;
        for (j = 0; j < 3; j++) {
            usun_sum = usun_sum + usun[3 * i + j] * usun[3 * i + j];
        }
        for (j = 0; j < 3; j++) {
            usun[3 * i + j] = usun[3 * i + j] / sqrt(usun_sum);
        }
    }

    /* define control point grid */
    if (want_verbose)
        printf("scansPerScene,ncol = %d %d \n", scansPerScene, ncol);
    for (i = 0; i < ncol; i++) {
        xctl[i] = (float) pxl[i];
    }
    for (i = 0; i < scansPerScene; i++) {
        row[i] = i * linesPerScan + det;
    }

    /* define output grid */
    eline = sline + nline - 1;
    nlon = npix; /* # of output pixs */
    if (nline < linesPerScan * scansPerScene) {
        nlat = nline; /* # of output lines */
    } else nlat = linesPerScan*scansPerScene;

    for (i = 0; i < nlon; i++) {
        ilon[i] = i + spix; /* pix #'s of output grid */
    }
    for (i = 0; i < nlat; i++) {
        ilat[i] = i + sline; /* line #'s of output grid */
    }

    /* create index of output lines to control point lines */
    k = 0;
    for (i = 0; i < nlat; i++) {
        for (j = 0; j < scansPerScene; j++) {
            if (ilat[i] == row[j]) {
                indx[k] = i;
                k++;
            }
        }
    }

    for (i = 0; i < nlat / linesPerScan; i++) {
        yctl[i] = (float) indx[i];
    }


    /* compute solar and sensor zenith and azimuth at input
    control points from info provided */
    if (want_verbose)
        printf("Computing view angles at input control points\n");
    for (i = 0; i < scansPerScene; i++) {
        for (j = 0; j < ncol; j++) {
            for (k = 0; k < 3; k++) {
                ipos[k] = pos[3 * i + k];
                iusun[k] = usun[3 * i + k];
            }
            CalcViewAngle(inlon[i * ncol + j], inlat[i * ncol + j], ipos, iusun);
            insolz[i * ncol + j] = solz1;
            insola[i * ncol + j] = sola1;
            insenz[i * ncol + j] = senz1;
            insena[i * ncol + j] = sena1;
        }
    }


    /* Compute unit vectors from lon/lat of control points */
    x_ctl_ll = (float *) calloc(ncol*scansPerScene, sizeof (float));
    y_ctl_ll = (float *) calloc(ncol*scansPerScene, sizeof (float));
    z_ctl_ll = (float *) calloc(ncol*scansPerScene, sizeof (float));

    x_ctl_sol = (float *) calloc(ncol*scansPerScene, sizeof (float));
    y_ctl_sol = (float *) calloc(ncol*scansPerScene, sizeof (float));
    z_ctl_sol = (float *) calloc(ncol*scansPerScene, sizeof (float));

    x_ctl_sen = (float *) calloc(ncol*scansPerScene, sizeof (float));
    y_ctl_sen = (float *) calloc(ncol*scansPerScene, sizeof (float));
    z_ctl_sen = (float *) calloc(ncol*scansPerScene, sizeof (float));


    for (i = 0; i < scansPerScene; i++) {
        for (j = 0; j < ncol; j++) {
            inlon[i * ncol + j] = inlon[i * ncol + j] / RADEG;
            inlat[i * ncol + j] = inlat[i * ncol + j] / RADEG;

            x_ctl_ll[i * ncol + j] = cos(inlat[i * ncol + j]) * cos(inlon[i * ncol + j]);
            y_ctl_ll[i * ncol + j] = cos(inlat[i * ncol + j]) * sin(inlon[i * ncol + j]);
            z_ctl_ll[i * ncol + j] = sin(inlat[i * ncol + j]);


            insola[i * ncol + j] = insola[i * ncol + j] / RADEG;
            insolz[i * ncol + j] = insolz[i * ncol + j] / RADEG;

            x_ctl_sol[i * ncol + j] = cos(insolz[i * ncol + j]) * cos(insola[i * ncol + j]);
            y_ctl_sol[i * ncol + j] = cos(insolz[i * ncol + j]) * sin(insola[i * ncol + j]);
            z_ctl_sol[i * ncol + j] = sin(insolz[i * ncol + j]);


            insena[i * ncol + j] = insena[i * ncol + j] / RADEG;
            insenz[i * ncol + j] = insenz[i * ncol + j] / RADEG;

            x_ctl_sen[i * ncol + j] = cos(insenz[i * ncol + j]) * cos(insena[i * ncol + j]);
            y_ctl_sen[i * ncol + j] = cos(insenz[i * ncol + j]) * sin(insena[i * ncol + j]);
            z_ctl_sen[i * ncol + j] = sin(insenz[i * ncol + j]);
        }
    }

    in_ptr[0][0] = x_ctl_ll;
    in_ptr[0][1] = y_ctl_ll;
    in_ptr[0][2] = z_ctl_ll;
    in_ptr[1][0] = x_ctl_sol;
    in_ptr[1][1] = y_ctl_sol;
    in_ptr[1][2] = z_ctl_sol;
    in_ptr[2][0] = x_ctl_sen;
    in_ptr[2][1] = y_ctl_sen;
    in_ptr[2][2] = z_ctl_sen;

    out_ptr[0][0] = lon;
    out_ptr[0][1] = lat;
    out_ptr[1][0] = sola;
    out_ptr[1][1] = solz;
    out_ptr[2][0] = sena;
    out_ptr[2][1] = senz;


    /* we now have all the info at each control point, so we
    can interpolate to all pixels, all lines */
    /* interpolate angles across each control point line */
    if (want_verbose)
        printf("Interpolating rows for longitude/azimuth\n");
    spl_aux = (float *) calloc(ncol, sizeof (float));

    for (i = 0; i < scansPerScene; i++) {
        jout = row[i] - sline;
        if ((row[i] >= sline) && (row[i] <= eline)) {
            for (l = 0; l < 3; l++) {
                spline(xctl, in_ptr[l][0] + i*ncol, ncol, 1e30, 1e30, spl_aux);
                for (j = 0; j < nlon; j++)
                    splint(xctl, in_ptr[l][0] + i * ncol, spl_aux, ncol,
                        (float) ilon[j], out_ptr[l][0] + jout * npix + j);

                spline(xctl, in_ptr[l][1] + i*ncol, ncol, 1e30, 1e30, spl_aux);
                for (j = 0; j < nlon; j++)
                    splint(xctl, in_ptr[l][1] + i * ncol, spl_aux, ncol,
                        (float) ilon[j], out_ptr[l][1] + jout * npix + j);
            }
        }
    }
    free(spl_aux);


    /* fill missing lines by interpolating columns */
    if (want_verbose)
        printf("Interpolating columns for longitude/azimuth\n");
    spl_aux = (float *) calloc(nlat / linesPerScan, sizeof (float));

    for (i = 0; i < nlon; i++) {
        for (l = 0; l < 3; l++) {
            for (k = 0; k < nlat / linesPerScan; k++) {
                in1[k] = *(out_ptr[l][0] + indx[k] * npix + i);
                in2[k] = *(out_ptr[l][1] + indx[k] * npix + i);
            }
            spline(yctl, in1, nlat / linesPerScan, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in1, spl_aux, nlat / linesPerScan, (float) j, (float *) &out1[j]);

            spline(yctl, in2, nlat / linesPerScan, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in2, spl_aux, nlat / linesPerScan, (float) j, (float *) &out2[j]);

            for (j = 0; j < nlat; j++) {
                *(out_ptr[l][0] + j * npix + i) = atan2(out2[j], out1[j]) * RADEG;
                if (l >= 1 && *(out_ptr[l][0] + j * npix + i) < 0) {
                    *(out_ptr[l][0] + j * npix + i) += 360;
                }
            }
        }
    }
    free(spl_aux);


    if (want_verbose)
        printf("Interpolating rows for latitude/zenith\n");
    spl_aux = (float *) calloc(ncol, sizeof (float));

    for (i = 0; i < scansPerScene; i++) {
        jout = row[i] - sline;
        if ((row[i] >= sline) && (row[i] <= eline)) {
            for (l = 0; l < 3; l++) {
                spline(xctl, in_ptr[l][2] + i*ncol, ncol, 1e30, 1e30, spl_aux);
                for (j = 0; j < nlon; j++)
                    splint(xctl, in_ptr[l][2] + i * ncol, spl_aux, ncol,
                        (float) ilon[j], out_ptr[l][1] + jout * npix + j);
            }
        }
    }
    free(spl_aux);


    /* fill missing lines by interpolating columns */
    if (want_verbose)
        printf("Interpolating columns for latitude/zenith\n");
    spl_aux = (float *) calloc(nlat / linesPerScan, sizeof (float));

    for (i = 0; i < nlon; i++) {
        for (l = 0; l < 3; l++) {
            for (k = 0; k < nlat / linesPerScan; k++) {
                in1[k] = *(out_ptr[l][1] + indx[k] * npix + i);
            }
            spline(yctl, in1, nlat / linesPerScan, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in1, spl_aux, nlat / linesPerScan, (float) j, (float *) &out1[j]);

            for (j = 0; j < nlat; j++) {
                *(out_ptr[l][1] + j * npix + i) = asin(out1[j]) * RADEG;
            }
        }
    }
    free(spl_aux);


    free(x_ctl_ll);
    free(y_ctl_ll);
    free(z_ctl_ll);
    free(x_ctl_sol);
    free(y_ctl_sol);
    free(z_ctl_sol);
    free(x_ctl_sen);
    free(y_ctl_sen);
    free(z_ctl_sen);

    free(inlon);
    free(inlat);
    free(insolz);
    free(insola);
    free(insenz);
    free(insena);
    free(usun);
    free(pos);
    free(row);
    free(pxl);
    free(xctl);
    free(indx);
    free(yctl);
    free(in1);
    free(in2);

    return (0);
}



/* ------------------------------------------------------ */
/* CalcViewAngle() - calculates the solar and sensor view */
/*             angles from one geodetic longitude and     */
/*             latitude, position vector, and solar unit  */
/*             vector.                                    */
/*                                                        */
/* INPUT       DESCRIPTION                                */
/* -----       -----------                                */
/*  lon        longitude of pixel in degrees              */
/*  lat        latitude of pixel in degrees               */
/*  pos[3]     vector from earth center to spacecraft in  */
/*             ECR, km                                    */
/*  usun[3]    unit vector from earth center to sun in    */
/*             ECR                                        */
/*                                                        */
/* OUTPUT      DESCRIPTION                                */
/* ------      -----------                                */
/*  solz       solar zenith angle of pixel in degrees     */
/*  sola       solar azimuth angle of pixel in degrees    */
/*  senz       sensor zenith angle of pixel in degrees    */
/*  sena       sensor azimuth angle of pixel in degrees   */
/*                                                        */

/* ------------------------------------------------------ */
int CalcViewAngle(float lon1, float lat1, float pos[3],
        float usun[3]) {
    float rmat[3][3]; /* rotation matrix */
    float rlon, rlat;
    float up[3], upxy, ea[3], no[3];
    float phi, R;
    float gvec[3], scvec[3], senl[3], sunl[3];
    int i, j;

    rlon = lon1 * PI / 180.0;
    rlat = lat1 * PI / 180.0;

    /* First, we must define the axes (in ECR space) of a
       pixel-local coordinate system, with the z-axis along
       the geodetic pixel normal, x-axis pointing east, and
       y-axis pointing north.  */
    up[0] = cos(rlat) * cos(rlon);
    up[1] = cos(rlat) * sin(rlon);
    up[2] = sin(rlat);
    upxy = sqrt(up[0] * up[0 ] + up[1] * up[1]);
    ea[0] = -up[1] / upxy;
    ea[1] = up[0] / upxy;
    ea[2] = 0.0;
    /* calculate cross product of up and ea */
    no[0] = up[1] * ea[2] - ea[1] * up[2];
    no[1] = up[2] * ea[0] - ea[2] * up[0];
    no[2] = up[0] * ea[1] - ea[0] * up[1];


    /* Now we need the pixel-to-spacecraft vector in ECR,
       Compute geocentric pixel location vector in km and subtract
       from spacecraft position. */
    /* geocentric latitude */
    phi = atan(tan(rlat)*(1 - FL)*(1 - FL));
    /* dist to Earth surface */
    R = RADIUS * (1 - FL) / sqrt(1 - (2 - FL) * FL * (cos(phi) * cos(phi)));
    gvec[0] = R * cos(phi) * cos(rlon);
    gvec[1] = R * cos(phi) * sin(rlon);
    gvec[2] = R * sin(phi);
    for (i = 0; i < 3; i++) {
        scvec[i] = pos[i] - gvec[i];
    }

    /* Now we can transform the pixel-to-spacecraft and Sun 
       vectors into the local frame. */
    for (i = 0; i < 3; i++) {
        rmat[0][i] = ea[i];
        rmat[1][i] = no[i];
        rmat[2][i] = up[i];
    }
    for (i = 0; i < 3; i++) {
        senl[i] = 0;
        sunl[i] = 0;
        for (j = 0; j < 3; j++) {
            senl[i] = senl[i] + rmat[i][j] * scvec[j];
            sunl[i] = sunl[i] + rmat[i][j] * usun[j];
        }
    }

    /* Compute the solar zenith and azimuth */
    solz1 = RADEG * atan(sqrt(sunl[0] * sunl[0] + sunl[1] * sunl[1]) / sunl[2]);
    if (solz1 < 0.0) solz1 += 180.0;

    if (solz1 > 0.05) {
        sola1 = RADEG * (atan2(sunl[0], sunl[1]));
    } else {
        sola1 = 0.0;
    }
    if (sola1 < 0.0) {
        sola1 = sola1 + 360.0;
    }

    /* Compute the sensor zenith and azimuth  */
    senz1 = RADEG * atan(sqrt(senl[0] * senl[0] + senl[1] * senl[1]) / senl[2]);
    if (senz1 < 0.0) senz1 += 180.0;
    if (senz1 > 0.05) {
        sena1 = RADEG * (atan2(senl[0], senl[1]));
    } else {
        sena1 = 0.0;
    }
    if (sena1 < 0.0) {
        sena1 = sena1 + 360.0;
    }
    return (0);
}


/* ------------------------------------------------------ */
/* LeapCheck() - checks if the year is a leap year:       */
/*             returns -1 if leap year                    */
/*             returns 0 if non-leap year                 */
/*                                                        */

/* ------------------------------------------------------ */
int LeapCheck(int yr) {
    int leap;

    leap = 0;
    if ((yr % 4) == 0) {
        leap = -1;
    }
    if ((yr % 100) == 0) {
        leap = 0;
    }
    if ((yr % 400) == 0) {
        leap = -1;
    }
    return (leap);
}

/* ------------------------------------------------------ */
/* closel1_octs() - closes the level 1 OCTS HDF file  */
/*                                                        */

/* ------------------------------------------------------ */
int closel1_octs(filehandle *l1file) {
    free(lon);
    free(lat);
    free(senz);
    free(sena);
    free(solz);
    free(sola);

    if (l1file->format == FT_OCTSL1A) {
        free(start_scan);
        free(tilt_flag);
        free(inst_temp);
    }

    /* End access to the HDF file */
    SDend(l1file->sd_id);

    return (0);
}

void reform_octs_time(char *time) {
    int year, month, day, hour, min;
    double sec;

    sscanf(time, "%4d%2d%2d %2d:%2d:%lf", &year, &month, &day, &hour, &min, &sec);
    sec = sec + 60. * min + 3600. * hour;
    strncpy(time, ydhmsf(ymds2unix(year, month, day, sec), 'G'), 16);
    time[16] = '\0';
}
