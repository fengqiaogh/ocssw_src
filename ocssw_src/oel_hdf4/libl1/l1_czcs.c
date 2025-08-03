/*
 *  W. Robinson, SAIC, 10 Dec 2004  I/O routines for CZCS
 */
#include "l1.h"
#include <hdf4utils.h>
#include "mfhdf.h"
#include "l1_czcs.h"
#include <math.h>
#include <libnav.h>

static int16 syear, sday; /* data start date                  */
static int32 smsec; /* data start time                  */
static int16 eyear, eday; /* data end date                    */
static int32 emsec; /* data end time                    */
static int32 nscan; /* number of scans                  */
static int32 npix; /* number pixels per scan           */

static char dtype[8];

static int32 nsta;
static int32 ninc;

static uint8 *counts, cz_band_present;
static int32 *msec;
static int16 *gain;
static float32 *tilt, *att_ang, *slope, *intercept;
static float32 *ctl_pt_lat, *ctl_pt_lon, *pos, *pos_err;
static int32 nctl_pt, *ctl_pt_cols;
static float *ctl_pt_vx, *ctl_pt_vy, *ctl_pt_vz, *y2_vx, *y2_vy, *y2_vz, *ctl_pt_x;
static float *lt750; /* internal 750 mn data source */
static char *ring_sat; /* set to 1 if 443, 520 or 550 are saturated for ringing 
                    mask computation */
static double fileEpochTime; // unix time of file start of day

int czcs_ring(int gain, float *lt750, char *ring_sat, l1str *l1rec)
/*******************************************************************

   czcs_ring

   purpose: use the Mueller ringing flagging algorithm to set the stray 
     light mask for CZCS

   Returns type: int - 0 if OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               gain             I      gain for the line of data
      float *           lt750            I      750 nm total rads
      char *            ring_sat         I      1 if bands 1,2,or 3 are 
                                                saturated = use the 750 
                                                excess rad in this case
      l1str *           l1rec           I/O     L1 struc including tot rads,
                                                flags

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       31 May 2005     Original development
      W. Robinson       31 Jul 2007     repaired problem that occurs in SGIs
                                        if bsum = 0, avoid log and set ring
                                        distance to 0

  Based on the algorithm of Mueller, J. L., Nimbus-7 CZCS: electronic
  overshoot due to cloud reflectance, Appl. Optics, Vol 27, No 3, 1 Feb, 
  1988, pp 438 - 440
  Modified to only include the excess radiance in 750 for pixels where 
  bands 1, 2, or 3 have saturated.  This is to reduce the masking from
  coastal areas where the vis is not as bright as for cloud

 *******************************************************************/ {
    static float gtrans[4] = {1.0, 1.25, 1.5, 2.1};
    /* note that the values below have +- of 9.7 and 3.3 in paper */
    /*  standard values followed by low, high confidence values
      static float mueller_int = 3.9, mueller_slp = 30.8;
      static float mueller_int = -6.2, mueller_slp = 27.5;
      static float mueller_int = 13.6, mueller_slp = 34.1;
     */
    static float mueller_int = 3.9, mueller_slp = 30.8;

    static float l0_750 = 2.45;

    float gfac, lthresh, gfac_ln, excess_brit[1968], bsub, bsum;
    int i, ipx, btarg, last_btarg, igrp, iavg, pxloc, n_ring;

    /*
     * set up gfac - gain factor for line, lthresh, brightness threshold
     * and flags for current and last pixel above threshold radiance
     */
    gfac = gtrans[ gain - 1];
    lthresh = l0_750 / gfac;
    gfac_ln = log(gfac);
    btarg = 0;
    last_btarg = 0;

    /*
     * go through pixels and set the mask for pixels downstream of bright target
     */
    for (i = 0; i < npix; i++) {
        if (*(lt750 + i) > lthresh) {
            btarg = 1;
            excess_brit[i] =
                    (*(ring_sat + i) == 1) ? *(lt750 + i) - lthresh : 0.;
            l1rec->stlight[i] = 1;
        } else
            excess_brit[i] = 0.;

        if ((btarg == 0) && (last_btarg == 1)) {
            /*
             *  sum up excess brightness in 5 groups of 10, weighted by f(dist)
             */
            bsum = 0.;
            for (igrp = 0; igrp < 5; igrp++) {
                bsub = 0.;
                for (iavg = 0; iavg < 10; iavg++) {
                    pxloc = i - 1 - (iavg + igrp * 10);
                    if (pxloc >= 0)
                        bsub += excess_brit[ pxloc ];
                }
                bsub = bsub / 10.;
                bsum += bsub * exp(-0.32 * (igrp + 1));
            }
            /*
             * compute # pixels to mask and mark the pixels in that range
             */
            bsum = (bsum > 450.) ? 450. : bsum;
            if (bsum > 0.) {
                n_ring = mueller_int + mueller_slp * (log(bsum) + gfac_ln);
                n_ring = (n_ring < 0) ? 0 : n_ring;
            } else
                n_ring = 0;

            for (ipx = 0; ipx < n_ring; ipx++) {
                pxloc = ipx + i;
                if (pxloc < npix)
                    l1rec->stlight[pxloc] = 1;
            }
        }
        /*
         *  lastly, reset to go to next pixel
         */
        last_btarg = btarg;
        btarg = 0;
    }
    return 0;
}

#define NBND 4
#define NGAIN 4
#define NEPOCH 5

int get_czcscal(char *file, int orbit, int16 year, int16 day, int32 msec, short l1acnt[], float slope750, float intercept750, int16 igain, float32 l1brads[])
/******************************************************************

   get_czcscal
   purpose: replicate the Evans & Gordon calibration

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      Sean Bailey       08 Feb 2006     Original development

  Based on the original czcscaleg.f
  Modified to use cal_czcs.hdf

 *******************************************************************/ {
    static int firstCall = 1;
    static float slope_coeff[NGAIN][NBND_CZCS];
    static float offsets[NGAIN][NBND_CZCS];
    static float rk[NGAIN][NBND_CZCS];
    static float tdecay[NEPOCH][NBND_CZCS];
    static int orbitepoch[NEPOCH];
    static float timeconst[2];
    static float time_dep[NBND_CZCS][NGAIN][3];

    int epochidx = 0;
    int i, j, iorbit;
    float rden, rnum;
    float tcorr, slp, S;
    static int orbepoch670 = 6750;
    float g[NBND];
    float64 jsec = yds2unix(year, day, ((double) (msec)) / 1000.0);
    float64 refjsec = yds2unix(1978, 278, (double) 0.0);
    float64 decayday = (jsec - refjsec) / 86400.;

    igain -= 1;

    if (firstCall) {
        char name [H4_MAX_NC_NAME] = "";
        char sdsname[H4_MAX_NC_NAME] = "";
        int32 sd_id;
        int32 sds_id;
        int32 rank;
        int32 nt;
        int32 dims[H4_MAX_VAR_DIMS];
        int32 nattrs;
        int32 start[3] = {0, 0, 0};
        int32 status;
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

        strcpy(sdsname, "slope_coeff");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) slope_coeff);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "offsets");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) offsets);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "rk");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) rk);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "dec");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) tdecay);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "iorbdec");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) orbitepoch);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "timeconst");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) timeconst);
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

        firstCall = 0;
    }

    iorbit = MIN(orbit, 39000);
    for (j = 0; j < NEPOCH - 1; j++) {
        if ((iorbit >= orbitepoch[j]) && (iorbit <= orbitepoch[j + 1]))
            epochidx = j;
    }

    rden = (float) orbitepoch[epochidx + 1] - (float) orbitepoch[epochidx];

    for (i = 0; i < NBND; i++) {
        /* Calculate time dependent correction... */
        /*  double exp	
                tcorr = time_dep[i][igain][0] - time_dep[i][igain][1] * (1.0 - exp(timeconst[0]*decayday)) - time_dep[i][igain][2] *(1.0 - exp(timeconst[1]*decayday));
         */
        /* single exp */
        tcorr = time_dep[i][igain][0] - time_dep[i][igain][1] * (1.0 - exp(-1 * time_dep[i][igain][2] * decayday));

        rnum = tdecay[epochidx + 1][i] - tdecay[epochidx][i];
        slp = rnum / rden;

        /* Antione modification for 670nm */
        if (i == 3 && iorbit > orbepoch670) {
            rden = (float) orbitepoch[4] - (float) orbepoch670;
            rnum = tdecay[4][i] - tdecay[3][i];
            slp = rnum / rden;
            S = (float) iorbit - (float) orbepoch670;
            g[i] = 1.0 / (S * slp + tdecay[3][i]);

        } else {
            S = (float) iorbit - (float) orbitepoch[epochidx];
            g[i] = 1.0 / (S * slp + tdecay[epochidx][i]);
        }
        /*fprintf(stderr,"i: %d g: %8.6f tcorr: %8.6f\n",i,g[i],tcorr);*/

        l1brads[i] = (tcorr * g[i] * rk[igain][i] * slope_coeff[igain][i] * (float) l1acnt[i]) + offsets[igain][i];

    }

    l1brads[4] = slope750 * (float) l1acnt[4] + intercept750;

    return 0;
}

/*
   W. Robinson, SAIC, 6 Jan 2006  add code to read the position and 
   position error SDSes
 */

int openl1_czcs(filehandle *file) {

    /*                                                                 */
    /* get_l1a_open interface                                          */
    /*                                                                 */
    int32 fileID;
    int32 status;
    int i;

    /* open the file and get the file ID */
    fileID = SDstart(file->name, DFACC_RDONLY);

    if (fileID < 0) {
        fprintf(stderr,
                "-E %s Line %d: Error opening %s for reading.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }
    status = getHDFattr(fileID, "Start Year", "", (VOIDP) & syear);
    status = getHDFattr(fileID, "Start Day", "", (VOIDP) & sday);
    status = getHDFattr(fileID, "Start Millisec", "", (VOIDP) & smsec);
    status = getHDFattr(fileID, "End Year", "", (VOIDP) & eyear);
    status = getHDFattr(fileID, "End Day", "", (VOIDP) & eday);
    status = getHDFattr(fileID, "End Millisec", "", (VOIDP) & emsec);
    status = getHDFattr(fileID, "Number of Scan Lines", "", (VOIDP) & nscan);
    status = getHDFattr(fileID, "Data Type", "", (VOIDP) & dtype);

    status = getHDFattr(fileID, "Pixels per Scan Line", "", (VOIDP) & npix);
    status = getHDFattr(fileID, "LAC Pixel Start Number", "", (VOIDP) & nsta);
    status = getHDFattr(fileID, "LAC Pixel Subsampling", "", (VOIDP) & ninc);

    status = getHDFattr(fileID, "Orbit Number", "", (VOIDP) & file->orbit_number);
    /*  WDR not in czcs l1 file, so comment for now
     status = getHDFattr(fileID, "Orbit Node Longitude", "", 
         (VOIDP)&file->orbit_node_lon);
     status = getHDFattr(fileID, "Node Crossing Time", "", 
         (VOIDP)&file->node_crossing_time[0]);
     */
    status = getHDFattr(fileID, "Number of Pixel Control Points", "",
            (VOIDP) & nctl_pt);
    status = getHDFattr(fileID, "Parameter Presence Code", "",
            (VOIDP) & cz_band_present);

    /* call cdata.f to initialize global FORTRAN common block data	*/
    cdata_();

    // set file epoch time
    fileEpochTime = yds2unix(syear, sday, 0.0);

    file->npix = npix;
    file->nscan = nscan;
    file->sensorID = CZCS;
    file->sd_id = fileID;
    strcpy(file->spatialResolution, "825 m");

    msec = (int32 *) calloc(nscan, sizeof (int32));
    status = rdSDS(file->sd_id, "msec", 0, 0, 0, 0, (VOIDP) msec);

    att_ang = (float32 *) calloc(nscan * 3, sizeof ( float32));
    status = rdSDS(file->sd_id, "att_ang", 0, 0, 0, 0, (VOIDP) att_ang);

    ctl_pt_cols = (int32 *) calloc(nctl_pt, sizeof ( int32));
    status = rdSDS(file->sd_id, "cntl_pt_cols", 0, 0, 0, 0,
            (VOIDP) ctl_pt_cols);
    ctl_pt_x = (float *) calloc(nctl_pt, sizeof ( float));
    for (i = 0; i < nctl_pt; i++)
        ctl_pt_x[i] = (float) ctl_pt_cols[i] - 1.;

    tilt = (float32 *) calloc(nscan, sizeof (float32));
    status = rdSDS(file->sd_id, "tilt", 0, 0, 0, 0, (VOIDP) tilt);

    slope = (float32 *) calloc(nscan * 6, sizeof (float32));
    status = rdSDS(file->sd_id, "slope", 0, 0, 0, 0, (VOIDP) slope);
    intercept = (float32 *) calloc(nscan * 6, sizeof (float32));
    status = rdSDS(file->sd_id, "intercept", 0, 0, 0, 0, (VOIDP) intercept);
    /* get nav orbit data  */
    pos = (float32 *) malloc(nscan * 3 * sizeof (float32));
    status = rdSDS(file->sd_id, "orb_vec", 0, 0, 0, 0, (VOIDP) pos);
    pos_err = (float32 *) malloc(nscan * sizeof (float32));
    status = rdSDS(file->sd_id, "pos_err", 0, 0, 0, 0, (VOIDP) pos_err);

    gain = (int16 *) malloc(nscan * sizeof ( int16));
    counts = (uint8 *) malloc(npix * NBND_CZCS * sizeof ( uint8));
    ctl_pt_lat = (float32 *) malloc(nctl_pt * sizeof ( float32));
    ctl_pt_lon = (float32 *) malloc(nctl_pt * sizeof ( float32));
    ctl_pt_vx = (float32 *) malloc(nctl_pt * sizeof ( float32));
    ctl_pt_vy = (float32 *) malloc(nctl_pt * sizeof ( float32));
    ctl_pt_vz = (float32 *) malloc(nctl_pt * sizeof ( float32));
    y2_vx = (float32 *) malloc(npix * sizeof ( float32));
    y2_vy = (float32 *) malloc(npix * sizeof ( float32));
    y2_vz = (float32 *) malloc(npix * sizeof ( float32));
    lt750 = (float *) malloc(npix * sizeof ( float));
    ring_sat = (char *) malloc(npix * sizeof ( char));

    return (status);
}

/*
 W. Robinson       31 May 2005     add ringing masking call
 W. Robinson, SAIC, 6 Jan 2005     add code to compute sat angles with
                                   cz_posll_2_satang if pos_err is acceptable
 J. Gales          23 Sep 2005     add kludge for subsetted pixel range in
                                   satang
   (note that satang may not give proper values in pixel subsets with this 
    but it may be academic with...)
 W. Robinson, SAIC, 10 Sep 2010    use vectors instead of lat, lon when
                                   doing interpolation of ctl pt lat, lon to
                                   all samples
 */

int readl1_czcs(filehandle *file, int32_t recnum, l1str *l1rec) {
    /*void czcscal_( int *, float[], float[], short *, int *, float * );*/
    void satang_(double *, double *, float *, float *, float *, float *,
            float *, float *, float *, float *);
    void sunangs_(int *, int *, float *, float *, float *, float *, float *);
    int ll2vec(float *, float *);
    int vec2ll(float *, float *);
    short cnt_vec[NBND_CZCS];
    float lt_lcl[NBND_CZCS], yout1, yout2, yout3, llvec[2], vec[3];
    int ipx, ibnd, orbit, i, tpix;
    uint8 cal_sum[5], cal_scan[6];
    int32 status, navbad, ltsat;
    double pi, radeg;
    int32_t ib;

    int32_t nwave = l1rec->l1file->nbands;
    int32_t *bindx = l1rec->l1file->bindx;
    char *cal_path = l1_input->calfile;
    float lonbuf[1968], latbuf[1968], senzbuf[1968], senabuf[1968];

    /*
     *  set scan time
     */
    if(msec[recnum] < msec[0]) {
        // adjust for day rollover (add 86400 sec in a day)
        l1rec->scantime = fileEpochTime + (double)msec[recnum] / 1000.0 + 86400.0;
    } else {
        l1rec->scantime = fileEpochTime + (double)msec[recnum] / 1000.0;
    }

    double scanSecOfDay, scanMsecOfDay;
    int16_t scanYear;
    int16_t scanDay;
    unix2yds(l1rec->scantime, &scanYear, &scanDay, &scanSecOfDay);
    scanMsecOfDay = scanSecOfDay * 1000.0;

    /*
     *  read in the line of counts, latitude, longitude
     */
    status = rdSDS(file->sd_id, "band1", recnum, 0, 1, npix, (VOIDP) counts);
    status = rdSDS(file->sd_id, "band2", recnum, 0, 1, npix,
            (VOIDP) (counts + npix));
    status = rdSDS(file->sd_id, "band3", recnum, 0, 1, npix,
            (VOIDP) (counts + 2 * npix));
    status = rdSDS(file->sd_id, "band4", recnum, 0, 1, npix,
            (VOIDP) (counts + 3 * npix));
    status = rdSDS(file->sd_id, "band5", recnum, 0, 1, npix,
            (VOIDP) (counts + 4 * npix));

    status = rdSDS(file->sd_id, "gain", recnum, 0, 1, 1, (VOIDP) gain);
    status = rdSDS(file->sd_id, "latitude", recnum, 0, 1, nctl_pt,
            (VOIDP) ctl_pt_lat);
    status = rdSDS(file->sd_id, "longitude", recnum, 0, 1, nctl_pt,
            (VOIDP) ctl_pt_lon);

    status = rdSDS(file->sd_id, "cal_sum", recnum, 0, 1, 5, (VOIDP) cal_sum);
    status = rdSDS(file->sd_id, "cal_scan", recnum, 0, 1, 6, (VOIDP) cal_scan);
    /*
     * flag setting for entire line: set navfail if cal_sum shows problems 
     * (suspect it is * bad ephemeris or atitude) and set hilt if bands 
     * missing
     */
    ltsat = 0;
    navbad = 0;
    if ((cz_band_present & 0XF8) != 0XF8)
        ltsat = 1;
    else {
        if ((cal_sum[3] != 0) || (cal_sum[4] != 0))
            navbad = 1;
        if ((cal_scan[0] != 0) || (cal_scan[1] != 0) ||
                (cal_scan[2] != 0) || (cal_scan[3] != 0) ||
                (cal_scan[4] != 0)) {
            ltsat = 1;
            navbad = 1;
        }
    }
    /*
     * calibrate the radiances and set flags for pixels
     */
    for (ipx = 0; ipx < npix; ipx++) {
        for (ibnd = 0; ibnd < NBND_CZCS; ibnd++) {
            cnt_vec[ibnd] = *(counts + ipx + ibnd * npix);
            if ((cnt_vec[ibnd] == 255) || (ltsat == 1)) l1rec->hilt[ipx] = 1;
        }
        *(ring_sat + ipx) = ((cnt_vec[0] == 255) || (cnt_vec[1] == 255)
                || (cnt_vec[2] == 255)) ? 1 : 0;
        /*
         *  call the fortran calibration routine 
         */
        orbit = (int) file->orbit_number;
        /*czcscal_( &orbit, slope + recnum * 6, 
          intercept + recnum * 6, cnt_vec, &lgain, lt_lcl );*/
        status = get_czcscal(cal_path, orbit, scanYear, scanDay, scanMsecOfDay,
                cnt_vec, slope[4], intercept[4], gain[0], lt_lcl);
        /*
         *  assign the calibrated radiances
         */
        for (ibnd = 0; ibnd < nwave; ibnd++) {
            ib = bindx[ ibnd ];
            l1rec->Lt[ ipx * nwave + ib ] = lt_lcl[ibnd];
        }
        /*
         *  save the 750 nm data here
         */
        *(lt750 + ipx) = lt_lcl[4];

        if (navbad == 1) l1rec->navfail[ipx] = 1;
    }
    /*
     *  set the ringing mask
     */
    czcs_ring(gain[0], lt750, ring_sat, l1rec);
    /*
     *  get navigation at all points from control point lat and lon values
     *  use spline interpolation to fill in, do in vector space to avoid date 
     *  line discontinuity
     */
    for (i = 0; i < nctl_pt; i++) {
        llvec[0] = ctl_pt_lat[i];
        llvec[1] = ctl_pt_lon[i];
        if (ll2vec(llvec, vec) == 0) {
            ctl_pt_vx[i] = *vec;
            ctl_pt_vy[i] = *(vec + 1);
            ctl_pt_vz[i] = *(vec + 2);
        } else {
            fprintf(stderr, "-E %s Line %d: ll2vec failure.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }
    spline(ctl_pt_x, ctl_pt_vx, nctl_pt, 1e30, 1e30, y2_vx);
    spline(ctl_pt_x, ctl_pt_vy, nctl_pt, 1e30, 1e30, y2_vy);
    spline(ctl_pt_x, ctl_pt_vz, nctl_pt, 1e30, 1e30, y2_vz);
    for (i = 0; i < npix; i++) {
        tpix = i * ninc /*+ spix*/;
        splint(ctl_pt_x, ctl_pt_vx, y2_vx, nctl_pt, tpix, &yout1);
        splint(ctl_pt_x, ctl_pt_vy, y2_vy, nctl_pt, tpix, &yout2);
        splint(ctl_pt_x, ctl_pt_vz, y2_vz, nctl_pt, tpix, &yout3);

        *vec = yout1;
        *(vec + 1) = yout2;
        *(vec + 2) = yout3;
        vec2ll(vec, llvec);

        l1rec->lon[i] = llvec[1];
        l1rec->lat[i] = llvec[0];
    }
    /*
     * calculate geometry.  For sensor angles, use the orbit info if the 
     * position error is available and within error threshold
     */
    if ((*(pos_err + recnum) >= 0.) &&
            (*(pos_err + recnum) < POS_ERR_THRESH)) {
        cz_posll_2_satang((pos + 3 * recnum), npix, l1rec->lat, l1rec->lon,
                l1rec->senz, l1rec->sena);
    } else {
        pi = OEL_PI;
        radeg = OEL_RADEG;
        for (i = 0; i < 1968; i++) {
            lonbuf[i] = 0.0;
            latbuf[i] = 0.0;
            senzbuf[i] = 0.0;
            senabuf[i] = 0.0;
        }

        for (i = 0; i < npix; i++) {
            lonbuf[i] = l1rec->lon[i];
            latbuf[i] = l1rec->lat[i];
        }

        satang_(&pi, &radeg, tilt + recnum, att_ang + 3 * recnum + 1,
                att_ang + 3 * recnum + 2, att_ang + 3 * recnum, lonbuf,
                latbuf, senzbuf, senabuf);

        for (i = 0; i < npix; i++) {
            l1rec->senz[i] = senzbuf[i];
            l1rec->sena[i] = senabuf[i];
        }
    }

    int intYear = scanYear;
    int intDay = scanDay;
    // gmt = time of day in hours
    float gmt = scanSecOfDay / 3600.0;
    for (i = 0; i < npix; i++) {
        sunangs_(&intYear, &intDay, &gmt, (l1rec->lon) + i, (l1rec->lat) + i,
                (l1rec->solz) + i, (l1rec->sola) + i);
    }

    return (status);
}

/*  W. Robinson, SAIC, 6 Jan 2005   include the pos, pos_err arrays in free */

int closel1_czcs(filehandle *file) {
    free(msec);
    free(tilt);
    free(att_ang);
    free(counts);
    free(ctl_pt_lat);
    free(ctl_pt_lon);
    free(ctl_pt_vx);
    free(ctl_pt_vy);
    free(ctl_pt_vz);
    free(y2_vx);
    free(y2_vy);
    free(y2_vz);
    free(ctl_pt_x);
    free(ctl_pt_cols);
    free(slope);
    free(intercept);
    free(lt750);
    free(ring_sat);
    free(pos);
    free(pos_err);

    SDend(file->sd_id);

    return (0);
}

#define re 6378.137
#define f 1. / 298.257
#define omf2 ( 1 - f ) * ( 1 - f )

int cz_posll_2_satang(float *pos, int npix, float *lat, float *lon,
        float *senz, float *sena)
/*******************************************************************

   cz_posll_2_satang

   purpose: convert satellite position, lat, lon into satellite zenith 
     and azimuth for a line of czcs data

   Returns type: int - 0 if no problems

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float *           pos              I      position (x, y, z) of 
                                                satellite in meters
      int               npix             I      # pixels of lat, lon to find
                                                sat angles for
      float *           lat              I      latitude array in degrees E
      float *           lon              I      longitude array in degrees E
      float *           senz             O      sensor zenith in degrees
      float *           sena             O      sensor azimuth in degrees

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC  5-Jan-2006     Original development

 *******************************************************************/ {
    double up[3], ea[3], no[3], radeg, upxy, xlmat[3][3], xlatg, gv[3];
    double r, rh[3], rl[3];
    int i, j;

    radeg = 90. / asin(1.);

    for (i = 0; i < npix; i++) {
        /*
         *  Compute the local vertical, East and North unit vectors
         */
        up[0] = cos(*(lat + i) / radeg) * cos(*(lon + i) / radeg);
        up[1] = cos(*(lat + i) / radeg) * sin(*(lon + i) / radeg);
        up[2] = sin(*(lat + i) / radeg);

        upxy = sqrt(up[0] * up[0] + up[1] * up[1]);

        ea[0] = -up[1] / upxy;
        ea[1] = up[0] / upxy;
        ea[2] = 0.;

        cross_prod(up, ea, no);

        /*
         *  Store vectors in 3x3 array
         */
        for (j = 0; j < 3; j++) {
            xlmat[0][j] = ea[j];
            xlmat[1][j] = no[j];
            xlmat[2][j] = up[j];
        }
        /*
         *  Compute geocentric position vector
         */
        xlatg = radeg * atan(tan(*(lat + i) / radeg) * omf2);
        gv[0] = cos(xlatg / radeg) * cos(*(lon + i) / radeg);
        gv[1] = cos(xlatg / radeg) * sin(*(lon + i) / radeg);
        gv[2] = sin(xlatg / radeg);

        r = re * (1. - f) /
                sqrt(1. - (2. - f) * f * pow(cos(xlatg / radeg), 2));

        /* note that pos is in m but rest of constants are in km, hence the 
           1000 factor  */
        for (j = 0; j < 3; j++)
            rh[j] = pos[j] / 1000. - r * gv[j];
        /*
         *  Transform the pixel-to-spacecraft vector into the local frame
         */
        matrix_mult(rh, xlmat, rl);
        /*
         *  Compute the sensor zenith and azimuth
         *  if zenith close to 0, set azimuth to 0 (and normalize azimuth)
         */
        *(senz + i) = radeg * atan2(sqrt(rl[0] * rl[0] + rl[1] * rl[1]),
                rl[2]);
        if (*(senz + i) > 0.05)
            *(sena + i) = radeg * atan2(rl[0], rl[1]);
        else
            *(sena + i) = 0.;

        if (*(sena + i) < 0.) *(sena + i) += 360.;
    }
    /*
     *  and end
     */
    return 0;
}

void matrix_mult(double vecin[3], double matrix[3][3], double vecout[3])
/*******************************************************************

   matrix_mult

   purpose: multiply a vector by a matrix and produce the output vector

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      double[3]         vecin            I      input vector
      double[3][3]      matrix           I      rotation matrix
      double[3]         vecout           O      rotated vector

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC  2-Nov-2005     Original development

 *******************************************************************/ {
    int i, j;
    for (i = 0; i < 3; i++) {
        vecout[i] = 0.;
        for (j = 0; j < 3; j++) {
            vecout[i] += matrix[i][j] * vecin[j];
        }
    }
}

void cross_prod(double *v1, double *v2, double *vout)
/*******************************************************************

   cross_prod

   purpose: produce the cross product of 2 vectors

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      double *          v1               I      vector 1
      double *          v2               I      vector 1
      double *          vout             O      v1 X v2

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC  4-Jan-2006     Original development

 *******************************************************************/ {
    *vout = v1[1] * v2[2] - v1[2] * v2[1];
    *(vout + 1) = v1[2] * v2[0] - v1[0] * v2[2];
    *(vout + 2) = v1[0] * v2[1] - v1[1] * v2[0];
}
