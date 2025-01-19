/*
 *  W. Robinson, SAIC, 10 Dec 2004  I/O routines for CZCS
 */
#include <netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#include <hdf4utils.h>
#include "mfhdf.h"
#include "l1_czcs_netcdf.h"
#include <math.h>
#include <libnav.h>

extern "C" {
    #include "l1_czcs.h"
}

#include <iostream>

static int32_t syear, sday; /* data start date                  */
static int32_t smsec; /* data start time                  */
static int32_t eyear, eday; /* data end date                    */
static int32_t emsec; /* data end time                    */
static int32_t nscan; /* number of scans                  */
static int32_t npix; /* number pixels per scan           */

static int32_t nsta;
static int32_t ninc;

static uint8_t *counts, cz_band_present;
static int32_t *msec;
static int16_t *gain;
static float *tilt, *att_ang, *slope, *intercept;
static float *ctl_pt_lat, *ctl_pt_lon, *pos, *pos_err;
static int32_t nctl_pt, *ctl_pt_cols;
static float *ctl_pt_vx, *ctl_pt_vy, *ctl_pt_vz, *y2_vx, *y2_vy, *y2_vz, *ctl_pt_x;
static float *lt750; /* internal 750 mn data source */
static char *ring_sat; /* set to 1 if 443, 520 or 550 are saturated for ringing 
                    mask computation */


#define NBND 4
#define NGAIN 4
#define NEPOCH 5

/*
   W. Robinson, SAIC, 6 Jan 2006  add code to read the position and 
   position error SDSes
 */

extern "C" int32_t openl1_czcs_netcdf(filehandle *file) {

    /*                                                                 */
    /* get_l1a_open interface                                          */
    /*                                                                 */
    int i;
    // Reading some Global Attributes and Dimentions
    try {
        NcFile dataFile(file->name, NcFile::read);

        char tempTimeCoverage[27];
        dataFile.getAtt("time_coverage_start").getValues((void*) &tempTimeCoverage);
        isodate2ydmsec(tempTimeCoverage, &syear, &sday, &smsec);

        dataFile.getAtt("time_coverage_end").getValues((void*) &tempTimeCoverage);
        isodate2ydmsec(tempTimeCoverage, &eyear, &eday, &emsec);

        npix = dataFile.getDim("pixels").getSize();
        nscan = dataFile.getDim("scans").getSize();
        nctl_pt = dataFile.getDim("control_points").getSize();
        dataFile.getAtt("orbit_number").getValues((void*) &file->orbit_number);
        dataFile.getAtt("parameter_presence_code").getValues((void*) &cz_band_present);
        dataFile.getAtt("LAC_pixel_start_number").getValues((void*) &nsta);
        dataFile.getAtt("LAC_pixel_subsampling").getValues((void*) &ninc);

        /* call cdata.f to initialize global FORTRAN common block data	*/
        cdata_();


        file->npix = npix;
        file->nscan = nscan;
        file->sensorID = CZCS;

        strcpy(file->spatialResolution, "825 m");
        // Allocate space for data and read them in
        msec = (int32_t *) calloc(nscan, sizeof (int32_t));
        dataFile.getVar("scan_time").getVar(msec);
        tilt = (float *) calloc(nscan, sizeof (float));
        dataFile.getVar("tilt").getVar(tilt);
        att_ang = (float *) calloc(nscan * 3, sizeof (float));
        dataFile.getVar("att_ang").getVar(att_ang);

        ctl_pt_cols = (int32_t *) calloc(nctl_pt, sizeof (int32_t));
        dataFile.getVar("cntl_pt_cols").getVar(ctl_pt_cols);

        ctl_pt_x = (float *) calloc(nctl_pt, sizeof ( float));
        for (i = 0; i < nctl_pt; i++)
            ctl_pt_x[i] = (float) ctl_pt_cols[i] - 1.;

        slope = (float *) calloc(nscan * 6, sizeof (float));
        dataFile.getVar("slope").getVar(slope);

        intercept = (float *) calloc(nscan * 6, sizeof (float));
        dataFile.getVar("intercept").getVar(intercept);

        /* get nav orbit data  */
        pos = (float *) malloc(nscan * 3 * sizeof (float));
        dataFile.getVar("orb_vec").getVar(pos);

        pos_err = (float *) malloc(nscan * sizeof (float));
        dataFile.getVar("pos_err").getVar(pos_err);

        gain = (int16_t *) malloc(nscan * sizeof ( int16_t));
        counts = (uint8_t *) malloc(npix * NBND_CZCS * sizeof ( uint8_t));
        ctl_pt_lat = (float *) malloc(nctl_pt * sizeof ( float));
        ctl_pt_lon = (float *) malloc(nctl_pt * sizeof ( float));
        ctl_pt_vx = (float *) malloc(nctl_pt * sizeof ( float));
        ctl_pt_vy = (float *) malloc(nctl_pt * sizeof ( float));
        ctl_pt_vz = (float *) malloc(nctl_pt * sizeof ( float));
        y2_vx = (float *) malloc(npix * sizeof ( float));
        y2_vy = (float *) malloc(npix * sizeof ( float));
        y2_vz = (float *) malloc(npix * sizeof ( float));
        lt750 = (float *) malloc(npix * sizeof ( float));
        ring_sat = (char *) malloc(npix * sizeof ( char));
    }
    catch (NcException& e) {
        cout << "-E- Error in reading input NetCDF: " << e.what() << endl;
        exit(1);
    }
    return (0);
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

extern "C" int32_t readl1_czcs_netcdf(filehandle *file, int32_t recnum, l1str *l1rec) {
    void satang_(double *, double *, float *, float *, float *, float *,
            float *, float *, float *, float *);
    int ll2vec(float *, float *);
    int vec2ll(float *, float *);
    short cnt_vec[NBND_CZCS];
    float lt_lcl[NBND_CZCS], yout1, yout2, yout3, llvec[2], vec[3], gmt;
    int ipx, ibnd, orbit, i, tpix;
    uint8_t cal_sum[5], cal_scan[6];
    int32_t status, navbad, ltsat;
    double pi, radeg;
    int32_t ib;

    int32_t nwave = l1rec->l1file->nbands;
    int32_t *bindx = l1rec->l1file->bindx;
    char *cal_path = l1_input->calfile;

    float lonbuf[1968], latbuf[1968], senzbuf[1968], senabuf[1968];

    try {
        NcFile dataFile(file->name, NcFile::read);

        /*
        *  read in the line of counts, latitude, longitude
        */

        NcVar l1aDataVar = dataFile.getVar("band1");

        vector<size_t> start(l1aDataVar.getDimCount());
        vector<size_t> edges(l1aDataVar.getDimCount());
        start[0] = recnum;
        start[1] = 0;
        edges[0] = 1;
        edges[1] = npix;

        l1aDataVar.getVar(start, edges, (VOIDP) counts);
        
        l1aDataVar = dataFile.getVar("band2");
        l1aDataVar.getVar(start, edges, (VOIDP) (counts + npix));
        l1aDataVar = dataFile.getVar("band3");
        l1aDataVar.getVar(start, edges, (VOIDP) (counts + 2 * npix));
        l1aDataVar = dataFile.getVar("band4");
        l1aDataVar.getVar(start, edges, (VOIDP) (counts + 3 * npix));
        l1aDataVar = dataFile.getVar("band5");
        l1aDataVar.getVar(start, edges, (VOIDP) (counts + 4 * npix));

        NcVar gainSDS = dataFile.getVar("gain");
        gainSDS.getVar(start, edges, gain);

        edges[1] = nctl_pt;
        NcVar latSDS = dataFile.getVar("latitude");
        latSDS.getVar(start, edges, ctl_pt_lat);       
        NcVar lonSDS = dataFile.getVar("longitude");
        lonSDS.getVar(start, edges, ctl_pt_lon);  

        edges[1] = 5;
        NcVar cal_sumSDS = dataFile.getVar("cal_sum");
        cal_sumSDS.getVar(start, edges, cal_sum);  

        edges[1] = 6;
        NcVar cal_scanSDS = dataFile.getVar("cal_scan");
        cal_scanSDS.getVar(start, edges, cal_scan); 
    }
    catch (NcException& e) {
        cout << "-E- Error in reading in readl1a_czcs_netcdf. NetCDF: " << e.what() << endl;
        exit(1);
    }
    /*
     * flag setting for entire line: set navfail if cal_sum shows problems 
     * (suspect it is * bad ephemeris or attitude) and set hilt if bands 
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
        status = get_czcscal(cal_path, orbit, syear, sday, msec[recnum],
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
        pi = PI;
        radeg = RADEG;
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


    /*  for( i = spix; i < epix; i = i + ninc )*/
    for (i = 0; i < npix; i++) {
        gmt = (float) msec[recnum] / (1000. * 3600);
        sunangs_(&syear, &sday, &gmt, (l1rec->lon) + i, (l1rec->lat) + i,
                (l1rec->solz) + i, (l1rec->sola) + i);
    }
    /*
     *  set scan time
     */
    double secs = (double) (msec[recnum] / 1000.0);
    int16_t year = syear;
    int16_t day = sday;

    if (recnum > 0 && msec[recnum] < msec[recnum - 1]) { /* adjust for day rollover */
        day += 1;
        if (day > (365 + (year % 4 == 0))) {
            year += 1;
            day = 1;
        }
    }
    l1rec->scantime = yds2unix(year, day, secs);

    return (status);
}

/*  W. Robinson, SAIC, 6 Jan 2005   include the pos, pos_err arrays in free */

extern "C" int32_t  closel1_czcs_netcdf(filehandle *file) {
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