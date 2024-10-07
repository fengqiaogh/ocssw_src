#include <netcdf>
#include <iostream>
#include "l1_seawifs_netcdf.h"
#include "l1.h"
#include "l1a.h"
#include "eng_qual.h"
#include "cal_l1a.h"
#include "call1a_proto.h"
#include "getcal_proto.h"
#include "l1a_proto.h"
#include "navigation.h"
#include "st_proto.h"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


#define   LAC_PIXEL_NUM           1285
#define   GAC_PIXEL_NUM           248
#define   NREC_IN_BUF             10
#define   STBUFSIZ                5
#define   NOTDONE                 0
#define   FIRST_KNEE              1
#define   MASK_HIGHLT1            16
#define   GENBUFSIZ               NREC_IN_BUF*sizeof(float)*40 /* size of inst_ana */

static int16_t syear, sday; /* data start date                  */
static int32_t smsec; /* data start time                  */
static int16_t eyear, eday; /* data end date                    */
static int32_t emsec; /* data end time                    */
static int32_t nscan; /* number of scans                  */
static int32_t npix; /* number pixels per scan           */
static int32_t spix; /* start pixel number (from 0)      */
static int32_t dpix; /* LAC pixel increment              */
static float dark_mean[8]; /* Mean of dark restore counts      */
static float dark_std[8]; /* Std dev  of dark restore counts  */

static char cal_path_tab[128];
static char dtype[8];

static float pcal_counts[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
static float pcal_rads[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
static short l2_flags_buffer[LAC_PIXEL_NUM];

static int32_t nsta;
static int32_t ninc;

static int16_t *l1a_data = NULL; /* raw radiances band-int-by-pix    */
static int16_t *l1a_back = NULL; /* raw radiances band-int-by-pix    */
static float *l1b_buffer = NULL; /* l1b radiances band-int-by-pix    */
static int32_t *msec;
static int16_t *side;
static int16_t *dark_rest;
static float *tilt;

static float ylat[LAC_PIXEL_NUM];
static float xlon[LAC_PIXEL_NUM];
static float solz[LAC_PIXEL_NUM];
static float sola[LAC_PIXEL_NUM];
static float senz[LAC_PIXEL_NUM];
static float sena[LAC_PIXEL_NUM];

static int16_t stray_light = -1; /* num pixels to apply sl correct.  */
static float Ltyp_frac = 0.25; /* fraction of B8 for sl threshold  */

typedef struct {
    int16_t gain[8];
    int16_t tdi[8];
    int16_t scan_temp[8];
    float inst_temp[40];
    float orb_vec[3];
    float sun_ref[3];
    float sen_mat[3 * 3];
    float scan_ell[6];
} inputBuffer;


static int32_t do_st = 1;


/**
 * Given a variable, get the data from the netcdf file and save it to genBuf
 * then copy what is in genBuf into the rdBuf
 * @param rdBuf - rdBuf referenced in get_l1a_rec
 * @param genBuf - arr referenced in get_l1a_rec
 * @param dataFile - netcdf file reference that is currently open
 * @param variable - what to read from the netcdf file
 * @param start<n> - start index of the netcdf
 * @param count<n> - defines how many values should be read relative to start
 * @param numBytes - when copying, how many bytes of data is being copied 
 * @param n_read
 * 
*/
extern "C" void getVariableAndCopyToRdBuffer(
    inputBuffer rdBuf[],
    byte genBuf[],
    NcFile &dataFile,
    string variable,
    size_t start0,
    size_t start1,
    size_t start2,
    size_t count0,
    size_t count1,
    size_t count2,
    int numBytes,
    int n_read
) {
    vector<size_t> start = {start0, start1, start2};
    vector<size_t> count = {count0, count1, count2};
    dataFile.getVar(variable).getVar(start, count, (void*)genBuf);
    for (int i = 0; i < n_read; i++)
        if (variable == "gain")
            memcpy(&rdBuf[i].gain, &genBuf[i * numBytes], numBytes);
        else if (variable == "tdi")
            memcpy(&rdBuf[i].tdi, &genBuf[i * numBytes], numBytes);
        else if (variable == "scan_temp")
            memcpy(&rdBuf[i].scan_temp, &genBuf[i * numBytes], numBytes);
        else if (variable == "inst_ana")
            memcpy(&rdBuf[i].inst_temp, &genBuf[i * numBytes], numBytes);
        else if (variable == "orb_vec")
            memcpy(&rdBuf[i].orb_vec, &genBuf[i * numBytes], numBytes);
        else if (variable == "sun_ref")
            memcpy(&rdBuf[i].sun_ref, &genBuf[i * numBytes], numBytes);
        else if (variable == "scan_ell")
            memcpy(&rdBuf[i].scan_ell, &genBuf[i * numBytes], numBytes);
        else if (variable == "sen_mat")
            memcpy(&rdBuf[i].sen_mat, &genBuf[i * numBytes], numBytes);
}   

extern "C" int32_t get_l1a_rec_netcdf(filehandle *file, int32_t recno, cal_mod_struc *cal_mod,
        int16_t *l1a_dum, float **l1b_data, int16_t **l2_flags) {
    /*                                                                 */
    /* local vars                                                      */
    /*                                                                 */
    int32_t ipix; /* pixel number                     */
    int32_t idet; /* detector number                  */
    int32_t irec;
    int16_t band_buf[LAC_PIXEL_NUM * 8];
    int32_t scan_no;
    int16_t gain8;
    int32_t AS_pixels;
    float Styp_frac = 0.9;
    int32_t sl_scan;
    int16_t *l1a_ptr;
    int16_t *gain_ptr;

    static inputBuffer rdBuf[NREC_IN_BUF];
    static inputBuffer bkBuf[STBUFSIZ - 2];
    byte genBuf[GENBUFSIZ];

    static int32_t crec = -1;
    static int16_t recursive_flag = 0;
    static int16_t n_read = 0;
    static int32_t max_rec_in_rdbuf;
    static int32_t offset;
    static int32_t j;
    static int32_t stray_light_scan_no;

    static float *st_l1b_data;
    static int16_t *st_l2_flags;

    static int32_t initial = 1;
    static byte first = 1;


    static int16_t hi_Lt[LAC_PIXEL_NUM]; /* whether radiance > knee */
    int pixvalue; /* l1a pixel value			*/
    int pixdark; /* dark_restore value for each pixel	*/
    float kneevalue; /* radiance at first knee		*/

   // Used to read variable data from a netcdf file
    vector<size_t> start = {0, 0, 0}; // Starting position index. 
    vector<size_t> count = {0, 0, 0}; // slice/poirtion of the data from start to read


    /*                                                                 */
    /* Read L1A data scan by scan, build L1B record, and write.        */

    /* If reading out-of-sequence, force restart from beginning to get stlight */
    /* baf, 05-oct-2001.                                                       */
    if (first || recno < crec) {
        max_rec_in_rdbuf = recno - 1;
        stray_light_scan_no = recno;
    }

    crec = recno;

    if (crec > nscan) return 0;

    if (crec > max_rec_in_rdbuf) {

        memcpy(&bkBuf, &rdBuf[NREC_IN_BUF - STBUFSIZ + 2],
                sizeof (inputBuffer) * (STBUFSIZ - 2));
        memcpy(l1a_back, &l1a_data[(NREC_IN_BUF - STBUFSIZ + 2) * npix * 8],
                npix * 8 * (STBUFSIZ - 2) * sizeof (int16_t));

        n_read = (nscan - max_rec_in_rdbuf >= NREC_IN_BUF) ?
                NREC_IN_BUF : nscan - max_rec_in_rdbuf;


        try {
            NcFile dataFile(file->name, NcFile::read);

            // starting position of the array to read. 3rd index only used for 3D arrays
            start[0] = 0;
            start[1] = 0;
            start[2] = 0;
            // how large the data to read from start. 3rd index only used for 3D arrays
            count[0] = 1;
            count[1] = (size_t) npix;
            count[2] = 8;

            for (irec = 0; irec < n_read; irec++) {
                start[0]= (size_t) (recno - 1 + irec);
                dataFile.getVar("l1a_data").getVar(start, count, (void*) band_buf);
                for (ipix = 0; ipix < npix; ipix++) {
                    for (idet = 0; idet < 8; idet++) {
                        l1a_data[irec * (8 * npix) + idet * npix + ipix] =
                            band_buf[ipix * 8 + idet];
                    }
                }
            }


            // Get variables and copy them to the reading buffer for calculation

            getVariableAndCopyToRdBuffer(
                rdBuf,      // reference the buffer and datafiles
                genBuf,
                dataFile,
                "gain",
                recno-1, 0, 0, // start
                n_read, 8, 0, // count
                8*2, // numBytes 8 columns, int_16 = 2 bytes 
                n_read
            );

            getVariableAndCopyToRdBuffer(
                rdBuf,
                genBuf,
                dataFile,
                "tdi",
                recno-1, 0, 0, // start
                n_read, 8, 0, // count
                8*2, // numBytes 8 columns, int_16 = 2 bytes 
                n_read
            );

            getVariableAndCopyToRdBuffer(
                rdBuf,
                genBuf,
                dataFile,
                "scan_temp",
                recno-1, 0, 0, // start
                n_read, 8, 0, // count
                8*2, // numBytes 8 columns, int_16 = 2 bytes 
                n_read
            );
            
            getVariableAndCopyToRdBuffer(
                rdBuf,
                genBuf,
                dataFile,
                "inst_ana",
                recno-1, 0, 0, // start
                n_read, 8, 0, // count
                40*4, // numBytes 40 columns, float = 4 bytes
                n_read
            );

            getVariableAndCopyToRdBuffer(
                rdBuf,
                genBuf,
                dataFile,
                "orb_vec",
                recno-1, 0, 0, // start
                n_read, 3, 0, // count
                3*4, // numBytes 3 columns, float = 4 bytes 
                n_read
            );

            getVariableAndCopyToRdBuffer(
                rdBuf,
                genBuf,
                dataFile,
                "sun_ref",
                recno-1, 0, 0, // start
                n_read, 3, 0, // count
                3*4, // numBytes 3 columns, float = 4 bytes 
                n_read
            );
            
            getVariableAndCopyToRdBuffer(
                rdBuf,
                genBuf,
                dataFile,
                "scan_ell",
                recno-1, 0, 0, // start
                n_read, 6, 0, // count
                6*4, // numBytes 6 columns, float = 4 bytes 
                n_read
            );

            getVariableAndCopyToRdBuffer(
                rdBuf,
                genBuf,
                dataFile,
                "sen_mat",
                recno-1, 0, 0, // start
                n_read, 3, 3, // count
                9*4, // numBytes 9 columns, float = 4 bytes 
                n_read
            );
            
            offset = max_rec_in_rdbuf + 1;
            max_rec_in_rdbuf += n_read;

        }
        catch (NcException& e) {
            cout << "-E- Error getting L1A Data in Seawifs l2gen reader: " << e.what() << endl;
            exit(1);
        }
    }

    if ((crec - offset) >= 0) {
        j = crec - offset;
        gain8 = rdBuf[j].gain[7];
        gain_ptr = rdBuf[j].gain;
        l1a_ptr = &l1a_data[j * 8 * npix];

        l1b_rad(syear, sday, smsec, msec[recno - 1],
                dtype, nsta, ninc, npix,
                dark_mean, rdBuf[j].gain, rdBuf[j].tdi,
                rdBuf[j].scan_temp, rdBuf[j].inst_temp, side[recno - 1],
                &l1a_data[j * 8 * npix], l1b_buffer, cal_mod);

        if (!recursive_flag) {
            geonav_(rdBuf[j].orb_vec, rdBuf[j].sen_mat, rdBuf[j].scan_ell,
                    rdBuf[j].sun_ref, &nsta, &ninc, &npix,
                    ylat, xlon, solz, sola, senz, sena);

        }
    } else {
        j = (STBUFSIZ - 2) + (crec - offset);
        gain8 = bkBuf[j].gain[7];
        gain_ptr = bkBuf[j].gain;
        l1a_ptr = &l1a_back[j * 8 * npix];

        if (!recursive_flag) {
            geonav_(bkBuf[j].orb_vec, bkBuf[j].sen_mat, bkBuf[j].scan_ell,
                    bkBuf[j].sun_ref, &nsta, &ninc, &npix,
                    ylat, xlon, solz, sola, senz, sena);
        }
    }
    if (first) {
        memcpy(pcal_counts, cal_counts, sizeof (cal_counts));
        memcpy(pcal_rads, cal_rads, sizeof (cal_rads));
        first = 0;
    }

    if (!recursive_flag) {
        for (ipix = 0; ipix < npix; ipix++) {
            hi_Lt[ipix] = FALSE;
            /* HILT check limited to bands 7 & 8, BAF, 23 Jan 2002 */
            /* Use dark_mean rather than dark_rest, BAF, 7 July 2003 */
            for (idet = 6; (idet < 8)&&(hi_Lt[ipix] == FALSE); idet++) {
                pixvalue = l1a_ptr[ipix + npix * idet];
                /*
                pixdark = dark_rest[8*(recno-1)+idet];
                 */
                pixdark = dark_mean[idet];
                kneevalue = pcal_counts[idet][gain_ptr[idet]][FIRST_KNEE];
                hi_Lt[ipix] = ((pixvalue - pixdark) >= kneevalue);
            }
        }
    }


    scan_no = crec;
    AS_pixels = stray_light;

    if (do_st == 0) {

        for (ipix = 0; ipix < npix; ipix++)
            l2_flags_buffer[ipix] = 0;

    } else if (!recursive_flag) {

        do {
            if (stray_light_scan_no != scan_no) {

                recursive_flag = 1;

                get_l1a_rec_netcdf(file, stray_light_scan_no, cal_mod,
                        l1a_dum, &st_l1b_data, &st_l2_flags);
            }
        } while (stray_light_corr(&initial, Ltyp_frac, Styp_frac,
                nscan, npix, stray_light_scan_no++, dtype,
                gain8, &pcal_rads[0][0][0], l1b_buffer, &sl_scan,
                l2_flags_buffer, &AS_pixels) == NOTDONE);
        recursive_flag = 0;
        crec = scan_no;

    }

    for (ipix = 0; ipix < npix; ipix++) {
        if (hi_Lt[ipix])
            l2_flags_buffer[ipix] |= MASK_HIGHLT1;
    }

    *l1b_data = l1b_buffer;
    *l2_flags = l2_flags_buffer;

    return (0);
}

extern "C" int openl1_seawifs_netcdf(filehandle *file) {

    /*                                                                 */
    /* input files                                                     */
    /*                                                                 */
    char *cal_path = l1_input->calfile;

    /*                                                                 */
    /* get_l1a_open interface                                          */
    /*                                                                 */
    char buf[32];

    // Reading some Global Attributes and Dimentions
    try {
        NcFile dataFile(file->name, NcFile::read);
       
        /*
            NetCDF contains only the isodate for the start and end time.
            The temp variable will store the isodate and then call the conversion function
            that will convert and set the start/end year, day and msec 
        */
        char tempTimeCoverage[27];

        // Start Time
        dataFile.getAtt("time_coverage_start").getValues((void*) &tempTimeCoverage);

        // temp to hold the value and then reassign it to the 16 bit int
        int32_t tempDay, tempYear;
        isodate2ydmsec(tempTimeCoverage, &tempYear, &tempDay, &smsec);
        sday = tempDay;
        syear = tempYear;

        // End Time
        dataFile.getAtt("time_coverage_end").getValues((void*) &tempTimeCoverage);
        isodate2ydmsec(tempTimeCoverage, &tempYear, &tempDay, &emsec);
        eday = tempDay;
        eyear = tempYear;


        npix = dataFile.getDim("pixels").getSize();
        nscan = dataFile.getDim("scans").getSize();

        dataFile.getAtt("data_type").getValues((void*) &dtype);
        dataFile.getAtt("LAC_pixel_start_number").getValues((void*) &nsta);
        dataFile.getAtt("LAC_pixel_subsampling").getValues((void*) &ninc);
        dataFile.getAtt("orbit_node_longitude").getValues((void*) &file->orbit_node_lon);
        dataFile.getAtt("orbit_number").getValues((void*) &file->orbit_number);
        dataFile.getAtt("node_crossing_time").getValues((void*) &buf);
        file->node_crossing_time = isodate2unix(buf);


        // Allocate space for data and read them in
        msec = (int32_t *) calloc(nscan, sizeof (int32_t));
        side = (int16_t *) calloc(nscan, sizeof (int16_t));
        dark_rest = (int16_t *) calloc(nscan * 8, sizeof (int16_t));
        tilt = (float *) calloc(nscan, sizeof (float));

        dataFile.getVar("scan_time").getVar(msec);
        dataFile.getVar("side").getVar(side);
        dataFile.getVar("dark_rest").getVar(dark_rest);
        dataFile.getVar("tilt").getVar(tilt);
    }
    catch (NcException& e) {
        cout << "-E- Error in reading in openl1a_seawifs_netcdf. NetCDF: " << e.what() << endl;
        exit(1);
    }

    Ltyp_frac = l1_input->sl_frac;
    stray_light = l1_input->sl_pixl;


    if (stray_light == 0)
        do_st = 0;
    else {
        do_st = 1;
        if (stray_light < 0) {
            if (strcmp(dtype, "GAC") == 0)
                stray_light = 4;
            else
                stray_light = 3;
        }
    }
    /* get the current calibration model solution */
    if (cal_path == NULL) {
        fprintf(stderr,
                "-E %s Line %d: No calibration file specified.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    read_caltable(cal_path);
    cdata_();

    file->npix = npix;
    file->nscan = nscan;
    file->sensorID = SEAWIFS;

    if (strcmp(dtype, "GAC") == 0)
        strcpy(file->spatialResolution, "4.5 km");
    else
        strcpy(file->spatialResolution, "1.1 km");

    strcpy(cal_path_tab, cal_path);

    l1a_data = (int16_t *) calloc(npix * 8 * NREC_IN_BUF, sizeof (int16_t));
    l1a_back = (int16_t *) calloc(npix * 8 * (STBUFSIZ - 2), sizeof (int16_t));
    l1b_buffer = (float *) calloc(npix * 8, sizeof (float));
   
    dark_rest_stat(dark_rest, nscan, dark_mean, dark_std);
    free(dark_rest);
    spix = nsta - 1;
    dpix = ninc;

    return (0);
}

extern "C" int readl1_seawifs_netcdf(filehandle *file, int32_t recnum, l1str *l1rec) {

    /*                                                                 */
    /* get_l1a_rec interface                                           */
    /*                                                                 */
    static cal_mod_struc cal_mod; /* cal modification structure       */
    static int16_t *l2_flags; /* radiance quality flags for L2    */
    static float *l1b_data;

    int32_t i;

    int32_t nwave = l1rec->l1file->nbands;
    int32_t *bindx = (int32_t*) l1rec->l1file->bindx;

    /*                                                                 */
    /* local vars                                                      */
    /*                                                                 */
    int32_t ipix; /* pixel number                     */
    int32_t iw, ib;

    int16_t *l1a_dum = NULL;

    static int32_t prev_recnum = -1;

    static int16_t first = 1;
    static int32_t ntilts;
    static int16_t tilt_flags[20];
    static int16_t tilt_ranges[2 * 20];

    int32_t nflag[8];


    /* If reading out-of-sequence, force restart from beginning to get stlight */
    /* baf, 05-oct-2001.                                                       */
    if (recnum < prev_recnum) {
        printf("Reading out-of-sequence %d %d\n", recnum, prev_recnum);
        prev_recnum = -1;
    }


    /* Get tilt info */
    /* ------------- */
    try {
        NcFile dataFile(file->name, NcFile::read);
        
        if (first) {
            ntilts = dataFile.getDim("tilts").getSize();
            dataFile.getVar("tilt_flags").getVar(tilt_flags);
            dataFile.getVar("tilt_ranges").getVar(tilt_ranges);
            first = 0;
        }

        /* Check for bad or changing tilt */
        /* ------------------------------ */
        for (i = 0; i < ntilts; i++) {
            if (tilt_ranges[2 * i] <= recnum + 1 && tilt_ranges[2 * i + 1] >= recnum)
                break;
        }

        // read navigation flags for each scan
        vector<size_t> start = {(size_t) recnum, 0};
        vector<size_t> count = {1, 8};
        dataFile.getVar("nflag").getVar(start, count, (void*)&nflag);
    }
    catch (NcException& e) {
        cout << "-E- Error in reading in readl1a_seawifs_netcdf. NetCDF: " << e.what() << endl;
        exit(1);
    }


    // /* Get l1a data */
    // /* ------------ */
    for (i = prev_recnum + 1; i <= recnum; i++)
        get_l1a_rec_netcdf(file, i + 1, &cal_mod, l1a_dum,
            &l1b_data, &l2_flags);

    // /*                                                              */
    // /* Copy scan geolocation and view geometry                      */
    // /*                                                              */
    memcpy(l1rec->lat, ylat, npix * sizeof (float));
    memcpy(l1rec->lon, xlon, npix * sizeof (float));
    memcpy(l1rec->solz, solz, npix * sizeof (float));
    memcpy(l1rec->sola, sola, npix * sizeof (float));
    memcpy(l1rec->senz, senz, npix * sizeof (float));
    memcpy(l1rec->sena, sena, npix * sizeof (float));

    // /*                                                              */
    // /* Copy L1B radiances, pixel interlaced by band. Add per-band   */
    // /* view angles.                                                 */
    // /*                                                              */
    for (ipix = 0; ipix < file->npix; ipix++) {
        if (l1rec->sena[ipix] > 180) {
            l1rec->sena[ipix] -= 360.0;
        }
        if (l1rec->sola[ipix] > 180) {
            l1rec->sola[ipix] -= 360.0;
        }
        l1rec->pixnum[ipix] = spix + ipix*dpix;
        for (iw = 0; iw < nwave; iw++) {
            ib = bindx[iw];
            l1rec->Lt [ipix * nwave + ib] = l1b_data[iw * npix + ipix];
        }
        l1rec->stlight[ipix] = ((l2_flags[ipix] & STRAYLIGHT) > 0);
        l1rec->hilt [ipix] = ((l2_flags[ipix] & HILT) > 0);
        l1rec->navwarn[ipix] = (nflag[7] & 1) | (nflag[0] & 1);
    }


    // /*                                                              */
    // /* Set scan time in L1B output record                           */
    // /*                                                              */
    double secs = (double) (msec[recnum] / 1.e3);
    int16_t year = syear;
    int16_t day = sday;

    if (msec[recnum] < smsec) { /* adjust for day rollover */
        day += 1;
        if (day > (365 + (year % 4 == 0))) {
            year += 1;
            day = 1;
        }
    }
    l1rec->scantime = yds2unix(year, day, secs);

    l1rec->tilt = tilt[recnum];
    l1rec->mside = side[recnum];
    prev_recnum = recnum;

    return (0); // no issues
}

extern "C" int readl1_lonlat_seawifs_netcdf(filehandle *file, int32_t recnum, l1str *l1rec) {
    inputBuffer rdBuf;

    if (recnum > nscan) return 0;

    try {   
        
        NcFile dataFile(file->name, NcFile::read);

        vector<size_t> start = {(size_t) recnum, 0};
        vector<size_t> count = {0, 0};


        count[0] = 1;
        count[1] = 3;
        dataFile.getVar("orb_vec").getVar(start, count, (void*)rdBuf.orb_vec);
        dataFile.getVar("sun_ref").getVar(start, count, (void*)rdBuf.sun_ref);

        count[1] = 6;
        dataFile.getVar("scan_ell").getVar(start, count, (void*)rdBuf.scan_ell);

        start = {(size_t)recnum, 0, 0};
        count = {1, 3, 3};
        dataFile.getVar("sen_mat").getVar(start, count, (void*)rdBuf.sen_mat);

    }
    catch (NcException& e) {
        cout << "-E- Error in reading seawifs lon and lat data netcdf." << endl;
        exit(1);
    }
    
    // Lon and lat
    geonav_(rdBuf.orb_vec, rdBuf.sen_mat, rdBuf.scan_ell,
            rdBuf.sun_ref, &nsta, &ninc, &npix,
            l1rec->lat, l1rec->lon, l1rec->solz, l1rec->sola, l1rec->senz, l1rec->sena);

    // Time
    double secs = (double) (msec[recnum] / 1.e3);
    int16_t year = syear;
    int16_t day = sday;

    if (msec[recnum] < smsec) { /* adjust for day rollover */
        day += 1;
        if (day > (365 + (year % 4 == 0))) {
            year += 1;
            day = 1;
        }
    }
    l1rec->scantime = yds2unix(year, day, secs);

    return (0);
}

extern "C" int closel1_seawifs_netcdf(filehandle *file) {
    free(l1a_data);
    free(l1a_back);
    free(l1b_buffer);
    free(msec);
    free(side);
    free(tilt);

    return (0);
}