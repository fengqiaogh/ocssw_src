#include <hdf4utils.h>
#include "l1a.h"
#include "eng_qual.h"
#include "l1.h"
#include "l1_seawifs.h"
#include "call1a_proto_osmi.h"
#include "getcal_proto_osmi.h"
#include "orbit.h"
#include "lunsol.h"

#include <hdf.h>
#include <mfhdf.h>

#define   LAC_PIXEL_NUM           1285
#define   NREC_IN_BUF             10
#define   NBANDS_OSMI             8

static int16_t syear, sday; /* data start date                  */
static int32_t smsec; /* data start time                  */
static int16_t eyear, eday; /* data end date                    */
static int32_t emsec; /* data end time                    */
static int32_t nscan; /* number of scans                  */
static int32_t npix; /* number pixels per scan           */

static int16_t cal_year;
static int16_t cal_day;
static int32_t cal_msec;

static float eoffset;
static float egain[8], temps[8], mirror[8];
static float slope[BANDS_DIMS_1A * 96];
static float dc[BANDS_DIMS_1A * 96];
static float sm[1044];

static float t_const[BANDS_DIMS_1A];
static float t_linear[BANDS_DIMS_1A];
static float t_quadratic[BANDS_DIMS_1A];

static char cal_path_tab[128];

static int32_t nsta;
static int32_t ninc;

static int16_t *l1a_data = NULL; /* raw radiances band-int-by-pix    */

static float *l1b_buffer = NULL; /* l1b radiances band-int-by-pix    */
static int32_t *msec;
static float *sc_ana;
static double *sc_ttag;
static double *sc_ttag_inc;

static char dtype[8];
static char cal_path_tab[128];

static short l2_flags_buffer[LAC_PIXEL_NUM];

static float ylat[LAC_PIXEL_NUM];
static float xlon[LAC_PIXEL_NUM];
static float solz[LAC_PIXEL_NUM];
static float sola[LAC_PIXEL_NUM];
static float senz[LAC_PIXEL_NUM];
static float sena[LAC_PIXEL_NUM];

int32_t get_l1a_rec_osmi(int32_t, int32_t, cal_mod_struc *, float **, int16_t **);

void intpos_(double *epoch1, double pos1[3], double vel1[3],
        double *epoch2, double pos2[3], double vel2[3],
        double *epoch, double pos [3], double vel [3]);

int openl1_osmi(filehandle *file) {

    /*                                                                 */
    /* input files                                                     */
    /*                                                                 */
    char *cal_path = l1_input->calfile;

    /*                                                                 */
    /* get_l1a_open interface                                          */
    /*                                                                 */
    int32_t fileID;
    int32_t status;

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


    /* get the current calibration model solution */
    if (cal_path == NULL) {
        fprintf(stderr,
                "-E %s Line %d: No calibration file specified.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    status = get_cal_osmi(cal_path, syear, sday, eday, smsec,
            &cal_year, &cal_day, &cal_msec,
            &eoffset, &egain[0], &temps[0],
            &mirror[0], &t_const[0], &t_linear[0],
            &t_quadratic[0], &slope[0], &dc[0],
            &sm[0]);
    if (status < 0) {
        fprintf(stderr,
                "-E- %s line %d: Error applying calibration table \"%s\".\n",
                __FILE__, __LINE__, cal_path);
        exit(status);
    }


    file->npix = npix;
    file->nscan = nscan;
    file->sensorID = OSMI;
    file->sd_id = fileID;
    strcpy(file->spatialResolution, "850 m");

    strcpy(cal_path_tab, cal_path);

    l1a_data = (int16_t *) calloc(npix * NBANDS_OSMI*NREC_IN_BUF, sizeof (int16_t));
    l1b_buffer = (float *) calloc(npix*NBANDS_OSMI, sizeof (float));

    msec = (int32_t *) calloc(nscan, sizeof (int32_t));
    status = rdSDS(file->sd_id, "msec", 0, 0, 0, 0, (VOIDP) msec);

    sc_ana = (float *) calloc(nscan * 40, sizeof (float));
    status = rdSDS(file->sd_id, "sc_ana", 0, 0, 0, 0, (VOIDP) sc_ana);

    sc_ttag = (double *) calloc(nscan, sizeof (double));
    status = rdSDS(file->sd_id, "sc_ttag", 0, 0, 0, 0, (VOIDP) sc_ttag);

    sc_ttag_inc = (double *) calloc(nscan, sizeof (double));
    status = rdSDS(file->sd_id, "sc_ttag_inc", 0, 0, 0, 0, (VOIDP) sc_ttag_inc);

    return (status);
}

int readl1_osmi(filehandle *file, int32_t recnum, l1str *l1rec) {
    /*recnum is 0-based */

    /*                                                                 */
    /* get_l1a_rec interface                                           */
    /*                                                                 */
    static cal_mod_struc cal_mod; /* cal modification structure       */
    static int16_t *l2_flags; /* radiance quality flags for L2    */
    static float *l1b_data;
    static int chindx[6] = {0, 1, 2, 4, 6, 7}; /* active bands in L1B data block*/

    int32_t status = 0;

    int32_t nwave = l1rec->l1file->nbands;
    int32_t *bindx = l1rec->l1file->bindx;

    /*                                                                 */
    /* local vars                                                      */
    /*                                                                 */
    int32_t ipix; /* pixel number                     */
    int32_t iw, ib;



    /* Get l1a data */
    /* ------------ */
    status = get_l1a_rec_osmi(file->sd_id, recnum + 1, &cal_mod,
            &l1b_data, &l2_flags);


    /*                                                              */
    /* Copy scan geolocation and view geometry                      */
    /*                                                              */
    memcpy(l1rec->lat, ylat, npix * sizeof (float));
    memcpy(l1rec->lon, xlon, npix * sizeof (float));
    memcpy(l1rec->solz, solz, npix * sizeof (float));
    memcpy(l1rec->sola, sola, npix * sizeof (float));
    memcpy(l1rec->senz, senz, npix * sizeof (float));
    memcpy(l1rec->sena, sena, npix * sizeof (float));

    /*                                                              */
    /* Copy L1B radiances, pixel interlaced by band. Add per-band   */
    /* view angles.                                                 */
    /*                                                              */
    for (ipix = 0; ipix < file->npix; ipix++) {

        if (l1rec->sena[ipix] < -180.0) l1rec->sena[ipix] += 360.0;

        for (iw = 0; iw < nwave; iw++) {
            ib = bindx[iw];
            l1rec->Lt [ipix * nwave + ib] = l1b_data[chindx[iw] * npix + ipix];
        }
        l1rec->stlight[ipix] = ((l2_flags[ipix] & STRAYLIGHT) > 0);
        l1rec->hilt [ipix] = ((l2_flags[ipix] & HILT) > 0);
        if (recnum < 96) l1rec->hilt[ipix] = 1; /* bad position interp */
        if (ipix < 15) l1rec->hilt[ipix] = 1; /* instrument shading? */
        if (ipix > 1028) l1rec->hilt[ipix] = 1; /* instrument shading? */
    }

    /*                                                              */
    /* Set scan time in L1B output record                           */
    /*                                                              */
    int16_t year = (int16_t) syear;
    int16_t day = (int16_t) sday;

    if ((msec[recnum] < msec[recnum + 1]) && (recnum < nscan - 1)) {
        /* adjust for day rollover */
        day += 1;
        if (day > (365 + (year % 4 == 0))) {
            year += 1;
            day = 1;
        }
    }

    l1rec->scantime = yds2unix(year, day, (double) (msec[recnum] / 1.e3));

    return (status);
}

int closel1_osmi(filehandle *file) {

    free(l1a_data);
    free(l1b_buffer);
    free(msec);
    free(sc_ana);
    free(sc_ttag);
    free(sc_ttag_inc);

    SDend(file->sd_id);

    return (0);
}

int32_t get_l1a_rec_osmi(int32_t sd_id, int32_t recno, cal_mod_struc *cal_mod,
        float **l1b_data, int16_t **l2_flags) {
    /* recno is 1-based */

    /*                                                                 */
    /* local vars                                                      */
    /*                                                                 */
    int32_t ipix; /* pixel number                     */
    int32_t idet; /* detector number                  */
    int32_t irec;
    int32_t start[3] = {0, 0, 0};
    int32_t edges[3];

    int16_t band_buf[LAC_PIXEL_NUM * NBANDS_OSMI];

    static int32_t crec;
    static int16_t n_read = 0;
    static int32_t max_rec_in_rdbuf;
    static int32_t offset;
    static int32_t i, j;

    static byte first = 1;

    int16_t i16dum;

    int32_t recno1, recno2;
    double uepoch, epoch, epoch1, epoch2;
    double pos1[3], vel1[3];
    double pos2[3], vel2[3];
    double pos [3], vel [3];

    struct GEODETIC obs;
    struct TOPOCENTRIC top;

    /*                                                                 */
    /* Read L1A data scan by scan, build L1B record, and write.        */

    if (first) {
        max_rec_in_rdbuf = recno - 1;
    }

    crec = recno;

    if (crec > nscan) return 0;


    if (crec > max_rec_in_rdbuf || crec < offset) {

        n_read = (nscan - recno >= NREC_IN_BUF) ?
                NREC_IN_BUF : nscan - recno + 1;


        edges[0] = 1;
        start[1] = 0;
        edges[1] = npix;
        edges[2] = 8;
        for (irec = 0; irec < n_read; irec++) {
            start[0] = recno - 1 + irec;
            SDreaddata(SDselect(sd_id, SDnametoindex(sd_id, "l1a_data")),
                    start, NULL, edges, (VOIDP) band_buf);
            for (ipix = 0; ipix < npix; ipix++)
                for (idet = 0; idet < 8; idet++)
                    l1a_data[irec * (NBANDS_OSMI * npix) + idet * npix + ipix] =
                        band_buf[ipix * NBANDS_OSMI + idet];
        }

        offset = recno;
        max_rec_in_rdbuf = offset + n_read - 1;

    }

    if ((crec - offset) >= 0) {
        j = crec - offset;

        calibrate_l1a_osmi(cal_path_tab, syear, sday, smsec, eday,
                msec[recno - 1], dtype, nsta, npix,
                (recno - 1) % 96, &i16dum, i16dum, 0,
                &l1a_data[j * NBANDS_OSMI * npix], l1b_buffer,
                cal_mod);


        rdSDS(sd_id, "lat", recno - 1, 0, 1, npix, (VOIDP) ylat);
        rdSDS(sd_id, "lon", recno - 1, 0, 1, npix, (VOIDP) xlon);

        /* Calculate sunz/a and senz/a */

        uepoch = itojul(700101, 0); /* unix time epoch   */

        if (recno <= 96) {
            recno1 = recno - 1;
            recno2 = recno - 1;
        } else {
            recno1 = recno - 1; /* current scan      */
            recno2 = recno - 96 - 1; /* next scan in time */
        }

        /* get time and ecef s/c pos/vel for bounding records */
        epoch1 = uepoch + sc_ttag[recno1] / 86400.0;
        pos1[0] = sc_ana[40 * (recno1) + 9] * 1000.0;
        pos1[1] = sc_ana[40 * (recno1) + 10] * 1000.0;
        pos1[2] = sc_ana[40 * (recno1) + 11] * 1000.0;
        vel1[0] = sc_ana[40 * (recno1) + 12] * 1000.0;
        vel1[1] = sc_ana[40 * (recno1) + 13] * 1000.0;
        vel1[2] = sc_ana[40 * (recno1) + 14] * 1000.0;

        epoch2 = uepoch + sc_ttag[recno2] / 86400.0;
        pos2[0] = sc_ana[40 * (recno2) + 9] * 1000.0;
        pos2[1] = sc_ana[40 * (recno2) + 10] * 1000.0;
        pos2[2] = sc_ana[40 * (recno2) + 11] * 1000.0;
        vel2[0] = sc_ana[40 * (recno2) + 12] * 1000.0;
        vel2[1] = sc_ana[40 * (recno2) + 13] * 1000.0;
        vel2[2] = sc_ana[40 * (recno2) + 14] * 1000.0;

        for (i = 0; i < npix; i++) {

            epoch = uepoch + (sc_ttag[recno - 1] +
                    /*
                    ((npix - 1 - i) * sc_ttag_inc[recno-1]))/86400.0;
                     */
                    ((i) * sc_ttag_inc[recno - 1])) / 86400.0;

            obs.lat = DTOR(ylat[i]);
            obs.lon = DTOR(xlon[i]);
            obs.height = 0;

            sunpos(tconv(epoch - AU / (VLIGHT * 86400.0), "utc:tdt", NULL), pos, NULL);
            ctotc(pos, &obs, gmha(epoch), &top, NULL);
            solz[i] = (float) (90.0 - RTOD(top.el));
            sola[i] = (float) RTOD(top.az);
            if (sola[i] < 0.0) sola[i] = sola[i] + 360.0;

            intpos_(&epoch1, pos1, vel1, &epoch2, pos2, vel2, &epoch, pos, vel);

            /*
            printf("%d %lf %lf %lf %lf\n",i,epoch1,pos1[0],pos1[1],pos1[2]);
            printf("%d %lf %lf %lf %lf\n",i,epoch ,pos [0],pos [1],pos [2]);
            printf("%d %lf %lf %lf %lf\n",i,epoch2,pos2[0],pos2[1],pos2[2]);
             */

            ctotc(pos, &obs, gmha(epoch), &top, NULL);
            senz[i] = (float) (90.0 - RTOD(top.el));
            sena[i] = (float) RTOD(top.az);
            if (sena[i] < 0.0) sena[i] = sena[i] + 360.0;

        }

    }

    *l1b_data = l1b_buffer;
    *l2_flags = l2_flags_buffer;
    first = 0;

    return (0);
}

