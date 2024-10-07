/* =========================================================== */
/* Module l1_hdf_generic_write.c                               */
/*                                                             */
/* Functions to open and write a multi-sensor (generic) l1     */
/* file in HDF format.                                         */
/*                                                             */
/* Written By:                                                 */
/*     B. A. Franz, SAIC GSC, SIMBIOS Project, January 1999.   */
/*                                                             */
/* Modified By:                                                */
/*     J. G. Gales, Futuretech, SIMBIOS Project, Sept. 1999.   */
/*           Generate standard L1B product                     */
/*     Gene Eplee, SAIC GSC, SeaWiFS Project, December 2000.   */
/*           Update time correction and mirror side            */
/*                       factors.                              */
/*     Gene Eplee, SAIC, SeaWiFS Project, March 2004.          */
/*           Convert time correction and mirror side           */
/*           factors to simultaneous exponentials.             */
/* =========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h> 
#include <time.h>
#include <math.h>

#include "l12_proto.h"

#include "hdf_utils.h"
#include <timeutils.h>
#include "mfhdf.h"

/* Global variables to facilitate communication
between some of the functions defined in this file */
static int32 sd_id;
static int32 numScans;
static int32 numPixels;
static int32 numBands;
static int32 spix;
static int32 cpix;
static int32 epix;
static int32 dpix;
static int32 cscan;
static int32 nctl;
static int32 *ictl;
static int32 *jctl;
static int32 sensorID;
static int32_t evalmask;
int32 *wavelength;
char rad_name[8];

int32 get_ctl_pts(int32 npix, int32 nscans, int32 ictl[], int32 jctl[]) {
    int32_t i;

    for (i = 0; i < nscans; i++)
        jctl[i] = i + 1;

    for (i = 0; i < npix; i++)
        ictl[i] = i + 1;

    return (npix);
}

/* -------------------------------------------------------- */
/* Assign SDSes to Vgroups                                  */

/* -------------------------------------------------------- */
int MakeVgroupsL1(filehandle *file) {

    int32 i;
    int32 h_id;
    int32 v_id;
    int32 sd_id = file->sd_id;

    /* Do we have the extra meta-data for SeaWiFS */
    int seawifs_meta = 0;
    if (sensorID == SEAWIFS) {
        int32 sds_id;
        if (sd_select(sd_id, "scan_ell", &sds_id) == 0)
            seawifs_meta = 1;
    }

    h_id = Hopen(file->name, DFACC_RDWR, 0);
    if (h_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: Hopen() failed for file, %s.\n",
                __FILE__, __LINE__, file->name);
        return (HDF_FUNCTION_ERROR);
    }
    Vstart(h_id);

    if (seawifs_meta) {

        /* Sensor Tilt */
        DPTB(v_attach(h_id, &v_id));
        Vsetclass(v_id, "Per File Data");
        Vsetname(v_id, "Sensor Tilt");
        DPTB(AddSdsToVgroup(sd_id, v_id, "ntilts"));
        DPTB(AddSdsToVgroup(sd_id, v_id, "tilt_flags"));
        DPTB(AddSdsToVgroup(sd_id, v_id, "tilt_ranges"));
        Vdetach(v_id);

    }

    /* Scan-Line Attributes */
    PTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Scan-Line Attributes");
    PTB(AddSdsToVgroup(sd_id, v_id, "year"));
    PTB(AddSdsToVgroup(sd_id, v_id, "day"));
    PTB(AddSdsToVgroup(sd_id, v_id, "msec"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "slon"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "clon"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "elon"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "slat"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "clat"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "elat"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "csol_z"));
    Vdetach(v_id);

    /* Image Data */
    PTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Geophysical Data");
    for (i = 0; i < numBands; i++) {
        sprintf(rad_name, "L_%3d", wavelength[i]);
        PTB(AddSdsToVgroup(sd_id, v_id, rad_name));
    }
    Vdetach(v_id);

    /* Navigation */
    PTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Navigation Data");
    PTB(AddSdsToVgroup(sd_id, v_id, "longitude"));
    PTB(AddSdsToVgroup(sd_id, v_id, "latitude"));
    PTB(AddSdsToVgroup(sd_id, v_id, "solz"));
    PTB(AddSdsToVgroup(sd_id, v_id, "sola"));
    PTB(AddSdsToVgroup(sd_id, v_id, "senz"));
    PTB(AddSdsToVgroup(sd_id, v_id, "sena"));
    PTB(AddSdsToVgroup(sd_id, v_id, "tilt"));

    if (seawifs_meta) {
        DPTB(AddSdsToVgroup(sd_id, v_id, "orb_vec"));
        DPTB(AddSdsToVgroup(sd_id, v_id, "sun_ref"));
        DPTB(AddSdsToVgroup(sd_id, v_id, "att_ang"));
        DPTB(AddSdsToVgroup(sd_id, v_id, "sen_mat"));
        DPTB(AddSdsToVgroup(sd_id, v_id, "scan_ell"));
        DPTB(AddSdsToVgroup(sd_id, v_id, "nflag"));
    }
    Vdetach(v_id);

    Vend(h_id);

    if (Hclose(h_id) != SUCCEED) {
        fprintf(stderr, "-E- %s line %d: Hclose(%d) failed for file, %s .\n",
                __FILE__, __LINE__, h_id, file->name);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}


/* -------------------------------------------------------- */
/* Create the level 2 SDSes for the scan-line data          */

/* -------------------------------------------------------- */
int32 CreateScanData(int32 np, int32 nb, int32 nscan,
        int32 sensorID, int32 bindx[]) {
    int32 i;
    int32 status;
    char longname[64];
    int32 *wavelen;

    /* Get sensor wavelengths */
    if (rdsensorinfo(sensorID, evalmask, "Lambda", (void **) &wavelen) != nb) {
        printf("-E- %s: Unable to determine sensor wavelengths\n", __FILE__);
        exit(FATAL_ERROR);
    }

    if ((wavelength = (int32 *) calloc(nb, sizeof (int32))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    for (i = 0; i < nb; i++)
        wavelength[i] = wavelen[bindx[i]];

    PTB(CreateSDS(
            sd_id, /* file id */
            "year", /* short name */
            "Scan year", /* long name */
            NULL, /* standard name   */
            "years", /* units */
            0, 0, /* valid range */
            0, 0, /* slope */
            DFNT_INT32, /* HDF number type */
            1, /* rank */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "day", /* short name */
            "Scan day of year", /* long name */
            NULL, /* standard name   */
            "days", /* units */
            0, 0, /* valid range */
            0, 0, /* slope */
            DFNT_INT32, /* HDF number type */
            1, /* rank */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "msec", /* short name */
            "Scan-line time, milliseconds of day", /* long name */
            NULL, /* standard name   */
            "milliseconds", /* units */
            0, 0, /* valid range */
            0, 0, /* slope */
            DFNT_INT32, /* HDF number type */
            1, /* rank */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "detnum", /* short name */
            "Detector Number", /* long name */
            NULL, /* standard name   */
            "", /* units */
            0, 0, /* valid range */
            0, 0, /* slope */
            DFNT_INT8, /* HDF number type */
            1, /* rank */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "mside", /* short name */
            "Mirror Side", /* long name */
            NULL, /* standard name   */
            "", /* units */
            0, 0, /* valid range */
            0, 0, /* slope */
            DFNT_INT8, /* HDF number type */
            1, /* rank */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "longitude", /* short name */
            "Longitude", /* long name */
            NULL, /* standard name   */
            "degree", /* units */
            -180.0, 180.0, /* valid range */
            0, 0, /* slope */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            nscan, np, 1, /* dimension sizes */
            "Number of Scan Lines", "Pixels per Scan Line", NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "latitude", /* short name */
            "Latitude", /* long name */
            NULL, /* standard name   */
            "degree", /* units */
            -90.0, 90.0, /* valid range */
            0, 0, /* slope */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            nscan, np, 1, /* dimension sizes */
            "Number of Scan Lines", "Pixels per Scan Line", NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "solz", /* short name */
            "Solar zenith angle", /* long name */
            NULL, /* standard name   */
            "degree", /* units */
            -90.0, 90.0, /* valid range */
            0, 0, /* slope */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            nscan, np, 1, /* dimension sizes */
            "Number of Scan Lines", "Pixels per Scan Line", NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "sola", /* short name */
            "Solar azimuth angle", /* long name */
            NULL, /* standard name   */
            "degree", /* units */
            -180.0, 180.0, /* valid range */
            0, 0, /* slope */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            nscan, np, 1, /* dimension sizes */
            "Number of Scan Lines", "Pixels per Scan Line", NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "senz", /* short name */
            "Sensor zenith angle", /* long name */
            NULL, /* standard name   */
            "degree", /* units */
            -90.0, 90.0, /* valid range */
            0, 0, /* slope */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            nscan, np, 1, /* dimension sizes */
            "Number of Scan Lines", "Pixels per Scan Line", NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "sena", /* short name */
            "Sensor azimuth angle", /* long name */
            NULL, /* standard name   */
            "degree", /* units */
            -180.0, 180.0, /* valid range */
            0, 0, /* slope */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            nscan, np, 1, /* dimension sizes */
            "Number of Scan Lines", "Pixels per Scan Line", NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id */
            "tilt", /* short name */
            "Sensor tilt angle", /* long name */
            NULL, /* standard name   */
            "degree", /* units */
            -20.1, 20.1, /* valid range */
            0, 0, /* slope */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));


    for (i = 0; i < nb; i++) {

        sprintf(rad_name, "L_%3d", wavelength[i]);
        sprintf(longname, "Top of atmosphere %3d nm radiance", wavelength[i]);

        PTB(CreateSDS(
                sd_id,
                rad_name, /* short name      */
                longname, /* long name       */
                NULL, /* standard name   */
                "mW cm^-2 um^-1 sr^-1", /* units           */
                0, 0, /* valid range     */
                0, 0, /* slope           */
                DFNT_FLOAT32, /* HDF number type */
                2, /* rank            */
                nscan, np, 1, /* dimension sizes */
                "Number of Scan Lines", "Pixels per Scan Line", NULL
                /* dimension names */
                ));
    }

    PTB(CreateSDS(
            sd_id, /* file id */
            "l2_flags", /* short name */
            "Bit masks and flags", /* long name */
            NULL, /* standard name   */
            "", /* units */
            0, 0, /* valid range */
            0, 0, /* slope */
            DFNT_INT32, /* HDF number type */
            2, /* rank */
            nscan, np, 1, /* dimension sizes */
            "Number of Scan Lines", "Pixels per Scan Line", NULL
            /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id         */
            "slon", /* short name      */
            "Starting Longitude", /* long name       */
            NULL, /* standard name   */
            "degree", /* units           */
            -180.0, 180.0, /* valid range     */
            0, 0, /* slope           */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank            */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id         */
            "clon", /* short name      */
            "Center Longitude", /* long name       */
            NULL, /* standard name   */
            "degree", /* units           */
            -180.0, 180.0, /* valid range     */
            0, 0, /* slope           */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank            */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id         */
            "elon", /* short name      */
            "Ending Longitude", /* long name       */
            NULL, /* standard name   */
            "degree", /* units           */
            -180.0, 180.0, /* valid range     */
            0, 0, /* slope           */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank            */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id         */
            "slat", /* short name      */
            "Starting Latitude", /* long name       */
            NULL, /* standard name   */
            "degree", /* units           */
            -90.0, 90.0, /* valid range     */
            0, 0, /* slope           */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank            */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id         */
            "clat", /* short name      */
            "Center Latitude", /* long name       */
            NULL, /* standard name   */
            "degree", /* units           */
            -90.0, 90.0, /* valid range     */
            0, 0, /* slope           */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank            */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id         */
            "elat", /* short name      */
            "Ending Latitude", /* long name       */
            NULL, /* standard name   */
            "degree", /* units           */
            -90.0, 90.0, /* valid range     */
            0, 0, /* slope           */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank            */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id         */
            "csol_z", /* short name      */
            "Center Solar Zenith Angle", /* long name       */
            NULL, /* standard name   */
            "degree", /* units           */
            -90.0, 90.0, /* valid range     */
            0, 0, /* slope           */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank            */
            nscan, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    free(wavelength);
    return (LIFE_IS_GOOD);
}


/* -------------------------------------------------------- */
/* Open an HDF file and store global attributes in it.      */

/* -------------------------------------------------------- */
int openl1_write_hdf(filehandle *l1file) {
    char *name = l1file->name;
    int32 nbands = (int32) l1file->nbands;
    int32 npix = (int32) l1file->npix;
    int32 nscans = (int32) l1file->nscan;
    int32 sensorID = (int32) l1file->sensorID;
    int32 format = (int32) l1file->format;
    int32 *bindx = (int32 *) l1file->bindx;
    char *pro_control = l1file->pro_control;

    char title [255];
    char soft_id[200]; /* software version info */

    /* set globals */
    numScans = nscans;
    numPixels = npix;
    numBands = nbands;
    spix = l1file->spix + 1;
    cpix = numPixels / 2;
    epix = l1file->epix + 1;
    dpix = l1file->dpix;
    cscan = numScans / 2;
    evalmask = l1file->input->evalmask;

    /* Get control-point array */
    if ((ictl = calloc(numPixels, sizeof (int32))) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate control-point array.\n",
                __FILE__, __LINE__);
        return (MEMORY_ALLOCATION_ERROR);
    }
    if ((jctl = calloc(numScans, sizeof (int32))) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate control-point array.\n",
                __FILE__, __LINE__);
        return (MEMORY_ALLOCATION_ERROR);
    }
    nctl = get_ctl_pts(numPixels, numScans, ictl, jctl);

    /* Create the HDF file */
    sd_id = SDstart(name, DFACC_CREATE);
    l1file->sd_id = sd_id;
    if (sd_id == FAIL) {
        fprintf(stderr,
                "-E- %s line %d: Could not create HDF file, %s .\n",
                __FILE__, __LINE__, name);
        return (HDF_FUNCTION_ERROR);
    }


    /* Write out some global attributes */
    strcpy(title, sensorName[sensorID]);
    strcat(title, " Level-1B");
    strcpy(soft_id, VERSION);

    idDS ds_id;
    ds_id.deflate = 0;
    ds_id.fid = sd_id;
    ds_id.fftype = DS_HDF;

    PTB(SetChrGA(ds_id, "product_name", basename(name)));
    PTB(SetChrGA(ds_id, "title", title));
    PTB(SetChrGA(ds_id, "software_name", "MSl1bgen"));
    PTB(SetChrGA(ds_id, "software_version", soft_id));
    PTB(SetChrGA(ds_id, "date_created", ydhmsf(now(), 'G')));
    PTB(SetChrGA(ds_id, "history", pro_control));
    //PTB( SetI32GA (ds_id, "Number of Bands", numBands)                 );
    //PTB( SetI32GA (ds_id, "Pixels per Scan Line", numPixels)           );
    //PTB( SetI32GA (ds_id, "Number of Scan Lines",numScans));
    //PTB( SetI32GA (ds_id, "Start Pixel",spix)                          );
    //PTB( SetI32GA (ds_id, "Pixel Subsampling Interval",dpix)           );
    //    PTB( SetI32GA (ds_id, "Scene Center Scan Line", cscan)             );
    //PTB( SetI32GA (ds_id, "Number of Scan Control Points",numScans)    );
    //PTB( SetI32GA (ds_id, "Number of Pixel Control Points",nctl)       );

    /* Write control-point arrays */
    PTB(CreateSDS(
            sd_id, /* file id         */
            "cntl_pt_cols", /* short name      */
            "Pixel control points", /* long name       */
            NULL, /* standard name   */
            NULL, /* units           */
            0, 0, /* range           */
            0, 0, /* slope           */
            DFNT_INT32, /* HDF number type */
            1, /* rank            */
            nctl, 1, 1, /* dimension sizes */
            "Number of Pixel Control Points", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(
            sd_id, /* file id         */
            "cntl_pt_rows", /* short name      */
            "Scan control points", /* long name       */
            NULL, /* standard name   */
            NULL, /* units           */
            0, 0, /* range           */
            0, 0, /* slope           */
            DFNT_INT32, /* HDF number type */
            1, /* rank            */
            numScans, 1, 1, /* dimension sizes */
            "Number of Scan Control Points", NULL, NULL /* dimension names */
            ));
    PTB(sd_writedata(sd_id, "cntl_pt_cols", ictl, 0, 0, 0, nctl, 0, 0));
    PTB(sd_writedata(sd_id, "cntl_pt_rows", jctl, 0, 0, 0, numScans, 0, 0));

    /* Create the scan-line SDSes */
    PTB(CreateScanData(numPixels, numBands, nscans, sensorID, bindx));

    return (LIFE_IS_GOOD);
}

/* -------------------------------------------------------- */
/* Create the SDSes for the scan-line data                  */

/* -------------------------------------------------------- */
int writel1_hdf(filehandle *l1file, int32_t recnum, l1str *l1rec) {
    int32 *year = (int32 *) l1rec->year;
    int32 *day = (int32 *) l1rec->day;
    int32 *msec = (int32 *) l1rec->msec;
    float32 *lon = (float32 *) l1rec->lon;
    float32 *lat = (float32 *) l1rec->lat;
    float32 *solz = (float32 *) l1rec->solz;
    float32 *sola = (float32 *) l1rec->sola;
    float32 *senz = (float32 *) l1rec->senz;
    float32 *sena = (float32 *) l1rec->sena;
    float32 *tilt = (float32 *)&(l1rec->tilt);
    int32 *l2_flags = (int32 *) l1rec->flags;
    int32 *bindx = (int32 *) l1file->bindx;

    static float32 *data = NULL;
    int32 i, j;

    idDS ds_id;
    ds_id.deflate = 0;
    ds_id.fid = sd_id;
    ds_id.fftype = DS_HDF;

    /* Write the scan-line data */
    PTB(sd_writedata(sd_id, "year", year, recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "day", day, recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "msec", msec, recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "longitude", lon, recnum, 0, 0, 1, numPixels, 1));
    PTB(sd_writedata(sd_id, "latitude", lat, recnum, 0, 0, 1, numPixels, 1));

    PTB(sd_writedata(sd_id, "detnum", &l1rec->detnum, recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "mside", &l1rec->mside, recnum, 0, 0, 1, 1, 1));

    PTB(sd_writedata(sd_id, "slon", &lon[0], recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "clon", &lon[cpix], recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "elon", &lon[numPixels - 1], recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "slat", &lat[0], recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "clat", &lat[cpix], recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "elat", &lat[numPixels - 1], recnum, 0, 0, 1, 1, 1));
    PTB(sd_writedata(sd_id, "csol_z", &solz[cpix], recnum, 0, 0, 1, 1, 1));

    PTB(sd_writedata(sd_id, "solz", solz, recnum, 0, 0, 1, numPixels, 1));
    PTB(sd_writedata(sd_id, "sola", sola, recnum, 0, 0, 1, numPixels, 1));
    PTB(sd_writedata(sd_id, "senz", senz, recnum, 0, 0, 1, numPixels, 1));
    PTB(sd_writedata(sd_id, "sena", sena, recnum, 0, 0, 1, numPixels, 1));
    PTB(sd_writedata(sd_id, "tilt", tilt, recnum, 0, 0, 1, 1, 1));

    /* Write the radiance data */
    if (data == NULL) {
        if ((data = (float32 *) malloc(numPixels * sizeof (float32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Memory allocation failure.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
    }

    for (i = 0; i < numBands; i++) {
        sprintf(rad_name, "L_%3d", wavelength[i]);
        for (j = 0; j < numPixels; j++) {
            data[j] = (float32) l1rec->Lt[j * numBands + bindx[i]];
        }
        PTB(sd_writedata(sd_id, rad_name, data,
                recnum, 0, 0, 1, numPixels, 0));
    }

    /* Write the l2_flag data */
    PTB(sd_writedata(sd_id, "l2_flags", l2_flags,
            recnum, 0, 0, 1, numPixels, 0));

    /* Write global attributes */
    if (recnum == (numScans - 1)) {
        scene_meta_write(ds_id);
    }

    return (LIFE_IS_GOOD);
}



/* -------------------------------------------------------- */
/* Finish access for the current file.                      */

/* -------------------------------------------------------- */
void closel1_hdf(filehandle *l1file) {
    /* Define Vgroups */
    MakeVgroupsL1(l1file);

    if (SDend(l1file->sd_id)) {
        fprintf(stderr, "-E- %s line %d: SDend(%d) failed.\n",
                __FILE__, __LINE__, sd_id);
    }
}


