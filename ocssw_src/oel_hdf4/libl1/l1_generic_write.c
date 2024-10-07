/* =========================================================== */
/* Module l1_generic_write.c                                   */
/*                                                             */
/* Functions to open and write a multi-sensor (generic) l1     */
/* file in HDF format.                                         */
/*                                                             */
/* Written By:                                                 */
/*     B. A. Franz, SAIC GSC, SIMBIOS Project, January 1999.   */
/*                                                             */
/* Modified By:                                                */
/*     J. M. Gales, Futuretech, SIMBIOS Project, Sept. 1999.   */
/*           Generate standard L1B product                     */
/*     Gene Eplee, SAIC GSC, SeaWiFS Project, December 2000.   */
/*           Update time correction and mirror side            */
/*                       factors.                              */
/*     Gene Eplee, SAIC, SeaWiFS Project, March 2004.          */
/*           Convert time correction and mirror side           */
/*           factors to simultaneous exponentials.             */
/*     J. M. Gales, Futuretech, Oct 2013.                      */
/*           Add support for NETCDF4 and CF-metadata           */
/* =========================================================== */

#include <netcdf.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h> 
#include <time.h>
#include <math.h>

#include "l1.h"
#include "l1_aci.h"

#include "dfutils.h"
#include <timeutils.h>
#include "mfhdf.h"

#include "scene_meta.h"

/* Global variables to facilitate communication
 between some of the functions defined in this file */
static int32_t sd_id;
static int32_t numScans;
static int32_t numPixels;
static int32_t numBands;
static int32_t numBandsIR;
static int32_t spix;
static int32_t cpix;
static int32_t epix;
static int32_t dpix;
static int32_t cscan;
static int32_t nctl;
static int32_t *ictl;
static int32_t *jctl;
static int32_t sensorID;
static int32_t evalmask;
int32_t *wavelength;
static int32_t *Lambda, *Lambda_p;
char rad_name[8];

int32_t get_ctl_pts(int32_t npix, int32_t nscans, int32_t ictl[], int32_t jctl[]) {
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

    int32_t i;
    int32_t h_id;
    int32_t v_id;
    int32_t sd_id = file->sd_id;

    /* Do we have the extra meta-data for SeaWiFS */
    int seawifs_meta = 0;
    if (sensorID == SEAWIFS) {
        int32_t sds_id;
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
        sprintf(rad_name, "Lt_%3d", wavelength[i]);
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
/* Open an HDF file and store global attributes in it.      */

/* -------------------------------------------------------- */
int openl1_write(filehandle *l1file) {
    char *name = l1file->name;
    int32_t nbands = (int32_t) l1file->nbands;
    int32_t npix = (int32_t) l1file->npix;
    int32_t nscans = (int32_t) l1file->nscan;
    int32_t sensorID = (int32_t) l1file->sensorID;
    int32_t *bindx = (int32_t *) l1file->bindx;
    char avhrrbird[10];

    char buf1[1024];
    char tmp_str[2048];
    char title[255];
    char soft_id[200]; /* software version info */
    char* bandNumberStr;
    char* totalBandNumberStr;
    float *Gain;
    float *Offset;
    float *Fonom;
    float *Fobar, *Fobar_p;
    float *Tau_r, *Tau_r_p;
    float *k_oz, *k_oz_p;
    float *k_no2, *k_no2_p;

    int i, n;
    static int firstCall = 1;
    int32_t dm[3];
    const char dm_name[3][80];
    int flagbits[32];

    /* set globals */
    numScans = nscans;
    numPixels = npix;
    numBands = nbands;
    numBandsIR = l1file->nbandsir;
    spix = l1file->spix + 1;
    cpix = numPixels / 2;
    epix = l1file->epix + 1;
    dpix = l1_input->dpixl;
    cscan = numScans / 2;
    evalmask = l1_input->evalmask;

    if (firstCall == 1) {
        if ((wavelength = (int32_t *) calloc(numBands, sizeof (int32_t))) == NULL) {
            printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((Lambda = (int32_t *) calloc((numBands + numBandsIR), sizeof (int32_t))) == NULL) {
            printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        firstCall = 0;
    }
    if ((Gain = (float *) calloc(numBands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((Offset = (float *) calloc(numBands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((Fonom = (float *) calloc(numBands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((Fobar = (float *) calloc(numBands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((Tau_r = (float *) calloc(numBands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((k_oz = (float *) calloc(numBands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((k_no2 = (float *) calloc(numBands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    /* Get control-point array */
    if ((ictl = calloc(numPixels, sizeof (int32_t))) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate control-point array.\n",
                __FILE__, __LINE__);
        return (MEMORY_ALLOCATION_ERROR);
    }
    if ((jctl = calloc(numScans, sizeof (int32_t))) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate control-point array.\n",
                __FILE__, __LINE__);
        return (MEMORY_ALLOCATION_ERROR);
    }
    nctl = get_ctl_pts(numPixels, numScans, ictl, jctl);

    /* Create the L1B file */
    idDS ds_id;
    int32_t nt_chr, nt_i32;

    if (l1file->format == FT_L1HDF) {
        ds_id.fftype = DS_HDF;
        nt_chr = DFNT_CHAR;
        nt_i32 = DFNT_INT32;
    } else if (l1file->format == FT_L1BNCDF) {
        ds_id.fftype = DS_NCDF;
        nt_chr = NC_CHAR;
        nt_i32 = NC_INT;
    }

    // Make sure startDS stores ds_id.fftype
    ds_id = startDS(name, ds_id.fftype, DS_WRITE, 0);
    l1file->sd_id = ds_id.fid;
    if (l1file->sd_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: Could not create L1B file, %s .\n",
                __FILE__, __LINE__, name);
        return (HDF_FUNCTION_ERROR);
    }

    totalBandNumberStr = "total band number";
    bandNumberStr = "band number";

    /* Get sensor-specific attributes */
    if ((n = rdsensorinfo(l1file->sensorID, l1_input->evalmask, NULL, NULL))
            != numBands) {
        fprintf(stderr, "-E- %s Line %d:  Error reading sensor table. %d %d\n",
                __FILE__, __LINE__, numBands, n);
        return (-1);
    }
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Lambda",
            (void **) &Lambda_p);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Fobar",
            (void **) &Fobar_p);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Tau_r",
            (void **) &Tau_r_p);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "k_oz",
            (void **) &k_oz_p);
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "k_no2",
            (void **) &k_no2_p);

    for (i = 0; i < numBands; i++) {
        Lambda[i] = Lambda_p[i];
        // multiply by 10 to put into W/m2/um, since internally all radiances are mW/cm2/um
        Fobar[i] = Fobar_p[i] * 10.0;
        Tau_r[i] = Tau_r_p[i];
        k_oz[i] = k_oz_p[i];
        k_no2[i] = k_no2_p[i];

        Gain[i] = l1_input->gain[i];
        Offset[i] = l1_input->offset[i];
        if (l1_input->outband_opt >= 2) {
            get_f0_thuillier_ext(Lambda[i], BANDW, &Fonom[i]);
            Fonom[i] *= 10.0;
        } else {
            Fonom[i] = Fobar[i];
        }
    }
    for (i = numBands; i < numBands + numBandsIR; i++)
        Lambda[i] = Lambda_p[i];

    if (l1file->format == FT_L1BNCDF) {
        int dumdim;
        totalBandNumberStr = "number_of_bands";
        bandNumberStr = "number_of_reflective_bands";
        if (nc_def_dim(ds_id.fid, "number_of_lines", numScans,
                &dumdim) != NC_NOERR)
            exit(1);
        if (nc_def_dim(ds_id.fid, "pixels_per_line", numPixels,
                &dumdim) != NC_NOERR)
            exit(1);
        if (nc_def_dim(ds_id.fid, "pixel_control_points", nctl,
                &dumdim) != NC_NOERR)
            exit(1);
        if (nc_def_dim(ds_id.fid, totalBandNumberStr, numBands + numBandsIR,
                &dumdim) != NC_NOERR)
            exit(1);
        if (nc_def_dim(ds_id.fid, bandNumberStr, numBands, &dumdim) != NC_NOERR)
            exit(1);

        nc_def_grp(ds_id.fid, "sensor_band_parameters", &l1file->grp_id[0]);
        nc_def_grp(ds_id.fid, "scan_line_attributes", &l1file->grp_id[2]);
        nc_def_grp(ds_id.fid, "geophysical_data", &l1file->grp_id[3]);
        nc_def_grp(ds_id.fid, "navigation_data", &l1file->grp_id[4]);
        nc_def_grp(ds_id.fid, "processing_control", &l1file->grp_id[5]);
    }

    /*                                                                  */
    /* Create the Sensor Band Parameters datasets                       */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if (l1file->format == FT_L1BNCDF)
        ds_id.fid = l1file->grp_id[0];

    dm[0] = numBands + numBandsIR;
    strcpy((char *) dm_name[0], totalBandNumberStr);
    PTB(createDS(ds_id, l1file->sensorID, "wavelength", dm, dm_name));

    dm[0] = numBands;
    strcpy((char *) dm_name[0], bandNumberStr);

    if (numBands > 0) {
        PTB(createDS(ds_id, l1file->sensorID, "vcal_gain", dm, dm_name));
        PTB(createDS(ds_id, l1file->sensorID, "vcal_offset", dm, dm_name));
        PTB(createDS(ds_id, l1file->sensorID, "F0", dm, dm_name));
        PTB(createDS(ds_id, l1file->sensorID, "k_oz", dm, dm_name));
        PTB(createDS(ds_id, l1file->sensorID, "k_no2", dm, dm_name));
        PTB(createDS(ds_id, l1file->sensorID, "Tau_r", dm, dm_name));
    }

    /* Write out some global attributes */
    strcpy(title, sensorId2SensorName(sensorID));
    strcat(title, " Level-1B");
    //strcpy(soft_id, VERSION);
    //sprintf(soft_id, "%d.%d.%d-%s", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GITSHA);

    if (l1file->format == FT_L1HDF) {
        PTB(SetChrGA(ds_id, "Sensor Name", sensorId2SensorName(sensorID)));
        PTB(SetChrGA(ds_id, "Product Name", basename(name)));
        PTB(SetChrGA(ds_id, "Title", title));
        PTB(SetChrGA(ds_id, "Software Name", "l1bgen"));
        PTB(SetChrGA(ds_id, "Software Version", "2.0"));
        PTB(SetChrGA(ds_id, "Processing Time", ydhmsf(now(), 'G')));
        //PTB(SetChrGA(ds_id, "Processing Control", pro_control));
        PTB(SetI32GA(ds_id, "Number of Bands", numBands));
        PTB(SetI32GA(ds_id, "Pixels per Scan Line", numPixels));
        PTB(SetI32GA(ds_id, "Number of Scan Lines", numScans));
        PTB(SetI32GA(ds_id, "Start Pixel", spix));
        PTB(SetI32GA(ds_id, "Pixel Subsampling Interval", dpix));
        PTB(SetI32GA(ds_id, "Scene Center Scan Line", cscan));
        PTB(SetI32GA(ds_id, "Number of Scan Control Points", numScans));
        PTB(SetI32GA(ds_id, "Number of Pixel Control Points", nctl));

        dm[0] = nctl;
        strcpy((char *) dm_name[0], "Number of Pixel Control Points");
        PTB(createDS(ds_id, sensorID, "cntl_pt_cols", dm, dm_name));
        PTB(writeDS(ds_id, "cntl_pt_cols", ictl, 0, 0, 0, nctl, 0, 0));

        dm[0] = numScans;
        strcpy((char *) dm_name[0], "Number of Scan Lines");
        PTB(createDS(ds_id, sensorID, "cntl_pt_rows", dm, dm_name));
        PTB(writeDS(ds_id, "cntl_pt_rows", jctl, 0, 0, 0, numScans, 0, 0));

    } else {

        // global attr
        ds_id.fid = l1file->sd_id;

        PTB(SetChrGA(ds_id, "product_name", basename(name)));
        PTB(SetChrGA(ds_id, "title", title));
        //PTB(SetChrGA(ds_id, "history", pro_control));

        PTB(SetChrGA(ds_id, "instrument", sensorId2InstrumentName(l1file->sensorID)));

        if (l1file->sensorID == AVHRR) {
            strcpy(avhrrbird, "NOAA-");
            strncat(avhrrbird, xsatid2name(l1file->subsensorID) + 2, 2);
            PTB(SetChrGA(ds_id, "platform", avhrrbird));
        } else {
            PTB(SetChrGA(ds_id, "platform", sensorId2PlatformName(l1file->sensorID)));
        }
        PTB(SetChrGA(ds_id, "Conventions", "CF-1.6 ACDD-1.3"));
        PTB(SetChrGA(ds_id, "license", LICENSE));
        PTB(SetChrGA(ds_id, "naming_authority", NAMING_AUTHORITY));
        // create the id -
        if (strcmp(l1_input->pversion, "Unspecified") != 0) {
            strcpy(buf1, l1_input->pversion);
            strcat(buf1, "/L1B/");
        } else {
            strcpy(buf1, "L1B/");
        }
        strcat(buf1, basename(l1file->name));
        PTB(SetChrGA(ds_id, "id", buf1));

        time_t tnow;
        time(&tnow);
        strcpy(buf1, unix2isodate(tnow, 'G'));
        PTB(SetChrGA(ds_id, "date_created", buf1));

        PTB(SetChrGA(ds_id, "standard_name_vocabulary", STDNAME_VOCABULARY));
        PTB(SetChrGA(ds_id, "institution", INSTITUTION));
        PTB(SetChrGA(ds_id, "creator_name", CREATOR_NAME));
        PTB(SetChrGA(ds_id, "creator_email", CREATOR_EMAIL));
        PTB(SetChrGA(ds_id, "creator_url", CREATOR_URL));
        PTB(SetChrGA(ds_id, "project", PROJECT));
        PTB(SetChrGA(ds_id, "publisher_name", PUBLISHER_NAME));
        PTB(SetChrGA(ds_id, "publisher_url", PUBLISHER_URL));
        PTB(SetChrGA(ds_id, "publisher_email", PUBLISHER_EMAIL));

        PTB(SetChrGA(ds_id, "processing_level", "L1B"));
        PTB(SetChrGA(ds_id, "cdm_data_type", "swath"));
        if (l1file->orbit_node_lon > -180.0 && l1file->orbit_node_lon < 180.0)
            PTB(
                SetF32GA(ds_id, "equatorCrossingLongitude",
                l1file->orbit_node_lon));
        if (l1file->orbit_number > 0)
            PTB(SetI32GA(ds_id, "orbit_number", l1file->orbit_number));
        if (l1file->node_crossing_time > 0) {
            strcpy(buf1, unix2isodate(l1file->node_crossing_time, 'G'));
            PTB(SetChrGA(ds_id, "equatorCrossingDateTime", buf1));
        }

        PTB(SetChrGA(ds_id, "spatialResolution", l1file->spatialResolution));

        dm[0] = nctl;
        strcpy((char *) dm_name[0], "pixel_control_points");
        PTB(createDS(ds_id, sensorID, "cntl_pt_cols", dm, dm_name));
        PTB(writeDS(ds_id, "cntl_pt_cols", ictl, 0, 0, 0, nctl, 0, 0));

        dm[0] = numScans;
        strcpy((char *) dm_name[0], "number_of_lines");
        PTB(createDS(ds_id, sensorID, "cntl_pt_rows", dm, dm_name));
        PTB(writeDS(ds_id, "cntl_pt_rows", jctl, 0, 0, 0, numScans, 0, 0));

        // Processing Control attibutes
        ds_id.fid = l1file->grp_id[5];
        PTB(SetChrGA(ds_id, "software_name", "l1bgen"));
        PTB(SetChrGA(ds_id, "software_version", soft_id));

        // Write input parameters metadata
        ds_id.fid = l1file->grp_id[5];
        int32_t grp_id_input_parms;
        nc_def_grp(l1file->grp_id[5], "input_parameters", &grp_id_input_parms);
        ds_id.fid = grp_id_input_parms;

        char *end_str;
        char *tmp_parms = replace_ocroots(l1_input->input_parms);
        char *token = strtok_r(tmp_parms, "\n", &end_str);
        while (token != NULL) {
            char *end_token;
            char *name = strtok_r(token, "=", &end_token);
            trimBlanks(name);
            char *value = strtok_r(NULL, ";", &end_token);
            trimBlanks(value);
            if (name[0] != '#') {
                PTB(SetChrGA(ds_id, name, value));
            }
            token = strtok_r(NULL, "\n", &end_str);
        }
        free(tmp_parms);

    }

    /* Create the scan-line SDSes */
    char longname[64];
    int32_t *wavelen;
    int32_t np = numPixels;
    int32_t nb = numBands;
    int32_t nscan = nscans;

    /* Get sensor wavelengths */
    if (rdsensorinfo(sensorID, evalmask, "Lambda", (void **) &wavelen) != nb) {
        printf("-E- %s: Unable to determine sensor wavelengths\n", __FILE__);
        exit(FATAL_ERROR);
    }

    for (i = 0; i < nb; i++)
        wavelength[i] = wavelen[bindx[i]];

    if (l1file->format == FT_L1BNCDF) {
        ds_id.fid = l1file->grp_id[2];
        strcpy((char *) dm_name[0], "number_of_lines");
    } else {
        strcpy((char *) dm_name[0], "Number of Scan Lines");
    }
    dm[0] = nscan;

    PTB(createDS(ds_id, sensorID, "year", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "day", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "msec", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "detnum", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "mside", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "slon", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "clon", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "elon", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "slat", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "clat", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "elat", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "csol_z", dm, dm_name));

    if (l1file->format == FT_L1BNCDF) {
        ds_id.fid = l1file->grp_id[4];
        strcpy((char *) dm_name[1], "pixels_per_line");
    } else {
        strcpy((char *) dm_name[1], "Pixels per Scan Line");
    }
    dm[1] = np;

    PTB(createDS(ds_id, sensorID, "longitude", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "latitude", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "solz", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "sola", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "senz", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "sena", dm, dm_name));
    PTB(createDS(ds_id, sensorID, "tilt", dm, dm_name));

    if (l1file->format == FT_L1BNCDF)
        ds_id.fid = l1file->grp_id[3];
    for (i = 0; i < nb; i++) {
        sprintf(rad_name, "Lt_%3d", wavelength[i]);
        sprintf(longname, "Top of atmosphere %3d nm radiance", wavelength[i]);

        PTB(createDS(ds_id, sensorID, rad_name, dm, dm_name));
    }

    PTB(createDS(ds_id, sensorID, "l2_flags", dm, dm_name));
    ds_id.sid = selectDS(ds_id, "l2_flags");

    tmp_str[0] = '\0';
    for (i = 0; i < L1_NFLAGS; i++) {
        strcat(tmp_str, l2_flag_lname[i]);
        flagbits[i] = pow(2, i);
        if (i < L1_NFLAGS - 1)
            strcat(tmp_str, " ");
        /*
         * Keep the old flag attribute set for HDF4 files
         */
        if (l1file->format == FT_L1HDF) {
            PTB(
                    setAttr(ds_id, l2_flag_sname[i], nt_chr,
                    strlen(l2_flag_lname[i]) + 1,
                    (VOIDP) l2_flag_lname[i]));
        }
    }
    PTB(setAttr(ds_id, "flag_masks", nt_i32, L1_NFLAGS, (VOIDP) flagbits));

    PTB(
            setAttr(ds_id, "flag_meanings", nt_chr, strlen(tmp_str) + 1,
            (VOIDP) tmp_str));
    endaccessDS(ds_id);

    if (l1file->format == FT_L1BNCDF)
        ds_id.fid = l1file->grp_id[0];
    PTB(
            writeDS(ds_id, "wavelength", Lambda, 0, 0, 0, numBands + numBandsIR,
            0, 0));

    if (numBands > 0) {
        PTB(writeDS(ds_id, "vcal_gain", Gain, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "vcal_offset", Offset, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "F0", Fobar, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "k_oz", k_oz, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "k_no2", k_no2, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "Tau_r", Tau_r, 0, 0, 0, numBands, 0, 0));
    }

    free(Gain);
    free(Offset);
    free(Fonom);
    free(Fobar);
    free(Tau_r);
    free(k_oz);
    free(k_no2);
    free(ictl);
    free(jctl);

    return (LIFE_IS_GOOD);
}

/* -------------------------------------------------------- */
/* Create the SDSes for the scan-line data                  */

/* -------------------------------------------------------- */
int writel1(filehandle *l1file, int32_t recnum, l1str *l1rec) {
    int32_t i, j;
    int16_t year, day;
    double dsec;
    unix2yds(l1rec->scantime, &year, &day, &dsec);
    int32_t msec = (int32_t) (dsec * 1000.0);

    float *lon = (float *) l1rec->lon;
    float *lat = (float *) l1rec->lat;

    void *pbuf;

    float *angledata[4];
    angledata[0] = (float *) l1rec->solz;
    angledata[1] = (float *) l1rec->senz;
    angledata[2] = (float *) l1rec->sola;
    angledata[3] = (float *) l1rec->sena;

    static char *anglenames[4] = {"solz", "senz", "sola", "sena"};
    static productInfo_t *p_info[4];
    if(!p_info[0]) {
        for (i = 0; i < 4; i++) {
            p_info[i] = allocateProductInfo();
            if (!p_info[i]) {
                printf("ERROR allocating Product Info memory for %s\n", anglenames[i]);
                exit(EXIT_FAILURE);
            }
            if (!findProductInfo(anglenames[i], l1rec->l1file->sensorID, p_info[i])) {
                printf("%s not found in XML product table\n", anglenames[i]);
                exit(EXIT_FAILURE);
            }
        }
    }


    float *tilt = (float *) &(l1rec->tilt);
    int32_t *l2_flags = (int32_t *) l1rec->flags;
    int32_t *bindx = (int32_t *) l1file->bindx;

    static float *data = NULL;

    idDS ds_id;
    ds_id.fid = l1file->sd_id;
    ds_id.deflate = 0;
    if (l1file->format == FT_L1BNCDF)
        ds_id.fftype = DS_NCDF;
    else
        ds_id.fftype = DS_HDF;

    /* Write the scan-line data */
    if (l1file->format == FT_L1BNCDF)
        ds_id.fid = l1file->grp_id[2];

    int32_t tmpI = year;
    PTB(writeDS(ds_id, "year", &tmpI, recnum, 0, 0, 1, 1, 1));
    tmpI = day;
    PTB(writeDS(ds_id, "day", &tmpI, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "msec", &msec, recnum, 0, 0, 1, 1, 1));

    PTB(writeDS(ds_id, "detnum", &l1rec->detnum, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "mside", &l1rec->mside, recnum, 0, 0, 1, 1, 1));

    PTB(writeDS(ds_id, "slon", &lon[0], recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "clon", &lon[cpix], recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "elon", &lon[numPixels - 1], recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "slat", &lat[0], recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "clat", &lat[cpix], recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "elat", &lat[numPixels - 1], recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "csol_z", &angledata[0][cpix], recnum, 0, 0, 1, 1, 1));

    /* Write navigation info */
    if (l1file->format == FT_L1BNCDF)
        ds_id.fid = l1file->grp_id[4];
    PTB(writeDS(ds_id, "longitude", lon, recnum, 0, 0, 1, numPixels, 1));
    PTB(writeDS(ds_id, "latitude", lat, recnum, 0, 0, 1, numPixels, 1));
    /* Write out the angles solz, senz, sola, sena */
    for (i = 0; i < 4; i++) {
        pbuf = scale_sds(angledata[i], p_info[i], numPixels);
        PTB(writeDS(ds_id, anglenames[i], pbuf, recnum, 0, 0, 1, numPixels, 1));
    }

    PTB(writeDS(ds_id, "tilt", tilt, recnum, 0, 0, 1, 1, 1));

    /* Write the radiance data */
    if (l1file->format == FT_L1BNCDF)
        ds_id.fid = l1file->grp_id[3];
    if (data == NULL) {
        if ((data = (float *) malloc(numPixels * sizeof (float))) == NULL) {
            fprintf(stderr, "-E- %s line %d: Memory allocation failure.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
    }

    for (i = 0; i < numBands; i++) {
        sprintf(rad_name, "Lt_%3d", wavelength[i]);
        for (j = 0; j < numPixels; j++) {
            data[j] = (float) l1rec->Lt[j * numBands + bindx[i]] * 10.0;
        }
        PTB(writeDS(ds_id, rad_name, data, recnum, 0, 0, 1, numPixels, 0));
    }

    /* Write the l2_flag data */
    PTB(writeDS(ds_id, "l2_flags", l2_flags, recnum, 0, 0, 1, numPixels, 0));

    /* Write global attributes */
    if (recnum == (numScans - 1)) {
        if (l1file->format == FT_L1BNCDF)
            ds_id.fid = l1file->sd_id;
        scene_meta_write(ds_id);
    }

    return (LIFE_IS_GOOD);
}

/* -------------------------------------------------------- */
/* Finish access for the current file.                      */

/* -------------------------------------------------------- */
void closel1_generic(filehandle *l1file) {
    idDS ds_id;
    ds_id.deflate = 0;
    ds_id.fid = l1file->sd_id;
    if (l1file->format == FT_L1BNCDF)
        ds_id.fftype = DS_NCDF;
    else
        ds_id.fftype = DS_HDF;

    /* Define Vgroups */
    if (l1file->format == FT_L1HDF)
        MakeVgroupsL1(l1file);

    if (endDS(ds_id)) {
        fprintf(stderr, "-E- %s line %d: SDend(%d) failed.\n",
                __FILE__, __LINE__, sd_id);
    }
}

