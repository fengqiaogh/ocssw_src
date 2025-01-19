/* =========================================================== */
/* Module l2_generic.c                                         */
/*                                                             */
/* Functions to open and write a multi-sensor (generic) l2     */
/* file in HDF/NCDF format.                                    */
/*                                                             */
/* Written By:                                                 */
/*     Bryan A. Franz, SAIC GSC, March 1998.                   */
/*     Gary Fu,        SAIC GSC, March 1999.                   */
/*     Joel M. Gales,  Futuretech, Sept. 1999.                 */
/*     Gene Eplee, SAIC GSC, SeaWiFS Project, December 2000.   */
/*           Update time correction and mirror side            */
/*                       factors.                              */
/*     Gene Eplee, SAIC, SeaWiFS Project, March 2004.          */
/*           Convert time correction and mirror side           */
/*           factors to simultaneous exponentials.             */
/*     Joel Gales, Futuretech, OBPG Project, Sept 2012.        */
/*           Add support for L2 NETCDF4 output                 */
/*     Joel Gales, Futuretech, OBPG Project, Feb 2013.         */
/*           Add support for INT8,UINT8 datatypes for NETCDF4  */
/*     Joel Gales, Futuretech, OBPG Project, Nov 2013.         */
/*           Add support for CF-compliant metadata             */
/* =========================================================== */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h> 
#include <time.h>
#include <math.h>
#include "l12_proto.h"
#include "l2_generic.h"
#include <timeutils.h>
#include "l2prod.h"
#include "l1_aci.h"
#include "flags_sst.h"
#include "flags_iop.h"
#include "mph_flags.h"
#include "version.h"
#include "l1_seabass.h"
#include <scene_meta.h>
#include "l1_hawkeye.h"

/* Global variables to facilitate communication */
static float bad_float = BAD_FLT;
//static float fill_float = FIL_FLT;
//static float fill_int   = FIL_INT;
//static float fill_byte  = FIL_BYT;
static int32 numScans;
static int32 numPixels;
static int32 numBands;
static int32 numLvlProf;
static int32 n_refl_loc;
static int32 n_cloud_phase;
static int32 numBandsIR;
static int32 spix;
static int32 cpix;
static int32 epix;
static int32 cscan;

static FILE *fp_meta = NULL;

#define   GEOBOX_INC 20.0

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
#define RFACTOR 100

/* -------------------------------------------------------- */
/* Assign SDSes to Vgroups                                  */

/* -------------------------------------------------------- */
int MakeVgroups(filehandle *l2file) {

    int32 h_id;
    int32 v_id;
    int32 sd_id = l2file->sd_id;
    int i;

    /* Do we have the extra meta-data for SeaWiFS */
    int seawifs_meta = 0;
    if (l2file->sensorID == SEAWIFS) {
        int32 sds_id;
        if (sd_select(sd_id, "scan_ell", &sds_id) == 0)
            seawifs_meta = 1;
    }

    h_id = Hopen(l2file->name, DFACC_RDWR, 0);
    if (h_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: Hopen() failed for file, %s.\n",
                __FILE__, __LINE__, l2file->name);
        return (HDF_FUNCTION_ERROR);
    }
    Vstart(h_id);

    /* Scan-Line Attributes */
    DPTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per File Data");
    Vsetname(v_id, "Sensor Band Parameters");
    DPTB(AddSdsToVgroup(sd_id, v_id, "wavelength"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "vcal_gain"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "vcal_offset"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "F0"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "k_oz"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "Tau_r"));
    Vdetach(v_id);

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
    DPTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Scan-Line Attributes");
    DPTB(AddSdsToVgroup(sd_id, v_id, "year"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "day"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "msec"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "slon"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "clon"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "elon"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "slat"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "clat"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "elat"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "csol_z"));
    Vdetach(v_id);

    /* Geophysical Data */
    DPTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Geophysical Data");
    for (i = 0; i < l2file->tot_prod; i++) {
        DPTB(AddSdsToVgroup(sd_id, v_id, l2file->l2_prod_names[i]));
    }
    Vdetach(v_id);

    /* Navigation */
    DPTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Navigation Data");
    DPTB(AddSdsToVgroup(sd_id, v_id, "longitude"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "latitude"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "cntl_pt_cols"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "cntl_pt_rows"));
    DPTB(AddSdsToVgroup(sd_id, v_id, "tilt"));
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
                __FILE__, __LINE__, h_id, l2file->name);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}


/*----------------------------------------------------------------- */
/* Open an L2 file and store global attributes in it.               */

/* ---------------------------------------------------------------- */
int openl2(filehandle *l2file) {
    char title[256];
    char soft_id[200]; /* software version info */
    int32_t tot_prod;
    VOIDP pbuf;

    int i, k;
    l2prodstr *p; 
    char tmp_str[2048];
    char avhrrbird[10];
    int flagbits[32];
    int16_t sstflagbits[16];
    uint8_t byteflagbits[8];
    char buf1[1024];

    idDS ds_id;

    int32 dm[3];
    const char dm_name[3][80];
    char *end_str;
    float tmpFloat;

    // strings for the different named global attributes
    char* calibrationDataStr;
    char* equatorCrossingLonStr;
    char* historyStr;
    char* inputFilesStr;
    char* maskNamesStr;
    char* orbitNumberStr;
    char* processingVersionStr;
    char* productNameStr;
    char* softwareNameStr;
    char* softwareVersionStr;
    char* titleStr;

    char* numberOfScanLinesStr;
    char* pixelsPerScanLineStr;
    char* bandsPerPixelStr;
    char* n_refl_loc_str;
    char* n_cloud_phase_str;
    char* profLvlPerPixelStr;
    char* totalBandNumberStr;
    char* bandNumberStr;

    char *wavelength_3d_str;
    int product_3d_exists = 0;

    int profileProductsExist = 0;

    if (strcmp(input->metafile, "") != 0) {
        fp_meta = fopen(input->metafile, "w");
        if (fp_meta == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to open specified meta-data file, %s .\n",
                    __FILE__, __LINE__, input->metafile);
            return (HDF_FUNCTION_ERROR);
        }
    }

    /* Create the L2 file */
    ds_format_t fileFormat = DS_NCDF;
    int32 nt_chr, nt_u8, nt_i16, nt_i32, nt_f32;
    if (l2file->format == FT_L2HDF) {
        fileFormat = DS_HDF;
        nt_chr = DFNT_CHAR;
        nt_u8  = DFNT_UINT8;
        nt_i16 = DFNT_INT16;
        nt_i32 = DFNT_INT32;
        nt_f32 = DFNT_FLOAT32;

        calibrationDataStr = "Calibration Data";
        equatorCrossingLonStr = "Orbit Node Longitude";
        historyStr = "Processing Control";
        inputFilesStr = "Input Files";
        maskNamesStr = "Mask Names";
        orbitNumberStr = "Orbit Number";
        processingVersionStr = "Processing Version";
        productNameStr = "Product Name";
        softwareNameStr = "Software Name";
        softwareVersionStr = "Software Version";
        titleStr = "Title";
        numberOfScanLinesStr = "Number of Scan Lines";
        pixelsPerScanLineStr = "Pixels per Scan Line";
        bandsPerPixelStr = "Bands per Pixel";
        n_refl_loc_str = "Number of Reflectance Location Values";
        n_cloud_phase_str = "Number of Cloud Phases";
        profLvlPerPixelStr = "Profile Levels per Pixel";
        totalBandNumberStr = "total band number";
        bandNumberStr = "band number";


    } else if (l2file->format == FT_L2NCDF) {
        fileFormat = DS_NCDF;
        nt_chr = NC_CHAR;
        nt_u8  = NC_UBYTE;
        nt_i16 = NC_SHORT;
        nt_i32 = NC_INT;
        nt_f32 = NC_FLOAT;

        calibrationDataStr = "calibration_data";
        equatorCrossingLonStr = "equatorCrossingLongitude";
        historyStr = "history";
        inputFilesStr = "input_sources";
        maskNamesStr = "mask_names";
        orbitNumberStr = "orbit_number";
        processingVersionStr = "processing_version";
        productNameStr = "product_name";
        softwareNameStr = "software_name";
        softwareVersionStr = "software_version";
        titleStr = "title";
        numberOfScanLinesStr = "number_of_lines";
        pixelsPerScanLineStr = "pixels_per_line";
        bandsPerPixelStr = "bands_per_pixel";
        profLvlPerPixelStr = "profile_levels_per_pixel";
        n_refl_loc_str = "number_of_reflectance_location_values";
        n_cloud_phase_str = "number_of_cloud_phases";
        totalBandNumberStr = "number_of_bands";
        bandNumberStr = "number_of_reflective_bands";
        wavelength_3d_str = "wavelength_3d"; 

    } else if ( l2file->format == FT_SEABASSRRS) {

      seabass *priv = (seabass*)l2file->private_data;
      priv->fp = fopen(l2file->name, "w");
      
      numPixels = l2file->npix;
      numScans = l2file->nscan;
      numBands = l2file->nbands;
      numBandsIR = l2file->nbandsir;
      spix     = 0;
      cpix     = numPixels/2;
      epix     = numPixels-1;
      cscan    = numScans/2;

      bindex_set(l2file->iwave,l2file->nbands+l2file->nbandsir,BANDW);

      tot_prod = prodlist(l2file->sensorID,l1_input->evalmask,l2file->l2prod, l2file->def_l2prod,l2file->l2_prod_names);

      printf("\n\nThe following products will be included in %s.\n",l2file->name);
      for (i=0; i<tot_prod; i++)
        printf("%d %s\n",i,l2file->l2_prod_names[i]);

      l2file->tot_prod = tot_prod;

      // cache the product and procuctInfo structures
      l2file->productInfos = (productInfo_t**) allocateMemory(tot_prod * sizeof(productInfo_t*), "l2file->productInfos");

      l2file->prodptr = (l2prodstr*) allocateMemory(tot_prod * sizeof(l2prodstr), "l2file->prodptr");
      for (i=0; i<tot_prod; i++) {
            l2file->productInfos[i] = allocateProductInfo();
            if(!findProductInfo(l2file->l2_prod_names[i], l2file->sensorID, l2file->productInfos[i])) {
              fprintf(stderr, "-E- %s line %d: product %s not found.\n", __FILE__,__LINE__, l2file->l2_prod_names[i]);
              return(1);            
            }
        }
    }
    ds_id = startDS(l2file->name, fileFormat, DS_WRITE, input->deflate);
    if (ds_id.fid == FAIL) return (HDF_FUNCTION_ERROR);
    l2file->sd_id = ds_id.fid;

    /* Get number of bands, pixels, scans */
    numPixels = l2file->npix;
    numScans = l2file->nscan;
    numBands = l2file->nbands;
    numBandsIR = l2file->nbandsir;
    numLvlProf = l2file->nlvl;
    n_refl_loc = l2file->n_refl_loc;
    n_cloud_phase = l2file->n_cloud_phase;
    spix = 0;
    cpix = numPixels / 2;
    epix = numPixels - 1;
    cscan = numScans / 2;

    /* Get number of bands and band indexing from sensor table */
    rdsensorinfo(l2file->sensorID, l1_input->evalmask, "Bindx", (void **) &l2file->bindx);

    /* set wavelength index */
    rdsensorinfo(l2file->sensorID, l1_input->evalmask, "Lambda", (void **) &l2file->iwave);
    rdsensorinfo(l2file->sensorID, l1_input->evalmask, "fwave", (void **) &l2file->fwave);
    rdsensorinfo(l2file->sensorID, l1_input->evalmask, "Fobar", (void **) &l2file->Fobar);
    rdsensorinfo(l2file->sensorID, l1_input->evalmask, "Tau_r", (void **) &l2file->Tau_r);
    rdsensorinfo(l2file->sensorID, l1_input->evalmask, "k_oz", (void **) &l2file->k_oz);
    rdsensorinfo(l2file->sensorID, l1_input->evalmask, "k_no2", (void **) &l2file->k_no2);

    bindex_set(l2file->iwave, l2file->nbands + l2file->nbandsir, BANDW);

    if (l2file->aw == NULL) {
        l2file->aw = (float *) allocateMemory(numBands * sizeof (float), "l2file->aw");
    }
    if (l2file->bbw == NULL) {
        l2file->bbw = (float *) allocateMemory(numBands * sizeof (float), "l2file->bbw");
    }
    if (l2file->Fonom == NULL) {
        l2file->Fonom = (float *) allocateMemory(numBands * sizeof (float), "l2file->Fonom");
    }

    for (i = 0; i < l2file->nbands; i++) {
        l2file->aw[i] = aw_spectra(l2file->iwave[i], BANDW);
        l2file->bbw[i] = bbw_spectra(l2file->iwave[i], BANDW);
        if (l1_input->outband_opt >= 2) {
            get_f0_thuillier_ext(l2file->iwave[i], BANDW, l2file->Fonom + i);
        } else {
            l2file->Fonom[i] = l2file->Fobar[i];
        }
    }

    /*                                                                  */
    /* Build list of dataset names from input parameters                */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    tot_prod = prodlist(l2file->sensorID, l1_input->evalmask, l2file->l2prod,
            l2file->def_l2prod, l2file->l2_prod_names);

    strcpy(l2file->l2_prod_names[tot_prod++], "l2_flags");

    printf("\n\nThe following products will be included in %s.\n", l2file->name);
    for (i = 0; i < tot_prod; i++)
        printf("%d %s\n", i, l2file->l2_prod_names[i]);

    l2file->tot_prod = tot_prod;

    // cache the product and procuctInfo structures
    l2file->productInfos = (productInfo_t**) allocateMemory(tot_prod * sizeof (productInfo_t*), "l2file->productInfos");

    l2file->prodptr = (l2prodstr*) allocateMemory(tot_prod * sizeof (l2prodstr), "l2file->prodptr");
    for (i = 0; i < tot_prod; i++) {
        l2file->productInfos[i] = allocateProductInfo();
        if (!findProductInfo(l2file->l2_prod_names[i], l2file->sensorID, l2file->productInfos[i])) {
            fprintf(stderr,
                    "-E- %s line %d: product %s not found.\n",
                    __FILE__, __LINE__, l2file->l2_prod_names[i]);
            return (1);
        }

        if ((p = get_l2prod_index(l2file->l2_prod_names[i], l2file->sensorID,
                numBands + numBandsIR, numPixels, numScans, l2file->iwave)) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: product index failure.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        l2file->prodptr[i] = *p; // need to actually copy all the memory

        // see if any 3D products are requested
        if(p->rank == 3) {
            product_3d_exists = 1;
        }

    }

    // see if there are any profile products so we don't add profile stuff that is not needed
    for (i = 0; i < tot_prod; i++) {
        if( strcmp( l2file->productInfos[i]->category, "Anc_profile" ) == 0 ) {
            profileProductsExist = 1;
            break;
        }
    }

    if (l2file->format == FT_L2NCDF) {
        int dumdim;
        if (nc_def_dim(ds_id.fid, numberOfScanLinesStr, numScans, &dumdim)
                != NC_NOERR) exit(1);
        if (nc_def_dim(ds_id.fid, pixelsPerScanLineStr, numPixels, &dumdim)
                != NC_NOERR) exit(1);
        if ( nc_def_dim(ds_id.fid, bandsPerPixelStr, numBands, &dumdim)
                != NC_NOERR) exit(1);
        if (profileProductsExist) {
            if ( nc_def_dim(ds_id.fid, profLvlPerPixelStr, numLvlProf, &dumdim)
                    != NC_NOERR) exit(1);
        }
        if (nc_def_dim(ds_id.fid, n_refl_loc_str, n_refl_loc, &dumdim)
                != NC_NOERR) exit(1);
        if (nc_def_dim(ds_id.fid, n_cloud_phase_str, n_cloud_phase, &dumdim)
                != NC_NOERR) exit(1);
        if (nc_def_dim(ds_id.fid, totalBandNumberStr, numBands + numBandsIR, &dumdim)
                != NC_NOERR) exit(1);
        if (nc_def_dim(ds_id.fid, bandNumberStr, numBands, &dumdim)
                != NC_NOERR) exit(1);
        if(product_3d_exists) {
            if (nc_def_dim(ds_id.fid, wavelength_3d_str, input->nwavelengths_3d, &dumdim)
                != NC_NOERR) exit(1);
        }

        nc_def_grp(ds_id.fid, "sensor_band_parameters", &l2file->grp_id[0]);
        nc_def_grp(ds_id.fid, "scan_line_attributes", &l2file->grp_id[2]);
        nc_def_grp(ds_id.fid, "geophysical_data", &l2file->grp_id[3]);
        nc_def_grp(ds_id.fid, "navigation_data", &l2file->grp_id[4]);
        nc_def_grp(ds_id.fid, "processing_control", &l2file->grp_id[5]);
    }


    /*                                                                  */
    /* Create the scan-line datasets                                    */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[2];

    dm[0] = numScans;
    strcpy((char *) dm_name[0], numberOfScanLinesStr);
    dm[1] = numPixels;
    strcpy((char *) dm_name[1], pixelsPerScanLineStr);

    PTB(createDS(ds_id, (int) l2file->sensorID, "year", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "day", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "msec", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "time", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "detnum", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "mside", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "slon", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "clon", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "elon", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "slat", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "clat", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "elat", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "csol_z", dm, dm_name));

    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[4];

    PTB(createDS(ds_id, (int) l2file->sensorID, "longitude", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "latitude", dm, dm_name));
    PTB(createDS(ds_id, (int) l2file->sensorID, "tilt", dm, dm_name));

    /*                                                                  */
    /* Create the geophysical datasets                                  */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[3];

    for (i = 0; i < tot_prod; i++) {
        // Skip parameters already included if user requested them specifically
        if (!strcmp(l2file->l2_prod_names[i], "detnum") ||
                !strcmp(l2file->l2_prod_names[i], "mside") ||
                !strcmp(l2file->l2_prod_names[i], "year") ||
                !strcmp(l2file->l2_prod_names[i], "day") ||
                !strcmp(l2file->l2_prod_names[i], "msec") ||
                !strcmp(l2file->l2_prod_names[i], "time") ||
                !strcmp(l2file->l2_prod_names[i], "slon") ||
                !strcmp(l2file->l2_prod_names[i], "clon") ||
                !strcmp(l2file->l2_prod_names[i], "elon") ||
                !strcmp(l2file->l2_prod_names[i], "slat") ||
                !strcmp(l2file->l2_prod_names[i], "clat") ||
                !strcmp(l2file->l2_prod_names[i], "elat") ||
                !strcmp(l2file->l2_prod_names[i], "csol_z") ||
                !strcmp(l2file->l2_prod_names[i], "latitude") ||
                !strcmp(l2file->l2_prod_names[i], "longitude") ||
                !strcmp(l2file->l2_prod_names[i], "tilt")
                )
            continue;

        /* the 3rd dimension can now change - so determine it per-product */
        if( strcmp( l2file->productInfos[i]->category, "Anc_profile" ) == 0 ) {
            dm[2] = numLvlProf;
            strcpy((char *) dm_name[2], profLvlPerPixelStr );
        }else if( strcmp( l2file->productInfos[i]->category, 
           "Reflectance_loc" ) == 0 ) {
            dm[2] = n_refl_loc;
            strcpy((char *) dm_name[2], n_refl_loc_str );
        }else if( strcmp( l2file->productInfos[i]->category,
           "CTH_parameters" ) == 0 ) {
            dm[2] = n_cloud_phase;
            strcpy((char *) dm_name[2], n_cloud_phase_str );
        }else {
            dm[2] = input->nwavelengths_3d;
            strcpy((char *) dm_name[2], wavelength_3d_str);
        }
  
        PTB(createDS2(ds_id, l2file->l2_prod_names[i], l2file->productInfos[i],
                dm, dm_name));

        ds_id.sid = selectDS(ds_id, l2file->l2_prod_names[i]);

        /*                                                              */
        /* Add flag-name attributes if this is a flag product           */
        /*                                                              */
        if ((strcmp(l2file->l2_prod_names[i], "l2_flags") == 0) ||
                (strcmp(l2file->l2_prod_names[i], "flags_sst") == 0) ||
                (strcmp(l2file->l2_prod_names[i], "flags_sst3") == 0) ||
                (strcmp(l2file->l2_prod_names[i], "flags_sst4") == 0) ||
                (strcmp(l2file->l2_prod_names[i], "qual_sst") == 0) ||
                (strcmp(l2file->l2_prod_names[i], "qual_sst3") == 0) ||
                (strcmp(l2file->l2_prod_names[i], "qual_sst4") == 0) ||
                (strcmp(l2file->l2_prod_names[i], "flags_habs") == 0) ||
                (strcmp(l2file->l2_prod_names[i], "flags_mph") == 0)) {

            if ((strcmp(l2file->l2_prod_names[i], "l2_flags") == 0)) {
                tmp_str[0] = '\0';
                int32_t val = 1;
                for (k = 0; k < L1_NFLAGS; k++) {
                    // in case if the geo_mask is set. GEOREGION  = SPARE1
                    if( k == 7 && input->georegionfile[0])
                        strcat(tmp_str, "GEOREGION");
                    else
                        strcat(tmp_str, l2_flag_lname[k]);

                    flagbits[k] = val;
                    val = val << 1;
                    if (k < L1_NFLAGS - 1)
                        strcat(tmp_str, " ");
                    /*
                     * Keep the old flag attribute set for HDF4 files
                     */
                    if (l2file->format == FT_L2HDF) {
                        PTB(setAttr(ds_id, l2_flag_sname[k], nt_chr,
                                strlen(l2_flag_lname[k]) + 1, (VOIDP) l2_flag_lname[k]));
                    }
                }
                PTB(setAttr(ds_id, "flag_masks", nt_i32, L1_NFLAGS, (VOIDP) flagbits));

            } else if ((strcmp(l2file->l2_prod_names[i], "flags_sst") == 0) ||
                    (strcmp(l2file->l2_prod_names[i], "flags_sst3") == 0) ||
                    (strcmp(l2file->l2_prod_names[i], "flags_sst4") == 0)) {
                tmp_str[0] = '\0';
                int16_t val = 1;
                for (k = 0; k < NSSTFLAGS; k++) {
                    sstflagbits[k] = val;
                    val = val << 1;
                    if (l2file->sensorID == AVHRR) {
                        strcat(tmp_str, avhrr_sst_flag_lname[k]);
                    } else if ((l2file->sensorID == VIIRSN) ||
			       (l2file->sensorID == VIIRSJ1) ||
			       (l2file->sensorID == VIIRSJ2)) {
                        strcat(tmp_str, viirs_sst_flag_lname[k]);
                    } else {
                        strcat(tmp_str, sst_flag_lname[k]);
                    }
                    if (k < NSSTFLAGS - 1)
                        strcat(tmp_str, " ");
                }
                PTB(setAttr(ds_id, "flag_masks", nt_i16, NSSTFLAGS, (VOIDP) sstflagbits));

            } else if ((strcmp(l2file->l2_prod_names[i], "qual_sst") == 0) ||
                    (strcmp(l2file->l2_prod_names[i], "qual_sst3") == 0) ||
                    (strcmp(l2file->l2_prod_names[i], "qual_sst4") == 0)) {
                tmp_str[0] = '\0';
                for (k = 0; k < NQSSTFLAGS; k++) {
                    sstflagbits[k] = k;
                    strcat(tmp_str, qual_sst_flag_lname[k]);
                    if (k < NQSSTFLAGS - 1)
                        strcat(tmp_str, " ");
                }
                PTB(setAttr(ds_id, "flag_masks", nt_i16, NQSSTFLAGS, (VOIDP) sstflagbits));

            } else if ((strcmp(l2file->l2_prod_names[i], "flags_mph") == 0)) {
                tmp_str[0] = '\0';
                uint8_t val = 1;
                for (k = 0; k < NMPHFLAGS; k++) {
                    byteflagbits[k] = val;
                    val = val << 1;
                    strcat(tmp_str, mph_flag_lname[k]);
                    if (k < NMPHFLAGS - 1)
                        strcat(tmp_str, " ");
                }
                PTB(setAttr(ds_id, "flag_masks", nt_u8, NMPHFLAGS, (VOIDP) byteflagbits));
            } else if ((strcmp(l2file->l2_prod_names[i], "flags_habs") == 0)) {
                tmp_str[0] = '\0';
                uint8_t val = 1;
                for (k = 0; k < NHABSFLAGS; k++) {
                    byteflagbits[k] = val;
                    val = val << 1;
                    strcat(tmp_str, habs_flag_lname[k]);
                    if (k < NHABSFLAGS - 1)
                        strcat(tmp_str, " ");
                }
                PTB(setAttr(ds_id, "flag_masks", nt_u8, NHABSFLAGS, (VOIDP) byteflagbits));            
            }
            PTB(setAttr(ds_id, "flag_meanings", nt_chr, strlen(tmp_str) + 1, (VOIDP) tmp_str));

        }
        /*                                                              */
        /* Add solar irradiance meta-data if this is the nLw or Rrs     */
        /*                                                              */
        p = l2file->prodptr + i; // get product structure from cache
        switch (p->cat_ix) {
        case CAT_nLw:
        case CAT_Rrs:
            if(p->prod_ix != -1) {
                tmpFloat = l2file->Fonom[p->prod_ix] * 10.0;
                PTB(setAttr(ds_id, "solar_irradiance", nt_f32, 1, &tmpFloat));
            }
        }
        /*                                                            */
        /* Add bad value attributes to HDF4 files                     */
        /*                                                            */
        if ((l2file->format == FT_L2HDF) && ((strcmp(l2file->l2_prod_names[i], "l2_flags") != 0) ||
                (strcmp(l2file->l2_prod_names[i], "flags_sst") != 0) ||
                (strcmp(l2file->l2_prod_names[i], "flags_sst3") != 0) ||
                (strcmp(l2file->l2_prod_names[i], "flags_sst4") != 0))) {

            switch (p->datatype) {
            case DFNT_UINT8:
                pbuf = (VOIDP) float2uint8((VOIDP) & bad_float, 0, 1, 1, p->slope, p->offset);
                setAttr(ds_id, "bad_value_scaled", nt_chr, 1, pbuf);
                pbuf = (VOIDP) unscale_sds(pbuf, l2file->productInfos[i], 0, 1, 1);
                setAttr(ds_id, "bad_value_unscaled", nt_f32, 1, pbuf);
                break;
            case DFNT_INT16:
                pbuf = (VOIDP) float2int16((VOIDP) & bad_float, 0, 1, 1, p->slope, p->offset);
                setAttr(ds_id, "bad_value_scaled", nt_i16, 1, pbuf);
                pbuf = (VOIDP) unscale_sds(pbuf, l2file->productInfos[i], 0, 1, 1);
                setAttr(ds_id, "bad_value_unscaled", nt_f32, 1, pbuf);
                break;
            case DFNT_INT32:
                break;
            case DFNT_FLOAT32:
                setAttr(ds_id, "bad_value_scaled", nt_f32, 1, &bad_float);
                setAttr(ds_id, "bad_value_unscaled", nt_f32, 1, &bad_float);
                break;
            default:
                break;
            }
        }

        endaccessDS(ds_id);
    }

    /*                                                                  */
    /* Create the Sensor Band Parameters datasets                       */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[0];

    dm[0] = numBands + numBandsIR;
    strcpy((char *) dm_name[0], totalBandNumberStr);
    PTB(createDS(ds_id, l2file->sensorID, "wavelength", dm, dm_name));

    if(product_3d_exists) {
        dm[0] = input->nwavelengths_3d;
        strcpy((char *) dm_name[0], wavelength_3d_str);
        PTB(createDS(ds_id, l2file->sensorID, "wavelength_3d", dm, dm_name));
    }

    dm[0] = numBands;
    strcpy((char *) dm_name[0], bandNumberStr);

    if (numBands > 0) {
        PTB(createDS(ds_id, l2file->sensorID, "vcal_gain", dm, dm_name));
        PTB(createDS(ds_id, l2file->sensorID, "vcal_offset", dm, dm_name));
        PTB(createDS(ds_id, l2file->sensorID, "F0", dm, dm_name));
        PTB(createDS(ds_id, l2file->sensorID, "aw", dm, dm_name));
        PTB(createDS(ds_id, l2file->sensorID, "bbw", dm, dm_name));
        PTB(createDS(ds_id, l2file->sensorID, "k_oz", dm, dm_name));
        PTB(createDS(ds_id, l2file->sensorID, "k_no2", dm, dm_name));
        PTB(createDS(ds_id, l2file->sensorID, "Tau_r", dm, dm_name));
    }

    /*                                                                  */
    /* Write out some global attributes                                 */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if (input->suite[0])
        sprintf(title, "%s Level-2 Data %s", sensorId2SensorName(l2file->sensorID), input->suite);
    else
        sprintf(title, "%s Level-2 Data", sensorId2SensorName(l2file->sensorID));
    sprintf(soft_id, "%d.%d.%d-%s", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GITSHA);

    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->sd_id;

    PTB(SetChrGA(ds_id, titleStr, title));
    PTB(SetChrGA(ds_id, productNameStr, basename(l2file->name)));
    PTB(SetChrGA(ds_id, processingVersionStr, l1_input->pversion));
    if (l2file->orbit_node_lon > -180.0 && l2file->orbit_node_lon < 180.0)
        PTB(SetF32GA(ds_id, equatorCrossingLonStr, l2file->orbit_node_lon));
    if (l2file->orbit_number > 0)
        PTB(SetI32GA(ds_id, orbitNumberStr, l2file->orbit_number));
    PTB(SetChrGA(ds_id, historyStr, input->pro_control));

    // Processing Control attibutes
    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[5];
    PTB(SetChrGA(ds_id, softwareNameStr, PROGRAM));
    PTB(SetChrGA(ds_id, softwareVersionStr, soft_id));
    PTB(SetChrGA(ds_id, inputFilesStr, l1_input->input_files));

    // cleanup and populate calfile metadata
    char *calfile_str;
    if ((calfile_str = malloc(strlen(l1_input->calfile) + 1)) == NULL) {
        fprintf(stderr, "-E- %s line %d: Unable to copy calfile string.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    strcpy(calfile_str, "\0");
    char *tmp_calfile;
    if ((tmp_calfile = strdup(l1_input->calfile)) == NULL) {
        fprintf(stderr, "-E- %s line %d: Unable to copy calfile string.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    char *token = strtok_r(tmp_calfile, ",", &end_str);
    while (token != NULL) {
        strcpy(tmp_str, token);
        trimBlanks(tmp_str);
        strcat(calfile_str, basename(tmp_str));
        token = strtok_r(NULL, ",", &end_str);
        if (token != NULL)
            strcat(calfile_str, ", ");
    }
    free(tmp_calfile);

    PTB(SetChrGA(ds_id, calibrationDataStr, calfile_str));

    PTB(SetChrGA(ds_id, maskNamesStr, input->mask_names));
    
    if (l2file->format == FT_L2NCDF) {
        // write out netCDF specific attributes

        // global attr
        ds_id.fid = l2file->sd_id;
        PTB(SetChrGA(ds_id, "instrument", sensorId2InstrumentName(l2file->sensorID)));

        if (l2file->sensorID == AVHRR) {
            strcpy(avhrrbird, "NOAA-");
#ifdef BUILD_HISTORICAL
            strncat(avhrrbird, xsatid2name(l2file->subsensorID) + 2, 2);
#endif
            PTB(SetChrGA(ds_id, "platform", avhrrbird));
        } else {
            PTB(SetChrGA(ds_id, "platform", sensorId2PlatformName(l2file->sensorID)));
        }
        PTB(SetChrGA(ds_id, "Conventions", "CF-1.8, ACDD-1.3"));
        PTB(SetChrGA(ds_id, "license", LICENSE));
        PTB(SetChrGA(ds_id, "naming_authority", NAMING_AUTHORITY));
        // create the id -
        if (strcmp(l1_input->pversion, "Unspecified") != 0) {
            strcpy(buf1, l1_input->pversion);
            strcat(buf1, "/L2/");
        } else {
            strcpy(buf1, "L2/");
        }
        strcat(buf1, basename(l2file->name));
        PTB(SetChrGA(ds_id, "id", buf1));


        time_t tnow;
        time(&tnow);
        strcpy(buf1, unix2isodate(tnow, 'G'));
        PTB(SetChrGA(ds_id, "date_created", buf1));

        const char* keywordStr = getGCMDKeywords(input->suite);
        if(keywordStr) {
            PTB(SetChrGA(ds_id, "keywords_vocabulary", KEYWORDS_VOCABULARY));
            PTB(SetChrGA(ds_id, "keywords", keywordStr));
        }

        PTB(SetChrGA(ds_id, "standard_name_vocabulary", STDNAME_VOCABULARY));
        PTB(SetChrGA(ds_id, "institution", INSTITUTION));
        PTB(SetChrGA(ds_id, "creator_name", CREATOR_NAME));
        PTB(SetChrGA(ds_id, "creator_email", CREATOR_EMAIL));
        PTB(SetChrGA(ds_id, "creator_url", CREATOR_URL));
        PTB(SetChrGA(ds_id, "project", PROJECT));
        PTB(SetChrGA(ds_id, "publisher_name", PUBLISHER_NAME));
        PTB(SetChrGA(ds_id, "publisher_url", PUBLISHER_URL));
        PTB(SetChrGA(ds_id, "publisher_email", PUBLISHER_EMAIL));

        //Some missions have DOIs
        if(strlen(input->doi) > 0) {
            PTB(SetChrGA(ds_id, "identifier_product_doi_authority", "http://dx.doi.org"));
            PTB(SetChrGA(ds_id, "identifier_product_doi", input->doi));
        }
        
        PTB(SetChrGA(ds_id, "processing_level", "L2"));
        PTB(SetChrGA(ds_id, "cdm_data_type", "swath"));
        //        if (strlen(l2file->node_crossing_time) > 1) {
        if (l2file->node_crossing_time > 0) {
            //            double eqcrosstm = zulu2unix(l2file->node_crossing_time);
            strcpy(buf1, unix2isodate(l2file->node_crossing_time, 'G'));
            PTB(SetChrGA(ds_id, "equatorCrossingDateTime", buf1));
        }

        PTB(SetChrGA(ds_id, "spatialResolution", l2file->spatialResolution));

        if (l2file->sensorID == HAWKEYE) {
            hawkeye_t *data = (hawkeye_t*)l2file->private_data;
            PTB(SetI32GA(ds_id, "exposureID",  (data->exposureID)));

            // change to navigation group
            ds_id.fid = l2file->grp_id[4];
            PTB(SetF32GA(ds_id, "roll_offset", (data->roll_offset)));
            PTB(SetF32GA(ds_id, "time_offset", (data->time_offset)));
        }

        // Write input parameters metadata
        ds_id.fid = l2file->grp_id[5];
        int32_t grp_id_input_parms;
        nc_def_grp(l2file->grp_id[5], "input_parameters", &grp_id_input_parms);
        ds_id.fid = grp_id_input_parms;

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

    } else {

        // write out HDF4 specific attributes
        PTB(SetChrGA(ds_id, "Sensor Name", sensorId2SensorName(l2file->sensorID)));
        PTB(SetChrGA(ds_id, "Mission", sensorId2PlatformName(l2file->sensorID)));
        PTB(SetI32GA(ds_id, "Number of Bands", numBands));
        if(profileProductsExist)
            PTB(SetI32GA(ds_id, "Number of Profile Levels", numLvlProf));
        PTB(SetI32GA(ds_id, "Number of Scan Lines", numScans));
        PTB(SetI32GA(ds_id, "Pixels per Scan Line", numPixels));
        PTB(SetI32GA(ds_id, "Scene Center Scan Line", cscan));

        PTB(SetChrGA(ds_id, "Node Crossing Time", ydhmsf(l2file->node_crossing_time, 'G')));
        PTB(SetChrGA(ds_id, "Processing Time", ydhmsf(now(), 'G')));
        PTB(SetChrGA(ds_id, "Input Parameters", l1_input->input_parms));

    }

    /*                                                                  */
    /* Write out some global datasets                                   */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[0];
    PTB(writeDS(ds_id, "wavelength", l2file->iwave, 0, 0, 0, numBands + numBandsIR, 0, 0));

    if(product_3d_exists) {
        PTB(writeDS(ds_id, "wavelength_3d", input->wavelength_3d, 0, 0, 0, input->nwavelengths_3d, 0, 0));
    }

    if (numBands > 0) {
        PTB(writeDS(ds_id, "vcal_gain", l1_input->gain, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "vcal_offset", l1_input->offset, 0, 0, 0, numBands, 0, 0));

        // multiply by 10 to put into W/m2/um, since internally all radiances are mW/cm2/um
        float tmpFobar[numBands];
        for (i = 0; i < numBands; i++)
            tmpFobar[i] = l2file->Fobar[i] * 10.0;
        PTB(writeDS(ds_id, "F0", tmpFobar, 0, 0, 0, numBands, 0, 0));

        PTB(writeDS(ds_id, "aw", l2file->aw, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "bbw", l2file->bbw, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "k_oz", l2file->k_oz, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "k_no2", l2file->k_no2, 0, 0, 0, numBands, 0, 0));
        PTB(writeDS(ds_id, "Tau_r", l2file->Tau_r, 0, 0, 0, numBands, 0, 0));
    }

    free(calfile_str);

    return (LIFE_IS_GOOD);
}


/*----------------------------------------------------------------- */
/* Update scan-line datasets for the specified scan.                */

/* ---------------------------------------------------------------- */

int writel2(filehandle *l2file, int32_t recnum, l2str *l2rec, int outfile_number) {
    static int32 *buf = NULL;
    static float fsol = -1.0;
    VOIDP pbuf;
    int32 i, j;
    l2prodstr *p;

    idDS ds_id;
    static int32_t l2_flag_cnt[L1_NFLAGS];
    static int32_t sst_flag_cnt[NSSTFLAGS];
    static int32_t sst3_flag_cnt[NSSTFLAGS];
    static int32_t sst4_flag_cnt[NSSTFLAGS];
    static int32_t giop_flag_cnt[NGIOPFLAGS];
    static int32_t qualsst_flag_cnt[NQSSTFLAGS];
    static int32_t qualsst3_flag_cnt[NQSSTFLAGS];
    static int32_t qualsst4_flag_cnt[NQSSTFLAGS];
    static int32_t qaa_flag_cnt[1] = {0};
    static int32_t carder_flag_cnt[1] = {0};
    static int32_t niwa_flag_cnt[1] = {0};
    static const char *flag_lname[1] = {"PRODFAIL"};

    static uint8_t first = 1;
    static float last_lat;
    static float geobox[4][100];
    static int32 geobox_cnt = 0;
    float gring_fval[100];
    int32 gring_ival[100];

    
    if (l2file->format == FT_SEABASSRRS) {
      seabass *priv = (seabass*)l2file->private_data;
      if (first) {
        char buffer[2048];
        priv->fp = fopen(l2rec->l1rec->l1file->name, "r");
        while (1) {
          fgets(buffer, 2047, priv->fp);

          int32_t fld_cnt = 0;
          int32_t count = 0;
          char value[128];

          // Write only Rrs fieldnames to output file
          if (strncmp(buffer, "/fields=", 8) == 0) {
            char *token = strtok(&buffer[8], ",");
            memcpy(value, buffer, 8);
            value[8] = 0;
            fputs(value, priv->fp);
    
            // Loop through comma-separated fieldnames
            while (token != NULL) {
              strcpy(value, token);

              // Get index of next Rrs field 
              char *ptr =
                (char *) (l2rec->l1rec->l1file->private_data+count*sizeof(int32_t));
              int32_t fld_idx;
              memcpy(&fld_idx, ptr, sizeof(int32_t));

              // If field in data string is desired Rrs field than write to output file
              if (fld_cnt == fld_idx) {
                if (count != (l2file->nbands-1)) {
                  strcat(value, ",");
                } else {
                  if (value[strlen(value)-1] == '\n') value[strlen(value)-1] = 0;
                  strcat(value, ",");
                }
                //                printf("%s\n", value);
                fputs(value, priv->fp);
                count++;
                if (count == l2file->nbands) {
                  strcpy(value, l2file->productInfos[0]->prefix);
                  strcat(value, l2file->productInfos[0]->suffix);
                  fputs(value, priv->fp);
                  fputs("\n", priv->fp);
                }
              }
              
              token = strtok(NULL, ",");
              fld_cnt++;
            } // end while
          } else if (strncmp(buffer, "/units=", 7) == 0) {
            char *token = strtok(&buffer[8], ",");
            memcpy(value, buffer, 7);
            value[7] = 0;
            fputs(value, priv->fp);
    
            // Loop through comma-separated fieldnames
            while (token != NULL) {
              strcpy(value, token);

              // Get index of next Rrs unit 
              char *ptr = (char *) (l2rec->l1rec->l1file->private_data+count*sizeof(int32_t));
              int32_t fld_idx;
              memcpy(&fld_idx, ptr, sizeof(int32_t));
              
              // If unit in data string is desired Rrs unit than write to output file
              if (fld_cnt == fld_idx) {
                if (count != (l2file->nbands-1)) {
                  strcat(value, ",");
                } else {
                  if (value[strlen(value)-1] == '\n') value[strlen(value)-1] = 0;
                  strcat(value, ",");
                }
                fputs(value, priv->fp);
                count++;
                if (count == l2file->nbands) {
                  strcpy(value, l2file->productInfos[0]->units);
                  fputs(value, priv->fp);
                  fputs("\n", priv->fp);
                }
              }
              
              token = strtok(NULL, ",");
              fld_cnt++;
            } // end while
          } else if (strncmp(buffer, "/end_header", 11) == 0){
            fputs(buffer, priv->fp);
            break;
          } else {
            fputs(buffer, priv->fp);
          }
          
        } // while loop
        first = 0;
      } // if first
      // fclose(priv->fp);
      
      for (i=0; i<l2file->tot_prod; i++) {
        // get the product index record
        p = l2file->prodptr + i; // get product structure from cache
        // extract the product and scale (if needed)
        pbuf = prodgen(p,l2rec);
        if (p->slope != 1.0 || p->offset != 0.0) {
          pbuf = scale_sds( pbuf, l2file->productInfos[i], l2file->npix);
        }

        char buffer[2048];
        
        // if chlor_a then write Rrs value to output file
        if (strcmp(l2file->productInfos[i]->productName, "chlor_a") == 0) {
          for (j=0; j<l2file->nbands; j++) {
            sprintf(buffer, "%15.6e", l2rec->Rrs[j]);
            fputs(buffer, priv->fp);
          }
        }
        
        // write pbuf to output text file
        float f;
        memcpy(&f, pbuf, sizeof(float));
        sprintf(buffer, "%15.6f\n", f);
        fputs(buffer, priv->fp);

        return 0;
      }
    } // SeaBASS
    ds_id.deflate = 0;
    ds_id.fid = l2file->sd_id;
    if (l2file->format == FT_L2NCDF)
        ds_id.fftype = DS_NCDF;
    else
        ds_id.fftype = DS_HDF;

    if (recnum >= numScans) {
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "attempt to write rec %d of %d\n", recnum, numScans);
        return (1);
    }

    if (recnum >= cscan && fsol < 0.0) {
        fsol = l2rec->l1rec->fsol;
    }

    /* allocate buffer space */
    if (buf == NULL) {
        if ((buf = calloc(numPixels, sizeof (int32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    /* Write the scan-line data */
    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[2];
    int16_t year, day;
    double sec;

    unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);
    int32_t year32 = (int32_t) year;
    int32_t day32 = (int32_t) day;
    int32_t msec32 = (int32_t) round(sec * 1.e3);

    // find good geolocation for spix and epix
    int32_t good_spix = spix;
    int32_t good_epix = epix;

    for(int i=spix; i<epix; i++) {
        if(!l2rec->l1rec->navfail[i]) {
            good_spix = i;
            break;
        }
    }
    for(int i=epix; i>=spix; i--) {
        if(!l2rec->l1rec->navfail[i]) {
            good_epix = i;
            break;
        }
    }

    PTB(writeDS(ds_id, "year", &year32, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "day", &day32, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "msec", &msec32, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "time", &l2rec->l1rec->scantime, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "mside", &l2rec->l1rec->mside, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "detnum", &l2rec->l1rec->detnum, recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "slon", &(l2rec->l1rec->lon[good_spix]), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "clon", &(l2rec->l1rec->lon[cpix]), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "elon", &(l2rec->l1rec->lon[good_epix]), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "slat", &(l2rec->l1rec->lat[good_spix]), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "clat", &(l2rec->l1rec->lat[cpix]), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "elat", &(l2rec->l1rec->lat[good_epix]), recnum, 0, 0, 1, 1, 1));
    PTB(writeDS(ds_id, "csol_z", &(l2rec->l1rec->solz[cpix]), recnum, 0, 0, 1, 1, 1));

    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[4];
    PTB(writeDS(ds_id, "longitude", l2rec->l1rec->lon, recnum, 0, 0, 1, numPixels, 1));
    PTB(writeDS(ds_id, "latitude", l2rec->l1rec->lat, recnum, 0, 0, 1, numPixels, 1));
    PTB(writeDS(ds_id, "tilt", &(l2rec->l1rec->tilt), recnum, 0, 0, 1, 1, 1));
    
    if (outfile_number == 0) {
        if ((first == 1) || (fabs(last_lat - l2rec->l1rec->lat[cpix]) > GEOBOX_INC)
                || (recnum == (numScans - 1))) {
            // make sure the points used in the gring are valid
            if (!(l2rec->l1rec->flags[cpix] & NAVFAIL)) {

                if ((!(l2rec->l1rec->flags[good_spix] & NAVFAIL)
                        && !(l2rec->l1rec->flags[good_epix] & NAVFAIL))) {

                    first = 0;
                    geobox[0][geobox_cnt] = l2rec->l1rec->lon[good_spix];
                    geobox[1][geobox_cnt] = l2rec->l1rec->lat[good_spix];
                    geobox[2][geobox_cnt] = l2rec->l1rec->lon[good_epix];
                    geobox[3][geobox_cnt] = l2rec->l1rec->lat[good_epix];
                    last_lat = l2rec->l1rec->lat[cpix];
                    geobox_cnt++;
                }
            }
        }
        // just in case the last line is buggered...
        if ((first == 0) && (recnum < (numScans - 1))
                && (l2rec->l1rec->lat[cpix] != last_lat)) { // && (geobox_cnt < 2)) {
            if (!(l2rec->l1rec->flags[cpix] & NAVFAIL)) {

                if ((!(l2rec->l1rec->flags[good_spix] & NAVFAIL)
                        && !(l2rec->l1rec->flags[good_epix] & NAVFAIL))) {
                    geobox[0][geobox_cnt] = l2rec->l1rec->lon[good_spix];
                    geobox[1][geobox_cnt] = l2rec->l1rec->lat[good_spix];
                    geobox[2][geobox_cnt] = l2rec->l1rec->lon[good_epix];
                    geobox[3][geobox_cnt] = l2rec->l1rec->lat[good_epix];
                }
            }
        }
        if (recnum == (numScans - 1) && geobox_cnt == 1) geobox_cnt++;
    }
    /*                                                                  */
    /* Write the geophysical data                                       */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[3];
    for (i = 0; i < l2file->tot_prod; i++) {
        // Skip parameters already included if user requested them specifically
        if (!strcmp(l2file->l2_prod_names[i], "detnum") ||
                !strcmp(l2file->l2_prod_names[i], "mside") ||
                !strcmp(l2file->l2_prod_names[i], "year") ||
                !strcmp(l2file->l2_prod_names[i], "day") ||
                !strcmp(l2file->l2_prod_names[i], "msec") ||
                !strcmp(l2file->l2_prod_names[i], "time") ||
                !strcmp(l2file->l2_prod_names[i], "slon") ||
                !strcmp(l2file->l2_prod_names[i], "clon") ||
                !strcmp(l2file->l2_prod_names[i], "elon") ||
                !strcmp(l2file->l2_prod_names[i], "slat") ||
                !strcmp(l2file->l2_prod_names[i], "clat") ||
                !strcmp(l2file->l2_prod_names[i], "elat") ||
                !strcmp(l2file->l2_prod_names[i], "csol_z") ||
                !strcmp(l2file->l2_prod_names[i], "latitude") ||
                !strcmp(l2file->l2_prod_names[i], "longitude") ||
                !strcmp(l2file->l2_prod_names[i], "tilt")
                )
            continue;

        // get the product index record
        p = l2file->prodptr + i; // get product structure from cache

        // find the correct 3rd dimension
        int32_t num_3d;
        if(p->rank == 3) {
            // use another 3rd dim for the anc profiles
            num_3d = ( strcmp( l2file->productInfos[i]->category, 
                "Anc_profile" ) == 0 ) ?  numLvlProf : input->nwavelengths_3d;
            // for the reflectance diagnostics
            num_3d = ( strcmp( l2file->productInfos[i]->category, 
                "Reflectance_loc" ) == 0 ) ?  n_refl_loc : num_3d;
            // for the cloud phase
            num_3d = ( strcmp( l2file->productInfos[i]->category,
                "CTH_parameters" ) == 0 ) ?  n_cloud_phase : num_3d;
        }

        // extract the product and scale (if needed)
        pbuf = prodgen(p, l2rec);
        if (p->slope != 1.0 || p->offset != 0.0) {
            if(p->rank == 3)
                pbuf = scale_sds(pbuf, l2file->productInfos[i], l2file->npix * num_3d);
            else
                pbuf = scale_sds(pbuf, l2file->productInfos[i], l2file->npix);
        }
        /* update flag counters when appropriate */
        if ((strcmp(l2file->l2_prod_names[i], "flags_sst") == 0)) {
            update_flag_cnts16(sst_flag_cnt, pbuf, NSSTFLAGS, l2rec->l1rec->npix, 1L);
        } else if ((strcmp(l2file->l2_prod_names[i], "flags_sst4") == 0)) {
            update_flag_cnts16(sst4_flag_cnt, pbuf, NSSTFLAGS, l2rec->l1rec->npix, 1L);
        } else if ((strcmp(l2file->l2_prod_names[i], "flags_sst3") == 0)) {
            update_flag_cnts16(sst3_flag_cnt, pbuf, NSSTFLAGS, l2rec->l1rec->npix, 1L);
        } else if ((strcmp(l2file->l2_prod_names[i], "flags_giop") == 0)) {
            update_flag_cnts(giop_flag_cnt, pbuf, NGIOPFLAGS, l2rec->l1rec->npix, 1L);
        } else if ((strcmp(l2file->l2_prod_names[i], "flags_qaa") == 0)) {
            update_flag_cnts(qaa_flag_cnt, pbuf, 1, l2rec->l1rec->npix, PRODFAIL);
        } else if ((strcmp(l2file->l2_prod_names[i], "flags_carder") == 0)) {
            update_flag_cnts(carder_flag_cnt, pbuf, 1, l2rec->l1rec->npix, PRODFAIL);
        } else if ((strcmp(l2file->l2_prod_names[i], "flags_niwa") == 0)) {
            update_flag_cnts(niwa_flag_cnt, pbuf, 1, l2rec->l1rec->npix, PRODFAIL);
        } else if ((strcmp(l2file->l2_prod_names[i], "qual_sst") == 0)) {
            update_qual_cnts(qualsst_flag_cnt, pbuf, NQSSTFLAGS, l2rec->l1rec->npix);
        } else if ((strcmp(l2file->l2_prod_names[i], "qual_sst3") == 0)) {
            update_qual_cnts(qualsst3_flag_cnt, pbuf, NQSSTFLAGS, l2rec->l1rec->npix);
        } else if ((strcmp(l2file->l2_prod_names[i], "qual_sst4") == 0)) {
            update_qual_cnts(qualsst4_flag_cnt, pbuf, NQSSTFLAGS, l2rec->l1rec->npix);
        }

        // write to L2 file
        if (p->rank == 3) {
            PTB( writeDS(ds_id,l2file->l2_prod_names[i],pbuf,
                    recnum,0,0,1,numPixels,num_3d) );
        } else if (p->rank == 2) {
            PTB(writeDS(ds_id, l2file->l2_prod_names[i], pbuf,
                    recnum, 0, 0, 1, numPixels, 1));
        } else {
            PTB(writeDS(ds_id, l2file->l2_prod_names[i], pbuf,
                    recnum, 0, 0, 1, 1, 1));
        }
    }

    /* Update global flag counter */
    update_flag_cnts(l2_flag_cnt, l2rec->l1rec->flags, L1_NFLAGS, l2rec->l1rec->npix, 1L);

    if ((l2file->format == FT_L2NCDF) && !(recnum % 200)) {
        nc_sync(l2file->sd_id);
    }

    /* Write global attributes */
    if (recnum == (numScans - 1)) {
        float flag_perc[L1_NFLAGS];
        if (l2file->format == FT_L2NCDF) {

            // write out to netCDF file
            ds_id.fid = l2file->sd_id;

            scene_meta_write(ds_id);

            PTB(SetF64GA(ds_id, "earth_sun_distance_correction", fsol));

            // Write flag percentages metadata
            int32_t grp_id_flag_percentages;
            ds_id.sid = NC_GLOBAL;

            /* Report flag percentages */
            /* determine if there are any flag products */
            for (i = 0; i < l2file->tot_prod; i++) {
                if ((strcmp(l2file->l2_prod_names[i], "flags_sst") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "sst_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nSST: Percentage of pixels flagged:\n");
                    if (l2file->sensorID == AVHRR)
                        write_flag_pcnts(ds_id, fp_meta, sst_flag_cnt, NSSTFLAGS, avhrr_sst_flag_lname, numScans, numPixels);
                    else if ((l2file->sensorID == VIIRSN) ||
			     (l2file->sensorID == VIIRSJ1) ||
			     (l2file->sensorID == VIIRSJ2))
                        write_flag_pcnts(ds_id, fp_meta, sst_flag_cnt, NSSTFLAGS, viirs_sst_flag_lname, numScans, numPixels);
                    else
                        write_flag_pcnts(ds_id, fp_meta, sst_flag_cnt, NSSTFLAGS, sst_flag_lname, numScans, numPixels);
                } else if ((strcmp(l2file->l2_prod_names[i], "flags_sst4") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "sst4_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nSST4: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, sst4_flag_cnt, NSSTFLAGS, sst_flag_lname, numScans, numPixels);
                } else if ((strcmp(l2file->l2_prod_names[i], "flags_sst3") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "sst3_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nSST3: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, sst3_flag_cnt, NSSTFLAGS, viirs_sst_flag_lname, numScans, numPixels);
                } else if ((strcmp(l2file->l2_prod_names[i], "flags_giop") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "giop_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nGIOP: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, giop_flag_cnt, NGIOPFLAGS, giop_flag_lname, numScans, numPixels);
                } else if ((strcmp(l2file->l2_prod_names[i], "flags_qaa") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "qaa_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nQAA: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, qaa_flag_cnt, 1, flag_lname, numScans, numPixels);
                } else if ((strcmp(l2file->l2_prod_names[i], "flags_carder") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "carder_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nCARDER: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, carder_flag_cnt, 1, flag_lname, numScans, numPixels);
                } else if ((strcmp(l2file->l2_prod_names[i], "flags_niwa") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "niwa_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nNIWA: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, niwa_flag_cnt, 1, flag_lname, numScans, numPixels);
                } else if ((strcmp(l2file->l2_prod_names[i], "qual_sst") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "qual_sst_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nQUAL_SST: Percentage of pixels flagged:\n");
                    write_qual_flag_pcnts(ds_id, fp_meta, qualsst_flag_cnt, NQSSTFLAGS, qual_sst_flag_lname);
                } else if ((strcmp(l2file->l2_prod_names[i], "qual_sst4") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "qual_sst4_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nQUAL_SST4: Percentage of pixels flagged:\n");
                    write_qual_flag_pcnts(ds_id, fp_meta, qualsst4_flag_cnt, NQSSTFLAGS, qual_sst_flag_lname);
                } else if ((strcmp(l2file->l2_prod_names[i], "qual_sst3") == 0)) {
                    nc_def_grp(l2file->grp_id[5], "qual_sst3_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nQUAL_SST3: Percentage of pixels flagged:\n");
                    write_qual_flag_pcnts(ds_id, fp_meta, qualsst3_flag_cnt, NQSSTFLAGS, qual_sst_flag_lname);
                }

            }
            nc_def_grp(l2file->grp_id[5], "flag_percentages", &grp_id_flag_percentages);
            ds_id.fid = grp_id_flag_percentages;
            printf("\nPercentage of pixels flagged:\n");
            write_flag_pcnts(ds_id, fp_meta, l2_flag_cnt, L1_NFLAGS, l2_flag_lname, numScans, numPixels);


            // Geobox attributes
            if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->sd_id;
            j = 1;
            gring_fval[0] = geobox[0][0];
            for (i = 0; i < geobox_cnt; i++) {
                gring_fval[j++] = geobox[2][i];
            }
            for (i = 0; i < geobox_cnt - 1; i++) {
                gring_fval[j++] = geobox[0][geobox_cnt - 1 - i];
            }
            if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[4];
            PTB(setAttr(ds_id, "gringpointlongitude", NC_FLOAT, j, (VOIDP) gring_fval));

            j = 1;
            gring_fval[0] = geobox[1][0];
            gring_ival[0] = j;
            for (i = 0; i < geobox_cnt; i++) {
                gring_ival[j] = j + 1;
                gring_fval[j++] = geobox[3][i];
            }
            for (i = 0; i < geobox_cnt - 1; i++) {
                gring_ival[j] = j + 1;
                gring_fval[j++] = geobox[1][geobox_cnt - 1 - i];
            }
            if (l2file->format == FT_L2NCDF) ds_id.fid = l2file->grp_id[4];
            PTB(setAttr(ds_id, "gringpointlatitude", NC_FLOAT, j, (VOIDP) gring_fval));
            PTB(setAttr(ds_id, "gringpointsequence", NC_INT, j, (VOIDP) gring_ival));

        } else {

            // write out to HDF4
            scene_meta_write(ds_id);

            PTB(SetF64GA(ds_id, "Earth-Sun Distance Correction", fsol));

            /* Report flag percentages */
            printf("\nPercentage of pixels flagged:\n");
            for (i = 0; i < L1_NFLAGS; i++) {
                flag_perc[i] = ((float) l2_flag_cnt[i]) / numScans / numPixels * 100.0;
                printf("Flag #%2d: %16s %10d %8.4f\n",
                        i + 1, l2_flag_lname[i], l2_flag_cnt[i], flag_perc[i]);
                if (fp_meta != NULL)
                    fprintf(fp_meta, "Flag #%2d: %16s %10d %8.4f\n",
                        i + 1, l2_flag_lname[i], l2_flag_cnt[i], flag_perc[i]);
            }
            PTB(sd_setattr(ds_id.fid, "Flag Percentages", DFNT_FLOAT32,
                    L1_NFLAGS, (VOIDP) flag_perc));

        }

    }

    return (LIFE_IS_GOOD);
}

/* -------------------------------------------------------- */
/* Finish access for the current file.                      */

/* -------------------------------------------------------- */
int closel2(filehandle *l2file) {
    idDS ds_id;

    if ( l2file->format == FT_SEABASSRRS) {
      seabass *priv = (seabass*)l2file->private_data;
      fclose(priv->fp);
      return(LIFE_IS_GOOD);
    }
    ds_id.deflate = 0;
    ds_id.fid = l2file->sd_id;
    if (l2file->format == FT_L2NCDF)
        ds_id.fftype = DS_NCDF;
    else
        ds_id.fftype = DS_HDF;

    if (l2file->format == FT_L2HDF) {
        PTB(MakeVgroups(l2file));
    }

    if (endDS(ds_id)) {
        fprintf(stderr, "-E- %s line %d: endDS(%d) failed for file, %s.\n",
                __FILE__, __LINE__, ds_id.fid, l2file->name);
        return (HDF_FUNCTION_ERROR);
    }

    free(l2file->aw);
    l2file->aw = NULL;
    free(l2file->bbw);
    l2file->bbw = NULL;
    free(l2file->Fonom);
    l2file->Fonom = NULL;
    
    int i;
    for(i=0; i<l2file->tot_prod; i++)
        freeProductInfo(l2file->productInfos[i]);
    free(l2file->productInfos);
    l2file->productInfos = NULL;
    free(l2file->prodptr);
    l2file->prodptr = NULL;
    
    if (fp_meta != NULL)
        fclose(fp_meta);
    fp_meta = NULL;
    
    return (LIFE_IS_GOOD);
}

void update_flag_cnts(int32_t *flag_cnt, int32_t*flags, int32_t nflags, int32_t npix, uint32_t init_mask) {
    int32_t i, ip;
    uint32_t mask;

    /* Update flag counter */
    for (ip = 0; ip < npix; ip++) {
        mask = init_mask;
        for (i = 0; i < nflags; i++) {
            flag_cnt[i] += ((flags[ip] & mask) > 0);
            mask *= 2L;
        }
    }
    return;
}

void update_flag_cnts16(int32_t *flag_cnt, int16_t *flags, int32_t nflags, int32_t npix, uint32_t init_mask) {
    int32_t i, ip;
    uint32_t mask;

    /* Update flag counter */
    for (ip = 0; ip < npix; ip++) {
        mask = init_mask;
        for (i = 0; i < nflags; i++) {
            flag_cnt[i] += ((flags[ip] & mask) > 0);
            mask *= 2L;
        }
    }
    return;
}

void update_qual_cnts(int32_t *flag_cnt, int8_t* flags, int32_t nflags, int32_t npix) {
    int32_t i, ip;
    uint32_t mask;

    /* Update flag counter */
    for (ip = 0; ip < npix; ip++) {
        mask = 0;
        for (i = 0; i < nflags; i++) {
            flag_cnt[i] += ((flags[ip] == mask));
            mask++;
        }
    }
    return;
}

int write_flag_pcnts(idDS ds_id, FILE *fpmeta, int32_t *flag_cnt, int32_t nflags, const char * const flag_lname[], int32_t numScans, int32_t numPixels) {
    int32_t i;
    float *flag_perc;

    flag_perc = (float *) calloc(nflags, sizeof (float));
    char tmp_str[256];
    /* Report flag percentages */
    for (i = 0; i < nflags; i++) {
        if( i == 7 && input->georegionfile[0])
            strcpy(tmp_str, "GEOREGION");
        else
            strcpy(tmp_str,flag_lname[i]);
        flag_perc[i] = ((float) flag_cnt[i]) / numScans / numPixels * 100.0;
        printf("Flag #%2d: %16s %10d %8.4f\n", i + 1, tmp_str, flag_cnt[i], flag_perc[i]);
        if (fp_meta != NULL)
            fprintf(fpmeta, "Flag #%2d: %16s %10d %8.4f\n", i + 1, tmp_str, flag_cnt[i], flag_perc[i]);

        PTB(setAttr(ds_id, tmp_str, NC_FLOAT, 1, (VOIDP) (flag_perc + i)));
    }

    free(flag_perc);
    return 0;
}

int write_qual_flag_pcnts(idDS ds_id, FILE *fpmeta, int32_t *flag_cnt, int32_t nflags, const char * const flag_lname[]) {
    int32_t i, sumflags = 0;
    float *flag_perc;

    flag_perc = (float *) calloc(nflags, sizeof (float));

    for (i = 0; i < nflags; i++) sumflags += flag_cnt[i];

    /* Report flag percentages */
    for (i = 0; i < nflags; i++) {
        flag_perc[i] = ((float) flag_cnt[i]) / sumflags * 100.0;
        printf("Flag #%2d: %16s %10d %8.4f\n", i + 1, flag_lname[i], flag_cnt[i], flag_perc[i]);
        if (fp_meta != NULL)
            fprintf(fpmeta, "Flag #%2d: %16s %10d %8.4f\n", i + 1, flag_lname[i], flag_cnt[i], flag_perc[i]);

        PTB(setAttr(ds_id, flag_lname[i], NC_FLOAT, 1, (VOIDP) (flag_perc + i)));
    }

    return 0;
}
