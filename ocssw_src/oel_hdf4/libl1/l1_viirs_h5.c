/* ========================================================================= */
/* module l1_viirs_h5.c - functions to read VIIRS L1B for MSL12              */
/*                                                                           */
/* Written By: W. Robinson, SAIC, Feb, 2009 based on l1_modis_hdf.c          */
/*   W. Robinson, SAIC, 13 Nov 2012  add reading of height from the GEO file */
/*   W. Robinson, SAIC, 9 Apr 2013  add capability to process                */
/*      time-aggregated files                                                */
/*   W. Robinson, SAIC, 5 Jun 2014  fix mirror normal used in polarization calc */
/*                                                                           */
/* ========================================================================= */

#include "libnav.h"
#include "l1_viirs_h5.h"
#include "l1.h"
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include <libgen.h>
#include <ctype.h>


#define MAXBANDS 16
#define STDMAXSCANS 48  /* for std granule - max # scans */
#define NDET     16
#define N_MS     2
#define NAGGPX   3200
#define NUMERATOR 0
#define DENOMINATOR 1
#define near( x, y ) ( fabsf( x - y ) < 0.001 )

/*  define fill condition values  */
#define NA_FLOAT32_FILL   -999.9f  /* Algorithm Exclusions  */
#define NA_UINT16_FILL  65535
#define MISS_FLOAT32_FILL  -999.8f  /* Missing at Time of Processing  */
#define MISS_UINT16_FILL  65534
#define ONBOARD_PT_FLOAT32_FILL  -999.7f  /*  Onboard Pixel Trimi (bow tie) */
#define ONBOARD_PT_UINT16_FILL  65533
#define ONGROUND_PT_FLOAT32_FILL  -999.6f  /* On-ground Pixel Trim */
#define ONGROUND_PT_UINT16_FILL  65532
#define ERR_FLOAT32_FILL  -999.5f  /* Cannot Calculate */
#define ERR_UINT16_FILL  65531
#define ELINT_FLOAT32_FILL  -999.4f  /* Ellipsoid Intersection Failed */
#define ELINT_UINT16_FILL  65530
#define VDNE_FLOAT32_FILL  -999.3f  /* Value Does Not Exist */
#define VDNE_UINT16_FILL  65529
#define SOUB_FLOAT32_FILL  -999.2f  /* Scaling Out Of Bounds */
#define SOUB_UINT16_FILL  65528

int gen_sdr_suite(char *);
int set_f_cal_corr(h5io_str *, filehandle *, int64_t);
int rd_vir_f_tbl(char *, int64_t, int);

/*  all the band file names, file IDs for geo and bands and needed 
 group and dataset ids  */
static char vfnames[MAXBANDS][FILENAME_MAX];
static h5io_str geo_fid, sdr_fid[MAXBANDS];
static h5io_str geo_dat_id[7], bnd_dat_id[2][MAXBANDS];
static char sdr_band_typ[MAXBANDS]; /* 0 if scaled, 1 if float */
static unsigned short *scl_rad; /*  storage for the scaled radiance values */
static float *scale, *offset; /* scaling factors */
static float *flt_rad; /* for float radiance read from file */
static int64_t *u58_scn_st, u58_gran_st, u58_gran_en;
static float *Fobar;
static unsigned char *viirs_qual2; /* storage for the SDR scan quality info */
/*
 *  the 2 ..._map items track the locations of scans that have real data in 
 *  them, to account for the scans in the datasets that don't have instrument
 *  data in them and get concatenated together in an time aggregate
 *  the gran_map tells what granule the scan is in.
 */
static int *gran_map; /* granule assoc to a scan */
static int *scan_map; /* real dataset scan location */
/* from the QF2_SCAN_SDR, get mirror side
 from here too */
static float *pos, *vel, *att; /* position, velocity, sensor attitude */
static int16_t scn_fmt; /* local scan format store */
static int16_t margin[2]; /* margin storage */

#define VSWIR  0
#define THERM  1
#define CIRRUS 2
static int32_t btype[] = {0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 1, 1, 1};
static int32_t ivswir = 0;
static int32_t itherm = 0;
static double f_cal_corr[MAXBANDS * NDET * N_MS]; /* f table correction */

void fix_viirs_filename(char* name) {
    char tmp[FILENAME_MAX];
    char dir[FILENAME_MAX];
    char base[FILENAME_MAX];
    
    // see if file exists
    if( access(name, R_OK ) == 0 ) {
        return;
    }
    
    // see if we can find a file that matches the first 44 char of the file name
    // This will happen if the creation time is the last part of the filename.
    strcpy(tmp, name);
    strcpy(base, basename(tmp));
    strcpy(tmp, name);
    strcpy(dir, dirname(tmp));
    
    DIR *dp;
    struct dirent *entry;
    if((dp = opendir(dir)) == NULL) {
        fprintf(stderr, "%s:%d - cannot open directory: %s\n", __FILE__, __LINE__, dir);
        exit(EXIT_FAILURE);
    }

    int fileFound = 0;
    while((entry = readdir(dp)) != NULL) {
        if(strncmp(entry->d_name, base, 44) == 0){
            fileFound = 1;
            sprintf(name, "%s/%s", dir, entry->d_name);
            break;
        }
    }
    closedir(dp);

    if(!fileFound) {
        fprintf(stderr, "%s:%d - valid VIIRS file not found: %s\n", __FILE__, __LINE__, name);
        exit(EXIT_FAILURE);
    }
    
}

int gen_sdr_suite(char *in_file)
/*******************************************************************

 gen_sdr_suite

 purpose: check the initial L1b sdr name to conform to VIIRS format and
 generate the full SDR name suite from that 1 name

 Returns type: int - description

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 char *            in_file          I      input VIIRS file with path

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson       9-Feb-2009      Original development
 W. Robinson, SAIC 6Aug2010        add ability to read suite from text file

 *******************************************************************/ {
    char *path, temp[FILENAME_MAX], base[FILENAME_MAX];
    int i, j;
    h5io_str h5fid;
    FILE *fp;

    /*
     *  if the file is a text file with names (not hdf 5) extract file names
     */
    if (h5io_openr(in_file, 0, &h5fid) != 0) {
        if ((fp = fopen(in_file, "r")) == NULL) {
            fprintf(stderr, "-E- %s, %d: input file %s open problem.\n",
                    __FILE__, __LINE__, in_file);
            return 1;
        }
        /*
         *  get the file names of needed bands
         */
        for (i = 0, j = 1; i < MAXBANDS; i++, j++) {
            if (fscanf(fp, "%s", vfnames[i]) != 1) {
                fprintf(stderr,
                        "-E- %s, %d: Failed to read VIIRS input file # %d from  list file: %s/n",
                        __FILE__, __LINE__, i, in_file);
                return 1;
            }
        }
        fclose(fp);
    } else {
        h5io_close(&h5fid);
        /*
         *  Although standard VIIRS file names have much more, just accept
         *  SDR files with the base starting with SVMxx (xx is band #)
         */
        strcpy(temp, in_file);
        strcpy(base, basename(temp));
        strcpy(temp, in_file);
        path = dirname(in_file);
        /*  */
        if (strncmp(base, "SVM", 3) != 0) {
            fprintf(stderr,
                    "-E- %s, %d: Improper VIIRS file format for file base:\n",
                    __FILE__, __LINE__);
            fprintf(stderr, "%s\n", base);
            return 1;
        }
        if (!(isdigit(*(base + 3))) || !(isdigit(*(base + 4)))) {
            fprintf(stderr,
                    "-E- %s, %d: Improper VIIRS file format for file base:\n",
                    __FILE__, __LINE__);
            fprintf(stderr, "%s\n", base);
            return 1;
        }
        for (i = 0, j = 1; i < MAXBANDS; i++, j++) {
            sprintf(vfnames[i], "%s/SVM%2.2d%s", path, j, base + 5);
            fix_viirs_filename(vfnames[i]);
            if (want_verbose)
                printf("VIIRS filename for M%2.2d: %s\n", j, vfnames[i]);
        }
    }
    return (LIFE_IS_GOOD);
}

/* ------------------------------------------------------------------------
 set_f_cal_corr

 purpose: create the array of f table corrections to apply to the radiances
 using the input cal table and the SDR indicated cal table

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 h5io_str *        g_id             I      group id of the area with
 the granule metadata
 (specifically the tables used)
 filehandle *      file             I      input L1 file info, including
 processing options
 int64_t           u58_time         I      start time of the granule in
 IET units

 (the array of f table corrections, f_cal_corr is shared in the
 l1_viirs_h5.c routine)

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson, SAIC 11 Sep 2012     add the ability to process pseudo-L1A
 SDR files

 ------------------------------------------------------------------------*/
int set_f_cal_corr(h5io_str *g_id, filehandle *file, int64_t u58_time) {
    char *in_calfile, sdr_cal_fil[FILENAME_MAX], *csave = NULL;
    char *aux_arr, *fildir;
    int ncorr, icorr, f_ndim, f_dim_siz[20], f_sto_len, i;
    H5T_class_t f_class;
    hid_t f_native_typ;
    /*
     *  the f_cal_corr is used always to add any radiance correction.
     *  These values are only non-unity if a cal table is specified
     *
     *  set up default unity table
     */
    ncorr = N_MS * NDET * MAXBANDS;
    for (icorr = 0; icorr < ncorr; icorr++)
        *(f_cal_corr + icorr) = 1.;

    in_calfile = l1_input->calfile;
    /*
     *  If no cal table is specified, it is done
     */
    if (want_verbose)
        printf("\nVIIRS calibration starting\n");
    if (*in_calfile == 0) {
        if (want_verbose)
            printf("No cal file supplied, calibration remains unchanged\n");
        return 0;
    }
    if (want_verbose)
        printf("Input calibration from file: %s\n", in_calfile);
    /*
     *  The numerator of the f table correction comes from the input cal file, but
     *  is one if the table is 'Unity'
     *  otherwise, make the numerator of the correction from the specified f table
     *  If any cal file is specified, it will back out the existing cal and use
     *  the new cal.  If the existing cal was 'Unity' it will not have to back
     *  out the cal.
     */
    if (strstr(in_calfile, "Unity") == NULL) {
        if (rd_vir_f_tbl(in_calfile, u58_time, NUMERATOR) != 0) {
            printf("%s, %d - E - failed to retrieve input F table values\n",
                    __FILE__, __LINE__);
            return 1;
        }
    } else if (want_verbose)
        printf("Input calibration is from the unity F table, Using unity\n");
    /*
     *  any cal other than unity in the SDR is removed by placing it into the
     *  denominator of the correction array
     */
    /*
     *  get some info on the attribute
     */
    if (h5io_info(g_id, "N_Aux_Filename", &f_class, &f_native_typ, &f_ndim,
            f_dim_siz, &f_sto_len) != 0) {
        printf("%s, %d - E - could not get info on N_Aux_Filename\n", __FILE__,
                __LINE__);
        return 1;
    }
    /*
     *  allocate space for the times
     */
    if ((aux_arr = (char *) malloc(f_sto_len * f_dim_siz[0] * sizeof (char)))
            == NULL) {
        printf("%s, %d - E - unable to allocate aux filename list array\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  read the attribute
     */
    if (h5io_rd_attr(g_id, "N_Aux_Filename", aux_arr) != 0) {
        printf("%s, %d - failed to read the N_Aux_Filename attribute\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  find the F LUT in the strings
     */
    for (i = 0; i < f_dim_siz[0]; i++) {
        /* printf( "# %3d: %s\n", i, ( aux_arr + i * f_sto_len ) ); */
        if (strstr((aux_arr + i * f_sto_len), "VIIRS-SDR-F-LUT_npp") != NULL)
            csave = (aux_arr + i * f_sto_len);
    }
    if (csave == NULL) {
        /*
         *  there can be another, dynamic F-LUT now.  If it is found, we cannot
         *  deal with that now.  Advise user to NOT use a F-LUT to calibrate with
         */
        printf("%s ,%d: I: No F-LUT found, Looking for F-PREDICTED-LUT\n",
                __FILE__, __LINE__);
        for (i = 0; i < f_dim_siz[0]; i++) {
            /*printf( "# %3d: %s\n", i, ( aux_arr + i * f_sto_len ) ); */
            if (strstr((aux_arr + i * f_sto_len),
                    "VIIRS-SDR-F-PREDICTED-LUT_npp") != NULL)
                csave = (aux_arr + i * f_sto_len);
        }
        if (csave == NULL) {
            printf(
                    "%s, %d E - failed to find the F LUT name in N_Aux_Filename\n",
                    __FILE__, __LINE__);
            return 1;
        } else {
            printf(
                    "%s, %d - A F-PREDICTED-LUT was found, but code is not set-up\n",
                    __FILE__, __LINE__);
            printf(
                    "   to handle this kind of LUT.  Set 'calfile=' in par file\n");
            printf("   to leave cal unchanged\n");
            return 1;
        }
    }
    if (want_verbose)
        printf("\n\nThe F LUT found in the SDR is: %s\n", csave);
    /*
     *  look for a Unity file
     */
    if (strstr(csave, "Unity") != NULL) {
        if (want_verbose)
            printf("Found a unity F table designation in the SDR\n");
    } else {
        /*
         *  the name of the file should be what was found with a .h5 after it
         *  in path $OCVARROOT/viirsn/cal/EVAL/
         */
        if ((fildir = getenv("OCVARROOT")) == NULL) {
            printf("-E- %s, %d: OCVARROOT env variable undefined.\n", __FILE__,
                    __LINE__);
            exit(1);
        }
        strcpy(sdr_cal_fil, fildir);
        strcat(sdr_cal_fil, "/viirsn/cal/EVAL/");
        strcat(sdr_cal_fil, csave);
        strcat(sdr_cal_fil, ".h5");
        printf("F Lut file name derived from the SDR is: %s\n", sdr_cal_fil);
        /*
         *  Only for non-unity table, read the f table values into the correction
         *  array denominator
         */
        if (rd_vir_f_tbl(sdr_cal_fil, u58_time, DENOMINATOR) != 0) {
            printf("-E- %s, %d: Failed to retrieve SDR F table values\n",
                    __FILE__, __LINE__);
            return 1;
        }
    }
    /*
     *  The corrections are all set to be applied
     */
    if (want_verbose)
        printf("Completed calibration using F LUT file(s)\n\n");
    free(aux_arr);
    return 0;
}

/* ------------------------------------------------------------------------

 rd_vir_f_tbl

 purpose: read F table lut high gain corrections and apply them to the
 correction array for the specified time of the viirs granule

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 char *            file             I      F LUT table file name
 int64_t           u58_time         I      Granule time
 int               corr_loc         I      either NUMERATOR or DENOMINATOR
 to apply f table values to the
 proper part of the correction
 array

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson, SAIC 11 Sep 2012     Original development
 W. Robinson, SAIC 16 Nov 2012     mod to extrapolate using end 2 points

 ------------------------------------------------------------------------*/
int rd_vir_f_tbl(char *file, int64_t u58_time, int corr_loc) {
    /* f_lut_sel [ # ms, # det, # bands, # times ] (ms fastest)  */
    double f_lut_sel[2 * 16 * 16 * 2];
    int64_t *tarr, del, dist, u58_t1, u58_t2;
    int start[5], count[5];
    H5T_class_t f_class;
    hid_t f_native_typ;
    h5io_str h5fid, dsid;
    int f_ndim, f_dim_siz[20], f_sto_len, ntime, itim1, ncorr, i;
    float fact;
    /*
     *  open the LUT file and get the times available
     */
    if (want_verbose)
        printf("Reading in F LUT: %s\n", file);
    if (h5io_openr(file, 0, &h5fid) != 0) {
        printf("- E - %s, %d: Unable to open the F LUT file: %s\n", __FILE__,
                __LINE__, file);
        return 1;
    }
    /*
     *  get the info on the times
     */
    if (h5io_set_ds(&h5fid, "Beginning_Time_IET", &dsid) != 0) {
        printf("- E - %s, %d: Unable to set dataset Beginning_Time_IET\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_info(&dsid, NULL, &f_class, &f_native_typ, &f_ndim, f_dim_siz,
            &f_sto_len) != 0) {
        printf("- E - %s, %d: Unable to get info on Beginning_Time_IET\n",
                __FILE__, __LINE__);
        return 1;
    }
    ntime = f_dim_siz[0];
    /*
     *  allocate time array and read it
     */
    if ((tarr = (int64_t *) malloc(ntime * sizeof (int64_t))) == NULL) {
        printf("- E - %s, %d: Unable to allocate time storage\n", __FILE__,
                __LINE__);
        return 1;
    }
    if (h5io_rd_ds(&dsid, (void *) tarr) != 0) {
        printf("- E - %s, %d: Unable to read Beginning_Time_IET\n", __FILE__,
                __LINE__);
        return 1;
    }
    h5io_close(&dsid);
    /*
     *  find the 2 times to interpolate between / extrapolate beyond
     */
    if (want_verbose)
        printf("Granule time is: %ld\n", (long) u58_time);
    if (u58_time < tarr[0]) {
        printf(
                "\n\n\n- W - %s, %d: WARNING, granule time below times in the F table\n",
                __FILE__, __LINE__);
        itim1 = 0;
        printf("granule time of %ld is below table start time of %ld\n",
                (long) u58_time, (long) tarr[0]);
    } else if (u58_time > *(tarr + ntime - 1)) {
        printf(
                "\n\n\n- W - %s, %d: WARNING, granule time higher than times in the F table\n",
                __FILE__, __LINE__);
        itim1 = ntime - 2;
        printf("granule time of %ld is above table end time of %ld\n", (long) u58_time,
                (long) tarr[ntime - 1]);
    } else {
        for (itim1 = 0; itim1 < ntime; itim1++)
            if (*(tarr + itim1 + 1) > u58_time)
                break;
    }
    /*
     *  get the 2 times
     */
    u58_t1 = *(tarr + itim1);
    u58_t2 = *(tarr + itim1 + 1);
    if (want_verbose)
        printf("itim1: %d, t1: %ld, t2: %ld\n", itim1, (long) u58_t1, (long) u58_t2);
    free(tarr);
    /*
     *  flut is 2 ms, 3 gain, 32 det, 22 bands, # times, check this
     */
    if (h5io_set_ds(&h5fid, "VIIRS-SDR-F-LUT", &dsid) != 0) {
        printf("- E - %s, %d: Unable to set dataset VIIRS-SDR-F-LUT\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_info(&dsid, NULL, &f_class, &f_native_typ, &f_ndim, f_dim_siz,
            &f_sto_len) != 0) {
        printf("- E - %s, %d: Unable to get info on dataset VIIRS-SDR-F-LUT\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (f_ndim != 5) {
        printf("- E - %s, %d: # dimensions of LUT dataset is not 5\n", __FILE__,
                __LINE__);
        return 1;
    }
    if ((f_dim_siz[0] != ntime) || (f_dim_siz[1] != 22) || (f_dim_siz[2] != 32)
            || (f_dim_siz[3] != 3) || (f_dim_siz[4] != 2)) {
        printf("- E - %s, %d: LUT dataset dimension sizes are unexpected\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  read the selected part of the whole lut: both mirror sides,
     *  1st 16 of 32 detectors, M bands at index 5 - 20, and the 2 times
     *  starting at itim1
     */
    start[0] = itim1;
    start[1] = 5;
    start[2] = 0;
    start[3] = 0;
    start[4] = 0;
    count[0] = 2;
    count[1] = 16;
    count[2] = 16;
    count[3] = 1;
    count[4] = 2;
    if (h5io_rd_ds_slice(&dsid, start, count, (void *) f_lut_sel) != 0) {
        printf("- E - %s, %d: Unable to read the VIIRS-SDR-F-LUT dataset\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  interpolate/extrapolate and apply
     */
    ncorr = N_MS * NDET * MAXBANDS;
    del = u58_t2 - u58_t1;
    dist = u58_time - u58_t1;
    fact = (float) dist / (float) del;
    if (want_verbose)
        printf("Interp factor from t2: %f\n", fact);

    if (corr_loc == NUMERATOR)
        for (i = 0; i < ncorr; i++)
            *(f_cal_corr + i) *= *(f_lut_sel + i) * (1. - fact)
            + *(f_lut_sel + i + ncorr) * fact;
    else
        for (i = 0; i < ncorr; i++)
            *(f_cal_corr + i) /= *(f_lut_sel + i) * (1. - fact)
            + *(f_lut_sel + i + ncorr) * fact;
    /*
     *  de-allocate space for arrays and close the LUT file
     */
    h5io_close(&dsid);
    h5io_close(&h5fid);

    return 0;
}

/**
 * This function reads a numeric attribute and converts it to a 32 bit int.
 * 
 * @param id HDF5 id for group or file
 * @param attr_name name of the attribute to read
 * @param val pointer to the 32 bit int where the value will be written
 * @return 0 for OK, else 1 for bad
 */
int h5ReadAttrInt32(h5io_str *id, char *attr_name, int32_t *val) {
    hid_t *base_grp, spaceid, attr_id;

    base_grp = (id->type == H5IO_TYP_DAT_ID) ?
            &(id->dat_id) : &(id->grp_id);
    if ((attr_id = H5Aopen_name(*base_grp, attr_name)) < 0) 
        return 1;

    if ((spaceid = H5Aget_space(attr_id)) < 0) 
        return 1;

    if(H5Sget_simple_extent_npoints(spaceid) != 1) {
        printf("%s: # of elements of attribute: %s has to be 1\n", __FILE__, attr_name);
        return 1;
    }

    if (H5Aread(attr_id, H5T_NATIVE_INT, val) < 0) 
        return 1;
    
    /*
     *  close the allocated spaces before leaving
     */
    if (H5Sclose(spaceid) < 0) 
        return 1;
    if (H5Aclose(attr_id) < 0) 
        return 1;
    return 0;
}

/* ------------------------------------------------------------------------
 openl1_viirs_h5

 purpose: opens a VIIRS L1B file for reading.

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 filehandle *      file             I      input file information

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson, SAIC, 10 Feb, 2009   Original development

 ------------------------------------------------------------------------*/
int openl1_viirs_h5(filehandle * file) {
    int npix, nscan, ibnd, jbnd, ids, i, igran;
    int *nscan_inst; /* # instrument scans per granule */
    int *nscan_stat; /* same as nscan_inst, but for dark bands (turned */
    /* off at night), valued at -993 */
    int nscan_sto; /* # scans in files storage, can be > nscan_inst_tot */
    int nscan_inst_tot; /* total # instrument scans in this file */
    int nscan_per_gran;
    h5io_str g_id, id_tmp;
    char g_path[100], geo_all_data_gnam[200], *band_dark, **grp_obj_nm;
    char *geo_name[] = {"Latitude", "Longitude", "SatelliteAzimuthAngle",
        "SatelliteZenithAngle", "SolarAzimuthAngle", "SolarZenithAngle",
        "Height"};
    float *rad_fact;
    H5T_class_t h5_class;
    hid_t h5_native_typ;
    int ndim, dim_siz[4], sto_len, n_obj, *grp_obj_typ;
    int terr_corr;
    unsigned char *cmp_vq2;
    int16_t ua_scn_fmt, ua_ndets, ua_margin[2];
    int32_t ua_npix, ua_nlin, attr_exist;
    int32_t qual2_found = 0, ngran, dgran, iscan, dk;
    /*  dgran will be used for checking the match in different files */
    int64_t u58_gran;
    /*
     *  get the names of the suite of VIIRS file names from the input name
     */
    if (gen_sdr_suite(file->name) != 0)
        return 1;

    // try a reasonable geo file if one is not given
    if(file->geofile[0] == 0) {
        strcpy(file->geofile, file->name);
        file->geofile[0] = 'G';
        file->geofile[1] = 'M';
        file->geofile[2] = 'T';
        file->geofile[3] = 'C';
        file->geofile[4] = 'O';
        fix_viirs_filename(file->geofile);
    }
    
    /* Open the HDF input file */
    /*
     *  for all the SDR band files, open and get, check size and set
     *  up for later dataset access
     */
    for (ibnd = 0, jbnd = 1; ibnd < MAXBANDS; ibnd++, jbnd++) {
        /* open */
        if (h5io_openr(vfnames[ibnd], 0, (sdr_fid + ibnd)) != 0) {
            fprintf(stderr, "-E- %s Line %d: Failure to open %s\n", __FILE__,
                    __LINE__, vfnames[ibnd]);
            return 1;
        }
        /*
         *  Determine the # of granules in this file
         */
        sprintf(g_path, "Data_Products/VIIRS-M%d-SDR/VIIRS-M%d-SDR_Aggr", jbnd,
                jbnd);
        if (h5io_set_ds((sdr_fid + ibnd), g_path, &g_id) != 0) {
            fprintf(stderr,
                    "-E- %s, %d: Failure to set aggregation dataset, band %d:\n",
                    __FILE__, __LINE__, ibnd);
            fprintf(stderr, "name: %s\n", vfnames[ibnd]);
            return 1;
        }
        if (h5ReadAttrInt32(&g_id, "AggregateNumberGranules", (void *) &dgran)
                != 0) {
            fprintf(stderr,
                    "-E- %s, %d: Unable to read the AggregateNumberGranules attribute\n",
                    __FILE__, __LINE__);
            fprintf(stderr, "band %d\n, name: %s\n", ibnd, vfnames[ibnd]);
            return 1;
        }
        if (h5io_close(&g_id) != 0) {
            fprintf(stderr, "-E- %s Line %d: Unable to close granule, ds %d\n",
                    __FILE__, __LINE__, ibnd);
            return 1;
        }
        /*
         *  set up the scan and granule maps and check # granules among band files
         */
        if (ibnd == 0) {
            ngran = dgran;
            /*
             *  With ngran determined, allocate granule-dependent arrays
             */
            if ((scale = (float *) malloc(MAXBANDS * ngran * sizeof (float)))
                    == NULL) {
                fprintf(stderr,
                        "-E- %s, %d: failure to allocate scale storage\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if ((offset = (float *) malloc(MAXBANDS * ngran * sizeof (float)))
                    == NULL) {
                fprintf(stderr,
                        "-E- %s, %d: failure to allocate offset storage\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if ((nscan_inst = (int *) malloc(ngran * sizeof (int))) == NULL) {
                fprintf(stderr,
                        "-E- %s, %d: failure to allocate nscan_inst storage\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if ((nscan_stat = (int *) malloc(ngran * sizeof (int))) == NULL) {
                fprintf(stderr,
                        "-E- %s, %d: failure to allocate nscan_stat storage\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if ((band_dark = (char *) calloc(ngran * MAXBANDS, sizeof (char)))
                    == NULL) {
                fprintf(stderr,
                        "-E- %s, %d: failure to allocate band_dark storage\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if ((rad_fact = (float *) malloc(ngran * 2 * sizeof (float))) == NULL) {
                fprintf(stderr,
                        "-E- %s, %d: failure to allocate rad_fact storage\n",
                        __FILE__, __LINE__);
                return 1;
            }
        } else {
            if (dgran != ngran) {
                fprintf(stderr,
                        "-E- %s, %d: Mismatch in AggregateNumberOfGranules\n",
                        __FILE__, __LINE__);
                fprintf(stderr, "band %d\n, name: %s\n", ibnd, vfnames[ibnd]);
                return 1;
            }
        }
        iscan = 0;
        /*
         *  Get the # scans from the granule area for all granules
         */
        nscan_inst_tot = 0;
        for (igran = 0; igran < ngran; igran++) {
            sprintf(g_path, "Data_Products/VIIRS-M%d-SDR/VIIRS-M%d-SDR_Gran_%d",
                    jbnd, jbnd, igran);
            if (h5io_set_ds((sdr_fid + ibnd), g_path, &g_id) != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Failure to set granule dataset, ds %d:\n",
                        __FILE__, __LINE__, ibnd);
                fprintf(stderr, "name: %s", vfnames[ibnd]);
                return 1;
            }

            if (h5io_rd_attr(&g_id, "N_Number_Of_Scans",
                    (void *) (nscan_inst + igran)) != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Unable to read the N_Number_Of_Scans ds attr, ds %d:\n",
                        __FILE__, __LINE__, ibnd);
                fprintf(stderr, "name: %s", vfnames[ibnd]);
                return 1;
            }
            if (nscan_inst[igran] <= 0) {
                fprintf(stderr,
                        "-E- %s, %d: N_Number_Of_Scans[%d]: %d invalid\n",
                        __FILE__, __LINE__, igran, nscan_inst[igran]);
                return 1;
            }
            /*  accumulate grand total of scans and set map info */
            nscan_inst_tot += nscan_inst[igran];

            /*  With the # scans and # granules, set up the gran_map and scan_map */
            if ((igran == 0) && (ibnd == 0)) {
                /* this will insure that simulated granules (with no # scan limit)
                 or std format granules ( with STDMAXSCANS at most scans)
                 will get sufficient storage allocated */
                nscan_per_gran =
                        (nscan_inst[0] > STDMAXSCANS) ?
                        nscan_inst[0] : STDMAXSCANS;

                if ((gran_map = (int *) malloc(
                        ngran * nscan_per_gran * sizeof (int))) == NULL) {
                    fprintf(stderr,
                            "-E- %s, %d: failure to allocate gran_map storage\n",
                            __FILE__, __LINE__);
                    return 1;
                }
                if ((scan_map = (int *) malloc(
                        ngran * nscan_per_gran * sizeof (int))) == NULL) {
                    fprintf(stderr,
                            "-E- %s, %d: failure to allocate scan_map storage\n",
                            __FILE__, __LINE__);
                    return 1;
                }
                for (i = 0; i < ngran * nscan_per_gran; i++) {
                    gran_map[i] = -1;
                    scan_map[i] = -1;
                }
            }
            for (i = 0; i < nscan_inst[igran]; i++) {
                scan_map[iscan] = igran * nscan_per_gran + i;
                gran_map[iscan++] = igran;
            }
            /*
             *   Get the start and end time (in IET to match scan times)
             */
            if (h5io_rd_attr(&g_id, "N_Beginning_Time_IET", (void *) &u58_gran)
                    != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Unable to read the N_Beginning_Time_IET ds attr, ds %d:\n",
                        __FILE__, __LINE__, ibnd);
                fprintf(stderr, "name: %s", vfnames[ibnd]);
                return 1;
            }
            if (igran == 0)
                u58_gran_st = u58_gran;
            else if ((u58_gran > 0) && (u58_gran < u58_gran_st))
                u58_gran_st = u58_gran;

            if (h5io_rd_attr(&g_id, "N_Ending_Time_IET", (void *) &u58_gran)
                    != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Unable to read the N_Ending_Time_IET ds attr, ds %d:\n",
                        __FILE__, __LINE__, ibnd);
                fprintf(stderr, "name: %s", vfnames[ibnd]);
                return 1;
            }
            if (igran == 0)
                u58_gran_en = u58_gran;
            else if (u58_gran > u58_gran_en)
                u58_gran_en = u58_gran;
            /*
             *  we need one of the 'GranX' datasets open to set up the
             *  calibration (calfile=<name>)
             *
             *  NOTE that a granule-specific or line-by-line cal may get switched
             *  to if needed - this makes cal based on just start time
             */
            if ((ibnd == 0) && (igran == ngran - 1))
                if (set_f_cal_corr(&g_id, file, u58_gran_st) != 0)
                    return 1;

            if (h5io_close(&g_id) != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Unable to close granule, ds %d\n",
                        __FILE__, __LINE__, ibnd);
                return 1;
            }
        }

        /*  Get the nscan_stat = NumberOfScans (from the All Data area) each time */
        sprintf(g_path, "All_Data/VIIRS-M%d-SDR_All/NumberOfScans", jbnd);
        if (h5io_grab_ds((sdr_fid + ibnd), g_path, (void *) nscan_stat) != 0) {
            fprintf(stderr,
                    "-E- %s, %d: Unable to read the # scans from NumberOfScans, band %d\n",
                    __FILE__, __LINE__, ibnd);
            return 1;
        }
        for (igran = 0; igran < ngran; igran++) {
            /*  note the dark bands  */
            if (nscan_stat[igran] == -993) {
                *(band_dark + ibnd + MAXBANDS * igran) = 1;
            } else {
                if (nscan_stat[igran] != nscan_inst[igran]) {
                    fprintf(stderr,
                            "-E- %s, %d: dataset #scans: %d not = to granule # scans: %d \n",
                            __FILE__, __LINE__, nscan_stat[igran],
                            nscan_inst[igran]);
                    fprintf(stderr, "band: %d, granule: %d\n", ibnd, igran);
                    return 1;
                }
            }
        }
        /*
         *  for the first band, get the count of scans in the actual
         *  datasets (we'll get that from the 'Radiance' dataset).
         *  This is informational only (the operational VIIRS SDRs
         *  have 48 scans of data space even if the # scans in the
         *  attributes is 47 - care must be taken due to this choice)
         *
         *  Also, set up storage for the quality 2
         */
        if (ibnd == 0) {
            sprintf(g_path, "All_Data/VIIRS-M%d-SDR_All/Radiance", jbnd);
            if (h5io_set_ds((sdr_fid + ibnd), g_path, &g_id) != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Getting nscan_sto, Unable to open path\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if (h5io_info(&g_id, NULL, &h5_class, &h5_native_typ, &ndim,
                    dim_siz, &sto_len) != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Getting nscan_sto, Unable to read dataset\n",
                        __FILE__, __LINE__);
                return 1;
            }
            nscan_sto = dim_siz[0] / NDET;
            if (h5io_close(&g_id) != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Getting nscan_sto, Unable to close path\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if (want_verbose)
                printf("%s: total # inst scans: %d, dataset scans: %d\n",
                    __FILE__, nscan_inst_tot, nscan_sto);
            if (nscan_sto < nscan_inst_tot) {
                fprintf(stderr, "-E- %s,%d: # aggregate scans < inst scans\n",
                        __FILE__, __LINE__);
            }
            /*
             */
            if ((viirs_qual2 = malloc(nscan_sto * sizeof (unsigned char)))
                    == NULL) {
                fprintf(stderr,
                        "-E- %s Line %d: mem allocation error for viirs_qual2\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if ((cmp_vq2 = malloc(nscan_sto * sizeof (unsigned char))) == NULL) {
                fprintf(stderr,
                        "-E- %s Line %d: mem allocation error for cmp_vq2\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }
        /*
         *  get the quality 2 data
         */
        sprintf(g_path, "All_Data/VIIRS-M%d-SDR_All/QF2_SCAN_SDR", jbnd);
        if (h5io_grab_ds((sdr_fid + ibnd), g_path, (void *) cmp_vq2) != 0) {
            fprintf(stderr,
                    "-E- %s, %d: Unable to get QF2 for path %s, band %d\n",
                    __FILE__, __LINE__, g_path, ibnd);
            return 1;
        }
        /*
         *  check # scan, ... consistency between datastes
         */
        if (ibnd == 0)
            nscan = nscan_inst_tot;
        else {
            if (nscan_inst_tot != nscan) {
                fprintf(stderr,
                        "-E- %s, %d:  # inst scans: %d in band file: %s",
                        __FILE__, __LINE__, nscan_inst_tot, vfnames[ibnd]);
                fprintf(stderr,
                        "is not the same as initially established (scans: %d)\n",
                        nscan);
                return 1;
            }
        }
        /*
         *  collect the quality 2 and check it for the non-dark bands
         *  and report if there is a lunar incursion
         */
        dk = 0;
        for (i = 0; i < ngran; i++) {
            if (*(band_dark + ibnd + MAXBANDS * i) == 1) {
                dk = 1;
                break;
            }
        }
        if (dk == 0) {
            if (qual2_found == 0) {
                memcpy(viirs_qual2, cmp_vq2, nscan_sto);
                qual2_found = 1;
                for (i = 0; i < nscan_inst_tot; i++) {
                    if ((*(viirs_qual2 + scan_map[i]) & 2) != 0) {
                        file->sv_with_moon = 1;
                        printf("VIIRS_SV_LUNAR_INCURSION detected\n");
                        break;
                    }
                }
            } else {
                /*  check the quality 2 for consistency */
                for (i = 0; i < nscan_inst_tot; i++) {
                    if (*(viirs_qual2 + scan_map[i])
                            != *(cmp_vq2 + scan_map[i])) {
                        fprintf(stderr,
                                "-E- %s Line %d: QF2_SCAN_SDR mismatch, val# %d, bnd %d\n",
                                __FILE__, __LINE__, i, ibnd);
                        return 1;
                    }
                }
            }
        }
        /*
         *  for non-aggregated files, see if extra info is there, and record
         *  unaggregated file info
         */
        if ((attr_exist = h5io_attr_exist((sdr_fid + ibnd), "Data Scan Format"))
                == -1) {
            fprintf(stderr,
                    "-E- %s Line %d: failure while checking for attr Data Scan Format, band: %d\n",
                    __FILE__, __LINE__, ibnd);
            return 1;
        } else {
            if (attr_exist == 1) {
                /*  a std aggregated file, no extra attributes */
                ua_npix = NAGGPX;
                ua_scn_fmt = 0;
                ua_nlin = nscan * NDET;
                ua_ndets = NDET;
                ua_margin[0] = 0;
                ua_margin[1] = 0;
            } else {
                /*  read the unaggregated information */
                if (h5io_rd_attr((sdr_fid + ibnd), "Data Scan Format",
                        (void *) &ua_scn_fmt) != 0) {
                    fprintf(stderr,
                            "-E- %s Line %d: Unable to read the Data Scan Format attr, ds %d: name: %s\n",
                            __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                    return 1;
                }
                if (h5io_rd_attr((sdr_fid + ibnd), "Scan Margin (track, scan)",
                        (void *) ua_margin) != 0) {
                    fprintf(stderr,
                            "-E- %s Line %d: Unable to read the Scan Margin attr, ds %d: name: %s\n",
                            __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                    return 1;
                }
                if (h5io_rd_attr((sdr_fid + ibnd), "Pixels per Scan Line",
                        (void *) &ua_npix) != 0) {
                    fprintf(stderr,
                            "-E- %s Line %d: Unable to read the Pixels per Scan Line attr, ds %d: name: %s\n",
                            __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                    return 1;
                }
                if (h5io_rd_attr((sdr_fid + ibnd), "Number of Scan Lines",
                        (void *) &ua_nlin) != 0) {
                    fprintf(stderr,
                            "-E- %s Line %d: Unable to read the Number of Scan Lines attr, ds %d: name: %s\n",
                            __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                    return 1;
                }
                if (h5io_rd_attr((sdr_fid + ibnd),
                        "Number of Detectors per Scan", (void *) &ua_ndets)
                        != 0) {
                    fprintf(stderr,
                            "-E- %s Line %d: Unable to read the Number of Detectors per Scan attr, ds %d: name: %s\n",
                            __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                    return 1;
                }
            }
        }
        /*
         *  check unaggregated file info between files
         */
        if (ibnd == 0) {
            margin[0] = ua_margin[0];
            margin[1] = ua_margin[1];
            file->npix = ua_npix;
            npix = ua_npix;
            file->ndets = ua_ndets;
            file->nscan = ua_nlin;
            scn_fmt = ua_scn_fmt;
            file->sd_id = 0; /* this appears to not be needed  */
        } else {
            if ((margin[0] != ua_margin[0]) || (margin[1] != ua_margin[1])) {
                fprintf(stderr,
                        "-E- %s Line %d: margin mismatch on band: %d name: %s\n",
                        __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                return 1;
            }
            if (file->npix != ua_npix) {
                fprintf(stderr,
                        "-E- %s Line %d: ua npix mismatch on band: %d name: %s\n",
                        __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                return 1;
            }
            if (file->ndets != ua_ndets) {
                fprintf(stderr,
                        "-E- %s Line %d: ua ndets mismatch on band: %d name: %s\n",
                        __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                return 1;
            }
            if (file->nscan != nscan * ua_ndets) {
                fprintf(stderr,
                        "-E- %s Line %d: ua # lines mismatch on band: %d name: %s\n",
                        __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                return 1;
            }
            if (scn_fmt != ua_scn_fmt) {
                fprintf(stderr,
                        "-E- %s Line %d: ua scan format mismatch on band: %d name: %s\n",
                        __FILE__, __LINE__, ibnd, vfnames[ibnd]);
                return 1;
            }
        }

        /* set up the datasets to read from this file, a line at a time */
        /* The RadianceFactors have been moved out of this, so for now, only
         1 dataset is set up this way, but in the future, this method may
         be of use so keep it  */
        for (ids = 0; ids < 1; ids++) {
            if (btype[ibnd] == VSWIR || btype[ibnd] == CIRRUS)
                sprintf(g_path, "All_Data/VIIRS-M%d-SDR_All/%s", jbnd,
                    "Reflectance");
            else
                sprintf(g_path, "All_Data/VIIRS-M%d-SDR_All/%s", jbnd,
                    "Radiance");

            if (h5io_set_ds((sdr_fid + ibnd), g_path, &(bnd_dat_id[ids][ibnd]))
                    != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Error setting dataset, band %d, ds %d\n",
                        __FILE__, __LINE__, ibnd, ids);
                return 1;
            }
            /*
             *  for the radiance, determine if it is scaled by finding data class
             */
            if (ids == 0) {
                if (h5io_info(&(bnd_dat_id[ids][ibnd]), NULL, &h5_class,
                        &h5_native_typ, &ndim, dim_siz, &sto_len) != 0) {
                    fprintf(stderr,
                            "-E- %s Line %d: Error accessing radiance dataset info, band %d\n",
                            __FILE__, __LINE__, ibnd);
                    return 1;
                }
                if (h5_class == H5T_INTEGER)
                    *(sdr_band_typ + ibnd) = 0;
                else
                    *(sdr_band_typ + ibnd) = 1;
            }
        }
        /*
         *  set scale and offset for band rads here, 0 if unscaled
         */
        if (*(sdr_band_typ + ibnd) == 0) {
            if (btype[ibnd] == VSWIR || btype[ibnd] == CIRRUS)
                sprintf(g_path, "All_Data/VIIRS-M%d-SDR_All/ReflectanceFactors",
                    jbnd);
            else
                sprintf(g_path, "All_Data/VIIRS-M%d-SDR_All/RadianceFactors",
                    jbnd);

            if (h5io_grab_ds((sdr_fid + ibnd), g_path, (void *) rad_fact)
                    != 0) {
                fprintf(stderr,
                        "-E- %s Line %d: Unable to read the rad fact, ds %d\n",
                        __FILE__, __LINE__, ibnd);
                return 1;
            }
            for (igran = 0; igran < ngran; igran++) {
                *(scale + ibnd + MAXBANDS * igran) = rad_fact[igran * 2];
                *(offset + ibnd + MAXBANDS * igran) = rad_fact[1 + igran * 2];
            }
        } else {
            for (igran = 0; igran < ngran; igran++) {
                *(scale + ibnd + MAXBANDS * igran) = 1.;
                *(offset + ibnd + MAXBANDS * igran) = 0.;
            }
        }
    } /***** end of band-specific initial set-up  *****/
    free(cmp_vq2);
    /*
     *  quality 2 needs to be collected to get the mirror sides.
     *  if not collected, there is trouble
     */
    if (qual2_found == 0) {
        fprintf(stderr,
                "-E- %s Line %d: All M bands contain no data (NumberOfScans undefined for all bands)\n",
                __FILE__, __LINE__);
        return 1;
    }

    /* Open the geolocation file */
    if (h5io_openr(file->geofile, 0, &geo_fid) != 0) {
        fprintf(stderr,
                "-E- %s Line %d:  Unable to open the geolocation file: %s\n",
                __FILE__, __LINE__, file->geofile);
        return 1;
    }
    /*
     *  Determine if this is a terrain-corrected or uncorrected geolocation
     *  VIIRS-MOD-GEO - un-corrected, VIIRS-MOD-GEO-TC - corrected
     */
    if (h5io_set_grp(&geo_fid, "Data_Products", &id_tmp) != 0) {
        fprintf(stderr,
                "-E- %s, %d: Failure to set geolocation Data_Products group\n",
                __FILE__, __LINE__);
        fprintf(stderr, "name: %s\n", file->geofile);
        return 1;
    }
    if (h5io_grp_contents(&id_tmp, &n_obj, &grp_obj_nm, &grp_obj_typ) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: failed to get group Data_Products contents\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (strcmp(grp_obj_nm[0], "VIIRS-MOD-GEO-TC") == 0) {
        terr_corr = 1;
        if (want_verbose)
            fprintf(stderr, "-I- %s , %d: Geolocation is Terrain Corrected\n",
                __FILE__, __LINE__);
    } else if (strcmp(grp_obj_nm[0], "VIIRS-MOD-GEO") == 0) {
        terr_corr = 0;
        if (want_verbose)
            fprintf(stderr, "-I- %s , %d: Geolocation is NOT Terrain Corrected\n",
                __FILE__, __LINE__);
    } else {
        fprintf(stderr, "-E- %s, %d: Geo Data_Products has unknown sub-group\n",
                __FILE__, __LINE__);
        fprintf(stderr, "= %s\n", grp_obj_nm[0]);
        return 1;
    }
    if (h5io_close(&id_tmp) != 0) {
        fprintf(stderr, "-E- %s, %d: failed to close group Data_Products\n",
                __FILE__, __LINE__);
        return 1;
    }
    file->terrain_corrected = terr_corr;
    /*
     *  Check geolocation file # granules to match the band files
     */
    if (terr_corr == 1)
        sprintf(g_path, "Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Aggr");
    else
        sprintf(g_path, "Data_Products/VIIRS-MOD-GEO/VIIRS-MOD-GEO_Aggr");

    if (h5io_set_ds(&geo_fid, g_path, &g_id) != 0) {
        fprintf(stderr,
                "-E- %s, %d: Failure to set geolocation aggregation dataset\n",
                __FILE__, __LINE__);
        fprintf(stderr, "name: %s\n", file->geofile);
        return 1;
    }
    if (h5ReadAttrInt32(&g_id, "AggregateNumberGranules", (void *) &dgran) != 0) {
        fprintf(stderr,
                "-E- %s, %d: Unable to read the AggregateNumberGranules attribute\n",
                __FILE__, __LINE__);
        fprintf(stderr, "name: %s\n", file->geofile);
        return 1;
    }
    if (h5io_close(&g_id) != 0) {
        fprintf(stderr, "-E- %s, %d: Unable to close geofile granule\n",
                __FILE__, __LINE__);
        return 1;
    }

    if (dgran != ngran) {
        fprintf(stderr,
                "-E- %s, %d: Mismatch found in # granules in geofile vs band files\n",
                __FILE__, __LINE__);
        fprintf(stderr, "    geofile: %d, band files: %d\n", dgran, ngran);
        return 1;
    }
    /*
     * get the group name under 'All_Data' from this geo file
     * for a geo file, it can be either VIIRS-MOD-GEO for elipsoid nav
     * or VIIRS-MOD-GEO-TC for terrain corrected
     */
    if (h5io_set_grp(&geo_fid, "All_Data", &id_tmp) != 0) {
        fprintf(stderr, "-E- %s Line %d: Unable to open group All_Data\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_grp_contents(&id_tmp, &n_obj, &grp_obj_nm, &grp_obj_typ) != 0) {
        fprintf(stderr,
                "-E- %s Line %d: failed to get group All_Data contents\n",
                __FILE__, __LINE__);
        return 1;
    }
    strcpy(geo_all_data_gnam, grp_obj_nm[0]);
    if (h5io_close(&id_tmp) != 0) {
        fprintf(stderr, "-E- %s Line %d: failed to close group All_Data\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  get the entire actual scan start times, position, velocity, and sensor
     *  attitude for all frames
     *  check/correct scan starts to be within gran t range (fixes anc
     *  problem downstream)
     */
    u58_scn_st = (int64_t *) malloc(nscan_sto * sizeof (int64_t));
    sprintf(g_path, "All_Data/%s/StartTime", geo_all_data_gnam);
    if (h5io_grab_ds(&geo_fid, g_path, (void *) u58_scn_st) != 0) {
        fprintf(stderr, "-E- %s, %d: Unable to read the Geo StartTime\n",
                __FILE__, __LINE__);
        return 1;
    }
    for (i = 0; i < nscan_inst_tot; i++) {
        iscan = scan_map[i];
        if ((u58_scn_st[iscan] < u58_gran_st)
                || (u58_scn_st[iscan] > u58_gran_en)) {
            fprintf(stderr,
                    "-W- %s Line %d: scan start time on scan %d was outside granule range - repaired.\n",
                    __FILE__, __LINE__, i);
            u58_scn_st[iscan] = u58_gran_st
                    + i * (u58_gran_en - u58_gran_st) / nscan_inst_tot;
        }
    }

    pos = (float *) malloc(nscan_sto * 3 * sizeof (float));

    sprintf(g_path, "All_Data/%s/SCPosition", geo_all_data_gnam);
    if (h5io_grab_ds(&geo_fid, g_path, (void *) pos) != 0) {
        fprintf(stderr, "-E- %s Line %d: Unable to read the SCPosition\n",
                __FILE__, __LINE__);
        return 1;
    }

    vel = (float *) malloc(nscan_sto * 3 * sizeof (float));

    sprintf(g_path, "All_Data/%s/SCVelocity", geo_all_data_gnam);
    if (h5io_grab_ds(&geo_fid, g_path, (void *) vel) != 0) {
        fprintf(stderr, "-E- %s Line %d: Unable to read the SCVelocity\n",
                __FILE__, __LINE__);
        return 1;
    }

    att = (float *) malloc(nscan_sto * 3 * sizeof (float));
    sprintf(g_path, "All_Data/%s/SCAttitude", geo_all_data_gnam);

    if (h5io_grab_ds(&geo_fid, g_path, (void *) att) != 0) {
        fprintf(stderr, "-E- %s Line %d: Unable to read the SCAttitude\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  set to the geolocation datasets for the lat, lon, and view angles
     */
    for (ids = 0; ids < 7; ids++) {
        sprintf(g_path, "All_Data/%s/%s", geo_all_data_gnam, geo_name[ids]);
        if (h5io_set_ds(&geo_fid, g_path, (geo_dat_id + ids)) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d:  Unable to set ds # %d in geolocation file\n",
                    __FILE__, __LINE__, ids);
            return 1;
        }
    }
    /*
     *  On to set-up of any initial storage or factors that can be done
     */
    if ((scl_rad = (unsigned short *) malloc(npix * sizeof (unsigned short)))
            == NULL) {
        fprintf(stderr, "-E- %s Line %d:  scl_rad allocate failed\n", __FILE__,
                __LINE__);
        return 1;
    }

    if ((flt_rad = (float *) malloc(npix * sizeof (float))) == NULL) {
        fprintf(stderr, "-E- %s Line %d:  flt_rad allocate failed\n", __FILE__,
                __LINE__);
        return 1;
    }
    /*
     *  get the Fobar here to set up Fo
     */
    rdsensorinfo(file->sensorID, l1_input->evalmask, "Fobar", (void **) &Fobar);
    /*
     *  finally, free locally used space
     */
    free(nscan_inst);
    free(nscan_stat);
    free(band_dark);
    free(rad_fact);

    return (LIFE_IS_GOOD);
}

/* ------------------------------------------------------------------------
 readl1_viirs_h5

 purpose: reads 1 line (scan) from a VIIRS L1B file and load l1rec.

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 filehandle *      file             I      input file information
 int32             dline            I      scan line to read, 0 origin
 l1str *           l1rec           I/O     data for that line
 int               lonlat           I      signal to read only lon and lat

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson, SAIC, Feb, 2009      Original development
 W. Robinson, SAIC, 10 sep 2012    use moon in space view to flag affected
 scans with nav failure

 ------------------------------------------------------------------------*/
int readl1_viirs_h5(filehandle * file, int32_t dline, l1str * l1rec, int lonlat) {
    static int32_t lastframe = -1;
    static double fsol;
    int start[2], count[2];
    int32_t i, ibnd, ipix, moon_affected;
    float *iptr, pos1[3], vel1[3], att1[3], sen_mat[3][3], coeff[10], rval;
    short year, day;
    double dsec, esdist, f_corr;
    static double mnorm[3];
    static int mside;

    int32_t detnum, igran, dscan, ascan, aline;
    int32_t nbands = (int32_t) file->nbands;

    l1rec->npix = file->npix;
    /*
     *  get proper line in data and detector #
     */
    dscan = dline / file->ndets;
    detnum = (dline % (NDET + 2 * margin[0])) - margin[0];
    ascan = scan_map[dscan]; /* scan and line in the datasets */
    aline = ascan * (NDET + 2 * margin[0]) + detnum;
    igran = gran_map[dscan];
    /*
     *  load scan start time
     */
    if (detnum < 0) /* can't have a detector outside instrument's existing ones */
        detnum = 0;
    if (detnum > NDET - 1)
        detnum = NDET - 1;
    double tai58 = u58_scn_st[ascan] / 1000000.0;
    l1rec->scantime = tai58_to_unix(tai58);
    unix2yds(l1rec->scantime, &year, &day, &dsec);

    /*
     *  get the location and view information from the geo file
     */
    for (ibnd = 0; ibnd < 7; ibnd++) {

        // when lonlat is on and ibnd is not 0, 1 and 5, skip the iteration
        // lonlat mode only needs lon, lat and solar zenith
        if (lonlat && (ibnd != 0) && (ibnd != 1) && (ibnd != 5))
            continue;

        switch (ibnd) {
        case 0:
            iptr = l1rec->lat;
            break;
        case 1:
            iptr = l1rec->lon;
            break;
        case 2:
            iptr = l1rec->sena;
            break;
        case 3:
            iptr = l1rec->senz;
            break;
        case 4:
            iptr = l1rec->sola;
            break;
        case 5:
            iptr = l1rec->solz;
            break;
        case 6:
            iptr = l1rec->height;
            break;
        }

        count[0] = 1;
        count[1] = l1rec->npix;
        start[0] = aline;
        start[1] = 0;
        if (h5io_rd_ds_slice(&(geo_dat_id[ibnd]), start, count, (void *) iptr)
                != 0) {
            fprintf(stderr,
                    "-E- %s, %d:  Failed read to geo line %d of band %d\n",
                    __FILE__, __LINE__, aline, ibnd);
            return 1;
        }
    }

    // stop if lonlat only
    if (lonlat)
        return LIFE_IS_GOOD;

    if (ascan != lastframe) {
        /*
         *  Get mirror side (0=A,1=B) from the low bit of the QF2_SCAN_SDR file
         */
        mside = ((*(viirs_qual2 + ascan) & 1) == 1) ? 1 : 0;
        /*
         *  convert the (pos, vel, att) from (m, m/s, arcsec) to (km, km/s, deg)
         *  and compute the sensor attitude in ECF
         */
        for (i = 0; i < 3; i++) {
            pos1[i] = *(pos + i + ascan * 3) / 1000.;
            vel1[i] = *(vel + i + ascan * 3) / 1000.;
            att1[i] = *(att + i + ascan * 3) / 3600.;
        }
        ocorient_(pos1, vel1, att1, sen_mat, coeff);
        /*
         *  extract the 'mirror normal' - WDR reverse sen_mat ixdexes used
         */
        for (i = 0; i < 3; i++)
            mnorm[i] = sen_mat[i][0];

        lastframe = ascan;
    }
    l1rec->mside = (int32_t) mside;

    /* get earth-sun distance correction for this frame */
    int32_t msec32 = (int32_t) (dsec * 1.e3);
    int32_t yr = (int32_t) year;
    int32_t dy = (int32_t) day;

    esdist = esdist_(&yr, &dy, &msec32);
    fsol = pow(1.0 / esdist, 2);

    /* Compute polarization frame rotation angles */
    compute_alpha(l1rec->lon, l1rec->lat, l1rec->senz, l1rec->sena, mnorm,
            l1rec->npix, l1rec->alpha);

    for (ipix = 0; ipix < l1rec->npix; ipix++)
        *(l1rec->pixnum + ipix) = ipix;
    /*
     *  it is possible that some flags can be set using the quality
     *  info in the files (maybe)  may be able to add this later
     *  hilt, navfail done in modis io
     */
    /*
     *  get the radiance information from band files and un-scale it if necessary
     *  check values against the fill values for VIIRS
     *  also make into Mw / cm^2 sr sec from W / m^2...
     *
     *  The refl bands go to nbands, after that, the IR get stored in Ltir
     */
    count[0] = 1;
    count[1] = l1rec->npix;
    start[1] = 0;
    ivswir = 0;
    itherm = 0;
    for (ibnd = 0; ibnd < MAXBANDS; ibnd++) {
        start[0] = aline;

        if (btype[ibnd] == VSWIR)
            l1rec->Fo[ivswir] = Fobar[ivswir] * fsol;
        /*  find if the line is moon affected  */
        moon_affected = ((*(viirs_qual2 + ascan) & 2) != 0) ? 1 : 0;
        /* get specific f table cal correction  */
        f_corr = *(f_cal_corr + mside + N_MS * (detnum + NDET * ibnd));

        if (sdr_band_typ[ibnd] == 0) {
            /*  re-scale */
            if (h5io_rd_ds_slice(&(bnd_dat_id[0][ibnd]), start, count,
                    (void *) scl_rad) != 0) {
                fprintf(stderr,
                        "-E- %s, %d:  Failed to read line %d of band %d\n",
                        __FILE__, __LINE__, aline, ibnd);
                return 1;
            }
            for (ipix = 0; ipix < l1rec->npix; ipix++) {
                /*
                 *  for vis bands in the dark, just set value = 0
                 */
                if ((btype[ibnd] == VSWIR)
                        && (*(l1rec->solz + ipix) > SOLZNIGHT)) {
                    *(l1rec->Lt + nbands * ipix + ivswir) = 0.;
                    continue;
                }
                /*
                 *  bad geo info will trigger nav fail flag
                 */
                if ((*(l1rec->lat + ipix) < -90.)
                        || (*(l1rec->lat + ipix) > 90.)
                        || (*(l1rec->lon + ipix) < -180.)
                        || (*(l1rec->lon + ipix) > 180.)) {
                    if (btype[ibnd] == VSWIR) {
                        *(l1rec->Lt + nbands * ipix + ivswir) = 0.;
                        l1rec->navfail[ipix] = 1;
                    } else if (btype[ibnd] == THERM)
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.;
                    else if (btype[ibnd] == CIRRUS)
                        l1rec->rho_cirrus[ipix] = BAD_FLT;
                    continue;
                }
                /*
                 *  moon in space view will set the nav warn flag
                 */
                if (moon_affected == 1)
                    l1rec->navwarn[ipix] = 1;

                /*
                 *  check data flag values
                 */
                switch (*(scl_rad + ipix)) {
                case SOUB_UINT16_FILL:
                    if (btype[ibnd] == VSWIR) {
                        *(l1rec->Lt + nbands * ipix + ivswir) = 1000.;
                        l1rec->hilt[ipix] = 1;
                    } else if (btype[ibnd] == THERM) {
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.;
                    } else if (btype[ibnd] == CIRRUS) {
                        l1rec->rho_cirrus[ipix] = BAD_FLT;
                    }
                    break;
                case NA_UINT16_FILL:
                case MISS_UINT16_FILL:
                case ONBOARD_PT_UINT16_FILL:
                case ONGROUND_PT_UINT16_FILL:
                case ERR_UINT16_FILL:
                case VDNE_UINT16_FILL:
                    if (btype[ibnd] == VSWIR) {
                        *(l1rec->Lt + nbands * ipix + ivswir) = BAD_FLT;
                    } else if (btype[ibnd] == THERM) {
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.;
                    } else if (btype[ibnd] == CIRRUS) {
                        l1rec->rho_cirrus[ipix] = BAD_FLT;
                    }
                    break;
                case ELINT_UINT16_FILL:
                    if (btype[ibnd] == VSWIR) {
                        *(l1rec->Lt + nbands * ipix + ivswir) = 0.;
                        l1rec->navfail[ipix] = 1;
                    } else if (btype[ibnd] == THERM)
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.;
                    else if (btype[ibnd] == CIRRUS)
                        l1rec->rho_cirrus[ipix] = BAD_FLT;
                    break;
                default:
                    if (btype[ibnd] == VSWIR)
                        *(l1rec->Lt + nbands * ipix + ivswir) = (*(scl_rad
                            + ipix) * *(scale + ibnd + MAXBANDS * igran)
                            + *(offset + ibnd + MAXBANDS * igran))
                        * l1rec->Fo[ivswir] * f_corr
                            * cos(l1rec->solz[ipix] / OEL_RADEG) / OEL_PI;
                    else if (btype[ibnd] == THERM)
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.1
                            * (*(scl_rad + ipix)
                            * *(scale + ibnd + MAXBANDS * igran)
                            + *(offset + ibnd + MAXBANDS * igran))
                        * f_corr;
                    else if (btype[ibnd] == CIRRUS)
                        l1rec->rho_cirrus[ipix] = (*(scl_rad + ipix)
                            * *(scale + ibnd + MAXBANDS * igran)
                            + *(offset + ibnd + MAXBANDS * igran)) * f_corr;
                    break;
                }
            }
        } else /* SDR band data in float form */ {
            if (h5io_rd_ds_slice(&(bnd_dat_id[0][ibnd]), start, count,
                    (void *) flt_rad) != 0) {
                fprintf(stderr,
                        "-E- %s, %d:  Failed to read scan %d of band %d\n",
                        __FILE__, __LINE__, aline, ibnd);
                return 1;
            }
            for (ipix = 0; ipix < l1rec->npix; ipix++) {
                if ((btype[ibnd] == VSWIR)
                        && (*(l1rec->solz + ipix) > SOLZNIGHT)) {
                    *(l1rec->Lt + nbands * ipix + ivswir) = 0.;
                    continue;
                }
                /*
                 *  bad geo info will trigger nav fail flag
                 */
                if ((*(l1rec->lat + ipix) < -90.)
                        || (*(l1rec->lat + ipix) > 90.)
                        || (*(l1rec->lon + ipix) < -180.)
                        || (*(l1rec->lon + ipix) > 180.)) {
                    if (btype[ibnd] == VSWIR) {
                        *(l1rec->Lt + nbands * ipix + ivswir) = 0.;
                        l1rec->navfail[ipix] = 1;
                    } else if (btype[ibnd] == THERM)
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.;
                    else if (btype[ibnd] == CIRRUS)
                        l1rec->rho_cirrus[ipix] = BAD_FLT;
                    continue;
                }
                /*
                 *  moon in space view will set the nav warn flag
                 */
                if (moon_affected == 1)
                    l1rec->navwarn[ipix] = 1;

                rval = *(flt_rad + ipix);
                if (near(rval, SOUB_FLOAT32_FILL)) {
                    if (btype[ibnd] == VSWIR) {
                        *(l1rec->Lt + nbands * ipix + ivswir) = 1000.;
                        l1rec->hilt[ipix] = 1;
                    } else if (btype[ibnd] == THERM) {
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.;
                    } else if (btype[ibnd] == CIRRUS) {
                        l1rec->rho_cirrus[ipix] = BAD_FLT;
                    }
                } else if (near(rval, (float)NA_FLOAT32_FILL)
                        || near(rval, MISS_FLOAT32_FILL)
                        || near(rval, ONBOARD_PT_FLOAT32_FILL)
                        || near(rval, ONGROUND_PT_FLOAT32_FILL)
                        || near(rval, ERR_FLOAT32_FILL)
                        || near(rval, VDNE_FLOAT32_FILL)) {
                    if (btype[ibnd] == VSWIR) {
                        *(l1rec->Lt + nbands * ipix + ivswir) = BAD_FLT;
                    } else if (btype[ibnd] == THERM) {
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.;
                    } else if (btype[ibnd] == CIRRUS) {
                        l1rec->rho_cirrus[ipix] = BAD_FLT;
                    }
                } else if (near(rval, ELINT_FLOAT32_FILL)) {
                    if (btype[ibnd] == VSWIR) {
                        *(l1rec->Lt + nbands * ipix + ivswir) = 0.;
                        l1rec->navfail[ipix] = 1;
                    } else if (btype[ibnd] == THERM)
                        *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.;
                    else if (btype[ibnd] == CIRRUS)
                        l1rec->rho_cirrus[ipix] = BAD_FLT;
                } else if (btype[ibnd] == VSWIR)
                    *(l1rec->Lt + nbands * ipix + ivswir) = rval * f_corr
                        * l1rec->Fo[ivswir]
                        * cos(l1rec->solz[ipix] / OEL_RADEG) / OEL_PI;
                else if (btype[ibnd] == THERM)
                    *(l1rec->Ltir + NBANDSIR * ipix + itherm) = 0.1 * rval
                        * f_corr;
                else if (btype[ibnd] == CIRRUS)
                    l1rec->rho_cirrus[ipix] = rval * f_corr;
            }
        }
        switch (btype[ibnd]) {
        case VSWIR:
            ivswir++;
            break;
        case THERM:
            itherm++;
            break;
        }
    }
    l1rec->detnum = (int32_t) detnum;
    l1rec->margin_s = margin[1];
    l1rec->scn_fmt = scn_fmt;

    /* Convert IR bands to brightness temperature */
    radiance2bt(l1rec, 1000);
    // flag bowtie deleted pixels
    for (ipix = 0; ipix < l1rec->npix; ipix++)
        flag_bowtie_deleted(l1rec, ipix, 0);

    return (LIFE_IS_GOOD);
}

int closel1_viirs_h5(filehandle * file)
/* ------------------------------------------------------------------------
 close1_viirs_h5

 purpose: Close all the ids for the band and geolocation file(s)

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 filehandle *      file             I      input file information

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson, SAIC, 11 Feb, 2009   Original development

 ------------------------------------------------------------------------*/ {
    int ict;
    /*
     *  first close the band dataset ids and file ids
     *  (may need extra loop for mult datasets in each band)
     */
    for (ict = 0; ict < MAXBANDS; ict++) {
        if (h5io_close(&(bnd_dat_id[0][ict])) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d:  Failed to close band data id #%d\n",
                    __FILE__, __LINE__, ict);
            return 1;
        }

        if (h5io_close(&(sdr_fid[ict])) != 0) {
            fprintf(stderr, "-E- %s Line %d:  Failed to close band file %s\n",
                    __FILE__, __LINE__, vfnames[ict]);
            return 1;
        }
    }
    /*
     *  close the geolocation data ids and file id
     */
    for (ict = 0; ict < 7; ict++) {
        if (h5io_close(&(geo_dat_id[ict])) != 0) {
            fprintf(stderr,
                    "-E- %s Line %d:  Failed to close geo data id #%d\n",
                    __FILE__, __LINE__, ict);
            return 1;
        }
    }

    if (h5io_close(&geo_fid) != 0) {
        fprintf(stderr, "-E- %s Line %d:  Failed to close geo data file %s\n",
                __FILE__, __LINE__, file->geofile);
        return 1;
    }
    /*
     *  Free allocated space for the work and exit
     */
    free(u58_scn_st);
    free(scl_rad);
    free(flt_rad);
    free(pos);
    free(vel);
    free(att);
    free(scale);
    free(offset);
    free(gran_map);
    free(scan_map);

    return (LIFE_IS_GOOD);
}
