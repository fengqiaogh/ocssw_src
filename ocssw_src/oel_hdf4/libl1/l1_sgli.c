/* ============================================================================ */
/* module l1_sgli.c - functions to read SGLI to MSL12
 * W. Robinson, SAIC 28 Sep 2016
 * ============================================================================ */

#include <math.h>

#include "l1.h"
#include <h5io.h>
#include "sgli.h"
#include "l1_sgli.h"

#define SGLI_NVNIR 9
#define SGLI_NVNIR_LG 2
#define SGLI_NVNIR_GEOM 2
#define SGLI_NSWIR 4

/* swirir_exist indicates (1) if the swir/ir file is available and
   resqkm indicates (1) if data is at quarter km resolution, else 1 km */
static int32_t *scan_times, swirir_exist, resqkm, resample, resamp_1km;
static int32_t n_vnir = SGLI_NVNIR, n_vnir_lg = SGLI_NVNIR_LG;
static int32_t n_vnir_geom = SGLI_NVNIR_GEOM, n_swir = SGLI_NSWIR;
static float *flt_buf, vnir_scale[SGLI_NVNIR], vnir_off[SGLI_NVNIR];
static float swir_scale[SGLI_NSWIR], swir_off[SGLI_NVNIR];
static float vnir_sat[SGLI_NVNIR], swir_sat[SGLI_NSWIR];
static float vnir_geom_scale[SGLI_NVNIR_GEOM], vnir_geom_off[SGLI_NVNIR_GEOM];
static float vnir_lg_scale[SGLI_NVNIR_LG], vnir_lg_off[SGLI_NVNIR_LG];
static float vnir_lg_sat[SGLI_NVNIR_LG];
static uint16_t vnir_mx_dn[SGLI_NVNIR], vnir_lg_mx_dn[SGLI_NVNIR_LG];
static uint16_t swir_mx_dn[SGLI_NSWIR];
static uint16_t *ui16_buf_1km;
static int32_t npix, npix_1km, nscan, npix_tie, nscan_tie;
static int32_t npix_tie_1km, nscan_tie_1km;
static h5io_str vnir_fid, swir_ir_fid, vnir_hg_dsid[SGLI_NVNIR];
static h5io_str vnir_lg_dsid[SGLI_NVNIR_LG], gen_dsid;
static h5io_str vnir_geom_dsid[SGLI_NVNIR_GEOM], swir_dsid[SGLI_NSWIR];
static double st_day_unix;
static char *geo_q, *swir_q, *sen_q, *sol_q;
static gsl_interp_accel *accel_x, *accel_y;
static gsl_spline2d *geo_int_id[3];
static double *xa, *ya, *xa_1km, *ya_1km;
static double *geo_x, *geo_y, *geo_z;
static float *tie_el, *tie_az;
static grid_res_str grid_res[2];
static band_geom_str band_geom_vnir[SGLI_NVNIR]; /* for sensor geom_per_band */
static band_geom_str band_geom_vnir_lg[SGLI_NVNIR_LG];
static band_geom_str band_geom_swir[SGLI_NSWIR];
static band_geom_str band_geom_nominal[2]; /* for nominal sen, sol angles */

sgli_t* createPrivateData_sgli() {

    sgli_t* data = (sgli_t*) calloc(1, sizeof (sgli_t));
    if (data == NULL) {
        fprintf(stderr, "-E- %s line %d: unable to allocate private data for sgli\n",
                __FILE__, __LINE__);
        exit(1);
    }

    return data;
}

/*
 *
 */
int sgli_file_ver(h5io_str *fid) {
    /*
     sgli_file_ver

     purpose: do a basic validity check on a SGLI file
  
     Returns 0 if all checks are OK
     Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     h5io_str *         fid             I      file id of file to check

     Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     Wayne Robinson    14 Oct 2016     Original development

     -----------------------------------------------------------------------*/
    int32_t stage_pass = 0;
    h5io_str g_id;
    char str_store[255];
    /*
     This will check (in group Global_attributes/attributes):
     Satellite: Global Change Observation Mission - Climate (GCOM-C)
     Sensor: Second-generation Global Imager (SGLI)
     Product_level: Level-1B
     Product_name:  Top of atmosphere radiance (reflectance)
     */

    if (h5io_set_grp(fid, "Global_attributes", &g_id) == 0) {
        if (h5io_attr_exist(&g_id, "Satellite") == 0) {
            if (h5io_rd_attr(&g_id, "Satellite", (void *) str_store) == 0) {
                if (strncmp(str_store, "Global Change Observation Mission - Climate (GCOM-C)", 52) == 0) {
                    stage_pass = 1;
                }
            }
        }
        if ((stage_pass == 1) && (h5io_attr_exist(&g_id, "Sensor") == 0)) {
            stage_pass = 0;
            if (h5io_rd_attr(&g_id, "Sensor", (void *) str_store) == 0) {
                if (strncmp(str_store, "Second-generation Global Imager (SGLI)", 38) == 0) {
                    stage_pass = 1;
                }
            }
        }
        if ((stage_pass == 1) && (h5io_attr_exist(&g_id, "Product_level") == 0)) {
            stage_pass = 0;
            if (h5io_rd_attr(&g_id, "Product_level", (void *) str_store) == 0) {
                if (strncmp(str_store, "Level-1B", 8) == 0) {
                    stage_pass = 1;
                }
            }
        }
        if ((stage_pass == 1) && (h5io_attr_exist(&g_id, "Product_name") == 0)) {
            stage_pass = 0;
            if (h5io_rd_attr(&g_id, "Product_name", (void *) str_store) == 0) {
                if (strncmp(str_store, "Top of atmosphere radiance (reflectance)", 40) == 0) {
                    stage_pass = 1;
                }
            }
        }
        h5io_close(&g_id);
    }
    return ( stage_pass == 0) ? 1 : 0;
}

int sgli_rad_info(h5io_str *fid, char *ds_nam, h5io_str *dsid,
        int *dim_siz, uint16_t *mx_dn, float *scale, float *off, float *sat)
/*-----------------------------------------------------------------------
 sgli_rad_info

 purpose: get information about the radiance arrays in SGLI files

 Returns 0 if all checks are OK
   Parameters: (in calling order)
   Type              Name            I/O     Description
   ----              ----            ---     -----------
   h5io_str *        fid              I      file id of file
   char *            ds_nam           I      dataset name in file
   h5io_str *        dsid             O      returned id of opened dataset
   int *             dim_siz          O      size of the dimensions
   uint16_t *        mx_dn            O      maximum DN for band counts
   float *           scale            O      scaling value for the data
   float *           off              O      offset for scaling
   float *           sat              O      saturation value of rad

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson       20 Oct 2016     Original development

 -----------------------------------------------------------------------*/ {
    H5T_class_t h5_class;
    hid_t h5_native_typ;
    int ndim, sto_len;

    if (h5io_set_ds(fid, ds_nam, dsid) != 0) {
        fprintf(stderr, "%s, %d, E - Cannot open dataset: %s\n", __FILE__,
                __LINE__, ds_nam);
        return 1;
    }
    /* get the sizes: npix, nscan and the resolution */
    if (h5io_info(dsid, NULL, &h5_class, &h5_native_typ, &ndim, dim_siz,
            &sto_len) != 0) {
        fprintf(stderr, "%s, %d, E - h5io_info failure, dataset: %s\n",
                __FILE__, __LINE__, ds_nam);
        return 1;
    }
    /*  get the max dn, scale, offset for these bands */
    if (h5io_rd_attr(dsid, "Maximum_valid_DN", (uint16_t *) mx_dn) != 0) {
        fprintf(stderr, "%s, %d, E - Cannot read Maximum_valid_DN for %s\n",
                __FILE__, __LINE__, ds_nam);
        return 1;
    }
    if (h5io_rd_attr(dsid, "Slope", (float *) scale) != 0) {
        fprintf(stderr, "%s, %d, E - Cannot read Slope for %s\n",
                __FILE__, __LINE__, ds_nam);
        return 1;
    }
    if (h5io_rd_attr(dsid, "Offset", (float *) off) != 0) {
        fprintf(stderr, "%s, %d, E - Cannot read Offset for %s\n",
                __FILE__, __LINE__, ds_nam);
        return 1;
    }
    if (h5io_rd_attr(dsid, "Saturation_radiance", (float *) sat) != 0) {
        fprintf(stderr, "%s, %d, E - Cannot read Saturation_radiance for %s\n",
                __FILE__, __LINE__, ds_nam);
        return 1;
    }
    return 0;
}

int view_ang_tie_init2(band_geom_str *band_geom, int32_t nbands )
/*-----------------------------------------------------------------------
 view_ang_tie_init2

 purpose: initialize tie point view angle data for interpolation

 Returns 0

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 band_geom_str *   band_geom       I/O     array of structures for tie 
                                           point data
 int32_t           nbands           I      # bands (above array size

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 Wayne Robinson    9 Dec 2016      Original development
 Wayne Robinson    6 Mar 2019      Update for file format change

 -----------------------------------------------------------------------*/ {
    int32_t ib, ip, il, ipl, npix_tie, nlin_tie;
    int16_t *tie_el_16, *tie_az_16;
    int start[2], count[2];
    double *geo[3], el_cvt, az_cvt, xy_mod, *xa, *ya;
    double deg2rad(double);
    /*
     *  initial set up
     */
    tie_el_16 = (int16_t *) tie_el;
    tie_az_16 = (int16_t *) tie_az;
    /*
     *  process each band, getting the interpolation object for a unit vector 
     *  version of the sensor angle
     */
    for (ib = 0; ib < nbands; ib++) {
        npix_tie = band_geom[ib].grd_desc->npix_tie;
        nlin_tie = band_geom[ib].grd_desc->nscn_sub_tie;
        xa = band_geom[ib].grd_desc->xa;
        ya = band_geom[ib].grd_desc->ya;
        if ((band_geom[ib].qual =
                (char *) calloc(npix_tie * nlin_tie, sizeof ( char))) == NULL) {
            fprintf(stderr,
                    "%s, %d, E - Failed to allocate vnir sensor, solar quality storage\n",
                    __FILE__, __LINE__);
            return 1;
        }
        /*
         *  read each zenith, azimuth
         */
        start[0] = band_geom[ib].grd_desc->tie_st_lin;
        start[1] = 0;
        count[0] = nlin_tie;
        count[1] = npix_tie;

        if ((h5io_rd_ds_slice(&band_geom[ib].dsid[1], start, count,
                (void *) tie_el_16) != 0) ||
                (h5io_rd_ds_slice(&band_geom[ib].dsid[0], start, count,
                (void *) tie_az_16) != 0)) {
            fprintf(stderr, "%s, %d, E Error reading vnir view angle tie points\n",
                    __FILE__, __LINE__);
            return 1;
        }

        for (il = 0; il < nlin_tie; il++) {
            for (ip = 0; ip < npix_tie; ip++) {
                ipl = ip + il * npix_tie;
                /*  scale to float */
                el_cvt = band_geom[ib].offset[1] +
                        band_geom[ib].scale[1] * tie_el_16[ipl];
                az_cvt = band_geom[ib].offset[0] +
                        band_geom[ib].scale[0] * tie_az_16[ipl];

                if ((el_cvt < 0.) || (el_cvt > 180.) ||
                        (az_cvt < -180.) || (az_cvt > 180.)) {
                    el_cvt = 0;
                    az_cvt = 0;
                    band_geom[ib].qual[ipl] = 1;
                }
                xy_mod = sin(deg2rad(el_cvt));
                geo_x[ipl] = xy_mod * cos(deg2rad(az_cvt));
                geo_y[ipl] = xy_mod * sin(deg2rad(az_cvt));
                geo_z[ipl] = cos(deg2rad(el_cvt));
            }
        }
        geo[0] = geo_x;
        geo[1] = geo_y;
        geo[2] = geo_z;
        /*  set up for the 2D spline interpolation and accelerator */
        for (ip = 0; ip < 3; ip++) {
            band_geom[ib].int_id_sen[ip] =
                    gsl_spline2d_alloc(gsl_interp2d_bicubic, npix_tie, nlin_tie);
            if ((gsl_spline2d_init(band_geom[ib].int_id_sen[ip], xa,
                    ya, geo[ip], npix_tie, nlin_tie)) != 0) {
                fprintf(stderr, "%s, %d, E vnir gsl_spline2d_init error\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }
    }
    /*
     *  interpolation coefficients are made
     */
    return 0;
}

int view_ang_tie_eval(band_geom_str *band_geom, int32_t ip, int32_t scan,
        float *azi_v, float *zen_v, char *qual)
/*-----------------------------------------------------------------------
 view_ang_tie_eval

 purpose: interpolate view angle to a pixel

 Returns 0

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 band_geom_str *   band_geom        I      structure for tie
                                           point data for that band
 int32_t           ip               I      specific pixel to evaluate
 int32_t           scan             I      specific scan to evaluate
 float *           azi_v            O      returned azimuth
 float *           zen_v            O      returned zenith
 char *            qual             O      quality 1 = bad

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 Wayne Robinson    12 Dec 2016     Original development

 -----------------------------------------------------------------------*/ {
    double xvec, yvec, zvec;
    double rad2deg(double);
    int32_t ipb, scan_tie, tie_en_line, npix_tie, nscan_tie;
    size_t ip_tie, il_tie;

    npix_tie = band_geom->grd_desc->npix_tie;
    nscan_tie = band_geom->grd_desc->nscn_sub_tie;
    tie_en_line = band_geom->grd_desc->tie_st_lin + nscan_tie - 1;
    scan_tie = ( ( scan >= band_geom->grd_desc->resamp * tie_en_line)) ?
        band_geom->grd_desc->resamp * tie_en_line - 1 : scan;

    if ((gsl_spline2d_eval_e(band_geom->int_id_sen[0], ip,
            scan, accel_x, accel_y, &xvec) != 0) ||
            (gsl_spline2d_eval_e(band_geom->int_id_sen[1], ip,
            scan, accel_x, accel_y, &yvec) != 0) ||
            (gsl_spline2d_eval_e(band_geom->int_id_sen[2], ip,
            scan, accel_x, accel_y, &zvec) != 0)) {
        fprintf(stderr, "%s, %d, E Error in vnir_lg gsl_spline2d_eval\n",
                __FILE__, __LINE__);
        return 1;
    }
    *zen_v = (float) rad2deg(acos(zvec));
    *azi_v = (float) rad2deg(atan2(yvec, xvec));
    /*  it is possible that tie point view geom is bad, so */
    ip_tie = gsl_interp_accel_find(accel_x, band_geom->grd_desc->xa, npix_tie, ip);
    il_tie = gsl_interp_accel_find(accel_y, band_geom->grd_desc->ya, nscan_tie, scan_tie);
    ipb = ip_tie + il_tie * npix_tie;
    *qual = 0;
    if ((band_geom->qual[ipb] != 0) ||
            (band_geom->qual[ ipb + 1 ] != 0) ||
            (band_geom->qual[ ipb + npix_tie ] != 0) ||
            (band_geom->qual[ ipb + 1 + npix_tie ] != 0))
        *qual = 1;
    return 0;
}

int openl1_sgli(filehandle * l1file) {
    char *cptr, *vis_nir_file, *swir_ir_file, ds_nam[FILENAME_MAX];
    char str_store[FILENAME_MAX];
    int32_t ibnd, ip;
    H5T_class_t h5_class;
    hid_t h5_native_typ;
    h5io_str g_id;
    int ndim, sto_len, dim_siz[4];
    float resolution;
    int st_year, st_mon, st_day;
    char *vnir_hg_bnam[] ={"01", "02", "03", "04", "05", "06", "07", "09", "10"};
    char *vnir_lg_bnam[] = {"08", "11"};
    char *swir_bnam[] = {"01", "02", "03", "04"};
    char *vnir_geom_nam[] = {"Latitude", "Longitude"};
    char *view_nam[] = {"Sensor", "Solar"};
    char *ang_nam[] = {"azimuth", "zenith"};
    char *sen_ang[2] = {"Sensor_azimuth", "Sensor_zenith"};

    //printf("%s, %d - I - SGLI vis/nir file name: %s\n", __FILE__, __LINE__,
    //        l1file->name);
    
    sgli_t *data = l1file->private_data = createPrivateData_sgli();

    /*
     *  Derive the SWIR/IR file name from the l1 file if not set
     */
    vis_nir_file = l1file->name;
    swir_ir_file = data->swir_ir_file;

    if (swir_ir_file[0] == 0) {
        strcpy(swir_ir_file, vis_nir_file);
        if ((cptr = strstr(swir_ir_file, "VNRD")) == NULL) {
            printf("%s, %d, I - VIS / NIR file has non-standard name - cannot create SWIR/IR file name\n", __FILE__, __LINE__);
            printf("VIS/NIR file name: %s\n", vis_nir_file);
            printf("Processing will only be done with the VIS/NIR data\n");
        } else {
            memcpy(cptr, "IRS", 3);
        }
    }
    /*
     *  open the vis/nir file and all the band datasets
     */
    if (h5io_openr(vis_nir_file, 0, &vnir_fid) != 0) {
        fprintf(stderr, "%s, %d, E - Failure to open %s\n", __FILE__,
                __LINE__, vis_nir_file);
        return 1;
    }
    /* set bands for reading, get info and check consistency */
    for (ibnd = 0; ibnd < n_vnir; ibnd++) {
        sprintf(ds_nam, "Image_data/Lt_VN%s", vnir_hg_bnam[ibnd]);
        if (sgli_rad_info(&vnir_fid, ds_nam, (vnir_hg_dsid + ibnd),
                dim_siz, (vnir_mx_dn + ibnd), (vnir_scale + ibnd),
                (vnir_off + ibnd), (vnir_sat + ibnd)) != 0)
            return 1;

        if (ibnd == 0) {
            l1file->npix = npix = dim_siz[1];
            l1file->nscan = nscan = dim_siz[0];

            /*  officially set the resolution */
            if (h5io_rd_attr((vnir_hg_dsid + ibnd), "Spatial_resolution",
                    (void *) &resolution) != 0) {
                fprintf(stderr, "%s, %d, E - Cannot read band resolution\n", __FILE__,
                        __LINE__);
                return 1;
            }
            resqkm = (resolution == 1000.) ? 0 : 1;

            /* set the read buffer for rad info */
            if ((flt_buf = (float *) malloc(npix * sizeof (float))) == NULL) {
                fprintf(stderr, "%s, %d, E - Cannot allocate the SGLI read buffer\n",
                        __FILE__, __LINE__);
                return 1;
            }
        } else {
            if ((dim_siz[1] != npix) || (dim_siz[0] != nscan)) {
                fprintf(stderr, "%s, %d, E - Vis/NIR band array size mismatch\n",
                        __FILE__, __LINE__);
                fprintf(stderr, "for band %s\n", ds_nam);
                return 1;
            }
        }
    }
    /* open the low gain 673 and 865 nm bands and check */
    for (ibnd = 0; ibnd < n_vnir_lg; ibnd++) {
        sprintf(ds_nam, "Image_data/Lt_VN%s", vnir_lg_bnam[ibnd]);

        if (sgli_rad_info(&vnir_fid, ds_nam, (vnir_lg_dsid + ibnd),
                dim_siz, (vnir_lg_mx_dn + ibnd), (vnir_lg_scale + ibnd),
                (vnir_lg_off + ibnd), (vnir_lg_sat + ibnd)) != 0)
            return 1;

        if ((dim_siz[1] != npix) || (dim_siz[0] != nscan)) {
            fprintf(stderr, "%s, %d, E - Vis/NIR band array size mismatch\n",
                    __FILE__, __LINE__);
            fprintf(stderr, "for band %s\n", ds_nam);
            return 1;
        }
    }
    /* get the start time */
    if (h5io_set_grp(&vnir_fid, "Global_attributes", &g_id) != 0) {
        fprintf(stderr, "%s, %d, E - Cannot set Global_attributes group\n",
                __FILE__, __LINE__);
        return 1;
    }
    if (h5io_rd_attr(&g_id, "Scene_start_time", (void *) str_store) != 0) {
        fprintf(stderr, "%s, %d, E - Cannot read Scene_start_time\n",
                __FILE__, __LINE__);
        return 1;
    }
    h5io_close(&g_id);
    sscanf(str_store, "%4d%2d%2d", &st_year, &st_mon, &st_day);
    st_day_unix = ymds2unix((int16_t) st_year, (int16_t) st_mon, (int16_t) st_day,
            0);

    /* get the time for each line - read it all here */
    strcpy(ds_nam, "Image_data/Line_msec");
    if (h5io_set_ds(&vnir_fid, ds_nam, &gen_dsid) != 0) {
        fprintf(stderr, "%s, %d, E - Cannot open dataset: %s\n", __FILE__,
                __LINE__, ds_nam);
        return 1;
    }
    scan_times = (int32_t *) malloc(nscan * sizeof (int32_t));
    if (h5io_rd_ds(&gen_dsid, (void *) scan_times) != 0) {
        fprintf(stderr, "%s, %d, E - Error reading dataset: %s\n", __FILE__,
                __LINE__, ds_nam);
        return 1;
    }
    if (h5io_close(&gen_dsid) != 0) {
        fprintf(stderr, "%s, %d, E - Error closing dataset: %s\n", __FILE__,
                __LINE__, ds_nam);
        return 1;
    }

    /* open the view geometry and geolocation */
    /* and get sizes */
    for (ibnd = 0; ibnd < n_vnir_geom; ibnd++) {
        sprintf(ds_nam, "Geometry_data/%s", vnir_geom_nam[ibnd]);
        if (h5io_set_ds(&vnir_fid, ds_nam, (vnir_geom_dsid + ibnd)) != 0) {
            fprintf(stderr, "%s, %d, E - Cannot open dataset: %s\n", __FILE__,
                    __LINE__, ds_nam);
            return 1;
        }
        /* for first band, get the sizes: npix, nscan and the resampling interval */
        if (ibnd == 0) {
            if (h5io_info((vnir_geom_dsid + ibnd), NULL, &h5_class, &h5_native_typ,
                    &ndim, dim_siz, &sto_len) != 0) {
                fprintf(stderr, "%s, %d, E - h5io_info failure\n", __FILE__,
                        __LINE__);
                return 1;
            }
            npix_tie = dim_siz[1];
            nscan_tie = dim_siz[0];
            if (h5io_rd_attr((vnir_geom_dsid + ibnd), "Resampling_interval",
                    (void *) &resample) != 0) {
                fprintf(stderr,
                        "%s, %d, E - Failed to read resampling for geo data\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }
        if ((h5io_rd_attr((vnir_geom_dsid + ibnd), "Slope",
                (void *) (vnir_geom_scale + ibnd)) != 0) ||
                (h5io_rd_attr((vnir_geom_dsid + ibnd), "Offset",
                (void *) (vnir_geom_off + ibnd)) != 0)) {
            fprintf(stderr,
                    "%s, %d, E - Failed to read scale or offset for geo data\n",
                    __FILE__, __LINE__);
            return 1;
        }
        
    }
    /*
     *  For the nominal sensor, solar angles
     */
    for (ibnd = 0; ibnd < 2; ibnd++) {
        for (ip = 0; ip < 2; ip++) {
            sprintf(ds_nam, "Geometry_data/%s_%s",
                    view_nam[ibnd], ang_nam[ip]);
            if (h5io_set_ds(&vnir_fid, ds_nam,
                    &band_geom_nominal[ibnd].dsid[ip]) != 0) {
                fprintf(stderr, "%s, %d, E - Cannot open dataset: %s\n", __FILE__,
                        __LINE__, ds_nam);
                return 1;
            }
            if ((h5io_rd_attr(&band_geom_nominal[ibnd].dsid[ip], "Slope",
                    (void *) (&band_geom_nominal[ibnd].scale[ip])) != 0) ||
                    (h5io_rd_attr(&band_geom_nominal[ibnd].dsid[ip], "Offset",
                    (void *) (&band_geom_nominal[ibnd].offset[ip])) != 0)) {
                fprintf(stderr,
                        "%s, %d, E - Failed to read scale or offset for nominal geo data\n",
                        __FILE__, __LINE__);
                return 1;
            }
            /*  as it is unsure that all geom_per_band senX sizes are same as
                the tie point size, do a check to be sure */
            if (h5io_info(&band_geom_nominal[ibnd].dsid[ip], NULL, &h5_class,
                    &h5_native_typ, &ndim, dim_siz, &sto_len) != 0) {
                fprintf(stderr, "%s, %d, E - h5io_info failure for band_geom_nominal\n",
                        __FILE__, __LINE__);
                return 1;
            }
            if ((dim_siz[1] != npix_tie) || (dim_siz[0] != nscan_tie)) {
                fprintf(stderr,
                        "%s, %d, E - Unexpected band_geom_nominal sizes found\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }
        band_geom_nominal[ibnd].grd_desc = grid_res;
    }
    /*
     *  if the geom per band is used, open all the sena, senz for vnir and
     *  low gain bands
     */
    if (l1_input->geom_per_band == 1) {
        for (ibnd = 0; ibnd < n_vnir; ibnd++) {
            for (ip = 0; ip < 2; ip++) {
                sprintf(ds_nam, "Geometry_data/%s_VN%s",
                        sen_ang[ip], vnir_hg_bnam[ibnd]);
                if (h5io_set_ds(&vnir_fid, ds_nam,
                        &band_geom_vnir[ibnd].dsid[ip]) != 0) {
                    fprintf(stderr, "%s, %d, E - Cannot open dataset: %s\n", __FILE__,
                            __LINE__, ds_nam);
                    return 1;
                }
                if ((h5io_rd_attr(&band_geom_vnir[ibnd].dsid[ip], "Slope",
                        (void *) (&band_geom_vnir[ibnd].scale[ip])) != 0) ||
                        (h5io_rd_attr(&band_geom_vnir[ibnd].dsid[ip], "Offset",
                        (void *) (&band_geom_vnir[ibnd].offset[ip])) != 0)) {
                    fprintf(stderr,
                            "%s, %d, E - Failed to read scale or offset for geo data\n",
                            __FILE__, __LINE__);
                    return 1;
                }
                /*  as it is unsure that all geom_per_band senX sizes are same as
                    the tie point size, do a check to e sure */
                if (h5io_info(&band_geom_vnir[ibnd].dsid[ip], NULL, &h5_class,
                        &h5_native_typ, &ndim, dim_siz, &sto_len) != 0) {
                    fprintf(stderr, "%s, %d, E - h5io_info failure for geom_per_band\n",
                            __FILE__, __LINE__);
                    return 1;
                }
                if ((dim_siz[1] != npix_tie) || (dim_siz[0] != nscan_tie)) {
                    fprintf(stderr, "%s, %d, E - Unexpected geom_per_band sizes found\n", __FILE__, __LINE__);
                    return 1;
                }
            }
        band_geom_vnir[ibnd].grd_desc = grid_res;
        }
        /*  low gain band-dependent geometry  */
        for (ibnd = 0; ibnd < n_vnir_lg; ibnd++) {
            for (ip = 0; ip < 2; ip++) {
                sprintf(ds_nam, "Geometry_data/%s_VN%s",
                        sen_ang[ip], vnir_lg_bnam[ibnd]);
                if (h5io_set_ds(&vnir_fid, ds_nam,
                        &band_geom_vnir_lg[ibnd].dsid[ip]) != 0) {
                    fprintf(stderr, "%s, %d, E - Cannot open dataset: %s\n", __FILE__,
                            __LINE__, ds_nam);
                    return 1;
                }
                if ((h5io_rd_attr(&band_geom_vnir_lg[ibnd].dsid[ip], "Slope",
                        (void *) (&band_geom_vnir_lg[ibnd].scale[ip])) != 0) ||
                        (h5io_rd_attr(&band_geom_vnir_lg[ibnd].dsid[ip], "Offset",
                        (void *) (&band_geom_vnir_lg[ibnd].offset[ip])) != 0)) {
                    fprintf(stderr,
                            "%s, %d, E - Failed to read scale or offset for geo data\n",
                            __FILE__, __LINE__);
                    return 1;
                }
                /*  as it is unsure that all geom_per_band senX sizes are same as
                    the tie point size, do a check to e sure */
                if (h5io_info(&band_geom_vnir_lg[ibnd].dsid[ip], NULL, &h5_class,
                        &h5_native_typ, &ndim, dim_siz, &sto_len) != 0) {
                    fprintf(stderr,
                            "%s, %d, E - h5io_info failure for vis lg geom_per_band\n",
                            __FILE__, __LINE__);
                    return 1;
                }
                if ((dim_siz[1] != npix_tie) || (dim_siz[0] != nscan_tie)) {
                    fprintf(stderr,
                            "%s, %d, E - Unexpected vis_lg geom_per_band sizes found\n",
                            __FILE__, __LINE__);
                    return 1;
                }
            }
        band_geom_vnir_lg[ibnd].grd_desc = grid_res;
        }
    }
    /*  now for the SWIR/IR file  */
    swirir_exist = 0;
    if (swir_ir_file[0] != 0) {
        /*  open the swir file and make sure it is right type */
        if (h5io_openr(swir_ir_file, 0, &swir_ir_fid) != 0) {
            fprintf(stderr, "%s, %d, I - The SWIR/IR file: %s\n",
                    __FILE__, __LINE__, swir_ir_file);
            fprintf(stderr, "does not exist.  Using VIS/NIR file only\n");
        } else {
            swirir_exist = 1;
            fprintf(stderr, "%s, %d, I - SGLI SWIR/IR file will be used\n",
                    __FILE__, __LINE__);
            fprintf(stderr, "    NAME: %s\n", swir_ir_file );
            if (sgli_file_ver(&swir_ir_fid) != 0) {
                fprintf(stderr,
                        "%s, %d, E - SWIR/IR file: %s is not a SGLI L1b file\n", __FILE__,
                        __LINE__, swir_ir_file);
                return 1;
            }
            /*  open the swir bands and check that the array sizes match */
            for (ibnd = 0; ibnd < n_swir; ibnd++) {
                sprintf(ds_nam, "Image_data/Lt_SW%s", swir_bnam[ibnd]);
                if (sgli_rad_info(&swir_ir_fid, ds_nam, (swir_dsid + ibnd),
                        dim_siz, (swir_mx_dn + ibnd), (swir_scale + ibnd),
                        (swir_off + ibnd), (swir_sat + ibnd)) != 0)
                    return 1;
                if ((resqkm == 1) && (ibnd != 2)) {
                    /*  the SW01, 02, 04 are at 1 km res for the qkm dataset
                        get # pixels for this */
                    if (ibnd == 0) {
                        npix_1km = dim_siz[1];
                    } else {
                        if (dim_siz[1] != npix_1km) {
                            fprintf(stderr, "%s, %d, E - image data size mismatch between\n",
                                    __FILE__, __LINE__);
                            fprintf(stderr, "VIS/NIR file: %s\nand SWIR/IR file: %s\n",
                                    vis_nir_file, swir_ir_file);
                            return 1;
                        }
                    }
                } else {
                    if ((dim_siz[1] != npix) || (dim_siz[0] != nscan)) {
                        fprintf(stderr,
                                "%s, %d, E - image data size mismatch between\n",
                                __FILE__, __LINE__);
                        fprintf(stderr, "VIS/NIR file: %s\nand SWIR/IR file: %s\n",
                                vis_nir_file, swir_ir_file);
                        return 1;
                    }
                }
            } /* end swir/ir radiance work */
            /*  also check the geoloc tie datasets to match */
            strcpy(ds_nam, "Geometry_data/Latitude");
            if (h5io_set_ds(&swir_ir_fid, ds_nam, &gen_dsid) != 0) {
                fprintf(stderr, "%s, %d, E - Cannot open dataset: %s\n", __FILE__,
                        __LINE__, ds_nam);
                return 1;
            }
            if (h5io_info(&gen_dsid, NULL, &h5_class, &h5_native_typ,
                    &ndim, dim_siz, &sto_len) != 0) {
                fprintf(stderr, "%s, %d, E - h5io_info failure\n", __FILE__,
                        __LINE__);
                return 1;
            }
            if ((npix_tie != dim_siz[1]) || (nscan_tie != dim_siz[0])) {
                fprintf(stderr,
                        "%s, %d, E - tie point data size mismatch between VIS/NIR file:\n",
                        __FILE__, __LINE__);
                fprintf(stderr, "VIS/NIR file: %s\nand SWIR/IR file: %s\n",
                        vis_nir_file, swir_ir_file);
                return 1;
            }
            h5io_close(&gen_dsid);

        /*  if the geom per band is used, open all the sena, senz for swir bands  */
        if (l1_input->geom_per_band == 1) {
            for (ibnd = 0; ibnd < n_swir; ibnd++) {
                for (ip = 0; ip < 2; ip++) {
                    sprintf(ds_nam, "Geometry_data/%s_SW%s",
                            sen_ang[ip], swir_bnam[ibnd]);
                    if (h5io_set_ds(&swir_ir_fid, ds_nam,
                            &band_geom_swir[ibnd].dsid[ip]) != 0) {
                        fprintf(stderr, "%s, %d, E - Cannot open dataset: %s\n", __FILE__,
                                __LINE__, ds_nam);
                        return 1;
                    }
                    if ((h5io_rd_attr(&band_geom_swir[ibnd].dsid[ip], "Slope",
                            (void *) (&band_geom_swir[ibnd].scale[ip])) != 0) ||
                            (h5io_rd_attr(&band_geom_swir[ibnd].dsid[ip], "Offset",
                            (void *) (&band_geom_swir[ibnd].offset[ip])) != 0)) {
                        fprintf(stderr,
                                "%s, %d, E - Failed to read scale or offset for geo data\n",
                                __FILE__, __LINE__);
                        return 1;
                    }
                    /*  check size of the tie point data */
                    if (h5io_info(&band_geom_swir[ibnd].dsid[ip], NULL, &h5_class,
                            &h5_native_typ, &ndim, dim_siz, &sto_len) != 0) {
                        fprintf(stderr,
                                "%s, %d, E - h5io_info failure for swir geom_per_band\n",
                                __FILE__, __LINE__);
                        return 1;
                    }
                    if ((resqkm == 0) || (ibnd == 2) ) {
                        if ((dim_siz[1] != npix_tie) || (dim_siz[0] != nscan_tie)) {
                            fprintf(stderr,
                                    "%s, %d, E - Unexpected swir geom_per_band sizes found\n",
                                    __FILE__, __LINE__);
                            return 1;
                        }
                    band_geom_swir[ibnd].grd_desc = grid_res;
                    }
                    else {
                        if( ( ibnd == 0 ) && ( ip == 0 ) ) {
                            npix_tie_1km = dim_siz[1];
                            nscan_tie_1km = dim_siz[0];
                            resamp_1km = resample * 4;
                        } else {
                            if ((dim_siz[1] != npix_tie_1km) || (dim_siz[0] != nscan_tie_1km)) {
                            fprintf(stderr,
                                    "%s, %d, E - Unexpected swir geom_per_band sizes found\n",                      
                                    __FILE__, __LINE__);
                            return 1; 
                            }
                        }
                    band_geom_swir[ibnd].grd_desc = ( grid_res + 1 );
                    }
                }
            }
        }
        }
    }
    return (0);
}

/*
  readl1_sgli - read a line of sgli data in
 */
int readl1_sgli(filehandle *file, int32_t scan, l1str *l1rec) {
    int32_t ip_tie, il_tie, sscan, escan, tie_st_line, tie_st_lin_1km;
    int32_t il, ipl, ip_1km, ip_off, scan_1km, isub, scan_tie, nscn_sub_tie;
    uint16_t *ui16_buf = (uint16_t *) flt_buf;
    int32_t start[2], count[2], ib, ip, ipb, ib_lg, ilr, iv;
    int32_t lt_lg_state, lt_hg_state;
    double lat_cvt, lon_cvt;
    double xvec, yvec, zvec, xy_mod, *geo[3];
    double rad2deg(double), deg2rad(double);
    char nav_bad;
    float zen_v, azi_v;
    static int32_t firstcall = 1, cur_scan_1km[3] = {-1, -1, -1};
    static int32_t tie_en_line, tie_en_lin_1km;

    enum lt_state {
        LT_BAD, LT_SAT, LT_NORM
    };
    int32_t nbands = l1rec->l1file->nbands;
    uint16_t lt_strip = pow(2, 14) - 1;
    float w_m2_to_mw_cm2 = 0.1; /* to store initial W m^-2 as mW cm^-2 */
    float lt_tmp;

    /*
     *  do some setup in the first call
     */
    /*
     *  if required for that record, set up the geom_per_band storage
     */
    if ((l1_input->geom_per_band == 1) && (l1rec->geom_per_band == NULL)) {
        init_geom_per_band(l1rec);
    }

    if (firstcall == 1) {
        /*  Set up for interpolating the geo and view information that is stored
            in tie point format at resample rate: resample 
            Only set up for the lines to be processed and always have bounding 
            tie point lines so as not to extrapolate.
            This is especially important so pad out the tie lines */
        sscan = l1_input->sline - 1;
        escan = l1_input->eline - 1;
        tie_st_line = sscan / resample - 1; /* pad out start tie line */
        tie_en_line = escan / resample + 2; /* pad out end tie line */

        tie_st_line = (tie_st_line < 0) ? 0 : tie_st_line;
        tie_en_line = (tie_en_line >= nscan_tie) ? nscan_tie - 1 : tie_en_line;
        /*
         *  there needs to be 4 lines of tie point data for the spline 
         *  to work - this will assure that there are 4
         */
        if (tie_en_line - tie_st_line < 3) {
            if (tie_en_line >= 3)
                tie_st_line = tie_en_line - 3;
            else
                tie_en_line = tie_st_line + 3;
        }
        if ((tie_st_line < 0) || (tie_en_line > nscan_tie - 1)) {
            fprintf(stderr, "%s, %d, E - L1 file has fewer than 4 tie point lines\n",
                    __FILE__, __LINE__);
            fprintf(stderr, "   Spline interpolation is not possible\n");
            /*  in the unlikely event of NEEDING to process such short files,
             *  a linear interpolation will need to be added as a fallback option
             */
            return 1;
        }
       /* set up the 1st grid description */
        nscn_sub_tie = tie_en_line - tie_st_line + 1;
        xa = (double *) malloc(npix_tie * sizeof ( double));
        ya = (double *) malloc(nscn_sub_tie * sizeof ( double));
        grid_res[0].nscn_sub_tie = nscn_sub_tie;
        grid_res[0].xa = xa;
        grid_res[0].ya = ya;
        grid_res[0].npix_tie = npix_tie;
        grid_res[0].tie_st_lin = tie_st_line;
        grid_res[0].resamp = resample;

        /*  First, latitude, longitude - read them in */
        start[0] = tie_st_line;
        start[1] = 0;
        count[0] = nscn_sub_tie;
        count[1] = npix_tie;

        if (((tie_el = (float *) malloc(npix_tie * nscn_sub_tie * sizeof ( float)))
                == NULL) ||
                ((tie_az = (float *) malloc(npix_tie * nscn_sub_tie * sizeof ( float)))
                == NULL) ||
                ((geo_x = (double *) malloc(npix_tie * nscn_sub_tie * sizeof ( double)))
                == NULL) ||
                ((geo_y = (double *) malloc(npix_tie * nscn_sub_tie * sizeof ( double)))
                == NULL) ||
                ((geo_z = (double *) malloc(npix_tie * nscn_sub_tie * sizeof ( double)))
                == NULL) ||
                ((geo_q = (char *) calloc(npix_tie * nscn_sub_tie, sizeof ( char)))
                == NULL)) {
            fprintf(stderr, "%s, %d, E - Failed to allocate lat, lon storage\n",
                    __FILE__, __LINE__);
            return 1;
        }

        if (h5io_rd_ds_slice(vnir_geom_dsid, start, count, (void *) tie_el)
                != 0) {
            fprintf(stderr, "%s, %d, E Error reading latitude tie points\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (h5io_rd_ds_slice((vnir_geom_dsid + 1), start, count, (void *) tie_az)
                != 0) {
            fprintf(stderr, "%s, %d, E Error reading longitude tie points\n",
                    __FILE__, __LINE__);
            return 1;
        }

        /* convert to x,y,z unit vectors */
        for (il = 0; il < nscn_sub_tie; il++) {
            ya[il] = (il + start[0]) * resample;
            for (ip = 0; ip < npix_tie; ip++) {
                if (il == 0) xa[ip] = ip * resample;
                ipl = ip + il * npix_tie;
                /*  lat, lon  not scaled
                      lat_cvt = vnir_geom_off[0] + vnir_geom_scale[0] * tie_el[ipl];
                      lon_cvt = vnir_geom_off[1] + vnir_geom_scale[1] * tie_az[ipl];
                 */
                lat_cvt = tie_el[ipl];
                lon_cvt = tie_az[ipl];
                if ((lat_cvt < -90.) || (lat_cvt > 90.) ||
                        (lon_cvt < -180.) || (lon_cvt > 180.)) {
                    lat_cvt = 0;
                    lon_cvt = 0;
                    geo_q[ipl] = 1;
                }
                xy_mod = cos(deg2rad(lat_cvt));
                geo_x[ipl] = xy_mod * cos(deg2rad(lon_cvt));
                geo_y[ipl] = xy_mod * sin(deg2rad(lon_cvt));
                geo_z[ipl] = sin(deg2rad(lat_cvt));
            }
        }
        geo[0] = geo_x;
        geo[1] = geo_y;
        geo[2] = geo_z;
        /*  set up for the 2D spline interpolation and accelerator */
        for (ip = 0; ip < 3; ip++) {
            geo_int_id[ip] =
                    gsl_spline2d_alloc(gsl_interp2d_bicubic, npix_tie, nscn_sub_tie);
            if ((gsl_spline2d_init(geo_int_id[ip], xa, ya, geo[ip], npix_tie,
                    nscn_sub_tie)) != 0) {
                fprintf(stderr, "%s, %d, E gsl_spline2d_init error\n", __FILE__,
                        __LINE__);
                return 1;
            }
        }
        accel_x = gsl_interp_accel_alloc();
        accel_y = gsl_interp_accel_alloc();
        /*
         *  set up sensor and solar zenith and azimuth
         */
        if (view_ang_tie_init2(band_geom_nominal, 2 ) != 0)
            return 1;
        /*
         *  include the band-dependent geometry, if requested
         */
        if (l1_input->geom_per_band == 1) {
            /* set up the spline coefficients for band dependent vnir */
            if (view_ang_tie_init2(band_geom_vnir, n_vnir) != 0)
                return 1;
            /*
             *  set up the spline coefficients for band dependent vnir low gain
             */
            if (view_ang_tie_init2(band_geom_vnir_lg, n_vnir_lg ) != 0)
                return 1;
        }

        /*  
         *  for the SWIR, set up various things
         *  for replication, we only need 1 line, but need to remember it
         *  for 3/4 swir bands
         */
        if (swirir_exist == 1) {
            /*  
             *  prepare for handling the 1km SWIR in a qkm file - 
             *  for replication, we only need 1 line, but need to remember it
             *  for 3/4 swir bands
             */
            if (resqkm == 1) {
                if (((ui16_buf_1km = (uint16_t *) malloc(npix_1km * 3 *
                        sizeof ( uint16_t))) == NULL) ||
                        ((swir_q = (char *) calloc(npix_1km * 3, sizeof ( char)))
                        == NULL)) {
                    fprintf(stderr, "%s, %d, E - Failed to allocate lat, lon storage\n",
                            __FILE__, __LINE__);
                    return 1;
                }
            }
            /*
             *  if band-dependent sensor angles needed, set up the interpolation 
             *  for this
             */
            if (l1_input->geom_per_band == 1) {
                /* WDR move to after the if( resqkm == 1 ...
                   and set-up of the grid_res[1] values???
                if (view_ang_tie_init2(band_geom_swir, n_swir ) != 0)
                    return 1;
*/
               /* set up the grid information for the tie point grids that
                  are used for the 1 km data */
                if( resqkm == 1 ) {
                    tie_st_lin_1km = sscan / resamp_1km - 1; /* pad out start tie line */
                    tie_en_lin_1km = escan / resamp_1km + 2; /* pad out end tie line */

                    tie_st_lin_1km = (tie_st_lin_1km < 0) ? 0 : tie_st_lin_1km;
                    tie_en_lin_1km = (tie_en_lin_1km >= nscan_tie_1km) ? nscan_tie_1km - 1 : tie_en_lin_1km;
                    /*
                     *  there needs to be 4 lines of tie point data for the spline 
                     *  to work - this will assure that there are 4
                     */
                    if (tie_en_lin_1km - tie_st_lin_1km < 3) {
                        if (tie_en_lin_1km >= 3)
                            tie_st_lin_1km = tie_en_lin_1km - 3;
                        else
                            tie_en_lin_1km = tie_st_lin_1km + 3;
                    }
                    if ((tie_st_lin_1km < 0) || (tie_en_lin_1km > nscan_tie_1km - 1)) {
                        fprintf(stderr, "%s, %d, E - L1 file has fewer than 4 tie point lines\n",
                                __FILE__, __LINE__);
                        fprintf(stderr, "   Spline interpolation is not possible\n");
                        /*  in the unlikely event of NEEDING to process such short files,
                         *  a linear interpolation will need to be added as a fallback option
                         */
                        return 1;
                    }
                   /* set up the 2nd grid description */
                    grid_res[1].nscn_sub_tie = tie_en_lin_1km - tie_st_lin_1km + 1;
                    xa_1km = (double *) malloc(npix_tie_1km * sizeof(double) );
                    for (ip = 0; ip < npix_tie_1km; ip++)
                        xa_1km[ip] = ip * resamp_1km;
                    ya_1km = (double *) malloc(grid_res[1].nscn_sub_tie 
                        * sizeof(double) );
                    for (il = 0; il < grid_res[1].nscn_sub_tie; il++)
                        ya_1km[il] = (il + tie_st_lin_1km) * resamp_1km;
                    grid_res[1].xa = xa_1km;
                    grid_res[1].ya = ya_1km;
                    grid_res[1].npix_tie = npix_tie_1km;
                    grid_res[1].tie_st_lin = tie_st_lin_1km;
                    grid_res[1].resamp = resamp_1km;
                }
                if (view_ang_tie_init2(band_geom_swir, n_swir ) != 0)
                    return 1;
            }
        }
        firstcall = 0;
    }
    /*
     *  the use of the last scan line in gsl_interp_accel_find causes problems,
     *  so, make this adjustment. scan_tie is modified scan so its not at 
     *  the end line of the tie point grid
     */
    scan_tie = (scan >= (resample * tie_en_line)) ?
            resample * tie_en_line - 1 : scan;

    /*  get the VIS/NIR read and scaled */
    for (ib = 0; ib < n_vnir; ib++) {
        /*  read in the scan line of Lt data */
        start[0] = scan;
        start[1] = 0;
        count[0] = 1;
        count[1] = npix;
        if (h5io_rd_ds_slice((vnir_hg_dsid + ib), start, count,
                (void *) ui16_buf) != 0) {
            fprintf(stderr, "%s, %d, E Error reading VIS/NIR hg, index %d\n",
                    __FILE__, __LINE__, ib);
            return 1;
        }
        /*  for data that is within the valid DN, scale it to mw m^-2...
            and if above the saturation radiance, set the hilt */
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + ib;
            if (*(ui16_buf + ip) < *(vnir_mx_dn + ib)) {
                /* They say the 1st 2 bits are for stray light and sign of corr
                   but examples have none of this.  For now, strip off top 2 bits
                   when applying the scale, offset */
                lt_tmp = *(vnir_off + ib) + *(vnir_scale + ib) *
                        (*(ui16_buf + ip) & lt_strip);
                if ((lt_tmp >= *(vnir_sat + ib)) &&
                        (ib != 6) && (ib != 8))
                    l1rec->hilt[ip] = 1;
                l1rec->Lt[ipb] = lt_tmp * w_m2_to_mw_cm2;
            }
        }
    }
    /*
     *  nominal view geometry
     */
    for (iv = 0; iv < 2; iv++) {
        for (ip = 0; ip < npix; ip++) {
            if (view_ang_tie_eval(&band_geom_nominal[iv], ip, scan, 
                    &azi_v, &zen_v, &nav_bad) != 0) return 1;
            if (iv == 0) {
                l1rec->senz[ip] = zen_v;
                l1rec->sena[ip] = azi_v;
            } else {
                l1rec->solz[ip] = zen_v;
                l1rec->sola[ip] = azi_v;
            }
            if (nav_bad != 0) l1rec->navfail[ip] = 1;
        }
    }
    /*  high gain band-dependent sena  */
    /*  and solar band-dependent values (just mirror the nominal) */
    if (l1_input->geom_per_band == 1) {
        for (ib = 0; ib < n_vnir; ib++) {
            /*  
             *  get interpolated view vectors and azimuth, zenith
             */
            for (ip = 0; ip < npix; ip++) {
                ipl = ip * nbands + ib;
                //printf( "ib = %d, ip = %d\n", ib, ip );
                if (view_ang_tie_eval(&band_geom_vnir[ib], ip, scan,
                        &azi_v, &zen_v, &nav_bad) != 0) return 1;
                l1rec->geom_per_band->senz[ipl] = zen_v;
                l1rec->geom_per_band->sena[ipl] = azi_v;
                if (nav_bad != 0) l1rec->navfail[ip] = 1;

                l1rec->geom_per_band->solz[ipl] = l1rec->solz[ip];
                l1rec->geom_per_band->sola[ipl] = l1rec->sola[ip];
            }
        }
    }
    /*  for the low gain bands, fill any saturated values with low gain part */
    for (ib = 0; ib < n_vnir_lg; ib++) {
        ib_lg = (ib == 0) ? 6 : 8;
        if (h5io_rd_ds_slice((vnir_lg_dsid + ib), start, count,
                (void *) ui16_buf) != 0) {
            fprintf(stderr, "%s, %d, E Error reading VIS/NIR hg, index %d\n",
                    __FILE__, __LINE__, ib);
            return 1;
        }
        /*  for data that is within the valid DN, scale it to mw m^-2...
            and if above the saturation radiance, set the hilt */

        /* logic to get the correct data between the low-and high-gain bands  */
        /* NOTE that due to observed non-linearities near saturation in the 
           high gain, the switch-over to low agin data use is set to .75 of 
           the high gain saturation - ...* 0.75... */
        for (ip = 0; ip < npix; ip++) {
            /* set status of the Lt for the high gain data */
            ipb = ip * nbands + ib_lg;
            if (l1rec->Lt[ipb] == BAD_FLT)
                lt_hg_state = LT_BAD;
            else {
                lt_hg_state = LT_NORM;
                if (l1rec->Lt[ipb] >= *(vnir_sat + ib_lg) * 0.75 * 
                    w_m2_to_mw_cm2)
                    lt_hg_state = LT_SAT;
            }
            /* if the hg is normal, nothing left, otherwise check the low gain */
            if (lt_hg_state != LT_NORM) {
                /*  get lg rad and lg rad state */
                if (*(ui16_buf + ip) < *(vnir_lg_mx_dn + ib)) {
                    lt_lg_state = LT_NORM;
                    lt_tmp = *(vnir_lg_off + ib) + *(vnir_lg_scale + ib) *
                            (*(ui16_buf + ip) & lt_strip);
                    if (lt_tmp >= *(vnir_lg_sat + ib))
                        lt_lg_state = LT_SAT;
                    lt_tmp *= w_m2_to_mw_cm2;
                    /*  get the low gain sensor zenith, azimuth - may be needed */
                    if (l1_input->geom_per_band == 1) {
                        if (view_ang_tie_eval(&band_geom_vnir_lg[ib], ip, scan, 
                                &azi_v, &zen_v, &nav_bad) != 0) return 1;
                    }
                } else
                    lt_lg_state = LT_BAD;
                /* go through hg, lg states to set Lt, detector used, and hilt */
                /* setting.  if we have band-dependent view angles, make sure 
                   the lg angles are used */
                if (lt_lg_state == LT_NORM) {
                    l1rec->Lt[ipb] = lt_tmp;
                    if (l1_input->geom_per_band == 1) {
                        l1rec->geom_per_band->sena[ipb] = azi_v;
                        l1rec->geom_per_band->senz[ipb] = zen_v;
                        l1rec->navfail[ip] = nav_bad;
                    }
                } else if (lt_lg_state == LT_SAT) {
                    l1rec->Lt[ipb] = lt_tmp;
                    l1rec->hilt[ip] = 1;
                    if (l1_input->geom_per_band == 1) {
                        l1rec->geom_per_band->sena[ipb] = azi_v;
                        l1rec->geom_per_band->senz[ipb] = zen_v;
                        l1rec->navfail[ip] = nav_bad;
                    }
                } else if (lt_hg_state == LT_SAT) {
                    /* lg is bad, so stick with hg value */
                    l1rec->hilt[ip] = 1;
                    /* use hg detector position - already there as default */
                } else {
                    /* both Lt bad, use bad (already there) and 
                       hg detector position - already there as default */
                }
            }
        }
        /* end logic for hg and lg bands */
    }
    /* set the time for the line */
    l1rec->scantime = st_day_unix + (double) *(scan_times + scan) / 1.e3;
    /* set the geolocation and view geometry */
    for (ip = 0; ip < npix; ip++) {
        if ((gsl_spline2d_eval_e(geo_int_id[0], ip, scan, accel_x,
                accel_y, &xvec) != 0) ||
                (gsl_spline2d_eval_e(geo_int_id[1], ip, scan, accel_x,
                accel_y, &yvec) != 0) ||
                (gsl_spline2d_eval_e(geo_int_id[2], ip, scan, accel_x,
                accel_y, &zvec) != 0)) {
            fprintf(stderr, "%s, %d, E Error in gsl_spline2d_eval\n",
                    __FILE__, __LINE__);
            return 1;
        }
        l1rec->lat[ip] = (float) rad2deg(asin(zvec));
        l1rec->lon[ip] = (float) rad2deg(atan2(yvec, xvec));
        /*  it is possible that tie point nav is bad, so */
        ip_tie = gsl_interp_accel_find(accel_x, xa, npix_tie, ip);
        il_tie = gsl_interp_accel_find(accel_y, ya, nscan_tie, scan_tie);
        ipb = ip_tie + il_tie * npix_tie;
        if ((geo_q[ipb] != 0) || (geo_q[ ipb + 1 ] != 0) ||
                (geo_q[ ipb + npix_tie ] != 0) ||
                (geo_q[ ipb + 1 + npix_tie ] != 0))
            l1rec->navfail[ip] = 1;
    }
    /*
     *  read the SWIR bands, if available
     *  this will do replication for the 1 km bands in a qkm file
     */
    if (swirir_exist == 1) {
        for (ib = 0; ib < n_swir; ib++) {
            if ((resqkm == 1) && (ib != 2)) {
                /*  replicate the 1km data to the qkm size */
                ilr = ib;
                if (ib == 3) ilr = 2;
                scan_1km = scan / 4;
                if (cur_scan_1km[ilr] != scan_1km) {
                    start[0] = scan_1km;
                    start[1] = 0;
                    count[0] = 1;
                    count[1] = npix_1km;
                    if ((h5io_rd_ds_slice((swir_dsid + ib), start, count,
                            (void *) (ui16_buf_1km + ilr * npix_1km))) != 0) {
                        fprintf(stderr, "%s, %d, E Error reading SWIR, band index %d\n",
                                __FILE__, __LINE__, ib);
                        return 1;
                    }
                    cur_scan_1km[ilr] = scan_1km;
                }
                for (ip_1km = 0; ip_1km < npix_1km; ip_1km++) {
                    ip = ip_1km * 4;
                    ipb = ip * nbands + (ib + n_vnir);
                    ip_off = ip_1km + npix_1km * ilr;
                    if (*(ui16_buf_1km + ip_off) < *(swir_mx_dn + ib)) {
                        lt_tmp = *(swir_off + ib) + *(swir_scale + ib) *
                                (*(ui16_buf_1km + ip_off) & lt_strip);
                        if (lt_tmp >= *(swir_sat + ib))
                            for (isub = 0; isub < 4; isub++)
                                l1rec->hilt[ ip + isub ] = 1;
                        lt_tmp *= w_m2_to_mw_cm2;
                        for (isub = 0; isub < 4; isub++)
                            l1rec->Lt[ipb + isub * nbands] = lt_tmp;
                    }
                }
            } else {
                /*  deal with bands at all the same resolution */
                start[0] = scan;
                start[1] = 0;
                count[0] = 1;
                count[1] = npix;
                if (h5io_rd_ds_slice((swir_dsid + ib), start, count,
                        (void *) ui16_buf) != 0) {
                    fprintf(stderr, "%s, %d, E Error reading VIS/NIR hg, index %d\n",
                            __FILE__, __LINE__, ib);
                    return 1;
                }
                /*  check for valid data in band, note initialized to BAD_FLT */
                for (ip = 0; ip < npix; ip++) {
                    ipb = ip * nbands + (ib + n_vnir);
                    if (*(ui16_buf + ip) <= *(swir_mx_dn + ib)) {
                        lt_tmp = *(swir_off + ib) + *(swir_scale + ib) *
                                (*(ui16_buf + ip) & lt_strip);
                        if (lt_tmp >= *(swir_sat + ib))
                            l1rec->hilt[ip] = 1;
                        l1rec->Lt[ipb] = lt_tmp * w_m2_to_mw_cm2;
                    }
                }
            }
        }
        /*
         *  get the geom_per_band for the swir, note that all sena/z arrays
         *  are relative to the 1/4 km resolution data and so will not need
         *  the same treatment as the Lt did above for 3/4 bands
         */
        if (l1_input->geom_per_band == 1) {
            for (ib = 0; ib < n_swir; ib++) {
                /*
                 *  get interpolated view vectors and azimuth, zenith
                 *  As the resolution of the grids changes for qkm SWIR 
                 *  in bands 0, 2, and 3 and then back after, do these 
                 *  calls to gsl_interp_accel_reset
                 */
                if ( ( resqkm == 1 ) && ( ib != 1 ) ) {
                    gsl_interp_accel_reset( accel_x );
                    gsl_interp_accel_reset( accel_y );
                }
                for (ip = 0; ip < npix; ip++) {
                    ipb = ip * nbands + (ib + n_vnir);
                    if (view_ang_tie_eval(&band_geom_swir[ib], ip, scan, 
                            &azi_v, &zen_v, &nav_bad) != 0) return 1;
                    l1rec->geom_per_band->senz[ipb] = zen_v;
                    l1rec->geom_per_band->sena[ipb] = azi_v;
                    if (nav_bad != 0) l1rec->navfail[ip] = 1;
                    /*  solar values, though constant are also per-band */
                    l1rec->geom_per_band->solz[ipb] = l1rec->solz[ip];
                    l1rec->geom_per_band->sola[ipb] = l1rec->sola[ip];
                }
            }
            if ( resqkm == 1 ) {
                gsl_interp_accel_reset( accel_x );
                gsl_interp_accel_reset( accel_y );
            }
        }
    } else {
       /* if geom per band is used but no swir is available, fill 
          view geom for swir bands with the nominal values */
        if (l1_input->geom_per_band == 1) {
            for (ib = 0; ib < n_swir; ib++) {
                for (ip = 0; ip < npix; ip++) {
                    ipb = ip * nbands + (ib + n_vnir);
                    l1rec->geom_per_band->senz[ipb] = l1rec->senz[ip];
                    l1rec->geom_per_band->sena[ipb] = l1rec->sena[ip];
                    l1rec->geom_per_band->solz[ipb] = l1rec->solz[ip];
                    l1rec->geom_per_band->sola[ipb] = l1rec->sola[ip];
                }
            }
        }
    }

    return (0);
}

int closel1_sgli(filehandle *file) {
    int32_t ibnd, iv;
    sgli_t *data = file->private_data;

    /*
     *  free allocated space, close all IDs and files
     */
    free(data);
    free(flt_buf);
    free(tie_el);
    free(tie_az);
    free(geo_x);
    free(geo_y);
    free(geo_z);
    free(geo_q);
    free(sen_q);
    free(sol_q);
    free(xa);
    free(ya);

    for (ibnd = 0; ibnd < n_vnir; ibnd++)
        h5io_close(vnir_hg_dsid + ibnd);
    for (ibnd = 0; ibnd < n_vnir_lg; ibnd++)
        h5io_close(vnir_lg_dsid + ibnd);
    for (ibnd = 0; ibnd < n_swir; ibnd++)
        h5io_close(swir_dsid + ibnd);
    for (ibnd = 0; ibnd < n_vnir_geom; ibnd++)
        h5io_close(vnir_geom_dsid + ibnd);

    h5io_close(&vnir_fid);

    for (ibnd = 0; ibnd < 3; ibnd++)
        gsl_spline2d_free(geo_int_id[ibnd]);
    /*
     *  free splines for nominal sen, sol angles
     */
    for (ibnd = 0; ibnd < 2; ibnd++)
        for (iv = 0; iv < 3; iv++)
            gsl_spline2d_free(band_geom_nominal[ibnd].int_id_sen[iv]);
    /*
     *  free splines for band-dependent sensor angles
     */
    if (l1_input->geom_per_band == 1) {
        for (ibnd = 0; ibnd < n_vnir; ibnd++)
            for (iv = 0; iv < 3; iv++)
                gsl_spline2d_free(band_geom_vnir[ibnd].int_id_sen[iv]);
        for (ibnd = 0; ibnd < n_vnir_lg; ibnd++)
            for (iv = 0; iv < 3; iv++)
                gsl_spline2d_free(band_geom_vnir_lg[ibnd].int_id_sen[iv]);
    }
    /*
     *  for the SWIR/IR file, free allocated space, close all IDs and files
     */
    if (swirir_exist) {
        if (resqkm == 1) {
            free(xa_1km);
            free(ya_1km);
            free(ui16_buf_1km);
            }

        if (l1_input->geom_per_band == 1)
            for (ibnd = 0; ibnd < n_swir; ibnd++)
                for (iv = 0; iv < 3; iv++)
                    gsl_spline2d_free(band_geom_swir[ibnd].int_id_sen[iv]);

        for (ibnd = 0; ibnd < n_swir; ibnd++)
            h5io_close(swir_dsid + ibnd);

        h5io_close(&swir_ir_fid);
    }
    return (0);
}

