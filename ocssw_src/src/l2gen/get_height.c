#include "l12_proto.h"
#include <terrain.h>
#define MAX_LAT 89.95

static int (*correct_terrain)(
        const char* demfile,
        float *lon, float *lat,
        float *senz, float *sena,
        float *height);
static int (*interp_height)(const char* demfile, float *xlon, float *xlat, float *height);

/**
 * Load DEM height for one pixel
 * @param[in] demfile
 * @param[in,out] l1str
 * @param[in] ip
 * @param[in] terrain_corrected
 * @return
 */
int get_height(l1str *l1rec, int32_t ip, int terrain_corrected) {
    static int firstCall = 1;
    int status = 1;
    float *xlon = &l1rec->lon [ip];
    float *xlat = &l1rec->lat [ip];
    float *senz = &l1rec->senz[ip];
    float *sena = &l1rec->sena[ip];
    float *height = &l1rec->height[ip];
    /*
     *  for band-dependent view angles, we currently do not have a 
     *  band-dependent height, which would be the exact treatment.  
     *  Seeing as the view angle change is small, the height should not 
     *  vary enough to be of concern.  However, the nominal sensor zenith 
     *  angle gets a correction and that amount of correction should be 
     *  applied to the band-dependent sensor zenith angles, if applicable
     *  (see code after call to correct_terrain)
     */
    float senz_sav;

    /* Initial file tests */
    if (firstCall) {
        int ncid;
        int netcdf_dem;
        firstCall = 0;

        /* input file defined? */
        if (input->demfile == NULL || input->demfile[0] == 0) {
            fprintf(stderr, "-E- %s line %d: "
                    "Elevation file is NULL.\n",
                    __FILE__, __LINE__);
            return 1;
        }
        printf("Loading DEM info from %s\n", input->demfile);

        /* test for NetCDF input file */
        if(Hishdf(input->demfile))
	        status = NC2_ERR;
	    else
	        status = nc_open(input->demfile, NC_NOWRITE, &ncid);
        
        netcdf_dem = (status == NC_NOERR);
        if (netcdf_dem) nc_close(ncid);

        /* set function pointers according to input file type */
        if (netcdf_dem) {
            interp_height = interp_nc_height;
            correct_terrain = get_nc_height;
        } else {
            interp_height = interp_dem_height;
            correct_terrain = get_dem_height;
        }
    }

    /* Interpolate DEM height if terrain correction already done,
       or target too close to poles */
    if (terrain_corrected
            || (fabs((double) *xlat) > MAX_LAT)) {
        status = interp_height(input->demfile, xlon, xlat, height);
        if (status) {
            fprintf(stderr, "-E- %s line %d: interp_height():\n",
                    __FILE__, __LINE__);
            fprintf(stderr,
                    "xlon=%f xlat=%f height=%f\n",
                    (double) *xlon, (double) *xlat, (double) *height);
        }
    }
        /* Otherwise, do terrain correction */
    else {
        senz_sav = *senz;
        status = correct_terrain(input->demfile, xlon, xlat, senz, sena, height);
        if (status) {
            fprintf(stderr, "-E- %s line %d: correct_terrain():\n",
                    __FILE__, __LINE__);
            fprintf(stderr,
                    "xlon=%f xlat=%f senz=%f sena=%f height=%f\n",
                    (double) *xlon, (double) *xlat,
                    (double) *senz, (double) *sena, (double) *height);
        }
        if (l1rec->geom_per_band != NULL) {
            float senz_corr = *senz - senz_sav;
            int32_t nwave = l1rec->l1file->nbands, iw;
            for (iw = 0; iw < nwave; iw++) {
                l1rec->geom_per_band->senz[ ip * nwave + iw ] += senz_corr;
            }
        }
    }
    return status;
}
