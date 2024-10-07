#ifndef TERRAIN_H
#define TERRAIN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Terrain correction function prototypes */
int get_nc_height(
        const char *demfile,
        float *lon, float *lat,
        float *senz, float *sena,
        float *height);
int get_dem_height(
        const char *demfile,
        float *lon, float *lat,
        float *senz, float *sena,
        float *height);
/* DEM interpolation function prototypes */
int interp_nc_height(const char* demfile, float *xlon, float *xlat, float *height);
int interp_dem_height(const char* demfile, float *xlon, float *xlat, float *height);

#ifdef __cplusplus
}
#endif

#endif // TERRAIN_H
