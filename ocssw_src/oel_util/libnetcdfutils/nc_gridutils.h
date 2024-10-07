#ifndef _NC_GRIDUTILS_
#define _NC_GRIDUTILS_

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <netcdf.h>

#define BAD_VALUE -32767

#ifdef __cplusplus
extern "C" {
#endif

/* Datatype Definitions */

typedef struct grid_info_struct {
    char file[NC_MAX_NAME]; // file path
    int ncid; // netCDF file ID
    int varid; // netCDF variable ID
    char varname[NC_MAX_NAME]; // netCDF variable name
    double FillValue;
    int nodeOffset; // 0=grid-registered, 1=cell-registered

    size_t numLat;
    double startLat;
    double deltaLat;
    double validLat[2];
    int globalLat; // 0=regional, 1=global

    size_t numLon;
    double startLon;
    double deltaLon;
    double validLon[2];
    int globalLon; // 0=regional, 1=global

} grid_info_t;

typedef struct var_info_struct {
    char name[NC_MAX_NAME]; /* variable name */
    int ndims; /* number of dimensions */
    int dimids[NC_MAX_VAR_DIMS]; /* dimension IDs */
    size_t dimlen[NC_MAX_VAR_DIMS]; /* dimension lengths */
    double FillValue; /* fill value */
} var_info_t;

typedef struct grid_area_t {
    grid_info_t *grid;
    size_t numLat;
    double startLat;
    double endLat;
    size_t numLon;
    double startLon;
    double endLon;
    double *values;
} grid_area_t;

/* Function Prototypes */
void handle_nc_error(int status, char *file, int line);

var_info_t* allocate_varinfo();
int find_varid(int ncid, const char *varnames[], int *varid);
var_info_t *load_varinfo(int ncid, int varid);

grid_info_t* allocate_gridinfo();
int init_gridinfo(const char *filename, const char *varnames[], grid_info_t *grid);
int fix_latlon(float *lat, float *lon);
int latlon_findex(grid_info_t *grid, float lat, float lon, double *findex);
int latlon_index(grid_info_t *grid, float lat, float lon, size_t *index);
int get_bylatlon(grid_info_t *grid, float lat, float lon, double *value);
int get_grid_area(grid_info_t *grid,
        float north, float south, float east, float west,
        grid_area_t *area);
int interp_gridvar(grid_info_t *grid, float lat, float lon, double *result);

void print_gridinfo(grid_info_t grid);

#ifdef __cplusplus
}
#endif

#endif /* ndef _NC_GRIDUTILS_ */
