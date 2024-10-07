/**
   @file nc_gridutils.c
   @brief Convenience utilities for handling NetCDF4 gridded variables.
   @author Gwyn Fireman
   @author Ocean Ecology Laboratory, NASA Goddard Space Flight Center
 */
#include "nc_gridutils.h"

#define NEW_CACHE_SIZE 240000000
#define NEW_CACHE_NELEMS 1001
#define NEW_CACHE_PREEMPTION .75
/* ------------------------------------------------------------------ */
/* Utilities */

/**
Print NetCDF error message and exit.
@param[in] status NetCDF error status.
@param[in] file Originating __FILE__.
@param[in] line Originating __LINE__.
 */
void handle_nc_error(int status, char *file, int line) {
    fprintf(stderr, "-E- %s:%d: %s\n",
            file, line, nc_strerror(status));
    exit(1);
}

/* ------------------------------------------------------------------ */
/* Variables */

/**
Allocate memory for structure describing NetCDF variable.
@retval var Empty structure for describing NetCDF variable.
 */
var_info_t* allocate_varinfo() {
    var_info_t* var;
    var = (var_info_t*) malloc(sizeof (var_info_t));
    if (var == NULL) {
        fprintf(stderr,
                "-E- %s:%d: Could not allocate variable info memory.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    return var;
}

/**
Read NetCDF variable attributes and load into structure.
@param[in] ncid NetCDF ID, from a previous call to nc_open(),
@param[in] varid NetCDF variable ID
@retval var Structure describing NetCDF variable.
 */
var_info_t *load_varinfo(int ncid, int varid) {
    var_info_t *var = {0};
    int idim; /* dimension ID index */

    /* initialize variable structure */
    var = allocate_varinfo();

    /* load variable info */
    nc_inq_var(ncid, varid, var->name, NULL,
            &var->ndims, var->dimids, NULL);

    /* load dimension info */
    for (idim = 0; idim < var->ndims; idim++)
        nc_inq_dimlen(ncid,
            var->dimids[idim],
            &var->dimlen[idim]);

    /* read fill value from variable attributes */
    if (nc_get_att_double(ncid, varid, "_FillValue", &var->FillValue))
        var->FillValue = BAD_VALUE; // otherwise set to default

    return var;
}

/**
Find ID of first variable found in input list.
@param[in] ncid NetCDF ID, from a previous call to nc_open(),
@param[in] varnames Null-terminated list of possible variable names.
@param[out] varid NetCDF variable ID
@return Error if no variable found.
 */
int find_varid(int ncid, const char *varnames[], int *varid) {
    int i;
    int status; /* error status */
    i = 0;
    while (varnames[i] != NULL) {
        status = nc_inq_varid(ncid, varnames[i], varid);
        if (status == NC_NOERR) break;
        i++;
    }
    return status;
}

/* ------------------------------------------------------------------ */
/* Grid */

/** Possible latitude variable names. */
static const char* latNames[] = {"y", "lat", "latitude", "Latitude", NULL};
/** Possible longitude variable names. */
static const char* lonNames[] = {"x", "lon", "longitude", "Longitude", NULL};

/**
Allocate memory for structure describing grid.
@retval grid Empty structure for describing gridded variable.
 */
grid_info_t* allocate_gridinfo() {
    grid_info_t* grid = {0};
    // allocate
    grid = (grid_info_t*) malloc(sizeof (grid_info_t));
    if (grid == NULL) {
        fprintf(stderr,
                "-E- %s:%d: Could not allocate grid info memory.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    // initialize
    grid->file[0] = '\0';
    grid->ncid = -1;
    grid->varid = -1;
    grid->varname[0] = '\0';
    return grid;
}

/**
Load information about gridded variable from a file.
@param[in] filename 
@param[in] varnames Null-terminated list of possible variable names.
@param[out] grid Structure describing grid.
@return Error status.
 */
int init_gridinfo(const char *filename, const char *varnames[], grid_info_t *grid) {
    int status;
    int ncid;
    var_info_t *var = {0};
    int varid;
    size_t index[2];
    double first, last;
    int rangeatt; // =1 if coordinate range read from attribute

    if (nc_set_chunk_cache(NEW_CACHE_SIZE, NEW_CACHE_NELEMS,
            NEW_CACHE_PREEMPTION) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_set_chunk_cache (%s) failed.\n",
                __FILE__, __LINE__, filename);
        return (1);
    }

    // open input file
    status = nc_open(filename, NC_NOWRITE, &ncid);
    if (status != NC_NOERR) return status;

    // initialize pre-allocated grid structure
    strcpy(grid->file, filename);
    grid->ncid = ncid;

    // find desired variable
    status = find_varid(ncid, varnames, &varid);
    if (status != NC_NOERR) return status;
    grid->varid = varid;

    // check dimensions
    var = load_varinfo(ncid, varid);
    if (var->ndims != 2) return 1;
    grid->numLat = var->dimlen[0];
    grid->numLon = var->dimlen[1];
    grid->FillValue = var->FillValue;
    strcpy(grid->varname, var->name);
    free(var);

    /* --- read grid registration type from global attributes --- */
    status = nc_get_att_int(ncid, NC_GLOBAL, "node_offset", &grid->nodeOffset);
    if (status != NC_NOERR)
        grid->nodeOffset = 0; // otherwise assume grid-registered

    // otherwise assume odd # of cells means grid-registered
    // ( see http://www.ngdc.noaa.gov/mgg/global/gridregistration.html )
    //        grid->nodeOffset = (grid->numLon % 2 == 0);  // 0=grid, 1=cell


    /* --- read LATITUDE limits from global attributes --- */
    rangeatt = 1;
    if (nc_get_att_double(ncid, NC_GLOBAL, "lower_lat", &first) ||
            nc_get_att_double(ncid, NC_GLOBAL, "upper_lat", &last)) {

        // otherwise read array values
        rangeatt = 0;
        if (find_varid(ncid, latNames, &varid)) {
            fprintf(stderr, "-E- %s:%d: "
                    "No latitude attributes or variables found in %s.",
                    __FILE__, __LINE__, grid->file);
            return 1;
        }
        // check dimensions
        var = load_varinfo(ncid, varid);
        if ((var->ndims != 1) || (var->dimlen[0] != grid->numLat)) {
            fprintf(stderr, "-E- %s:%d: "
                    "Latitude array %s has wrong dimensions.",
                    __FILE__, __LINE__, var->name);
            return 1;
        }
        // read first and last value of array
        index[0] = 0;
        index[1] = var->dimlen[0] - 1;
        status = nc_get_var1_double(ncid, varid, &index[0], &first);
        status = nc_get_var1_double(ncid, varid, &index[1], &last);
        free(var);
    }

    /* calculate step size */
    if (rangeatt && grid->nodeOffset)
        grid->deltaLat = (last - first) / (double) (grid->numLat);
    else
        grid->deltaLat = (last - first) / (double) (grid->numLat - 1);

    /* extend range to cover outer half-pixels */
    first -= grid->deltaLat / 2.0;
    last  += grid->deltaLat / 2.0;

    /* save parameters */
    grid->startLat = first;
    grid->validLat[0] = (first < last) ? first : last;
    grid->validLat[1] = (first < last) ? last : first;
    grid->globalLat = fabs(last - first) >= 180.0;

    /* --- read LONGITUDE limits from global attributes --- */
    rangeatt = 1;
    if (nc_get_att_double(ncid, NC_GLOBAL, "left_lon", &first) ||
            nc_get_att_double(ncid, NC_GLOBAL, "right_lon", &last)) {

        // otherwise read array values
        rangeatt = 0;
        if (find_varid(ncid, lonNames, &varid)) {
            fprintf(stderr, "-E- %s:%d: "
                    "No longitude attributes or variables found in %s.",
                    __FILE__, __LINE__, grid->file);
            return 1;
        }
        // check dimensions
        var = load_varinfo(ncid, varid);
        if ((var->ndims != 1) || (var->dimlen[0] != grid->numLon)) {
            fprintf(stderr, "-E- %s:%d: "
                    "Longitude array %s has wrong dimensions.",
                    __FILE__, __LINE__, var->name);
            return 1;
        }
        // read first and last value of array
        index[0] = 0;
        index[1] = var->dimlen[0] - 1;
        status = nc_get_var1_double(ncid, varid, &index[0], &first);
        status = nc_get_var1_double(ncid, varid, &index[1], &last);
        free(var);
    }

    /* calculate step size */
    if (first > last) last += 360.0;
    if (rangeatt && grid->nodeOffset)
        grid->deltaLon = (last - first) / (double) (grid->numLon);
    else
        grid->deltaLon = (last - first) / (double) (grid->numLon - 1);


    /* extend range to cover outer half-pixels */
    first -= grid->deltaLon / 2.0;
    last  += grid->deltaLon / 2.0;

    /* save parameters */
    grid->startLon = first;
    grid->validLon[0] = first;
    grid->validLon[1] = last;
    grid->globalLon = (last - first) >= 360.0;

    return 0;
}

/**
Force coordinates into range:
  {-90 <= lat < 90} and {-180 <= lon < 190}
@param[in,out] lat, lon Desired geographic coordinates.
@return Error if either coordinate is NaN.
 */
int fix_latlon(float *lat, float *lon) {
    int status;

    /* check for valid input */
    status = (isnan(*lat) || isnan(*lon));
    if (status) return status;

    /* latitude changes direction as you cross the poles. */
    while (*lat < -90.0) {
        *lat = -180.0 - *lat;
        //*lon += 180.0;  // and you cross to the opposite longitude
    }
    while (*lat > 90.0) {
        *lat = 180.0 - *lat;
        //*lon += 180.0;
    }

    /* longitude keeps wrapping around the globe. */
    while (*lon < -180.0) *lon += 360.0;
    while (*lon > 180.0) *lon -= 360.0;

    return status;
}

/**
Find floating point index of specified coordinates within gridded variable.
@param[in] grid Structure describing grid.
@param[in] lat, lon Desired geographic coordinates.
@param[out] findex Indices of [y,x] corresponding to [lat,lon].
@return Error if coordinates outside of grid range.
 */
int latlon_findex(grid_info_t *grid, float lat, float lon, double *findex) {
    int status; // 0=good, 1=bad
    double normLat, normLon;

    /* enforce valid coordinate range */
    findex[0] = findex[1] = -999.0;
    status = fix_latlon(&lat, &lon);
    if (status) return status;

    /* make sure final lat/lon are in grid range. */
    status = ((grid->validLat[0] > lat) || (lat > grid->validLat[1]) ||
            (grid->validLon[0] > lon) || (lon > grid->validLon[1]));
    if (status)
        return status;

    /* calculate array index from normalized lat/lon */
    normLat = lat - grid->startLat;
    normLon = lon - grid->startLon;
    findex[0] = normLat / grid->deltaLat;
    findex[1] = normLon / grid->deltaLon;
    return status;
}

/**
Find integer index of specified coordinates within gridded variable.
@param[in] grid Structure describing grid.
@param[in] lat, lon Desired geographic coordinates.
@param[out] index Indices of [y,x] corresponding to [lat,lon].
@return Error if coordinates outside of grid range. (0=good)
 */
int latlon_index(grid_info_t *grid, float lat, float lon, size_t *index) {
    int status; // 0=good, 1=bad
    double findex[2];

    /* calculate float array index */
    status = latlon_findex(grid, lat, lon, findex);
    if (status) return status;

    /* make sure max value gives valid index */
    if (lat == grid->validLat[1]) findex[0] -= FLT_EPSILON;
    if (lon == grid->validLon[1]) findex[1] -= FLT_EPSILON;

    /* truncate to integer */
    index[0] = (size_t) findex[0];
    index[1] = (size_t) findex[1];

    return status; // inherited from latlon_findex
}

/**
Read gridded variable value nearest to specified coordinates.
@param[in] grid Structure describing grid.
@param[in] lat, lon Desired geographic coordinates.
@param[out] value Pointer to value, or FillValue if not found.
@return Error if coordinates outside of grid range. (0=good)
 */
int get_bylatlon(grid_info_t *grid, float lat, float lon, double *value) {
    int status = 1;
    size_t index[2];
    if (grid != NULL) {
        *value = grid->FillValue;
        status = latlon_index(grid, lat, lon, index);
        if (!status) {
            status = nc_get_var1_double(grid->ncid,
                    grid->varid,
                    index, value);
            if (status != NC_NOERR)
                handle_nc_error(status, __FILE__, __LINE__);
            // test for fill value
            status = fabs(*value - grid->FillValue) < DBL_EPSILON;
        }
    }
    return status;
}

/**
Determine whether grid has extra longitude column.
Needed when extracting area crossing +/- 180 degrees.
@param[in] grid Structure describing grid.
@return TRUE: has extra longitude column.
@return FALSE: no extra longitude column.
 */
int has_extra_column(grid_info_t *grid) {
    return ( (grid->validLon[1] - grid->validLon[0])
            > (360.0 + grid->deltaLon / 2.0));
}

/**
Read gridded variable values covering specified boundaries.
@param[in] grid Structure describing grid.
@param[in] north,south,east,west Geographic boundaries of desired area.
@param[out] area Structure contains and describes extracted area.
@return Error if any boundary set outside of grid range.
 */
int get_grid_area(grid_info_t *grid,
        float north, float south,
        float east, float west,
        grid_area_t *area) {
    int status = 0;
    size_t SWindex[2], NEindex[2];
    size_t start[] = {0, 0};
    size_t count[] = {2, 2};
    double *values;
    size_t iy, ix, ny, nx;
    size_t nxWest, nxEast;
    double *westvals, *eastvals;
    double *outrow, *inrow;

    /* make sure north > south */
    if (north < south) { // never extract > half the globe
        float coord = north;
        north = south;
        south = coord;
    }

    /* pad north & south a little */
    if (south < (grid->validLat[0] + grid->deltaLat)) south = grid->validLat[0];
    if (north > (grid->validLat[1] - grid->deltaLat)) north = grid->validLat[1];

    /* get subset indices */
    status = latlon_index(grid,
            south - FLT_EPSILON,
            west - FLT_EPSILON,
            SWindex);
    if (status) return status;
    status = latlon_index(grid,
            north + FLT_EPSILON,
            east + FLT_EPSILON,
            NEindex);
    if (status) return status;

    start[0] = SWindex[0]; // set latitude values
    ny = count[0] = NEindex[0] - SWindex[0] + 1;

    /* simplest case: east > west */
    if (NEindex[1] > SWindex[1]) {
        nx = NEindex[1] - SWindex[1] + 1;
        if ((nx == grid->numLon) && (has_extra_column(grid)))
            nx -= 1; // skip duplicate end column

        start[1] = SWindex[1];
        count[1] = nx;
        values = (double*) malloc(ny * nx * sizeof (double));
        status = nc_get_vara_double(grid->ncid, grid->varid,
                start, count, values);
        if (status != NC_NOERR) handle_nc_error(status, __FILE__, __LINE__);
    }
        /* crossing the dateline: west > east */
    else {
        nxWest = grid->numLon - SWindex[1];
        if (has_extra_column(grid))
            nxWest -= 1; // skip duplicate end column
        nxEast = NEindex[1] + 1;
        nx = nxWest + nxEast;

        /* west side */
        start[1] = SWindex[1];
        count[1] = nxWest;
        westvals = (double*) malloc(ny * nxWest * sizeof (double));
        status = nc_get_vara_double(grid->ncid, grid->varid,
                start, count, westvals);
        if (status != NC_NOERR) handle_nc_error(status, __FILE__, __LINE__);

        /* east side */
        start[1] = 0;
        count[1] = nxEast;
        eastvals = (double*) malloc(ny * nxEast * sizeof (double));
        status = nc_get_vara_double(grid->ncid, grid->varid,
                start, count, eastvals);
        if (status != NC_NOERR) handle_nc_error(status, __FILE__, __LINE__);

        /* merge west and east */
        values = (double*) malloc(ny * nx * sizeof (double));
        for (iy = 0; iy < ny; iy++) {
            /* copy west values into first nxWest columns */
            outrow = values + iy * nx;
            inrow = westvals + iy * nxWest;
            for (ix = 0; ix < nxWest; ix++) outrow[ix] = inrow[ix];

            /* copy east values into last nxEast columns */
            outrow += nxWest;
            inrow = eastvals + iy * nxEast;
            for (ix = 0; ix < nxEast; ix++) outrow[ix] = inrow[ix];
        }
        free(westvals);
        free(eastvals);
    }

    /* calculate actual boundaries of extracted area */
    south = SWindex[0] * grid->deltaLat + grid->startLat;
    north = south + ny * grid->deltaLat;
    west = SWindex[1] * grid->deltaLon + grid->startLon;
    east = west + nx * grid->deltaLon;

    /* fill output structure */
    area->grid = grid;
    area->numLat = ny;
    area->startLat = south;
    area->endLat = north;
    area->numLon = nx;
    area->startLon = west;
    area->endLon = east;
    area->values = values;

    return status;
}

/**
Get bilinear interpolation of any gridded variable.
@param[in] grid Structure describing grid.
@param[in] lat, lon Desired geographic coordinates.
@param[out] result Pointer to interpolated value, or FillValue if not found.
@return Error if coordinates outside of grid range.
@remark
  nc_get_var1_double will automatically convert native type to double - see
  http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Type-Conversion.html
 */
int interp_gridvar(grid_info_t *grid, float lat, float lon, double *result) {
    int status = 1;
    double findex[2];
    size_t index[] = {0, 0};
    size_t x0, x1, y0, y1;
    double values[4];
    double x, y;

    *result = grid->FillValue;

    /* find integer array index for current point */
    status = latlon_index(grid, lat, lon, index);
    if (status) return status;
    y0 = index[0];
    x0 = index[1];

    /* find integer array index for adjacent points */
    status = latlon_index(grid, lat + grid->deltaLat, lon + grid->deltaLon, index);
    y1 = (index[0] < grid->numLat) ? index[0] : y0;
    x1 = (index[1] < grid->numLon) ? index[1] : x0;

    /* retrieve nearest 4 data points */
    index[0] = y0;
    index[1] = x0; // x=0, y=0
    status = nc_get_var1_double(grid->ncid, grid->varid, index, &values[0]);
    if (status != NC_NOERR) handle_nc_error(status, __FILE__, __LINE__);

    index[0] = y0;
    index[1] = x1; // x=1, y=0
    status = nc_get_var1_double(grid->ncid, grid->varid, index, &values[1]);
    if (status != NC_NOERR) handle_nc_error(status, __FILE__, __LINE__);

    index[0] = y1;
    index[1] = x0; // x=0, y=1
    status = nc_get_var1_double(grid->ncid, grid->varid, index, &values[2]);
    if (status != NC_NOERR) handle_nc_error(status, __FILE__, __LINE__);

    index[0] = y1;
    index[1] = x1; // x=1, y=1
    status = nc_get_var1_double(grid->ncid, grid->varid, index, &values[3]);
    if (status != NC_NOERR) handle_nc_error(status, __FILE__, __LINE__);

    /* find float array index */
    /* shift to unit coordinate system */
    status = latlon_findex(grid, lat, lon, findex);
    y = findex[0] - y0;
    x = findex[1] - x0;

    /* calculate interpolated value */
    /* NOTE: need to fix interp_gridvar to not include FillValue in avg! */
    *result
            = values[0] * (1 - x)*(1 - y) // x=0, y=0
            + values[1] * x * (1 - y) // x=1, y=0
            + values[2] * (1 - x) * y // x=0, y=1
            + values[3] * x * y; // x=1, y=1

    return status;
}

/**
Print grid description structure.
@param[in] grid Structure describing grid.
 */
void print_gridinfo(grid_info_t grid) {
    fprintf(stdout, "file\t= %s\n", grid.file);
    fprintf(stdout, "ncid\t= %d\n", grid.ncid);
    fprintf(stdout, "varid\t= %d\n", grid.varid);
    fprintf(stdout, "varname\t= %s\n", grid.varname);
    fprintf(stdout, "FillValue\t= %f\n", grid.FillValue);
    fprintf(stdout, "nodeOffset\t= %d\n", grid.nodeOffset);

    fprintf(stdout, "numLat\t= %d\n", (int) grid.numLat);
    fprintf(stdout, "startLat\t= %f\n", grid.startLat);
    fprintf(stdout, "deltaLat\t= %.12f\n", grid.deltaLat);
    fprintf(stdout, "validLat\t= [%f,%f]\n", grid.validLat[0], grid.validLat[1]);
    fprintf(stdout, "globalLat\t= %d\n", grid.globalLat);

    fprintf(stdout, "numLon\t= %d\n", (int) grid.numLon);
    fprintf(stdout, "startLon\t= %f\n", grid.startLon);
    fprintf(stdout, "deltaLon\t= %.12f\n", grid.deltaLon);
    fprintf(stdout, "validLon\t= [%f,%f]\n", grid.validLon[0], grid.validLon[1]);
    fprintf(stdout, "globalLon\t= %d\n", grid.globalLon);
}
