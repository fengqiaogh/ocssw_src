/**
 * @file get_nc_height.c
 * @brief Reads DEM from NetCDF4 file, and performs terrain correction.
 * @author Gwyn Fireman
 * @author Ocean Ecology Laboratory, NASA Goddard Space Flight Center
 */
#include <stdio.h>
#include <hdfi.h>
#include "nc_gridutils.h"
#include "terrain.h"
#include <genutils.h>

/** Number of tile (latitude) rows, adjustable for performance. */
#define TILE_ROWS 120
/** Maximum sensor zenith angle; avoid limb. */
#define MAXSENZEN 75
/** Degrees latitude around each tile; calculated from MAXSENZEN. */
#define BORDER 0.30
/** Earth mean radius (m). */
#define REMM 6371000.000

// only init the dem file once
static int firstCall = 1;

/* global variables */
static grid_info_t* dem_grid = {0};
/** Possible variable names for height. */
static const char* demNames[] ={"height", "z", "depth", "elevation", "altitude", NULL};
static grid_info_t* surfGrid = {0};
/** Possible variable names for non-land surface. */
static const char* surfNames[] = {"water_surface_height", NULL};

/**
 * Geometry of a single terrain (DEM) tile.
 */
typedef struct {
    int number;
    double NSEW[4]; /**< ideal boundaries of tile, including border */
    size_t numLat;
    double startLat;
    double endLat;
    double deltaLat;
    size_t numLon;
    double startLon;
    double endLon;
    double deltaLon;
    double minval;
    double maxval;
    double *height; /**< values filled as needed */
} tile_struct;

/**
 * Contains all terrain (DEM) tiles.
 */
static struct {
    size_t ntiles;
    size_t ncols[TILE_ROWS];
    size_t start_tile[TILE_ROWS];
    tile_struct *tiles;
} dem;

/**
 * Print info for a rectangular area.
 * @param[in] area
 */
void print_area(grid_area_t area) {
    fprintf(stdout, "\tnumLat   = %d\n", (int) area.numLat);
    fprintf(stdout, "\tstartLat = %8.3f\n", area.startLat);
    fprintf(stdout, "\tendLat   = %8.3f\n", area.endLat);
    fprintf(stdout, "\tnumLon   = %d\n", (int) area.numLon);
    fprintf(stdout, "\tstartLon = %8.3f\n", area.startLon);
    fprintf(stdout, "\tendLon   = %8.3f\n", area.endLon);
}

/**
 * Print info for a single tile.
 * @param[in] tile
 */
void print_tile(tile_struct tile) {
    fprintf(stdout, "tile # %d\n", tile.number);
    fprintf(stdout, "\tstartLat = %8.3f\n", tile.startLat);
    fprintf(stdout, "\tSouth    = %8.3f\n", tile.NSEW[1]);
    fprintf(stdout, "\tNorth    = %8.3f\n", tile.NSEW[0]);
    fprintf(stdout, "\tendLat   = %8.3f\n", tile.endLat);
    fprintf(stdout, "\tnumLat   = %d\n", (int) tile.numLat);
    fprintf(stdout, "\tdeltaLat = %f\n", tile.deltaLat);
    fprintf(stdout, "\tstartLon = %8.3f\n", tile.startLon);
    fprintf(stdout, "\tWest     = %8.3f\n", tile.NSEW[3]);
    fprintf(stdout, "\tEast     = %8.3f\n", tile.NSEW[2]);
    fprintf(stdout, "\tendLon   = %8.3f\n", tile.endLon);
    fprintf(stdout, "\tnumLon   = %d\n", (int) tile.numLon);
    fprintf(stdout, "\tdeltaLon = %f\n", tile.deltaLon);
    fprintf(stdout, "\tvalue (%p) range = {%d,%d}\n",
            tile.height, (int) tile.minval, (int) tile.maxval);
}

/**
 * Divide the globe into approximately equal-area chunks.
 */
void define_tile_geometry() {
    int itile = 0;
    int irow, nrows = TILE_ROWS;
    int icol, ncols, maxcols;
    double lonborder, latborder = BORDER;
    double minlat, maxlat, dlat;
    double minlon, maxlon, dlon;
    double latfrac;

    dlat = 180.0 / nrows; // latitude height for each row
    maxcols = 2 * nrows; // number of columns at equator

    /* Calculate number of columns in each row */
    for (irow = 0; irow < nrows; irow++) {

        // first tile in this row
        dem.start_tile[irow] = itile;

        // number of columns varies with latitude
        latfrac = (double) irow / (nrows - 1);
        ncols = (int) round(maxcols * sin(OEL_PI * latfrac));
        if (ncols == 0) ncols = 1;
        dem.ncols[irow] = ncols;

        // running total number of tiles
        itile += ncols;
    }
    // now itile holds total number of tiles.
    dem.ntiles = itile;

    /* Allocate space for tile array */
    dem.tiles = (tile_struct *) malloc(itile * sizeof (tile_struct));
    if (dem.tiles == NULL) {
        fprintf(stderr,
                "-E- %s:%d: Could not allocate memory for %d tiles\n",
                __FILE__, __LINE__, itile);
        exit(1);
    }

    /* Calculate extent of each tile */

    /* for each row */
    for (irow = 0; irow < nrows; irow++) {

        // nominal edge latitudes for this row
        minlat = dlat * irow - 90.0;
        maxlat = minlat + dlat;

        // adjust for border
        minlat -= latborder;
        if (minlat < -90.0) minlat = -90.0;
        maxlat += latborder;
        if (maxlat > 90.0) maxlat = 90.0;

        // longitude width for each column of this row
        ncols = dem.ncols[irow];
        dlon = 360.0 / ncols;

        // calculate longitude border for this row
        if (ncols == 1) lonborder = 0; // unless at poles
        else if (minlat < 0)
            lonborder = latborder / cos(minlat / OEL_RADEG);
        else lonborder = latborder / cos(maxlat / OEL_RADEG);

        /* for each column */
        for (icol = 0; icol < ncols; icol++) {

            // nominal edge longitudes for this column
            minlon = dlon * icol - 180.0;
            maxlon = minlon + dlon;

            // adjust for border
            minlon -= lonborder;
            maxlon += lonborder;

            // store results
            itile = dem.start_tile[irow] + icol;
            dem.tiles[itile].NSEW[0] = maxlat;
            dem.tiles[itile].NSEW[1] = minlat;
            dem.tiles[itile].NSEW[2] = maxlon;
            dem.tiles[itile].NSEW[3] = minlon;
            dem.tiles[itile].number = -1; // flag for tile loaded
            dem.tiles[itile].height = NULL; // also init values to null

        } // icol
    } // irow
}

/**
 *
 * @param[in] demfile (via input struct)
 * @return
 */
int initialize_deminfo(const char* demfile) {
    int status;
    dem.ntiles = 0;

    /* allocate and load primary elevation info */
    dem_grid = allocate_gridinfo();
    status = init_gridinfo(demfile, demNames, dem_grid);
    if (status != NC_NOERR) {
        free(dem_grid);
        dem_grid = NULL;
        fprintf(stderr, "-E- %s line %d: "
                "Unable to initialize grid for digital elevation \"%s\".\n",
                __FILE__, __LINE__, demfile);
        return status;
    }

    /* allocate and load secondary elevation info */
    surfGrid = allocate_gridinfo();
    status = init_gridinfo(demfile, surfNames, surfGrid);
    if (status != NC_NOERR) {
        free(surfGrid);
        surfGrid = NULL;
        fprintf(stderr, "-W- %s line %d: "
                "Could not read secondary grid info from \"%s\"; continuing.\n",
                __FILE__, __LINE__, demfile);
        status = NC_NOERR;
    }

    /* Initialize tile geometry */
    define_tile_geometry();

    return status;
}

/**
 *
 */
void free_deminfo() {
    int i;
    if (dem_grid != NULL) {

        /* free space allocated for tile values */
        for (i = 0; i < dem.ntiles; i++)
            if (dem.tiles[i].height != NULL) {
                free(dem.tiles[i].height);
                dem.tiles[i].height = NULL;
            }

        /* free space allocated for tile array */
        free(dem.tiles);
        dem.tiles = NULL;

        /* free surfGrid and dem_grid */
        if (surfGrid != NULL) {
            free(surfGrid);
            surfGrid = NULL;
        }
        free(dem_grid);
        dem_grid = NULL;
    } 
}

void unload_tile(int itile) {
    tile_struct *tile;
    
    tile = &dem.tiles[itile];
    free(tile->height);
    tile->height = NULL;
    tile->number = -1;
}

void load_tile_cache(int itile) {
    const int list_size = 30;
    static int tile_list[30];
    static int num_in_list = 0;

    //printf("***** load tile = %d\n", itile);

    tile_list[num_in_list++] = itile;
    if(num_in_list == list_size) {
        unload_tile(tile_list[0]);
        num_in_list--;
        for(int i=0; i<num_in_list; i++) {
            tile_list[i] = tile_list[i+1];
        }
    }
}

/**
 *
 * @param[in] dem_tile_t
 * @param[in] lat
 * @param[in] lon
 * @return
 */
int load_tile(float lat, float lon) {
    int status;
    double dlat, dlon;
    int irow, icol, itile;
    tile_struct *tile;
    grid_area_t elevArea = {0};
    grid_area_t surfArea = {0};
    int i, nvalues;
    double minlat, minlon, maxlon;
    double minval, maxval;

    /* enforce valid coordinate range */
    status = fix_latlon(&lat, &lon);
    if (status) return -1;

    /* get tile number for this lat, lon */
    dlat = 180.0 / TILE_ROWS; // latitude height for each row
    irow = (int) floor((lat + 90.) / dlat); // row number
    // make sure we don't run off the end
    if (irow == TILE_ROWS)
        irow -= 1;
    dlon = 360.0 / dem.ncols[irow]; // longitude width for this row
    if (lon == 180.0)
        icol = 0;
    else
        icol = (int) floor((lon + 180.) / dlon); // column number

    if (icol > dem.ncols[irow] - 1) {
        fprintf(stderr, "\nERROR in load_tile():\n");
        fprintf(stderr, "\ticol = %d ncols = %d\n",
                icol, (int) dem.ncols[irow]);
        return -1;
    }
    itile = dem.start_tile[irow] + icol; // tile number
    tile = &dem.tiles[itile];

    /* load only if not already in memory */
    if (tile->number != itile) {

        load_tile_cache(itile);

        // read values for tile area
        status = get_grid_area(dem_grid,
                tile->NSEW[0], tile->NSEW[1],
                tile->NSEW[2], tile->NSEW[3],
                &elevArea);
        if (status) return -1;

        if (surfGrid != NULL) {
            status = get_grid_area(surfGrid,
                    tile->NSEW[0], tile->NSEW[1],
                    tile->NSEW[2], tile->NSEW[3],
                    &surfArea);
            if (status) return -1;
        }

        // adjust boundaries as needed
        minlat = elevArea.startLat;
        minlon = elevArea.startLon;
        maxlon = elevArea.endLon;
        if (maxlon < tile->NSEW[2]) minlon += 360.0;
        if (minlon > tile->NSEW[3]) minlon -= 360.0;

        // fill tile structure
        tile->number = itile;
        tile->height = elevArea.values;

        tile->numLat = elevArea.numLat;
        tile->startLat = minlat;
        tile->deltaLat = elevArea.grid->deltaLat;
        tile->endLat = tile->startLat + tile->numLat * tile->deltaLat;
        if(irow+1 == TILE_ROWS)
            tile->endLat = 90.0;

        tile->numLon = elevArea.numLon;
        tile->startLon = minlon;
        tile->deltaLon = elevArea.grid->deltaLon;
        tile->endLon = tile->startLon + tile->numLon * tile->deltaLon;

        // step through extracted area(s)
        nvalues = tile->numLat * tile->numLon;
        minval = fabs(dem_grid->FillValue);
        maxval = -1 * minval;
        for (i = 0; i < nvalues; i++) {

            // substitute values as appropriate
            if ((surfGrid != NULL) &&
                    (fabs(surfArea.values[i] - surfGrid->FillValue) > DBL_EPSILON))
                tile->height[i] = surfArea.values[i];

            // find range
            if (tile->height[i] < minval) minval = tile->height[i];
            if (tile->height[i] > maxval) maxval = tile->height[i];
        }
        tile->minval = minval;
        tile->maxval = maxval;
        free(surfArea.values);
    }

    status = ((tile->startLat > lat) || (lat > tile->endLat) ||
            (tile->startLon > lon)); //|| (lon > tile->endLon));
    if (status) {
        fprintf(stderr, "\nERROR in load_tile():\n");
        fprintf(stderr, "lat=%7.3f  lon=%8.3f\n", lat, lon);
        fprintf(stderr, "dlat=%f irow=%d dlon=%f icol=%d itile=%d\n",
                dlat, irow, dlon, icol, itile);
        fprintf(stderr, "lat diffs: %f %f\n", lat - tile->startLat, tile->endLat - lat);
        fprintf(stderr, "lon diffs: %f %f\n", lon - tile->startLon, tile->endLon - lon);
        fprintf(stderr, "elevArea:\n");
        print_area(elevArea);
        print_tile(*tile);
        fprintf(stderr, "\n");
        //exit(1);
    }

    return itile;
}

/**
 * Floating-point index of lat, lon within current tile.
 * @param[in] tile
 * @param[in] lat,lon
 * @param[out] findex
 * @return
 */
int tile_findex(tile_struct tile, float lat, float lon, double *findex) {
    int status;
    double normLat, normLon;

    /* enforce valid coordinate range */
    findex[0] = findex[1] = -999.0;
    status = fix_latlon(&lat, &lon);
    if (status) return status; // flag NaN inputs

    /* calculate array index from normalized lat/lon */
    normLat = lat - tile.startLat;
    normLon = lon - tile.startLon;
    if (normLon < 0.0) normLon += 360.0; // adjust for dateline
    if (normLon > 360.0) normLon -= 360.0; // crossing
    findex[0] = normLat / tile.deltaLat;
    findex[1] = normLon / tile.deltaLon;

    /* cap latitude at North pole */
    if ((findex[0] >= tile.numLat) &&
            ((tile.endLat + FLT_EPSILON) > 90.0))
        findex[0] = (double) tile.numLat - FLT_EPSILON;

    /* longitudes wrap for global tile only */
    if ((findex[1] >= tile.numLon) &&
            ((tile.endLon - tile.startLon + FLT_EPSILON) > 360.0))
        findex[1] -= tile.numLon;

    /* make sure final index is inside current tile */
    status = ((0 > findex[0]) || (findex[0] >= tile.numLat) ||
            (0 > findex[1]) || (findex[1] >= tile.numLon));
    if (status) {
        fprintf(stderr, "-E- %s:%d: tile_findex():\n",
                __FILE__, __LINE__);
        fprintf(stderr, "input    lat = %7.3f   lon = %8.3f\n", lat, lon);
        fprintf(stderr, "norm     lat = %7.3f   lon = %8.3f\n", normLat, normLon);
        fprintf(stderr, "output  flat = %7.3f  flon = %8.3f\n",
                findex[0], findex[1]);
        fprintf(stderr, "      numLat = %7d   Lon = %8d\n",
                (int) tile.numLat, (int) tile.numLon);
        print_tile(tile);
        fprintf(stderr, "\n");
    }
    return status;
}

/**
 * Integer index of lat, lon within current tile.
 * @param[in] tile
 * @param[in] lat,lon
 * @param[out] index
 * @return
 */
int tile_index(tile_struct tile, float lat, float lon, size_t *index) {
    int status;
    double findex[2];

    /* calculate float array index */
    status = tile_findex(tile, lat, lon, findex);

    /* make sure max value gives valid index */
    if (findex[0] > tile.numLat - 1) findex[0] -= FLT_EPSILON;
    if (findex[1] > tile.numLon - 1) findex[1] -= FLT_EPSILON;

    /* truncate to integer */
    index[0] = (int) findex[0];
    index[1] = (int) findex[1];

    return status; // inherited from tile_findex
}

/**
 *
 * @param[in] tile
 * @param[in] lat,lon
 * @param[out] value
 * @return
 */
int get_tilevalue(tile_struct tile, float lat, float lon, double *value) {
    int status = 1;
    size_t index[2];
    if (tile.height != NULL) {
        status = tile_index(tile, lat, lon, index);
        if (!status)
            *value = tile.height[index[0] * tile.numLon + index[1]];
    }
    return status;
}

/**
 *
 * @param[in] tile
 * @param[in] lat,lon
 * @param[out] result
 * @return
 */
int interp_tilevalue(tile_struct tile, float lat, float lon, double *result) {
    int status = 1;
    size_t index[2];
    double findex[2];
    size_t x0, x1, y0, y1;
    size_t nx;
    double x, y;

    if (tile.height == NULL) return 1;

    /* find integer array index for current point */
    status = tile_index(tile, lat, lon, index);
    if (status) return status;

    /* find adjacent points */
    y0 = index[0];
    y1 = y0 + 1;
    x0 = index[1];
    x1 = x0 + 1;

    /* cap latitude at North pole */
    if (y1 == tile.numLat)
        y1 = tile.numLat - 1;

    /* longitudes wrap for global tile only */
    if ((x1 == tile.numLon) &&
            ((tile.endLon - tile.startLon + FLT_EPSILON) > 360.0))
        x1 = 0;

    /* make sure all four points within tile bounds */
    status = ((y1 > tile.numLat - 1) || (x1 > tile.numLon - 1));
    if (status) {
        tile_findex(tile, lat, lon, findex);
        fprintf(stderr, "-E- %s:%d: interp_tilevalue():\n",
                __FILE__, __LINE__);
        fprintf(stderr, "   lat = %7.3f   lon = %8.3f\n", lat, lon);
        fprintf(stderr, "  flat = %7.3f  flon = %8.3f\n",
                findex[0], findex[1]);
        fprintf(stderr, "  ilat = %7d  ilon = %8d\n", (int) y0, (int) x0);
        fprintf(stderr, "numLat = %7d   Lon = %8d\n",
                (int) tile.numLat, (int) tile.numLon);
        print_tile(tile);
        fprintf(stderr, "\n");
        return status;
    }

    /* find float array index */
    /* shift to unit coordinate system */
    tile_findex(tile, lat, lon, findex);
    y = findex[0] - index[0];
    x = findex[1] - index[1];

    /* calculate interpolated value */
    nx = tile.numLon;
    *result
            = tile.height[y0 * nx + x0] * (1 - x)*(1 - y) // x=0, y=0
            + tile.height[y0 * nx + x1] * x * (1 - y) // x=1, y=0
            + tile.height[y1 * nx + x0] * (1 - x) * y // x=0, y=1
            + tile.height[y1 * nx + x1] * x * y; // x=1, y=1

    return status;
}

/**
 * Perform terrain correction; return DEM value
 * @param[in,out] xlon,xlat
 * @param[in,out] senz
 * @param[in] sena
 * @param[out] height
 * @return
 */
int get_nc_height(const char* demfile, float *xlon, float *xlat,
        float *senz, float *sena,
        float *height) {
    int status;
    static double ddist;
    double step = 0.5;
    double tszmin = 0.01;
    double tansnz, sinsna, cossna, coslat;
    double lat, lon, dlat, dlon;
    double dd, dist;
    double dem_hgt, dem_old, los_hgt, los_old;
    int stepnum, maxsteps;
    int itile;
    tile_struct tile;

    if (firstCall) {
        /* Initialize DEM info */
        if (initialize_deminfo(demfile)) return 1;
        firstCall = 0;
    }
 
    /* Set step size in meters */
    if(ddist == 0.0)
        ddist = step * REMM * 180.0 / (dem_grid->numLat * OEL_RADEG);

    /* Load tile for this lat,lon */
    itile = load_tile(*xlat, *xlon);
    if ((itile < 0) || (itile > dem.ntiles - 1)) {
            fprintf(stderr, "\nERROR from load_tile():");
            fprintf(stderr, "\txlat=%8.3f  xlon=%8.3f  itile=%d\n",
                    *xlat, *xlon, itile);
            return 1;
        }
    tile = dem.tiles[itile];

    /* Default for sea-level tile */
    if (((int) tile.minval == 0) &&
            ((int) tile.maxval == 0)) {
        *height = 0.0f;
        return 0;
    }

    /* Skip if close to limb */
    if (*senz > MAXSENZEN) {
        status = interp_tilevalue(tile, *xlat, *xlon, &dem_hgt);
        *height = (float) dem_hgt;
        return status;
    }

    /* Pre-compute trig functions */
    tansnz = tan(*senz / OEL_RADEG);
    sinsna = sin(*sena / OEL_RADEG);
    cossna = cos(*sena / OEL_RADEG);
    coslat = cos(*xlat / OEL_RADEG);
    dlat = cossna * OEL_RADEG / REMM;
    dlon = sinsna * OEL_RADEG / (REMM * coslat);
    if(fabs(dlon) > 360.0)
        dlon = 0.0;

    /* If flat tile */
    if ((tile.maxval - tile.minval) < 1.0) {
        *height = (float) tile.maxval;
        dist = *height * tansnz * REMM / (REMM + *height);
    }
        /* Check for sensor zenith close to zero */
    else if (tansnz < tszmin) {
        status = interp_tilevalue(tile, *xlat, *xlon, &dem_hgt);
        *height = (float) dem_hgt;
        dist = *height * tansnz;
    }
        /* Else step downward along line-of-sight */
    else {

        /* Determine maximum number of steps above ellipsoid */
        maxsteps = tile.maxval * tansnz / ddist + 2;

        /* Initialize search parameters */
        dem_hgt = 0.0;
        los_hgt = 0.0;
        stepnum = maxsteps;

        /* Step downward along line-of-sight until terrain intersection */
        do {
            los_old = los_hgt;
            dem_old = dem_hgt;
            dist = stepnum * ddist;

            /* Interpolated height from DEM */
            lat = *xlat + dist * dlat;
            lon = *xlon + dist * dlon;
            status = interp_tilevalue(tile, lat, lon, &dem_hgt);
            if (status) {
                fprintf(stderr, "-E- %s:%d: get_nc_height():\n",
                        __FILE__, __LINE__);
                fprintf(stderr, "\titer = %6d\n", maxsteps - stepnum);
                return status;
            }

            /* LOS height with correction for Earth curvature effect */
            los_hgt = dist
                    * (dist * tansnz / 2.0 + REMM)
                    / (REMM * tansnz - dist);
            stepnum--;
        } while (los_hgt > dem_hgt);

        /* Interpolate to find final distance along line of sight */
        dd = (dem_hgt - los_hgt) / (dem_hgt - los_hgt + los_old - dem_old);
        dist += dd * ddist;
    }

    /* Compute adjustments to lon, lat, solar zenith and azimuth */
    /* (need to look at corrections for polar regions) */
    *xlat += (float) dist * dlat;
    *xlon += (float) dist * dlon;
    *senz -= (float) dist * OEL_RADEG / REMM;
    status = interp_tilevalue(tile, *xlat, *xlon, &dem_hgt);
    *height = (float) dem_hgt;

    /* Check for longitude transition over date line */
    if (*xlon > 180.0 || *xlon < -180.0) {
        *xlon = fmodf(*xlon, 360.0);
        if(*xlon > 180.0)
            *xlon -= 360.0;
        if(*xlon < -180.0)
            *xlon += 360.f;
    }

    // if lat goes over the pole and back down
    if(*xlat > 90.0 || *xlat < -90.0) {
        *xlat = fmodf(*xlat, 180.0);
        if(*xlat > 90.0)
            *xlat = 180.0 - *xlat;
        if(*xlat < -90.0)
            *xlat = -180.0 - *xlat;

        if(*xlon > 0.0)
            *xlon -= 180.0;
        else
            *xlon += 180.0;
    }

    return status;
}


/**
 * Return DEM value obtained by interpolating between four nearest points.
 * @param[in] xlon,xlat
 * @param[out] height
 * @return
 */
int interp_nc_height(const char* demfile, float *xlon, float *xlat, float *height) {
    int status;
    double dem_hgt;
    int itile;
    tile_struct tile;
    if (firstCall) {
        /* Initialize DEM info */
        if (initialize_deminfo(demfile)) return 1;
        firstCall = 0;
    }
    /* Load tile for this lat,lon */
    itile = load_tile(*xlat, *xlon);
    if ((itile < 0) || (itile > dem.ntiles - 1)) {
        fprintf(stderr, "\nERROR from load_tile():");
        fprintf(stderr, "\txlat=%8.3f  xlon=%8.3f  itile=%d\n",
                *xlat, *xlon, itile);
        return 1;
    }
    tile = dem.tiles[itile];
    status = interp_tilevalue(tile, *xlat, *xlon, &dem_hgt);
    *height = (float) dem_hgt;
    return status;
}
