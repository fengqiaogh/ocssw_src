#include <stdio.h>
#include <math.h>
#include <hdfi.h>
#include <mfhdf.h>
#include <genutils.h>
#include "terrain.h"

#define RADEG 57.29577951

/*
  Exactly replicate the functionality of get_dem_height.f
  Originally written by: Frederick S. Patt, SAIC GSC, June 1, 2000
  Translated by: Gwyn F. Fireman, SAIC, July 2013 (with help from f2c)
 */

/* DEM structure (from dem_s.fin) */
struct {
    int32 nrows; // Number of rows in grid
    int32 neq; // Number of tiles in equatorial row
    int32 ntiles; // Number of tiles in grid
    int32 nftiles; // Number of filled tiles
    int32 isize; // Tile array size
    int32 itsize; // Array size of nominal area
    int32 ioff; // Offset of nominal tile in array
    int16 ncols[121]; // Number of tiles in each row
    int16 nstart[121]; // Start tile # for each row (0-base)
    int16 irecno[18338]; // Record no. for tile (0-base; -1 if not filled)
    int16 iflag[18338]; // Tile flag (0=sea, 1=flag, 2=filled)
    int16 ihmin[18338]; // Minimum terrain height in tile
    int16 ihmax[18338]; // Maximum terrain height in tile
} dem;

/* Function prototypes */
float bilin(float *, float *, int32 *, int16 *, int32 *);
int initialize_dem(const char *, int32 *);

/* get_dem_height */
int get_dem_height(const char *demfile,
        float *xlon, float *xlat,
        float *senz, float *sena,
        float *height) {
    /* Local variables */
    int32 status = 0;
    static float tilat, tilon, tslon, tslat, tansnz, sinsna, cossna;
    static float x0, y0, x, y, dx, dy, dd, dist, distn, ddist;
    static float coslat, hgt1, hgt2, los1, los2;
    static int16 npts, iflag;
    static int16 i, irow, icol, itile;
    static int16* tile;
    static int32 tileid;

    /* Initialized data */
    static int firstCall = 1;
    static int32 idims[3] = {1, 220, 220};
    static int32 istart[3] = {0, 0, 0};
    static int32 itilem = -999;
    static float step = .5f;
    static float tszmin = .01f;
    static float snzmax = 56.f;
    static float remm = 6371000.f; // Earth mean radius in meters

    /* Check if file is open */
    if (firstCall) {
        firstCall = 0;

        tile = (int16*) allocateMemory(48400 * sizeof (int16), "get_dem_height - tile");

        /* Open DEM file and fill structure */
        if (initialize_dem(demfile, &tileid) != 0) {
            return 1;
        }

        /* Compute constants used for calculations */
        /* Tile size in latitude */
        tslat = 180.f / (dem.nrows - 1);

        /* Nominal tile grid point spacing in meters */
        distn = remm * tslat / (dem.itsize * RADEG);

        /* Set step size in meters */
        ddist = step * distn;
    }

    /* Compute tile # */
    irow = (*xlat - 0.0001f + 90.f) / tslat + .5f;
    icol = (*xlon - 0.0001f + 180.f) / 360.f * dem.ncols[irow];
    itile = dem.nstart[irow] + icol;
    iflag = dem.iflag[itile];

    /* Set nominal values (sea level) */
    *height = 0.f;
    dist = 0.f;

    /* If sea level, return input values */
    if (iflag == 0) {
        return 0;
    }

    /* Compute needed trig functions of angles */
    tansnz = tan(*senz / RADEG);
    sinsna = sin(*sena / RADEG);
    cossna = cos(*sena / RADEG);
    coslat = cos(*xlat / RADEG);

    /* If flat tile */
    if (iflag == 1) {
        npts = 2;
        *height = (float) dem.ihmax[itile];
        dist = *height * tansnz * remm / (remm + *height);

        /* Else terrain tile */
    } else {

        /* Check for tile already in memory */
        if (itile != itilem) {

            /* Read tile array from file */
            istart[0] = dem.irecno[itile];
            //            printf("Tile %d record %d\n",itile,istart[0]);
            itilem = itile;
            if (SDreaddata(tileid, istart, NULL, idims, tile) != 0) {
                printf("Error reading tile %d record %d\n",
                        itile, dem.irecno[itile]);
                return 1;
            }

            /* Compute longitude width of tile */
            tslon = 360.f / dem.ncols[irow];

            /* Compute coordinates of lower-left corner of nominal tile */
            tilat = (irow - .5f) * tslat - 90.f;
            tilon = icol * tslon - 180.f;
        }

        /* Get terrain heights under line-of-sight */

        /* Determine maximum number of points needed */
        npts = dem.ihmax[itile] * tansnz / ddist + 2;

        /* Determine location of input lon/lat on tile grid */
        /*  (recall that points are bin centers) */
        x0 = (*xlon - tilon) / tslon * dem.itsize + dem.ioff + .5f;
        y0 = (*xlat - tilat) / tslat * dem.itsize + dem.ioff + .5f;

        /* Get interpolated height at input lon/lat */
        *height = bilin(&x0, &y0, &dem.isize, tile, &status);
        dist = 0.f;

        /* Check for zero height and low-slope tile */
        if (*height == 0.f && iflag == 2) {
            return 0;
        }
        /* Compute grid step size */
        dx = step * sinsna;
        dy = step * cossna;

        /* Correct horizontal step size for latitude and tile size */
        dx = dx * tslat / (coslat * tslon);

        /* Find intersection point */

        /*  Check for sensor zenith close to zero */
        if (tansnz < tszmin) {
            dist = *height * tansnz;

            /*  Else step downward along line-of-sight */
        } else {

            /*  Initialize search parameters */
            hgt1 = 0.f;
            hgt2 = 0.f;
            i = npts;
            dist = i * ddist;

            /*  Compute LOS height with correction for Earth curvature effect */
            los1 = dist * (dist * tansnz / 2.f + remm)
                    / (remm * tansnz - dist);

            /*  Continue until line-of-sight crosses terrain */
            while (los1 > hgt1) {
                --i;
                los2 = los1;
                hgt2 = hgt1;
                dist = i * ddist;
                x = x0 + i * dx;
                y = y0 + i * dy;
                hgt1 = bilin(&x, &y, &dem.isize, tile, &status);

                /*  Check for current point outside tile area and print warning message */
                if (status != 0 && i < 1) {
                    printf("Warning: %d %d outside tile range\n", (int) x, (int) y);
                }
                los1 = dist * (dist * tansnz / 2.f + remm)
                        / (remm * tansnz - dist);
            }

            /*  Interpolate to find intersection point */
            dd = (hgt1 - los1) / (hgt1 - los1 + los2 - hgt2);
            dist += dd * ddist;
            x += dd * dx;
            y += dd * dy;
            hgt1 = bilin(&x, &y, &dem.isize, tile, &status);

            /*  Check for final height outside tile area */
            if (status != 0) {

                /*    If within GAC swath, return with error */
                if (*senz <= snzmax) {
                    printf("Bilinear interpolation error %d %d %d\n",
                            (int) x, (int) y, (int) dem.isize);
                    return 0;
                } else {
                    status = 0;
                }

            } else {
                *height = hgt1;
            }
        }

        /* Compute adjustments to lon, lat, solar zenith and azimuth */
        /*  (need to look at corrections for polar regions) */
        *xlon += dist * sinsna * RADEG / (remm * coslat);
        *xlat += dist * cossna * RADEG / remm;
        *senz -= dist * RADEG / remm;

        /* Check for longitude transition over date line */
        if (*xlon > 180.f) {
            *xlon -= 360.f;
        }
        if (*xlon < -180.f) {
            *xlon += 360.f;
        }
    }

    return 0;
} /* get_dem_height */


#define tile_ref(a_1,a_2) tile[(a_2)*tile_dim1 + a_1]

float bilin(float *x, float *y, int32 *isize, int16 *tile, int32 *status) {
    /* System generated locals */
    int tile_dim1, tile_offset;
    float ret_val;

    /* Local variables */
    static float dx, dy;
    static int32 ix, iy;

    /* Parameter adjustments */
    tile_dim1 = *isize;
    tile_offset = 1 + tile_dim1;
    tile -= (int32) tile_offset;

    /* Compute indices */
    ix = *x;
    iy = *y;

    /* Check indices against array size */
    if ((ix > 0) && (ix < *isize) && (iy > 0) && (iy < *isize)) {

        dx = *x - ix;
        dy = *y - iy;

        ret_val = (1 - dx) * (1 - dy) * tile_ref(ix, iy)
                + dx * (1 - dy) * tile_ref(ix + 1, iy)
                + (1 - dx) * dy * tile_ref(ix, iy + 1)
                + dx * dy * tile_ref(ix + 1, iy + 1);
        status[0] = 0;

    } else {
        ret_val = 0.f;
        status[0] = -1;
    }

    return ret_val;
} /* bilin */
#undef tile_ref

int initialize_dem(const char *demfile, int32 *tileid) {
    /* Initialized data */
    int32 status = 0;
    static int32 istart = 0;
    static int32 istr = 1;

    /* Local variables */
    static int32 ind;
    static int32 at_id, sd_id, idims, fileid;

    status = 0;
    tileid[0] = -1;

    /* Open file */
    fileid = SDstart(demfile, DFACC_RDONLY);
    if (fileid == -1) {
        printf("Error opening DEM file: %s\n", demfile);
        status = -1;
        return 0;
    }

    /* Read file attributes */

    at_id = SDfindattr(fileid, "Number of rows");
    if (at_id == -1) {
        printf("Error getting index for Number of rows\n");
        return 0;
    }
    status = SDreadattr(fileid, at_id, &dem.nrows);
    if (status == -1) {
        printf("Error getting value for Number of rows\n");
        return 0;
    }

    at_id = SDfindattr(fileid, "Equator tiles");
    if (at_id == -1) {
        printf("Error getting index for Equator tiles\n");
        return 0;
    }
    status = SDreadattr(fileid, at_id, &dem.neq);
    if (status == -1) {
        printf("Error getting value for Equator tiles\n");
        return 0;
    }

    at_id = SDfindattr(fileid, "Array size");
    if (at_id == -1) {
        printf("Error getting index for Array size\n");
        return 0;
    }
    status = SDreadattr(fileid, at_id, &dem.isize);
    if (status == -1) {
        printf("Error getting value for Array size\n");
        return 0;
    }

    at_id = SDfindattr(fileid, "Tile size");
    if (at_id == -1) {
        printf("Error getting index for Tile size\n");
        return 0;
    }
    status = SDreadattr(fileid, at_id, &dem.itsize);
    if (status == -1) {
        printf("Error getting value for Tile size\n");
        return 0;
    }

    at_id = SDfindattr(fileid, "Array tile offset");
    if (at_id == -1) {
        printf("Error getting index for Array tile offset\n");
        return 0;
    }
    status = SDreadattr(fileid, at_id, &dem.ioff);
    if (status == -1) {
        printf("Error getting value for Array tile offset\n");
        return 0;
    }

    at_id = SDfindattr(fileid, "Number of tiles");
    if (at_id == -1) {
        printf("Error getting index for Number of tiles\n");
        return 0;
    }
    status = SDreadattr(fileid, at_id, &dem.ntiles);
    if (status == -1) {
        printf("Error getting value for Number of tiles\n");
        return 0;
    }

    at_id = SDfindattr(fileid, "Number of filled tiles");
    if (at_id == -1) {
        printf("Error getting index for Number of filled tiles\n");
        return 0;
    }
    status = SDreadattr(fileid, at_id, &dem.nftiles);
    if (status == -1) {
        printf("Error getting value for Number of filled tiles\n");
        return 0;
    }

    /* Read index arrays */

    idims = dem.nrows;

    ind = SDnametoindex(fileid, "Row_number_of_tiles");
    if (ind == -1) {
        printf("Error getting index for Row_number_of_tiles\n");
        return 0;
    }
    sd_id = SDselect(fileid, ind);
    if (sd_id == -1) {
        printf("Error selecting Row_number_of_tiles\n");
        return 0;
    }
    status = SDreaddata(sd_id, &istart, &istr, &idims, dem.ncols);
    if (status == -1) {
        printf("Error reading Row_number_of_tiles\n");
        return 0;
    }
    status = SDendaccess(sd_id);

    ind = SDnametoindex(fileid, "Row_start_tile");
    if (ind == -1) {
        printf("Error getting index for Row_start_tile\n");
        return 0;
    }
    sd_id = SDselect(fileid, ind);
    if (sd_id == -1) {
        printf("Error selecting Row_start_tile\n");
        return 0;
    }
    status = SDreaddata(sd_id, &istart, &istr, &idims, dem.nstart);
    if (status == -1) {
        printf("Error reading Row_start_tile\n");
        return 0;
    }
    status = SDendaccess(sd_id);

    idims = dem.ntiles;

    ind = SDnametoindex(fileid, "DEM_tile_record");
    if (ind == -1) {
        printf("Error getting index for DEM_tile_record\n");
        return 0;
    }
    sd_id = SDselect(fileid, ind);
    if (sd_id == -1) {
        printf("Error selecting DEM_tile_record\n");
        return 0;
    }
    status = SDreaddata(sd_id, &istart, &istr, &idims, dem.irecno);
    if (status == -1) {
        printf("Error reading DEM_tile_record\n");
        return 0;
    }
    status = SDendaccess(sd_id);

    ind = SDnametoindex(fileid, "DEM_tile_flag");
    if (ind == -1) {
        printf("Error getting index for DEM_tile_flag\n");
        return 0;
    }
    sd_id = SDselect(fileid, ind);
    if (sd_id == -1) {
        printf("Error selecting DEM_tile_flag\n");
        return 0;
    }
    status = SDreaddata(sd_id, &istart, &istr, &idims, dem.iflag);
    if (status == -1) {
        printf("Error reading DEM_tile_flag\n");
        return 0;
    }
    status = SDendaccess(sd_id);

    ind = SDnametoindex(fileid, "Tile_minimum_height");
    if (ind == -1) {
        printf("Error getting index for Tile_minimum_height\n");
        return 0;
    }
    sd_id = SDselect(fileid, ind);
    if (sd_id == -1) {
        printf("Error selecting Tile_minimum_height\n");
        return 0;
    }
    status = SDreaddata(sd_id, &istart, &istr, &idims, dem.ihmin);
    if (status == -1) {
        printf("Error reading Tile_minimum_height\n");
        return 0;
    }
    status = SDendaccess(sd_id);

    ind = SDnametoindex(fileid, "Tile_maximum_height");
    if (ind == -1) {
        printf("Error getting index for Tile_maximum_height\n");
        return 0;
    }
    sd_id = SDselect(fileid, ind);
    if (sd_id == -1) {
        printf("Error selecting Tile_maximum_height\n");
        return 0;
    }
    status = SDreaddata(sd_id, &istart, &istr, &idims, dem.ihmax);
    if (status == -1) {
        printf("Error reading Tile_maximum_height\n");
        return 0;
    }
    status = SDendaccess(sd_id);

    /* Get SDS ID for tile array */
    ind = SDnametoindex(fileid, "DEM_tile_data");
    if (ind == -1) {
        printf("Error getting index for DEM_tile_data\n");
        return 0;
    }
    tileid[0] = SDselect(fileid, ind);
    if (tileid[0] == -1) {
        printf("Error selecting DEM_tile_flag\n");
    }
    return 0;
} /* initialize_dem */

/**************************************************************************/

int interp_dem_height(const char *demfile,
        float *xlon, float *xlat,
        float *height) {
    /* Local variables */
    int32 status = 0;
    static float tilat, tilon, tslon, tslat;
    static float x0, y0;
    static int16 iflag;
    static int16 irow, icol, itile;
    static int16* tile;
    static int32 tileid;

    /* Initialized data */
    static int firstCall = 1;
    static int32 idims[3] = {1, 220, 220};
    static int32 istart[3] = {0, 0, 0};
    static int32 itilem = -999;

    /* Check if file is open */
    if (firstCall) {
        firstCall = 0;

        tile = (int16*) allocateMemory(48400 * sizeof (int16), "interp_dem_height - tile");

        /* Open DEM file and fill structure */
        if (initialize_dem(demfile, &tileid) != 0) {
            return 1;
        }

        /* Compute constants used for calculations */
        /* Tile size in latitude */
        tslat = 180.f / (dem.nrows - 1);
    }

    /* Compute tile # */
    irow = (*xlat - 0.0001f + 90.f) / tslat + .5f;
    icol = (*xlon - 0.0001f + 180.f) / 360.f * dem.ncols[irow];
    itile = dem.nstart[irow] + icol;
    iflag = dem.iflag[itile];

    /* Set nominal values (sea level) */
    *height = 0.f;

    /* If sea level, return input values */
    if (iflag == 0) {
        return 0;
    }

    /* If flat tile */
    if (iflag == 1) {
        *height = (float) dem.ihmax[itile];

        /* Else terrain tile */
    } else {

        /* Check for tile already in memory */
        if (itile != itilem) {

            /* Read tile array from file */
            istart[0] = dem.irecno[itile];
            //            printf("Tile %d record %d\n",itile,istart[0]);
            itilem = itile;
            if (SDreaddata(tileid, istart, NULL, idims, tile) != 0) {
                printf("Error reading tile %d record %d\n",
                        itile, dem.irecno[itile]);
                return 1;
            }

            /* Compute longitude width of tile */
            tslon = 360.f / dem.ncols[irow];

            /* Compute coordinates of lower-left corner of nominal tile */
            tilat = (irow - .5f) * tslat - 90.f;
            tilon = icol * tslon - 180.f;
        }

        /* Determine location of input lon/lat on tile grid */
        /*  (recall that points are bin centers) */
        x0 = (*xlon - tilon) / tslon * dem.itsize + dem.ioff + .5f;
        y0 = (*xlat - tilat) / tslat * dem.itsize + dem.ioff + .5f;

        /* Get interpolated height at input lon/lat */
        *height = bilin(&x0, &y0, &dem.isize, tile, &status);
    }

    return 0;
}
