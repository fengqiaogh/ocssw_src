/*

Norman Kuring		28-Jul-1999
I modified map_octs_to_box.c to produce this code.

Norman Kuring		15-Feb-2000
I modified swmapdepth.c to produce this code.

Norman Kuring		04-Oct-2002
Change in input file format.  Depth is now represented as a floating point.

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mfhdf.h>
#include <genutils.h>

#define USAGE "\
Usage: %s filename\n\
where: filename is the full path to a SeaWiFS level-2 file that contains\n\
                computed depth data.\n\
This program maps SeaWiFS level-2 depth data to bins of a Plate Carree\n\
projection that extends from 35 North to 35 South and from 180 West\n\
to 180 East.  The results are written to standard output as a series\n\
of ASCII text records.  Each output record contains a bin number and\n\
pixel value separated by whitespace.  Each bin represents one 0.01 by\n\
0.01 degree subregion of the aforementioned area.  Bins are numbered\n\
from 0 to 251,999,999.  Depths are expressed in meters.\n"

#define FILENAME argv[1]
#define LAC_PIXEL_INSET 200
#define SCALE  100

#define READ_GLBL_ATTR(nam,ptr) {         \
  if(SDreadattr(sd_id,SDfindattr(sd_id,(nam)),(VOIDP)(ptr))){      \
    fprintf(stderr,           \
    "-E- %s line %d: Could not get global attribute, %s, from file, %s.\n", \
    __FILE__,__LINE__,(nam),FILENAME);         \
    return(EXIT_FAILURE);          \
  }             \
}

#define MALLOC(ptr,typ,num) {      \
  (ptr) = (typ *)malloc((num) * sizeof(typ));    \
  if((ptr) == NULL){       \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n", \
    __FILE__,__LINE__);       \
    exit(EXIT_FAILURE);       \
  }         \
}

#define READ_SDS(nam,ptr,s0,s1,s2,e0,e1,e2) {   \
  int32 start[3],edge[3];     \
  edge[0]=(e0); edge[1]=(e1); edge[2]=(e2);   \
  start[0]=(s0); start[1]=(s1); start[2]=(s2);   \
  if(SDreaddata(SDselect(sd_id, SDnametoindex(sd_id, (nam))), \
  start, NULL, edge, (VOIDP)(ptr)) == FAIL){   \
    fprintf(stderr,"-E- %s line %d: Could not read SDS, %s.\n", \
    __FILE__,__LINE__,(nam));     \
    exit(EXIT_FAILURE);      \
  }        \
}

void get_scan_coords(int32, int32, int32, float32 **, float32 **);

static int32 sd_id;

int main(int argc, char *argv[]) {

    int32 npix, nscans, num_cntl_pts, lacstart, lacsub;
    int32 pix_0, pix_n;
    float32 *lat, *lon;
    float32 *data;
    int32 nflag[8];
    int s, bin, x, y;

    if (argc < 2) {
        fprintf(stderr, USAGE, argv[0]);
        return (EXIT_FAILURE);
    }

    /* Open the HDF file */
    sd_id = SDstart(FILENAME, DFACC_RDONLY);
    if (sd_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart failed for file, %s.\n",
                __FILE__, __LINE__, FILENAME);
        return (EXIT_FAILURE);
    }

    /* Read the file dimensions. */
    READ_GLBL_ATTR("Pixels per Scan Line", &npix);
    READ_GLBL_ATTR("Number of Scan Lines", &nscans);
    READ_GLBL_ATTR("Number of Pixel Control Points", &num_cntl_pts);
    READ_GLBL_ATTR("LAC Pixel Start Number", &lacstart);
    READ_GLBL_ATTR("LAC Pixel Subsampling", &lacsub);
    if (lacsub == 1) {
        if (lacstart - 1 + npix < LAC_PIXEL_INSET
                || lacstart - 1 > 1284 - LAC_PIXEL_INSET) {
            fprintf(stderr,
                    "-E- %s line %d: Input pixels are too close to the edge of the swath.\n",
                    __FILE__, __LINE__);
            return (EXIT_FAILURE);
        }
        pix_0 = LAC_PIXEL_INSET - lacstart + 1;
        pix_n = 1285 - lacstart - LAC_PIXEL_INSET;
    } else if (lacsub == 4) {
        if (lacstart - 1 + 4 * npix < LAC_PIXEL_INSET
                || lacstart - 1 > 1284 - LAC_PIXEL_INSET) {
            fprintf(stderr,
                    "-E- %s line %d: Input pixels are too close to the edge of the swath.\n",
                    __FILE__, __LINE__);
            return (EXIT_FAILURE);
        }
        pix_0 = (LAC_PIXEL_INSET - lacstart + 1) / 4 + 0.5;
        pix_n = 247 - (int) ((lacstart - 1 + LAC_PIXEL_INSET) / 4 + 0.5);
    } else {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr,
                "Unexpected number of pixels per scanline (%d) in file, %s.\n",
                npix, FILENAME);
        return (EXIT_FAILURE);
    }
    if (pix_0 < 0) pix_0 = 0;
    if (pix_n >= npix) pix_n = npix - 1;
    if (nscans <= 0) {
        fprintf(stderr,
                "-E- %s line %d: Bad scene dimension: nscans=%d in file, %s.\n",
                __FILE__, __LINE__, nscans, FILENAME);
        return (EXIT_FAILURE);
    }

    /* Allocate some memory for the per-scanline data. */
    MALLOC(data, float32, npix);

    /* For each scan line ... */
    for (s = 0; s < nscans; s++) {

        int32 p;

        /*
        Read the navigation flags and skip scans having invalid navigation.
        Also skip scans that were collected while the sensor tilt was changing.
         */
        READ_SDS("nflag", nflag, s, 0, 0, 1, 8, 1);
        if (nflag[0] == 1) continue; /* Invalid navigation */
        if (nflag[6] == 2) continue; /* Tilt is changing */

        /* Read the data and their coordinates. */
        READ_SDS("depth", data, s, 0, 0, 1, npix, 1);
        get_scan_coords(s, npix, num_cntl_pts, &lat, &lon);

        /* For each pixel excluding the limbs ... */
        for (p = pix_0; p <= pix_n; p++) {

            /* Skip pixels outside of the 35 North to 35 South range. */
            if (lat[p] > 35 || lat[p] <= -35) continue;

            /* Plate Carree projection */
            x = (int) ((lon[p] + 180) * SCALE);
            y = (int) ((35 - lat[p]) * SCALE);

            /* Compute the bin number for this pixel. */
            bin = y * 36000 + x;

            /* Print out the bin number and pixel value. */
            printf("%d %g\n", bin, data[p]);

        }
    }

    free(lat);
    free(lon);
    free(data);

    if (SDend(sd_id) == FAIL) {
        fprintf(stderr, "-W- %s line %d: SDend failed for file, %s.\n",
                __FILE__, __LINE__, FILENAME);
    }
    return (EXIT_SUCCESS);
}

void
get_scan_coords(
        int32 s, /* (in) scan line number */
        int32 npix, /* (in) number of pixels per scan */
        int32 nctl, /* (in) number of control points per scan */
        float32 **lat, /* (out) latitude  for each pixel in the scan */
        float32 **lon /* (out) longitude for each pixel in the scan */
        ) {
    static float32 *lati = NULL, *longi = NULL, *latitude = NULL, *longitude = NULL;
    static float32 *cntl_pts, *spl2nderiv;
    static int32 np, nc;

    if (latitude == NULL) {
        /* Execute this block of code only the 1st time this function is called. */
        np = npix;
        nc = nctl;
        MALLOC(latitude, float32, nc);
        MALLOC(longitude, float32, nc);
        if (nc != np) {
            /*
            If there is not a lat/lon pair for each pixel we must interpolate.
            Read the pixel control point indices, and then convert them to
            floats because that is what the spline and splint functions expect.
             */
            int32 p, *cntl_pt_cols;
            MALLOC(cntl_pt_cols, int32, nc);
            READ_SDS("cntl_pt_cols", cntl_pt_cols, 0, 0, 0, nc, 0, 0);
            MALLOC(cntl_pts, float32, nc);
            for (p = 0; p < nc; p++) {
                cntl_pts[p] = (float32) cntl_pt_cols[p];
            }
            free(cntl_pt_cols);
            MALLOC(spl2nderiv, float32, nc);
            MALLOC(lati, float32, np);
            MALLOC(longi, float32, np);
        } else {
            /* No interpolation needed. */
            lati = latitude;
            longi = longitude;
        }
    }

    /* Read the coordinates for the specified scan. */
    READ_SDS("latitude", latitude, s, 0, 0, 1, nc, 1);
    READ_SDS("longitude", longitude, s, 0, 0, 1, nc, 1);

    if (nc != np) { /* Interpolation required */
        int p;

        /* Remove any dateline discontinuity in the longitudes. */
        for (p = 1; p < nc; p++) {
            float delta;
            delta = longitude[p] - longitude[p - 1];
            if (delta < -180) longitude[p] += 360;
            else if (delta > 180) longitude[p] -= 360;
        }

        /*
        Interpolate to provide a latitude and longitude for each pixel
        in the scan.  I assume here that the pixel control points stored
        in the HDF file are one-based.
         */
        spline(cntl_pts, latitude, nc, 1e30, 1e30, spl2nderiv);
        for (p = 0; p < np; p++) {
            splint(cntl_pts, latitude, spl2nderiv, nc, p + 1, &lati[p]);
        }
        spline(cntl_pts, longitude, nc, 1e30, 1e30, spl2nderiv);
        for (p = 0; p < np; p++) {
            splint(cntl_pts, longitude, spl2nderiv, nc, p + 1, &longi[p]);
            /* Put the longitudes back in the [-180,180] range. */
            while (longi[p] >= 180) longi[p] -= 360;
            while (longi[p] < -180) longi[p] += 360;
        }
    }

    *lat = lati;
    *lon = longi;
}
