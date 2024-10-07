#include "proto.h"
#include "mfhdf.h"

int get_attributes(int32 sd_id, char *prefix, int *Nlines, int *Npixels, int *startpix, int *subsampling, int *year, int *doy, int *orbno, char *title) {
    int32 int32tmp;

    /* get product name */
    if (read_attr(sd_id, "Product Name", prefix) == -1) return -1;

    /* get number of scans */
    if (read_attr(sd_id, "Number of Scan Lines", &int32tmp) == -1) return -1;
    *Nlines = int32tmp;

    /* get number of pixels per scan : should be 1285 for LAC and 248 for GAC */
    if (read_attr(sd_id, "Pixels per Scan Line", &int32tmp) == -1) return -1;
    *Npixels = int32tmp;

    /* get LAC number of 1st archived pixel : should be 1 for LAC and 147 for GAC */
    if (read_attr(sd_id, "LAC Pixel Start Number", &int32tmp) == -1) return -1;
    *startpix = int32tmp;

    /* get LAC pixel subsampling : should be 1 for LAC and 4 for GAC */
    if (read_attr(sd_id, "LAC Pixel Subsampling", &int32tmp) == -1) return -1;
    *subsampling = int32tmp;

    /* get orbit number */
    if (read_attr(sd_id, "Orbit Number", (int16 *) orbno) == -1) return -1;

    /* get product title */
    if (read_attr(sd_id, "Title", title) == -1) return -1;

    return 0;
}
