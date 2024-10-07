/* 
 * File:   lonlat2pixline.h
 * Author: dshea
 *
 * Created on March 8, 2010, 8:53 AM
 */

#ifndef _LONLAT2PIXLINE_H
#define _LONLAT2PIXLINE_H

#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Structure used to control lonlat2pixline and return the results.
 *
 * @see lonlat2pixline()
 */
typedef struct lonlat2pixline_struct {
    /* --------- input variables --------- */
    int want_pixbox; /**< use SWlon, SWlat as center and xbox, ybox as
                              number of pixels on each side */
    int pix_srch; /**< search for a single pixel (use SW lon & lat)*/
    int32_t xbox; /**< min number of x pixels */
    int32_t ybox; /**< min number of y lines */
    float SWlon; /**< SW longitude corner of box */
    float SWlat; /**< SW latitude corner of box */
    float NElon; /**< NE longitude corner of box */
    float NElat; /**< NE latitude corner of box */
    char input_filename[FILENAME_MAX]; /**< input filename */
    char geo_filename[FILENAME_MAX]; /**< MODIS GEO file */
    int32_t resolution; /**< output resolution for MODIS (-1,1000,500,250) */

    /* --------- output variables --------- */
    int box_failed; /**< 1 if the box was not fully extracted */
    float pixLon; /**< actual longitude of the pixel found */
    float pixLat; /**< actual latitude of the pixel found */
    int32_t spixl; /**< start pixel */
    int32_t epixl; /**< end pixel */
    int32_t sline; /**< start line */
    int32_t eline; /**< end line */
} lonlat2pixline_t;



/**
 * Function that returns the pixel, line of a lon, lat or lon, lat box for
 * the given L1, L2 or MODIS GEO file.
 *
 * Note that if a MODIS GEO file is given with a resolution of 500 or
 * 250 the pixl and line could be off by 1 or 3 respectively.
 *
 * box usage:<BR>
 *      params->SWlon          - SW longitude corner of box
 *      params->SWlat          - SW latitude corner of box
 *      params->NElon          - NE longitude corner of box
 *      params->NElat          - NE latitude corner of box
 *      params->want_pixbox    - 0<BR>
 *      params->pix_srch       - 0<BR>
 *      params->input_filename - L1, L2 or GEO file<BR>
 *      params->geo_filename   - MODIS GEO file<BR>
 *      params->resolution     - output resolution if MODIS (-1,1000,500,250) 
 * 
 * pix usage:<BR>
 *      params->SWlon          - longitude value of location<BR>
 *      params->SWlat          - latitude value of location<BR>
 *      params->want_pixbox    - 0 for single point<BR>
 *                               1 for box defined by xbox, ybox (pixels)<BR>
 *      params->pix_srch       - 1<BR>
 *      params->xbox           - num pixels on either side of location (want_pixbox=1)<BR>
 *      params->ybox           - num scan lines on either side of location<BR>
 *      params->input_filename - L1 or L2 file<BR>
 *      params->geo_filename   - MODIS GEO file<BR>
 *      params->resolution     - output resolution if MODIS (-1,1000,500,250) 
 *
 * other inputs:<BR>
 *      want_verbose(gloabal)  - 0 for no output to stdout<BR>
 *                             - 1 for status messages<BR>
 *
 * outputs:<BR>
 *      params->box_failed - 1 if the box was not fully extracted, else 0<BR>
 *      params->pixLon     - actual longitude of the pixel found<BR>
 *      params->pixLat     - actual latitude of the pixel found<BR>
 *      params->spix       - start pixel<BR>
 *      params->epix       - end pixel<BR>
 *      params->sline      - start line<BR>
 *      params->eline      - end line<BR>
 *
 * @param params a structure that controls the function.
 * @return 0 = all OK<BR>
 * 1=error<BR>
 * 110=full box not extracted<BR>
 * 120=whole file extracted
 */
int lonlat2pixline(lonlat2pixline_t *params);

int lonlat2pixline1(char *input_filename, char *geo_filename,
        int32_t resolution, float SWlon, float SWlat, float NElon, float NElat,
        int32_t *spixl, int32_t *epixl, int32_t *sline, int32_t *eline);

int lonlat2pixline2(char *input_filename, char *geo_filename,
        int32_t resolution, float lon, float lat, int32_t dx, int32_t dy,
        int32_t *spixl, int32_t *epixl, int32_t *sline, int32_t * eline);

int lonlat2pixline3(char *input_filename, char *geo_filename,
        int32_t resolution, float lon, float lat,
        int32_t *pixl, int32_t *line);

#ifdef __cplusplus
}
#endif

#endif /* _LONLAT2PIXLINE_H */

