/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/mops/navigation.h,v 4.11 1995/04/03 14:44:36 seawifsd Exp seawifsd $
$Log: navigation.h,v $
Revision 4.11  1995/04/03 14:44:36  seawifsd
changed references of 'long' to 'int' to support OSF/1.

Revision 4.10  1995/01/17 19:58:30  seawifsd
Jan. 17, 1994, V4.10

Revision 4.2  1995/01/17 14:51:28  seawifsd
added prototyping for FORTRAN cdata_().

Revision 4.1  1995/01/17 14:14:43  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:46:47  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:45  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.1  1994/10/04 16:41:26  frank
added prototyping for geonav_().

Revision 1.2  1994/05/10 18:50:18  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:36:12  seawifst
Initial revision


 */


#ifndef NAVIGATION_H
#define NAVIGATION_H

//#ifndef LAC_PIXEL_NUM
//#define LAC_PIXEL_NUM 1285
//#endif /* !LAC_PIXEL_NUM */

#ifdef __cplusplus
extern "C" {
#endif
void cdata_(void);

//extern void geonav_(float *orb_vec,float *sen_mat, float *scan_ell,
//                    float *sun_ref, int *nsta, int *ninc, int *npix,
//                    float ylat[LAC_PIXEL_NUM], float xlon[LAC_PIXEL_NUM],
//                    float solz[LAC_PIXEL_NUM], float sola[LAC_PIXEL_NUM],
//                    float senz[LAC_PIXEL_NUM], float sena[LAC_PIXEL_NUM]);

void geonav_(float *orb_vec, float *sen_mat, float *scan_ell,
        float *sun_ref, int32_t *nsta, int32_t *ninc, int32_t *npix,
        float ylat[], float xlon[],
        float solz[], float sola[],
        float senz[], float sena[]);

void geonav_lonlat_(float *orb_vec, float *sen_mat, float *scan_ell,
        float *sun_ref, int32_t *nsta, int32_t *ninc, int32_t *npix,
        float ylat[], float xlon[]);


#ifdef __cplusplus
}
#endif


#ifndef NAVBLOCK_STRUCT
#define NAVBLOCK_STRUCT

typedef struct navblockStruct {
    float *orb_vec;
    float *l_vert;
    float *sun_ref;
    float *att_ang;
    float *sen_mat;
    float *scan_ell;
    int *nflag;
} navblockType;
#endif /* NAVBLOCK_STRUCT */

#ifndef GEOLOC_STRUCT
#define GEOLOC_STRUCT

typedef struct geolocStruct {
    float *ylat; /* pixel geodetic latitudes		*/
    float *xlon; /* pixel geodetic longitude 		*/
    float *solz; /* pixel solar zenith angle		*/
    float *sola; /* pixel solar azimuth angle		*/
    float *senz; /* pixel sensor zenith angle		*/
    float *sena; /* pixel sensor azimuth angle		*/
} geolocType;
#endif /* GEOLOC_STRUCT */

#endif /* NAVIGATION_H */
