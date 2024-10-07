#ifndef _SWL1_STRUC_H
#define _SWL1_STRUC_H

#include "swl0_parms.h"
#include "swl0_types.h"

typedef struct swl1_struct {
    INT16 type; /* 0=LAC,1=LUN, ..., 15=GAC                */
    INT16 npix; /* pixels per scan (GAC: 248, LAC: 1285)   */

    /* The fields that begin here are mimic the L1A hdf file spec   */
    INT32 msec; /* frame time millisecs of day             */
    BYTE eng_qual[4]; /* Engineering data out-of-range flags     */
    BYTE s_flags[4]; /* scan-line quality flags                 */
    INT16 s_satp[8]; /* number of saturated pixels per band     */
    INT16 s_zerop[8]; /* number of zero pixels per band          */
    FLOAT32 slat; /* starting latitude of scan               */
    FLOAT32 slon; /* starting longitude of scan              */
    FLOAT32 clat; /* central latitude of scan                */
    FLOAT32 clon; /* central longitude of scan               */
    FLOAT32 elat; /* ending latitude of scan                 */
    FLOAT32 elon; /* ending longitude of scan                */
    FLOAT32 csol_z; /* central solar zenith angle of scan      */
    FLOAT32 tilt; /* tilt angle (deg)                        */
    INT16 scid[2]; /* raw S/C ID field                        */
    INT16 ttag[4]; /* raw frame timetag                       */
    BYTE soh[775]; /* raw SOH block                           */
    INT16 inst[44]; /* raw instrument or ancillary tlm         */
    INT16 data[1285][8]; /* raw image data (counts)               */
    INT16 startpix[8]; /* raw start pixel per band                */
    INT16 stoppix[8]; /* raw stop pixel per band                 */
    INT16 darkpix[8]; /* raw dark restore pixel per band         */
    INT16 gain[8]; /* gain setting per band                   */
    INT16 tdi[8]; /* time delay and integration set per band */
    FLOAT32 inst_ana[40]; /* converted instrument analog telemetry */
    BYTE inst_dis[32]; /* converted instrument discrete telem   */
    FLOAT32 sc_ana[40]; /* converted spacecraft analog telemetry */
    BYTE sc_dis[40]; /* converted spacecraft discrete telem   */
    INT16 scan_temp[8]; /* focal-plane temp per band (counts)    */
    INT16 side; /* mirror side                             */
    FLOAT32 orb_vec[3]; /* orbit position vector at scan time      */
    FLOAT32 l_vert[3]; /* local vertical vector in ECEF frame     */
    FLOAT32 sun_ref[3]; /* solar reference vector in ECEF frame    */
    FLOAT32 att_ang[3]; /* attitude angles (yaw,roll,pitch) (deg)  */
    FLOAT32 sen_mat[3][3]; /* ECEF-to-sensor-frame matrix           */
    FLOAT32 scan_ell[6]; /* scan-track ellipse coefficients         */
    INT32 nflag[8]; /* Navigation flags                        */

    /* These are bonus fields                                       */
    FLOAT32 lon[1285]; /* scan longitude (deg)                    */
    FLOAT32 lat[1285]; /* scan latitude (deg)                     */
    FLOAT32 solz[1285]; /* scan solar zenith angle (deg)           */
    FLOAT32 sola[1285]; /* scan solar azimuth angle (deg)          */
    FLOAT32 senz[1285]; /* scan sensor zenith angle (deg)          */
    FLOAT32 sena[1285]; /* scan sensor azimuth angle (deg)         */


} swl1rec;

typedef struct l1met_struct {
    INT32 gain1_satpix[NBANDS]; /*                            */
    INT32 gain2_satpix[NBANDS]; /*                            */
    INT32 gain1_nonsatpix[NBANDS]; /*                            */
    INT32 gain2_nonsatpix[NBANDS]; /*                            */
    INT32 zeropix[NBANDS]; /*                            */
    FLOAT32 gain1_mean_rad[NBANDS]; /*                            */
    FLOAT32 gain2_mean_rad[NBANDS]; /*                            */

} l1met;


#endif



