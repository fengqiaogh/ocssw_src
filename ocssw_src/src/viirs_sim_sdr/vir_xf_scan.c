#include "viirs_sim_sdr.h"

int vir_xf_scan(float scn_in, int scn_in_typ, int scn_out_typ, float *scn_out)
/*-----------------------------------------------------------------------------
    Routine:   vir_xf_scan

    Description:  utility to transform VIIRS scan values between
      pixel, scan angle, and AOI on HAM
      
    Returns type: int - 0 if good, else I/O error

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        float     scn_in        I    incoming scan quantity
        int       scn_in_typ    I    designation of the type of input scan
                                     quantity:
        int       scn_out_typ   I    designation of the type of output scan
                                     quantity
        float *   scn_out       O    output scan quantity

    The types of scan quantity can be
      VIR_SCAN_UASMP  unaggregated samples
      VIR_SCAN_AGSMP  aggregated samples
      VIR_SCAN_AOI    AOI on HAM mirror, degrees
      VIR_SCAN_ANG    scan angle, degrees

    Implemented transforms:
      VIR_SCAN_UASMP -> VIR_SCAN_AOI

    Modification history:

    W. Robinson, SAIC  16 Sep 2010  Original development

----------------------------------------------------------------------------*/ {
    float sind = 0.0003104; /*  # radians / sample  */
    float rad2deg = 180. / M_PI;
    float theta_oop = 28.6 / rad2deg;
    int iret, idir, in_chain;
    float ua_smp, scan, theta_ham, aoi;
    /*
     *  The types of scan quantity can be
     *  mnemonic        value     description
     *  VIR_SCAN_AGSMP      0     aggregated sample #
     *  VIR_SCAN_UASMP      1     unaggregated sample #
     *  VIR_SCAN_ANG        2     scan angle, degrees
     *  VIR_SCAN_AOI        3     AOI on HAM mirror, degrees
     *
     *  The transforms will fall in 1 of 2 main directions:
     *  1: VIR_SCAN_AGSMP -> VIR_SCAN_UASMP -> VIR_SCAN_ANG -> VIR_SCAN_AOI
     *  2: VIR_SCAN_AOI ->  VIR_SCAN_ANG -> VIR_SCAN_UASMP -> VIR_SCAN_AGSMP
     *  as the values of the mnemonic in (1) rise from 0 -> 3 the direction
     *  to choose is determined simply by seeing if scn_out_typ - scn_in_typ
     *  is > 0 (direction (1) ) or < 0 (direction (2) )
     */
    iret = 0;
    if (scn_out_typ == scn_in_typ) {
        *scn_out = scn_in;
    } else {
        idir = (scn_out_typ > scn_in_typ) ? 1 : 2;
        if (idir == 1) {
            /*
             *  set internal variable for the type and progress through transform 
             *  chain, doing what is needed
             */
            switch (scn_in_typ) {
            case VIR_SCAN_AGSMP:
                //ag_samp = scn_in;
                break;
            case VIR_SCAN_UASMP:
                ua_smp = scn_in;
                break;
            case VIR_SCAN_ANG:
                scan = scn_in;
                break;
            }

            in_chain = 0;
            if (scn_in_typ == VIR_SCAN_AGSMP) {
                /*
                 *  aggregated to unaggregated transform
                 */
                printf("%s, %d: Error, aggregated to unaggregated scan transform not supported yet\n", __FILE__, __LINE__);
                iret = 1;
                return iret;
            }
            if (scn_out_typ == VIR_SCAN_UASMP) {
                in_chain = 0;
                *scn_out = ua_smp;
            }
            if ((scn_in_typ == VIR_SCAN_UASMP) || (in_chain == 1)) {
                /*
                 *  unaggregated sample to scan angle
                 */
                scan = (ua_smp - 3151.5) * sind * rad2deg;
                in_chain = 1;
            }
            if (scn_out_typ == VIR_SCAN_ANG) {
                in_chain = 0;
                *scn_out = scan;
            }
            if (in_chain == 1) {
                /*
                 *  scan angle to AOI angle
                 */
                theta_ham = (scan - 46.) / 2.;
                aoi = acos(cos(theta_ham / rad2deg) * cos(theta_oop));
                *scn_out = aoi * rad2deg;
            }
        } else {
            /*
             *  nothing implemented yet
             */
            printf("%s, %d: Error, this scan transform is not supported yet\n",
                    __FILE__, __LINE__);
            printf("\t Path AOI on HAM -> scan angle -> unagg sample -> agg\n");
        }
    }
    return iret;
}
