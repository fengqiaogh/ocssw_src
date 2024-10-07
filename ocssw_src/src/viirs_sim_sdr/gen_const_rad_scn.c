#include "viirs_sim_sdr.h"

int gen_const_rad_scn(sdr_info_struc *sdr_info, int st_bnd,
        in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Program:   gen_const_rad_scn.c

    Description:  Make a constant value for the Lt lines in a scan
      possible elaboration for future

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        sdr_info_struc * sdr_info I  general sdr information
        int       st_bnd        I    start band to fill
        in_rec_struc * in_rec   I/O  input record controls

    Modification history:

    W. Robinson, SAIC  12 Dec 2008  Original development
    W. Robinson, SAIC  12 Mar 2010  adapt to scan processing
    W. Robinson, SAIC  7 Jul 2010   control the start band for filling

----------------------------------------------------------------------------*/ {
    int ibnd, ismp, loc, idet;
    /* Ltyp for Lt(410, 445, 488,555, 672, 751, 865, 1240, 1378, 1610, 2250 */
    /* low radiance for al M bands from Y4413_Refl_calib_SD_dual_gain.doc */
    /* (W m^-2 s^-1 sr^-1  */
    /* Also, include the Ttyp for the remaining BBT bands  */
    float l_typ[] = {44.9, 40., 32., 21., 10., 9.6, 6.4, 5.4, 6.0, 7.3, 0.12,
        270., 300., 270., 300., 300.};

    static int entry = 0; /* just set so that we don't re-make the constants 
                        every time  */
    /*
     *  loop through and make the values, probably Ltyp
     */
    if (entry == 0) {
        printf("%s-temp: initial set-up of the band radiances (I hope)\n", __FILE__);
        entry = 1;
        for (idet = 0; idet < in_rec->ndet_scan; idet++) {
            for (ibnd = st_bnd; ibnd < in_rec->nbnd; ibnd++) {
                for (ismp = 0; ismp < in_rec->npix; ismp++) {
                    loc = idet * in_rec->npix + ismp;
                    *(in_rec->bnd_lt[ibnd] + loc) = l_typ[ibnd];
                    *(in_rec->bnd_q[ibnd] + loc) = 0;
                }
            }
        }
    }
    return 0;
}
