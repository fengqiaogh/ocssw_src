#include "viirs_sim_sdr.h"
#include <math.h>

int mod_artifact(ctl_struc *ctl, in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Program:   mod_artifact.c

    Description:  control all the artifact addition in this routine

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    input controls
        in_rec_struc * in_rec    I/O  controls for input record reading

    Modification history:

    W. Robinson, SAIC  29 Apr 2010  Original development

----------------------------------------------------------------------------*/ {
    int ibnd, ilin, ipix;
    float *val_ptr;
    /*
     *  Near Field Response, stray light artifact section
     */
    if (ctl->stray_opt == 1)
        if (viirs_straylt(ctl, in_rec, 0) != 0)
            return 1;
    /*
     *  Optical Cross Talk section
     */
    if (ctl->oxt_mode == 1)
        if (viirs_oxt(ctl, in_rec) != 0)
            return 1;
    /*
     *  Noise artifact is added here
     */
    if (ctl->noise_mode == 1)
        if (viirs_noise(ctl, in_rec, 0) != 0)
            return 1;
    /*
     *  Electronic crosstalk needs to be applied to the dark-corrected counts,
     *  so the total rads must be de-calibrated first to apply this artifact
     *  and then calibrated again.  Also, if the effects of different cal and
     *  de-cal or integerization errors of counts id needed, this branch can 
     *  be run through.  Possibly in the future, dn values could be returned
     */
    if (ctl->count_cal_opt > 0) {
        /*
         *  NOTE that earlier, the de-cal and cal gain and rvs files 
         *  would be set based on their existance in the inputs
         *
         *  call the de-calibration step 
         */
        if (viirs_decal(ctl, in_rec) != 0) return 1;

        /*
         *  The dn's are either the dark subtracted (dn) or non-subtracted (DN)
         *  based on input controls.  
         *  this step will go through and integerize the dn values
         *
         *  I take the best case - where the integer dn is nearest to the actual
         *  value
         */
        if (ctl->count_cal_opt == 2) {
            for (ibnd = 0; ibnd < in_rec->nbnd; ibnd++)
                for (ilin = 0; ilin < in_rec->ndet_scan; ilin++)
                    for (ipix = 0; ipix < in_rec->npix; ipix++) {
                        val_ptr = in_rec->dn[ibnd] + ipix + in_rec->npix * ilin;
                        *val_ptr = floorf(*val_ptr + 0.5);
                    }
        }
        /*
         *  The electronic crosstalk will be applied to the dn, gain information
         */
        if (ctl->ext_opt == 1)
            if (viirs_ext(ctl, in_rec) != 0)
                return 1;
        /*
         *  Lastly, the dn values will get calibrated again, not necessarily 
         *  with the same calibration values
         */
        if (viirs_cal(ctl, in_rec) != 0) return 1;
    }
    return 0;
}
