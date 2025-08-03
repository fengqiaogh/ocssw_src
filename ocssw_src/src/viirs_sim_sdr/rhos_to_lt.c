#include "viirs_sim_sdr.h"
#include "l12_parms.h"

/*#include "l12_proto.h"  - Again, if we get the make more like ocssw 
   and maybe also need to move formally to that area */

int rhos_to_lt(int rhos_opt, float **rhos, in_rec_struc *in_rec, int idet,
        sdr_info_struc *sdr_info)
/*-----------------------------------------------------------------------------
    Program:   rhos_to_lt

    Description:  convert the reflectance into Lt value and replace
      the Lt conditionally for a line in a scan of lines

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       rhos_opt      I    reflectance use option: 0 - replace 
                                     where Lt from ocean color is missing
                                     1 - replace everywhere
        float **  rhos          I    a scan's worth of reflectance
        in_rec_struc *  in_rec  I/O  input record structure with input values
                                     of atmospheric values and output of 
                                     final lt for simulated scan
        int       idet          I    line in scan to process
        sdr_info_struc * sdr_info I  general SDR information

    Modification history:

    W. Robinson, SAIC  18 Nov 2009  Original development
    W. Robinson, SAIC  12 Apr 2011  for a reflectance that is missing -999.
                or just < 0., set the bnd_q = 2

----------------------------------------------------------------------------*/ {
    int ipx, npix, ibnd, repl, loc;
    float lt[N_VNIR_BND], solz, mu0;
    float tg_sol, tg_sen, t_sol, t_sen, t_o2, lr, refl;
    /*  positions of parameter fields in the l2 record, see rd_sim_init */
    int loc_tg_sol = N_VNIR_BND * 2, loc_tg_sen = N_VNIR_BND * 3;
    int loc_t_sol = N_VNIR_BND * 4, loc_t_sen = N_VNIR_BND * 5;
    int loc_t_o2 = N_VNIR_BND * 7;
    int loc_lr = N_VNIR_BND;

    npix = in_rec->npix;
    /*
     *  get lt from rhos for all pixels and lines in the scan
     */
    for (ipx = 0; ipx < npix; ipx++) {
        loc = ipx + idet * npix;
        solz = *(in_rec->solz + loc);
        mu0 = cos(solz / OEL_RADEG);

        /*
         *  derive the Lt for each band
         */
        for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
            refl = *(*rhos + ipx + (idet + ibnd * in_rec->ndet_scan) * npix);
            if (refl < 0.)
                lt[ibnd] = -999.;
            else {
                lr = in_rec->l2_str.l2_data[loc_lr + ibnd][ipx];
                t_sen = in_rec->l2_str.l2_data[loc_t_sen + ibnd][ipx];
                t_sol = in_rec->l2_str.l2_data[loc_t_sol + ibnd][ipx];
                t_o2 = in_rec->l2_str.l2_data[loc_t_o2 + ibnd][ipx];
                tg_sen = in_rec->l2_str.l2_data[loc_tg_sen + ibnd][ipx];
                tg_sol = in_rec->l2_str.l2_data[loc_tg_sol + ibnd][ipx];

                lt[ibnd] = tg_sol * tg_sen *
                        (refl * in_rec->f0[ibnd] * mu0 * t_sol * t_sen * t_o2 / PI +
                        RAD_CGS_2_MKS * lr);
            }
        }
        /*
         *  examine all the ocean lt and decide on replacing it with the lt(rhos)
         *  Options are to replace if Lt < 0 or a bad value (rhos_opt=0)
         *  or replace al points (rhos_opt=1)
         */
        if (rhos_opt == 1)
            repl = 1;
        else {
            repl = 0;
            for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
                if ((*(in_rec->bnd_lt[ibnd] + loc) < 0.) ||
                        (*(in_rec->bnd_q[ibnd] + loc) == 2)) {
                    repl = 1;
                    break;
                }
            }
            /*
             *  now, assume all reflectance is good and fill
             */
            if (repl == 1)
                for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
                    if (lt[ibnd] < -0.5) {
                        *(in_rec->bnd_lt[ibnd] + loc) = 0.;
                        *(in_rec->bnd_q[ibnd] + loc) = 2;
                    } else {
                        *(in_rec->bnd_lt[ibnd] + loc) = lt[ibnd];
                        *(in_rec->bnd_q[ibnd] + loc) = 0;
                    }
                }
        }
    }
    return 0;
}
