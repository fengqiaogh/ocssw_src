#include "viirs_sim_sdr.h"
#include "l12_parms.h"

int rd_sim_init(ctl_struc *ctl, sdr_info_struc *sdr_info,
        in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Routine:  rd_sim_init

    Description:  prepare the simulator input files for use
      It will prepare the input files, including the geolocation
      file and the level-2 file containing the TOA radiances

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    program argument controls
        sdr_info_struc * sdr_info I/O  general SDR information
        in_rec_struc * in_rec  I/O   input file information

    Modification history:

    W. Robinson, SAIC  20 Nov 2008  Original development
    W. Robinson, SAIC  16 Nov 2009  if reflectance file is used, set up to 
                                    read extra atmospheric parameters from
                                    the L2 file

----------------------------------------------------------------------------*/ {
    char parm_list[2000];
    char *atm_prm_root[] = {"Lr_", "tg_sol_", "tg_sen_", "t_sol_", "t_sen_",
        "t_h2o_", "t_o2_"};
    int iprm, ibnd, nchr, sptr;
    /*
     *  set up the geolocation file
     */
    if (rd_geo_init(ctl, sdr_info, in_rec)
            != 0) return 1;
    /*
     *  make space for the half angle mirror side array with the # of scans
     */
    if ((sdr_info->ham_side = (unsigned char *)
            malloc(in_rec->nscan * sizeof ( unsigned char))) == NULL) {
        printf("%s, %d: Error, failed to allocate ham_side storage\n",
                __FILE__, __LINE__);
        return 1;
    }
    /*
     *  the # bands to use is determined by the make_m
     */
    in_rec->nbnd = (ctl->make_m == 0) ? N_VNIR_BND : MAX_BND;
    /*
     *  Also, set up the input L2
     */
    if (ctl->l2_use == 1) {
        /*
         *  set the parm list to be the VIIRS bands and if rhos_use is set,
         *  also get required extra atmospheric parameters
         */
        sptr = 0;
        for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
            nchr = sprintf((parm_list + sptr), "vLt_%d:",
                    *(in_rec->lam_band + ibnd));
            sptr += nchr;
        }
        if (ctl->rhos_use == 1) {
            for (iprm = 0; iprm < 7; iprm++) {
                for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
                    nchr = sprintf((parm_list + sptr), "%s%d:",
                            *(atm_prm_root + iprm), *(in_rec->lam_band + ibnd));
                    sptr += nchr;
                }
            }
        }
        *(parm_list + sptr - 1) = 0;

        printf("%s: initializing input from L2 file: %s\n", __FILE__,
                ctl->l2_file);
        if (openL2(ctl->l2_file, parm_list, &in_rec->l2_str) != 0) {
            printf("%s, %d: Failure to open L2 TOA dataset: %s\n", __FILE__,
                    __LINE__, ctl->l2_file);
            return 1;
        }
    } else {
        printf("%s, %d: Dummy L2 data initialization (none currently)\n",
                __FILE__, __LINE__);
        /* call blnk_l2_init if ever needed or do in line */
    }
    return 0;
    /*
     *  and exit
     */
}
