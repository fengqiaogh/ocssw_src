#include "viirs_sim_sdr.h"

int gen_sdr_fname(int isdr, char *path, sdr_info_struc *sdr_info,
        int fname_opt, char *fname)
/*-----------------------------------------------------------------------------
    Program:   gen_sdr_fname

    Description:  generate the standard VIIRS file name from the sdr 
      information

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       isdr          I    sdr index: 0 is for geo, 1 - 13 is for 
                                     the band data
        char *    path          I    location to place the file
        sdr_info_struc * sdr_info I  general sdr information
        int       fname_opt     I    file name create option, 0 - std VIIRS,
                                     1 - omit create time
        char *    fname        I/O   final file name

    Modification history:

    W. Robinson, SAIC  3 Feb, 2009  Original development

----------------------------------------------------------------------------*/ {
    int mode = 0;
    char prefix[350], cre_t_str[13];
    /*
     *  make a simple name using the old method
     */
    if (mode == 1) {
        if (isdr == 0)
            sprintf(fname, "V%s_GEO.h5", sdr_info->ofile_base);
        else
            sprintf(fname, "V%s_SDRM%d.h5", sdr_info->ofile_base,
                isdr);
    } else {
        /*
         *  create the name based on the standard viirs naming convention
         */
        if (isdr == 0) {
            sprintf(prefix, "%s/GMTCO_npp_", path);
        } else {
            sprintf(prefix, "%s/SVM%2.2d_npp_", path, isdr);
        }
        /*
         *  use all of create time to microsecs but remove the '.'
         */
        if (fname_opt == 0) {
            strncpy(cre_t_str, sdr_info->cre_time, 6);
            strncpy(cre_t_str + 6, sdr_info->cre_time + 7, 6);
            sprintf(fname, "%sd%s_t%6.6s0_e%6.6s0_b00001_c%s%12.12s_%s_%s.h5",
                    prefix, sdr_info->st_date, sdr_info->st_time, sdr_info->en_time,
                    sdr_info->cre_date, cre_t_str, sdr_info->origin,
                    sdr_info->domain);
        } else {
            sprintf(fname, "%sd%s_t%6.6s0_e%6.6s0_b00001_%s_%s.h5",
                    prefix, sdr_info->st_date, sdr_info->st_time, sdr_info->en_time,
                    sdr_info->origin, sdr_info->domain);
        }
    }
    printf("%s, %d: created file name: %s\n", __FILE__, __LINE__, fname);
    return 0;
}
