#include "viirs_sim_sdr.h"
#include "l12_parms.h"

int bnd_ix_2_sen_info(char *pname, void *pval)
/*-----------------------------------------------------------------------------
    Program:   bnd_ix_2_sen_info

    Description: augment rdsensorinfo for the F0 and lambda lists to add
     VIIRS M9 band info into the list

    Arguments:
      Type      Name         I/O   Description
      ----      ----         ---   -----------
      char *    pname         I    parameter name to get data for
      void *    pval         I/O   pointer to allocated storage location to 
                                   hold the data

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development

----------------------------------------------------------------------------*/ {
    static int32_t m_bnd_ix[] = {0, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11,
        12, 13, 14};
    static int n_gap = 1;
    static int32_t lam_cirrus = 1378;
    static float f0_cirrus = 36.463;
    int iprm, ix_in, ix_out, ibnd;
    void *lcl_prm;
    int *lcl_int, *out_int;
    float *lcl_flt, *out_flt;
    /*
     *  depending on the parm name, set an integer value
     */
    if (strcmp(pname, "Lambda") == 0)
        iprm = 0;
    else if (strcmp(pname, "Fobar") == 0)
        iprm = 1;
    else {
        printf("%s, %d: parameter %s is outside values this routine can handle\n",
                __FILE__, __LINE__, pname);
        return 1;
    }
    /*
     *  get the parameter list with rdsensorinfo
     */
    if (rdsensorinfo(VIIRSN, 0, pname, &lcl_prm) < 0) {
        printf("%s, %d: failure to read sensor parameter %s\n",
                __FILE__, __LINE__, pname);
    }

    switch (iprm) {
    case 0: lcl_int = (int *) lcl_prm;
        out_int = (int *) pval;
        break;
    case 1: lcl_flt = (float *) lcl_prm;
        out_flt = (float *) pval;
        break;
    }
    /*
     *  insert the cirrus value, moving following values up 1
     */
    ix_in = MAX_BND - 1 - n_gap;
    for (ibnd = 0, ix_out = MAX_BND - 1; ibnd < MAX_BND; ibnd++, ix_out--) {
        if (m_bnd_ix[ix_out] >= 0) {
            switch (iprm) {
            case 0: *(out_int + ix_out) = *(lcl_int + ix_in);
                break;
            case 1: *(out_flt + ix_out) = *(lcl_flt + ix_in);
                break;
            }
            ix_in--;
        } else {
            switch (iprm) {
            case 0: *(out_int + ix_out) = lam_cirrus;
                break;
            case 1: *(out_flt + ix_out) = f0_cirrus;
                break;
            }
        }
    }
    return 0;
}
