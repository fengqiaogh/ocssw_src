#include "viirs_sim_sdr.h"

int get_cal_gain(char *file, unsigned char ham, int det, int *nf_bands,
        int *nf_gains, c0_sub, c1_sub, c2_sub, dark, lt_max_hg)
/*-----------------------------------------------------------------------------
    Routine:   get_cal_gain.c

    Description:  provide the gain coefficients, the dark count, and the 
      maximum Lt value in high gain mode for a specific ham side and line in 
      the scan.  Gain coefficients and dark count are for band, gain 
      setting and lt_max_hg is for the band
      
    Returns type: int - 0 if good, else I/O error

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    file          I    calibration file name
        unsigned char   ham     I    the ham side for this scan
        int       idet          I    detector number to access
        int *     nf_bands      O    # bands in gain file
        int *     nf_gains      O    # gains in gain file (should be 2 max)
        float *   c0_sub        O    0 order coefficients for gain, low gain
                                     followed by high gain values
        float *   c1_sub        O    1st order coefficients for gain, low gain
                                     followed by high gain values
        float *   c2_sub        O    2nd order coefficients for gain, low gain
                                     followed by high gain values
        float *   dark_sub      O    dark count for bands low gain
                                     followed by high gain values
        float *   lt_max_hg_sub O    Lt at max high gain state for bands

    Note that the c0_sub, c1_sub, c2_sub, dark_sub, lt_max_hg_sub
      are pointers passed in, the 
      routine will alocate the space needed for these

    Modification history:

    W. Robinson, SAIC  16 Sep 2010  Original development

----------------------------------------------------------------------------*/ {
    int nf_ham, nf_det;
    static float *c0, *c1, *c2, *dark, *lt_max_hg;
    /*
     *  read in all the info from the hdf file
     */
    if (first_time) {
        if (setup_gain1(nf_bands, nf_gains, &nf_ham, &nf_det,
                c0, c1, c2, dark, lt_max_hg) != 0) return 1;
        /*
         if( setup_gain2( file, nf_bands, nf_gains, &nf_ham, &nf_det,
           c0, c1, c2, dark, lt_max_hg ) != 0 )
         */
        /*
         *  
         }
        /*
         *  do mallocs for the space to read in the coefficients
         */
        c0

        c1

        c2

        dark_lcl

        lt_max_hg_lcl

        /*
         *  finish up with file
         */
        close
    }
    /*
     *  transfer the correct values for specific ham and detector
     */
    for (ibnd = 0; ibnd <#bands; ibnd++ )
        {
            lt_max_hg[ibnd] = lt_max_hg_lcl[ibnd];

            for (igain = 0; igain < 2; igain++) {
                *(c0_sub + ibnd + nbnd * igain) = *(c0 + INDEXING ? ? ?)

                        about same

                                   for c1, c2, dark
                }
        }
}
int setup_gain1(int *nf_bands, int *nf_gains, int *nf_ham, int *nf_det,
        float *c0, float *c1, float *c2, float *dark, float *lt_max_hg) {
    /*
     *  set the sizes
     *  allocate the storage
      /*
     *  open the hdf 5 file with the gain info and get the data
     */
    open

    read c0

    read c1

    read c2

    read dark_lcl

    read lt_max_hg_lcl

            * move all locally stored coeffs to storage
            * finished
            */
    return 0;
    /*
      for the gain2 which will read the hdf 5 file:
     */
}
