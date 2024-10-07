#include "viirs_sim_sdr.h"

int get_cal_rvs(char *file, int ham, int line, int ipix, rvs)
/*-----------------------------------------------------------------------------
    Routine:   get_cal_rvs.c

    Description:  provide the rvs value for all bands and the 2 gain states 
      for a particular ham side and unaggregated pixel along scan
      
    Returns type: int - 0 if good, else I/O error

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        char *    file          I    calibration file name
        int       ham           I    the ham side for this scan
        int       line          I    detector number to access
        int       ipix          I    unaggregated pixel number, can be negative 
                                       for a along scan margin
        float *   rvs           O    RVS for the pixel, detector and ham side,
                                     for all bands 

    Modification history:

    W. Robinson, SAIC  16 Sep 2010  Original development

----------------------------------------------------------------------------*/ {
    /*
     *  read in all the info from the hdf file
     */
    if (first_time) {
        /*
         *  do mallocs for the space to read in the coefficients
         */
        rvs_a0

        rvs_a1

        rvs_a2
        /*
         *  open the hdf 5 file with the gain info and get the data
         */
        open

        read rvs_a0

        read rvs_a1

        read rvs_a2
        /*
         *  finish up with file
         */
        close
    }
    /*
     *  get the AOI from the pixel number
     */
    vir_xf_scan(ipix, VIR_SCAN_UASMP, VIR_SCAN_AOI, &aoi);
    /*
     *  transfer the correct values for specific ham, detector, and AOI
     */
    for (ibnd = 0; ibnd <#bands; ibnd++ )
        {
            *(rvs + ibnd) = *(a0 + INDEX ?) + *(a1 + INDEX ?) * aoi +
                    *(a2 + INDEX ?) * aoi * aoi;
        }
    return 0;
}
