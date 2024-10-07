#include "l1czcs.h"

int cztimqual(char *file, timqual_struc *timqual, int *orb_st_day)
/*******************************************************************

   cztimqual

   purpose: get a time and quality record fro a CZCS L1 file

   Returns type: int - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      name of file to open
      timqual_struc *   timqual          O      structure with time, quality
                                                summary
      int *             orb_st_day      I/O     orbit start day - determined 
                                                from first file's start day 
                                                and < 0 if uninitialized

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 6 Aug 2004      Original development

 *******************************************************************/
 {
    l1_data_struc cz_dat;
    gattr_struc cz_attr;
    int i;
    int32_t fil_st_msec;
    /*
     *  read in the basic L1 data needed to generate the time and quality
     */
    if (cz_l1_read(file, 1, &cz_attr, &cz_dat) != 0) return -1;
    timqual->nscan = cz_attr.scan_lines;
    timqual->n_ctl_pt = cz_attr.n_ctl_pt;
    if (*orb_st_day < 0) *orb_st_day = cz_attr.start_day;
    /*
     *  set up the areas in the timqual struct to receive data
     */
    timqual->qual = (int *) malloc(cz_attr.scan_lines * sizeof ( int));
    timqual->msec = (int32_t *) malloc(cz_attr.scan_lines * sizeof ( int32_t));
    /*
     *  transfer the time and create the quality
     */
    fil_st_msec = cz_dat.msec[0];

    for (i = 0; i < cz_attr.scan_lines; i++) {
        /*
         *  if either the current file's start day is different (1 more) than
         *  the orbit start day or the msec is < the file's start msec,
         *  incriment the msec value in the timqual struct by 1 day (86400000 msec)
         *  this assures steadily rising time in the timqual struct
         */
        timqual->msec[i] = ((cz_dat.msec[i] < fil_st_msec) ||
                (cz_attr.start_day != *orb_st_day)) ?
                cz_dat.msec[i] + 86400000 : cz_dat.msec[i];

        if ((cz_attr.parm_presence & 0XF8) != 0XF8)
            timqual->qual[i] = 1; /* the needed channel missing from file */
        else {
            /*  if the nav may be bad or missing data in band 1-5 for this scan */
            if ((cz_dat.cal_sum[ i * 5 + 3 ] != 0) ||
                    (cz_dat.cal_sum[ i * 5 + 4 ] != 0) ||
                    (cz_dat.cal_scan[ i * 6 ] != 0) ||
                    (cz_dat.cal_scan[ i * 6 + 1 ] != 0) ||
                    (cz_dat.cal_scan[ i * 6 + 2 ] != 0) ||
                    (cz_dat.cal_scan[ i * 6 + 3 ] != 0) ||
                    (cz_dat.cal_scan[ i * 6 + 4 ] != 0))
                timqual->qual[i] = 1;
            else
                timqual->qual[i] = 0;
        }
    }
    /*
     *  remove allocated arrays in data struct
     */
    cz_dat_free(&cz_dat, 1);
    return 0;
}
