#include "l1czcs.h"

void fill_mstr(int *mstr_last, mstr_struc *mstr, timqual_struc *init_info,
        int ds_num, int st_scan, int en_scan)
/*******************************************************************

   fill_mstr

   purpose: fill the master list of czcs sczn info

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int *             mstr_last       I/O     last filled line of the master 
                                                list
      mstr_struc *      mstr            I/O     structure of the master list
      timqual_struc *   init_info        I      structure with time, quality
                                                summary to add to master list
      int               ds_num           I      current dataset in use
      int               st_scan          I      start scan of new data to add
                                                to master list
      int               en_scan          I      and end scan (see above)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 6 Aug 2004      Original development

 *******************************************************************/
 {
    int in_st_scan, nscan, i, j;
    int32_t step;
    /*
     *  this will preserve a 124 msec scan line step in the master list -
     *  any gaps will be marked as having non-existant data (exist = 0)
     *
     *  (if we start with scan 0 as input, don't worry about 124 msec spacing
     */
    printf("fill_mstr: filling master struct after index %d with input scans %d to %d\n",
            *mstr_last, st_scan, en_scan);
    if (st_scan == 0) {
        in_st_scan = 1;
        (*mstr_last)++;
        mstr->msec[*mstr_last] = init_info->msec[0];
        mstr->exist[*mstr_last] = 1;
        mstr->qual[*mstr_last] = init_info->qual[0];
        mstr->ds_num[*mstr_last] = ds_num;
        mstr->scan[*mstr_last] = (int32_t) st_scan;
    } else
        in_st_scan = st_scan;

    nscan = en_scan - in_st_scan;
    if (nscan > 0) {
        /*
         *  transfer the msec and quality over to the master list, expanding
         *  any time gaps 
         */
        for (i = 0; i < nscan; i++) {
            step = (init_info->msec[in_st_scan] - mstr->msec[*mstr_last] + 60) /
                    124;
            /*
             *  note that the msec step can be off by as much as +-30 msec, so the
             *  60 msec will generously account for it
             *  very large or < 1 time steps will cause skipping of the line 
             *  during fill
             */
            if ((step > 900) || (step <= 0)) {
                printf("fill_mstr: Unrealistic time gap found at scan %d\n",
                        i);
                printf("           Gap of %d lines.  Skipping this line\n",
                        step);
            } else {
                if (step > 1) {
                    printf(
                            "fill_mstr: gap of %d lines found in filling master at scan # %d\n",
                            step, i);
                    for (j = 0; j < step - 1; j++) {
                        mstr->exist[ *mstr_last + 1 ] = 0;
                        mstr->msec[ *mstr_last + 1 ] = mstr->msec[ *mstr_last ] + 124;
                        (*mstr_last)++;
                    }
                }
                /*
                 *  fill new entry
                 */
                (*mstr_last)++;
                mstr->msec[*mstr_last] = init_info->msec[in_st_scan];
                mstr->exist[*mstr_last] = 1;
                mstr->qual[*mstr_last] = init_info->qual[in_st_scan];
                mstr->ds_num[*mstr_last] = ds_num;
                mstr->scan[*mstr_last] = (int32_t) in_st_scan;
            }
            in_st_scan++;
        }
    }
    return;
}
