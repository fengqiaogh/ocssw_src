#include "l1czcs.h"

void olap_resolve(mstr_struc *mstr, int mstr_last, int new_start, int new_end,
        int ds_num)
/*******************************************************************

   olap_resolve

   purpose: resolve the overlap of new data over current data in master list
     Use the quality  / existance values to pick the best replacement scans

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      mstr_struc *      mstr            I/O     master structure of scan info
      int               mstr_last        I      last line in master list
      int               new_start        I      start of new overlapping data
      int               new_end          I      end of new overlapping data
      int               ds_num           I      dataset id of new data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 6 Aug 2004      Original developmento

  Note that if the overlap is inside the master range, it will either 
  completely replace it or not at all.  If the overlap is at the end, 
  the best scan line will be found in the overlap to begin inserting scans 
  from the new data.  This will preserve the scan continuity while getting
  the best quality possible.

 *******************************************************************/
 {
    int cum_cur, cum_new, *ca_cur, *ca_new, i, olap_st, olap_en, n_olap,
            break_scan, sum;
    /*
     *  first, decide if the new data range is inside (or at the end of) 
     *  the current data
     */
    cum_cur = 0;
    cum_new = 0;

    if (new_end < mstr_last) {
        /*
         *  overlap is inside current range, get the cumulative # scans with 
         *  bad status, non-existent scans for current and new scans
         */
        printf("olap_resolve: overlap is inside current range\n");

        for (i = new_start; i < new_end; i++) {
            if (mstr->exist[i] == 0)
                cum_cur++;
            else if (mstr->qual[i] != 0)
                cum_cur++;

            if (mstr->in_exist[i] == 0)
                cum_new++;
            else if (mstr->in_qual[i] != 0)
                cum_new++;
        }
        printf("olap_resolve: cum_new = %d, cum_cur = %d\n", cum_new, cum_cur);
        if (cum_new < cum_cur) {
            /*
             *  new data has better overall quality in the overlap than current data.
             *  so, replace the current data
             */
            printf("olap_resolve: new data better, cum_new = %d, cum_cur = %d\n",
                    cum_new, cum_cur);
            for (i = new_start; i < (new_end + 1); i++) {
                mstr->exist[i] = mstr->in_exist[i];
                mstr->qual[i] = mstr->in_qual[i];
                mstr->ds_num[i] = ds_num;
                mstr->msec[i] = mstr->in_msec[i];
                mstr->scan[i] = mstr->in_scan[i];
            }
        } else
            printf("olap_resolve: current data has equal of better quality\n");
    } else {
        /*
         *  overlap is at / beyond the current end of data.  get overlap range
         */
        printf("olap_resolve: overlap is at / beyond the current end\n");
        olap_st = new_start;
        olap_en = mstr_last;
        n_olap = olap_en - olap_st + 1;
        /*
         *  allocate cumulative quality arrays for the current and new overlaps
         */
        ca_cur = malloc(n_olap * sizeof ( int));
        ca_new = malloc(n_olap * sizeof ( int));
        /*
         * loop through the overlap range, forward for current data and backward for
         * new data.  accumulate # quality problems for each
         */
        for (i = 0; i < n_olap; i++) {
            if (mstr->exist[ i + olap_st ] == 0)
                cum_cur++;
            else if (mstr->qual[ i + olap_st ] != 0)
                cum_cur++;

            ca_cur[i] = cum_cur;

            if (mstr->in_exist[ olap_en - i ] == 0)
                cum_new++;
            else if (mstr->in_qual[ olap_en - i ] != 0)
                cum_new++;

            ca_new[ n_olap - 1 - i ] = cum_new;
        }
        printf("olap_resolve: current total bad: %d, new total bad: %d\n",
                cum_cur, cum_new);
        /*
         *  Cumulative arrays set.  go thru and find scan where total problems
         *  are least.  Note that this favors the new file in equally good data
         */
        sum = n_olap * 3; /* start at a high valus */
        break_scan = -1;
        for (i = 0; i < n_olap; i++) {
            if ((ca_new[i] + ca_cur[i]) < sum) {
                sum = ca_new[i] + ca_cur[i];
                break_scan = i;
            }
        }
        printf("olap_resolve: break at line %d of overlap, total # olap = %d\n",
                break_scan, n_olap);
        printf("olap_resolve: lowest sum is %d\n", sum);
        /*
         *  the best place to break the overlap is found.  replace current data
         *  with new data beyond this line
         */
        if (break_scan < (n_olap - 1)) {
            for (i = break_scan; i < n_olap; i++) {
                mstr->exist[ i + olap_st ] = mstr->in_exist[ i + olap_st ];
                mstr->qual[ i + olap_st ] = mstr->in_qual[ i + olap_st ];
                mstr->ds_num[ i + olap_st ] = ds_num;
                mstr->msec[ i + olap_st ] = mstr->in_msec[ i + olap_st ];
                mstr->scan[ i + olap_st ] = mstr->in_scan[ i + olap_st ];
            }
        }
        /*
         *  free the ca_cur, ca_new
         */
        free(ca_cur);
        free(ca_new);
    }
    return;
}
