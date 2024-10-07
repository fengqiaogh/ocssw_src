//#ifndef LINUX
//#include <ieeefp.h>
//#endif

#include <math.h>

#include "fmt_check.h"
extern int verbose; /* 0 - don't print info on each vgroup, 1 do  */
extern int fmt_status; /* status of format check, see fmt_check */

int chk_sds(int32 fid, fmt_str *fmt, int isds)
/*******************************************************************

   chk_sds

   purpose: for a particular science dataset, check the existance of it
            in the dataset and make sure it has the specified 
            attributes and that the type and ranges of the dataset are 
            correct. 

   Returns type: int - 0 if all went well, 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             fid              I      file id or sds id
      fmt_str *         fmt              I      format description struct
      int               isds             I      index of the current SDS 
                                                  being checked

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       18-Feb-1995     Original development
      W. Robinson       30-Jun-1995     set fmt_status if problems
      W. Robinson       4-Oct-1995      add a return of the actual count of 
                                        attr items read
      W. Robinson       1-Oct-1996      add optional checking of SDS contents 
                                        within min / max range
      W. Robinson       31-Oct-1996     upgrade to handle float64 data values
      W. Robinson       31-May-2000     check for NaN in floats and not finite
                                        in doubles
      W. Robinson       15-Mar-2001     update with sone linux changes libraries
      W. Robinson, SAIC 31 Mar 2005     changes for the table driven version

 *******************************************************************/
 {
    int32 ind, i, erng, dimsiz, a_rank, a_type, a_range[4], a_nattr,
            sds_id, index, irem, icount, tot_pts, start[4], *ptr_i32;
    char range_str[50], s_range_str[100];
    void *ptr;
    float32 *ptr_f32;
    float64 *ptr_f64;
    int16 *ptr_i16;
    int8 *ptr_i8;
    uint8 *ptr_ui8;
    int a_count;
    u_data value;
    s_data s_range; /* actual sds range found (if looking) */
    sds_info_str cur_sds; /* current SDS being looked at  */

    /*
     *  get the current SDS for convenience
     */
    cur_sds = fmt->sds_info[isds];

    /*
     *  get the index for the science dataset
     */
    if ((index = SDnametoindex(fid, cur_sds.name)) < 0) {
        printf("**** SDS: '%s' cannot be located in the file\n",
                cur_sds.name);
        fmt_status = fmt_status | 2;
        return -1;
    }
    /*
     *  next, get the sds id
     */
    if ((sds_id = SDselect(fid, index)) < 0) {
        printf("**** SDS: '%s', cannot get the ID from index\n",
                cur_sds.name);
        fmt_status = fmt_status | 2;
        return -1;
    }
    /*
     *  next, find out about the dataset and check with expectations
     */
    if (SDgetinfo(sds_id, NULL, &a_rank, a_range, &a_type, &a_nattr)
            != 0) {
        printf("**** SDS: '%s', Unable to read the info on SDS\n",
                cur_sds.name);
        fmt_status = fmt_status | 2;
    }

    if (a_rank != cur_sds.rank) {
        printf("**** SDS: '%s', Mismatch in the rank\n", cur_sds.name);
        printf("expected: %d, read: %d\n", cur_sds.rank, a_rank);
        fmt_status = fmt_status | 2;
    }

    for (i = 0; i < a_rank; i++) {
        erng = cur_sds.e_ranges[i];
        if (erng != 0) /* if we need to check dimensions... */ {
            /*
             *  get the dimension size: erng if > 0, dimension struct entry if < 0
             */
            dimsiz = (erng > 0) ? erng : fmt->dim_id[ (-erng) - 1 ].dim_size;

            if (a_range[i] != dimsiz) {
                printf("**** SDS: '%s', Mismatch in length of dim # %d\n",
                        cur_sds.name, i);
                printf("expected: %d, read: %d\n", dimsiz, a_range[i]);
                fmt_status = fmt_status | 2;
            }
        }
    }
    if (a_type != cur_sds.type) {
        printf("**** SDS: '%s', Mismatch data type\n", cur_sds.name);
        printf("expected: %d, read: %d\n", cur_sds.type,
                a_type);
        fmt_status = fmt_status | 2;
    }
    if ((cur_sds.n_attr >= 0) && (a_nattr != cur_sds.n_attr)) {
        printf("**** SDS: '%s', Mismatch in # attributes for sds\n",
                cur_sds.name);
        printf("expected: %d, read: %d\n", cur_sds.n_attr, a_nattr);
        fmt_status = fmt_status | 2;
    }

    /*
     *  if signaled, check the SDS values within a range
     */
    if (cur_sds.flg_rng == 1) {
        /*
         *  read the SDS
         */
        tot_pts = 1;
        for (i = 0; i < a_rank; i++) {
            start[i] = 0;
            tot_pts *= a_range[i];
        }
        ptr = (void *) malloc(tot_pts * cur_sds.byt_per_val);

        if (SDreaddata(sds_id, start, NULL, a_range, ptr) < 0) {
            printf("**** SDS: '%s', SDreaddata failed reading %d values\n",
                    cur_sds.name, tot_pts);
            fmt_status = fmt_status | 2;
        } else {
            /*
             *  for each type that is processed, check each value within the min
             *  / max range or is not a number (not finite for doubles)
             */
            icount = 0;
            irem = -1;
            switch (a_type) {
            case DFNT_FLOAT32:
                sprintf(range_str, "%f and %f", cur_sds.sds_rng.f32[0],
                        cur_sds.sds_rng.f32[1]);
                ptr_f32 = (float32 *) ptr;

                s_range.f32[0] = *ptr_f32;
                s_range.f32[1] = *ptr_f32;

                for (i = 0; i < tot_pts; i++) {
                    if (*(ptr_f32 + i) < s_range.f32[0])
                        s_range.f32[0] = *(ptr_f32 + i);
                    if (*(ptr_f32 + i) > s_range.f32[1])
                        s_range.f32[1] = *(ptr_f32 + i);

                    if (*(ptr_f32 + i) < cur_sds.sds_rng.f32[0] ||
                            *(ptr_f32 + i) > cur_sds.sds_rng.f32[1] ||
                            isnan(*(ptr_f32 + i))) {
                        irem = (irem < 0) ? i : irem; /* remember first occurence */
                        icount++; /* and total # of occurences */
                    }
                }
                sprintf(s_range_str, "%f to %f", s_range.f32[0], s_range.f32[1]);
                break;
            case DFNT_FLOAT64:
                sprintf(range_str, "%f and %f", cur_sds.sds_rng.f64[0],
                        cur_sds.sds_rng.f64[1]);
                ptr_f64 = (float64 *) ptr;

                s_range.f64[0] = *ptr_f64;
                s_range.f64[1] = *ptr_f64;

                for (i = 0; i < tot_pts; i++) {
                    if (*(ptr_f64 + i) < s_range.f64[0])
                        s_range.f64[0] = *(ptr_f64 + i);
                    if (*(ptr_f64 + i) > s_range.f64[1])
                        s_range.f64[1] = *(ptr_f64 + i);

                    if (*(ptr_f64 + i) < cur_sds.sds_rng.f64[0] ||
                            *(ptr_f64 + i) > cur_sds.sds_rng.f64[1] ||
                            !isfinite(*(ptr_f64 + i))) {
                        irem = (irem < 0) ? i : irem; /* remember first occurence */
                        icount++;
                        /* and total # of occurences */                    }
                }
                sprintf(s_range_str, "%f to %f", s_range.f64[0], s_range.f64[1]);
                break;
            case DFNT_INT32:
                sprintf(range_str, "%d and %d", cur_sds.sds_rng.i32[0],
                        cur_sds.sds_rng.i32[1]);
                ptr_i32 = (int32 *) ptr;

                s_range.i32[0] = *ptr_i32;
                s_range.i32[1] = *ptr_i32;

                for (i = 0; i < tot_pts; i++) {
                    if (*(ptr_i32 + i) < s_range.i32[0])
                        s_range.i32[0] = *(ptr_i32 + i);
                    if (*(ptr_i32 + i) > s_range.i32[1])
                        s_range.i32[1] = *(ptr_i32 + i);

                    if (*(ptr_i32 + i) < cur_sds.sds_rng.i32[0] ||
                            *(ptr_i32 + i) > cur_sds.sds_rng.i32[1]) {
                        irem = (irem < 0) ? i : irem; /* remember first occurence */
                        icount++;
                        /* and total # of occurences */                    }
                }
                sprintf(s_range_str, "%d to %d", s_range.i32[0], s_range.i32[1]);
                break;
            case DFNT_INT16:
                sprintf(range_str, "%d and %d", cur_sds.sds_rng.i16[0],
                        cur_sds.sds_rng.i16[1]);
                ptr_i16 = (int16 *) ptr;

                s_range.i16[0] = *ptr_i16;
                s_range.i16[1] = *ptr_i16;

                for (i = 0; i < tot_pts; i++) {
                    if (*(ptr_i16 + i) < s_range.i16[0])
                        s_range.i16[0] = *(ptr_i16 + i);
                    if (*(ptr_i16 + i) > s_range.i16[1])
                        s_range.i16[1] = *(ptr_i16 + i);

                    if (*(ptr_i16 + i) < cur_sds.sds_rng.i16[0] ||
                            *(ptr_i16 + i) > cur_sds.sds_rng.i16[1]) {
                        irem = (irem < 0) ? i : irem; /* remember first occurence */
                        icount++;
                        /* and total # of occurences */                    }
                }
                sprintf(s_range_str, "%d to %d", s_range.i16[0], s_range.i16[1]);
                break;
            case DFNT_INT8:
                sprintf(range_str, "%d and %d", cur_sds.sds_rng.i8[0],
                        cur_sds.sds_rng.i8[1]);
                ptr_i8 = (int8 *) ptr;

                s_range.i8[0] = *ptr_i8;
                s_range.i8[1] = *ptr_i8;

                for (i = 0; i < tot_pts; i++) {
                    if (*(ptr_i8 + i) < s_range.i8[0])
                        s_range.i8[0] = *(ptr_i8 + i);
                    if (*(ptr_i8 + i) > s_range.i8[1])
                        s_range.i8[1] = *(ptr_i8 + i);

                    if (*(ptr_i8 + i) < cur_sds.sds_rng.i8[0] ||
                            *(ptr_i8 + i) > cur_sds.sds_rng.i8[1]) {
                        irem = (irem < 0) ? i : irem; /* remember first occurence */
                        icount++;
                        /* and total # of occurences */                    }
                }
                sprintf(s_range_str, "%d to %d", s_range.i8[0], s_range.i8[1]);
                break;
            case DFNT_UINT8:
                sprintf(range_str, "%d and %d", cur_sds.sds_rng.ui8[0],
                        cur_sds.sds_rng.ui8[1]);
                ptr_ui8 = (uint8 *) ptr;

                s_range.ui8[0] = *ptr_ui8;
                s_range.ui8[1] = *ptr_ui8;

                for (i = 0; i < tot_pts; i++) {

                    if (*(ptr_ui8 + i) < s_range.ui8[0])
                        s_range.ui8[0] = *(ptr_ui8 + i);
                    if (*(ptr_ui8 + i) > s_range.ui8[1])
                        s_range.ui8[1] = *(ptr_ui8 + i);

                    if (*(ptr_ui8 + i) < cur_sds.sds_rng.ui8[0] ||
                            *(ptr_ui8 + i) > cur_sds.sds_rng.ui8[1]) {
                        irem = (irem < 0) ? i : irem; /* remember first occurence */
                        icount++;
                        /* and total # of occurences */                    }
                }
                sprintf(s_range_str, "%d to %d", s_range.ui8[0], s_range.ui8[1]);
                break;
            }
            /* 
             *  report any problem sdses Also, print dimension in '[]'
             */
            if (icount > 0) {
                printf("**** SDS: '%s', has a value outside valid range of %s.\n",
                        cur_sds.name, range_str);
                printf("     (or is not a number)\n");
                printf("     first point at linear loc. %d, # invalid: %d, dims: [",
                        irem, icount);
                for (i = 0; i < a_rank; i++)
                    printf(" %d ", a_range[i]);
                printf("]\n");
                fmt_status = fmt_status | 2;
            }

            /* 
             *  report sds range if desired 
             */
            if (verbose == 1) {
                printf("sds: '%s' has a min max range of: %s\n", cur_sds.name,
                        s_range_str);
            }
        }
        free(ptr);
    } else {
        if (verbose == 1) {
            printf("sds: '%s' was checked\n", cur_sds.name);
        }
    }
    /*
     *  Now, check the attributes for this sds
     */
    if (cur_sds.n_attr > 0) {
        for (i = 0; i < cur_sds.n_attr; i++) {
            ind = cur_sds.attr_ind[i];
            if (get_attr(sds_id, fmt->att[ ind ], &value, &a_count)
                    == 0) {
                /*
                 *  code to print attribute info
                 */
                attr_disp(fmt->att[ ind ], value, a_count);
            }
        }
    }
    return 0;
}
