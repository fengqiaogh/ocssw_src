#include "fmt_check.h"
extern int verbose;
extern int fmt_status; /* status of check see fmt_check */

void attr_disp(attr_str attr, u_data value, int a_count)
/*******************************************************************

   attr_disp   

   purpose: ~utility to print the info for hdf attributes

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      attr_str          attr            I       attribute info for this 
                                                  attribute
      u_data            value           I       actual value of the attribute
      int               a_count         I       actual count of data items
                                                read

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       1-May-1995     Original development
      W. Robinson       7-Jun-1995     add UINT8 case
      W. Robinson       30-jun-1995    add verbose flag
      W. Robinson       30-Jun-1995     set fmt_status if problems
      W. Robinson       4-Oct-1995     for char data, print the # bytes 
                                       in field: a_count
      W. Robinson       18-Sep-1996    for char data, use a temp string for 
                                       output and insert '?' in slots with 
                                       nulls so they fully print
      W. Robinson       31-Oct-1996    upgrade to handle float64 data type
      W. Robinson       18-Jul-1997    for copying character strings to tstr,
                                       use memcpy

 *******************************************************************/ {
    int32 j;
    char *tstr;

    /*
     *  if vervose is on, write the attribute information
     */
    if (verbose == 0) return;

    /*
     *  depending on whether it is global or part of an SDS, write out slightly
     *  different information
     */
    if (strcmp(attr.obj_nm, "gbl") == 0) {
        switch (attr.type) {
        case DFNT_CHAR:
            tstr = malloc(a_count + 1);
            memcpy(tstr, value.chr, a_count);
            *(tstr + a_count) = 0;
            for (j = 0; j < (a_count - 1); j++)
                if (*(tstr + j) == 0)
                    *(tstr + j) = '?';

            printf("'%s'(%d) is '%s'\n",
                    attr.access_nm, a_count, tstr);
            free(tstr);
            break;
        case DFNT_FLOAT32:
            for (j = 0; j < attr.count; j++) {
                printf("'%s'[%d] is %f\n", attr.access_nm,
                        j, value.f32[j]);
            }
            break;
        case DFNT_FLOAT64:
            for (j = 0; j < attr.count; j++) {
                printf("'%s'[%d] is %f\n", attr.access_nm,
                        j, value.f64[j]);
            }
            break;
        case DFNT_INT16:
            for (j = 0; j < attr.count; j++) {
                printf("'%s'[%d] is %d\n", attr.access_nm,
                        j, value.i16[j]);
            }
            break;
        case DFNT_INT32:
            for (j = 0; j < attr.count; j++) {
                printf("'%s'[%d] is %d\n", attr.access_nm,
                        j, value.i32[j]);
            }
            break;
        case DFNT_INT8:
            for (j = 0; j < attr.count; j++) {
                printf("'%s'[%d] is %d\n", attr.access_nm,
                        j, value.i8[j]);
            }
            break;
        case DFNT_UINT8:
            for (j = 0; j < attr.count; j++) {
                printf("'%s'[%d] is %d\n", attr.access_nm,
                        j, value.ui8[j]);
            }
            break;
        default:
            printf("************* Program Problem DEFAULT CASE OF VALUE CHECKING\n");
            fmt_status = fmt_status | 1;
            break;
        }
    } else {
        switch (attr.type) {
        case DFNT_CHAR:
            tstr = malloc(a_count + 1);
            strncpy(tstr, value.chr, a_count);
            *(tstr + a_count) = 0;
            for (j = 0; j < (a_count - 1); j++)
                if (*(tstr + j) == 0)
                    *(tstr + j) = '?';

            printf("sds: '%s', attr '%s'(%d) is '%s'\n", attr.obj_nm,
                    attr.access_nm, a_count, tstr);
            free(tstr);

            break;
        case DFNT_FLOAT32:
            for (j = 0; j < attr.count; j++) {
                printf("sds: '%s', attr '%s'[%d] is %f\n", attr.obj_nm,
                        attr.access_nm, j, value.f32[j]);
            }
            break;
        case DFNT_FLOAT64:
            for (j = 0; j < attr.count; j++) {
                printf("sds: '%s', attr '%s'[%d] is %f\n", attr.obj_nm,
                        attr.access_nm, j, value.f64[j]);
            }
            break;
        case DFNT_INT16:
            for (j = 0; j < attr.count; j++) {
                printf("sds: '%s', attr '%s'[%d] is %d\n", attr.obj_nm,
                        attr.access_nm, j, value.i16[j]);
            }
            break;
        case DFNT_INT32:
            for (j = 0; j < attr.count; j++) {
                printf("sds: '%s', attr '%s'[%d] is %d\n", attr.obj_nm,
                        attr.access_nm, j, value.i32[j]);
            }
            break;
        default:
            printf("************* Program Problem DEFAULT CASE OF VALUE CHECKING\n");
            fmt_status = fmt_status | 1;
            break;
        }
    }
}
