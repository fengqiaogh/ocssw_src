#ifdef IRIX
#include <ieeefp.h>
#endif

#include <math.h>
#include "fmt_check.h"

extern int fmt_status; /* format checking status, see fmt_check */

int get_attr(int32 fid, attr_str attr, u_data *value, int *d_count)
/*******************************************************************

   get_attr

   purpose: for a particular attribute, check the existance of it
            in the dataset and make sure it has the specified 
            attribute name, number type and count.  If so,
            read the attribute into the data area

   Returns type: int - 0 if all went well, 
                       -1 can't find the attribute
                       -2 bad match with name, number type or count
                       -3 unable to read the data

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             fid              I      file id or sds id
      attr_str          attr             I      attribute description
                                                structure. 
      u_data *          value            O      returned data structure
      int*              d_count          O      data item count read

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       27-Jan-1995     Original development
      W. Robinson       30-Jun-1995     assure that only 3000 bytes of
                                        attribute data are read
      W. Robinson       30-Jun-1995     set fmt_status if problems
      W. Robinson       4-Oct-1995      save actual count (a_count) for
                                        printing
      W. Robinson       20-Nov-1995     add check for strings so they have 
                                        only 1 null at the end (str_chk)
      L. Kumar		19-Dec-1995     after gets back from str_chk, if the
                                        string lack a null, place a null at
                                        the end of string
      W. Robinson       17-Sep-1996     increase the string handeling 
                                        capability to 66000 bytes (enough 
                                        for HDF 4.0)
      W. Robinson       26-Sep-1996     enable attribute range checking
      W. Robinson       31-Oct-1996     upgrade to handle float64 data type
      W. Robinson       31-May-2000     check the floats for NaN if read and
                                        the doubles to be finite
      W. Robinson       15-Mar-2001     update for a linux array size change
                                        (hope it always works)
      W. Robinson       3 Feb 2010      allow a small range in float attrib
                                        value checking

 *******************************************************************/
 {
    int i;
    int32 attr_index;
    char a_name[200], *temp;
    int32 a_n_type, a_count;

    void chk_str(attr_str, char *, int32);

    if ((attr_index = SDfindattr(fid, attr.access_nm))
            == -1) {
        printf("**** object: '%s', Error trying to get attribute '%s'\n",
                attr.obj_nm, attr.access_nm);
        fmt_status = fmt_status | 2;
        return -1;
    }
    /*
     *  Call SDattrinfo to get the info on this attribute
     */
    if (SDattrinfo(fid, attr_index, a_name, &a_n_type, &a_count)
            == -1) {
        printf("**** Error in SDattrinfo for object: '%s', attribute: '%s'\n",
                attr.obj_nm, attr.access_nm);
        fmt_status = fmt_status | 2;
        return -4;
    }
    /*
     *  put actual count in output variable
     */
    *d_count = a_count;

    /*
     *  Check the name, number type and count
     */
    i = 0;
    if (strcmp(a_name, attr.int_nm) != 0) {
        printf(
                "**** object: '%s', attribute: '%s', difference in attribute name found\n",
                attr.obj_nm, attr.access_nm);
        printf("expected: '%s', read: '%s'\n", attr.int_nm, a_name);
        fmt_status = fmt_status | 2;
        i = 1;
    }
    if (a_n_type != attr.type) {
        printf(
                "**** object: '%s', attribute: '%s', difference in number type found\n",
                attr.obj_nm, attr.access_nm);
        printf("expected: %d, read: %d\n", attr.type, a_n_type);
        fmt_status = fmt_status | 2;
        i = 1;
    }
    if (attr.count > 0 && a_count != attr.count) {
        printf("**** object: '%s', attribute: '%s', difference in count found\n",
                attr.obj_nm, attr.access_nm);
        printf("expected: %d, read: %d\n", attr.count, a_count);
        fmt_status = fmt_status | 2;
        i = 1;
    }
    if (i != 0) return -2;
    /*
     *  read in the attribute
     */
    if (attr.read != ATT_RD_NOREAD) {
        /*
         *  in the event the attribute gets larger than the # bytes
         *  in value (a generous 66000 bytes!), allocate a temp read area
         *  and transfer the first 66000 bytes
         *   ** Use variable FMT_8BLEN for the 66000 **
         */
        if (a_count * DFKNTsize(a_n_type) >= FMT_8BLEN * 8) {
            printf("**** object: '%s', attribute: '%s' has excessive count of \n",
                    attr.obj_nm, attr.access_nm);
            printf("%d items. Truncated list of bytes will be used\n",
                    a_count);
            /*   remove this as a problem
             fmt_status = fmt_status | 2;
             */

            temp = (char *) malloc(a_count * DFKNTsize(a_n_type));
            if (SDreadattr(fid, attr_index, temp) == -1) {
                printf("**** object: '%s', attribute: '%s', Error reading data\n",
                        attr.obj_nm, attr.access_nm);
                return -3;
                fmt_status = fmt_status | 2;
            }
            /*
             *  check for null at end only
             */
            chk_str(attr, temp, a_count);

            memcpy(value->chr, temp, FMT_8BLEN);
            free(temp);
        } else {
            if (SDreadattr(fid, attr_index, value->chr) == -1) {
                printf("**** object: '%s', attribute: '%s', Error reading data\n",
                        attr.obj_nm, attr.access_nm);
                fmt_status = fmt_status | 2;
                return -3;
            }
            /*
             *  check for null at end only
             */
            chk_str(attr, value->chr, a_count);
        }
        /*
         *  any read-in floats can be checked for NaN status and the doubles can
         *  be checked for finite status
         */
        switch (attr.type) {
        case DFNT_FLOAT32:
            for (i = 0; i < attr.count; i++) {
                if (isnan(value->f32[i])) {
                    printf(
                            "**** object: '%s', attribute '%s'[%d]\n",
                            attr.obj_nm, attr.access_nm, i);
                    printf("is not a number (NaN)\n");
                    fmt_status = fmt_status | 2;
                }
            }
            break;
        case DFNT_FLOAT64:
            for (i = 0; i < attr.count; i++) {
                if (!isfinite(attr.data.f64[i])) {
                    printf(
                            "**** object: '%s', attribute '%s'[%d]\n",
                            attr.obj_nm, attr.access_nm, i);
                    printf("is not a number (NaN) or not finite\n");
                    fmt_status = fmt_status | 2;
                }
            }
            break;
        }
        /*
         *  Check the constant attribute values (read = READ_ONE_VAL) 
         *  and handle variables in a range to calling routine
         */
        if (attr.read == ATT_RD_READ_ONE_VAL) {
            switch (attr.type) {
            case DFNT_CHAR:
                if (value->chr[a_count - 1] != 0) /* LK */
                    value->chr[a_count] = 0; /* LK */
                if (strcmp(value->chr, attr.data.chr) != 0) {
                    printf("**** object: '%s', difference in attribute '%s'\n",
                            attr.obj_nm, attr.access_nm);
                    printf("expected value: '%s'\nread value    : '%s'\n",
                            attr.data.chr, value->chr);
                    fmt_status = fmt_status | 2;
                }
                break;
            case DFNT_FLOAT32:
                for (i = 0; i < attr.count; i++) {
                    if ((value->f32[i] < attr.data.f32[i] -
                            fabsf(attr.data.f32[i]) * ERR_FRAC_F32) ||
                            (value->f32[i] > attr.data.f32[i] +
                            fabsf(attr.data.f32[i]) * ERR_FRAC_F32)) {
                        printf(
                                "**** object: '%s', difference in attribute '%s'[%d]\n",
                                attr.obj_nm, attr.access_nm, i);
                        printf("expected value: %f\nread value: %f\n",
                                attr.data.f32[i], value->f32[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_FLOAT64:
                for (i = 0; i < attr.count; i++) {
                    if ((value->f64[i] < attr.data.f64[i] -
                            fabs(attr.data.f64[i]) * ERR_FRAC_F64) ||
                            (value->f64[i] > attr.data.f64[i] +
                            fabs(attr.data.f64[i]) * ERR_FRAC_F64)) {
                        printf(
                                "**** object: '%s', difference in attribute '%s'[%d]\n", attr.obj_nm, attr.access_nm, i);
                        printf("expected value: %f\nread value: %f\n",
                                attr.data.f64[i], value->f64[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_INT8:
                for (i = 0; i < attr.count; i++) {
                    if (value->i8[i] != attr.data.i8[i]) {
                        printf(
                                "**** object: '%s', difference in attribute '%s'[%d]\n", attr.obj_nm, attr.access_nm, i);
                        printf("expected value: %d\nread value: %d\n",
                                attr.data.i8[i], value->i8[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_UINT8:
                for (i = 0; i < attr.count; i++) {
                    if (value->ui8[i] != attr.data.ui8[i]) {
                        printf(
                                "**** object: '%s', difference in attribute '%s'[%d]\n", attr.obj_nm, attr.access_nm, i);
                        printf("expected value: %d\nread value: %d\n",
                                attr.data.ui8[i], value->ui8[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_INT16:
                for (i = 0; i < attr.count; i++) {
                    if (value->i16[i] != attr.data.i16[i]) {
                        printf(
                                "**** object: '%s', difference in attribute '%s'[%d]\n",
                                attr.obj_nm, attr.access_nm, i);
                        printf("expected value: %d\nread value: %d\n",
                                attr.data.i16[i], value->i16[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_INT32:
                for (i = 0; i < attr.count; i++) {
                    if (value->i32[i] != attr.data.i32[i]) {
                        printf(
                                "**** object: '%s', difference in attribute '%s'[%d]\n",
                                attr.obj_nm, attr.access_nm, i);
                        printf("expected value: %d\nread value: %d\n",
                                attr.data.i32[i], value->i32[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            default:
                printf("************* DEFAULT CASE OF VALUE CHECKING\n");
                printf("Program error\n");
                fmt_status = fmt_status | 1;
                break;
            }
        } else if (attr.read == ATT_RD_READ_INCLUSIVE) {
            switch (attr.type) {
            case DFNT_CHAR:
                if (value->chr[a_count - 1] != 0) /* LK */
                    value->chr[a_count] = 0; /* LK */
                if (strcmp(value->chr, attr.data.chr) != 0) {
                    printf("**** object: '%s', difference in attribute '%s'\n",
                            attr.obj_nm, attr.access_nm);
                    printf("expected value: '%s'\nread value    : '%s'\n",
                            attr.data.chr, value->chr);
                    fmt_status = fmt_status | 2;
                }
                break;
            case DFNT_FLOAT32:
                for (i = 0; i < attr.count; i++) {
                    if (value->f32[i] < attr.data.f32[0] ||
                            value->f32[i] > attr.data.f32[1]) {
                        printf(
                                "**** object: '%s', outside range in attribute '%s'[%d]\n",
                                attr.obj_nm, attr.access_nm, i);
                        printf("valid range: %f - %f\nread value: %f\n",
                                attr.data.f32[0], attr.data.f32[1], value->f32[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_FLOAT64:
                for (i = 0; i < attr.count; i++) {
                    if (value->f64[i] < attr.data.f64[0] ||
                            value->f64[i] > attr.data.f64[1]) {
                        printf(
                                "**** object: '%s', outside range in attribute '%s'[%d]\n", attr.obj_nm, attr.access_nm, i);
                        printf("valid range: %f - %f\nread value: %f\n",
                                attr.data.f64[0], attr.data.f64[1], value->f64[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_INT8:
                for (i = 0; i < attr.count; i++) {
                    if (value->i8[i] < attr.data.i8[0] ||
                            value->i8[i] > attr.data.i8[1]) {
                        printf(
                                "**** object: '%s', outside range in attribute '%s'[%d]\n", attr.obj_nm, attr.access_nm, i);
                        printf("valid range: %d - %d\nread value: %d\n",
                                attr.data.i8[0], attr.data.i8[1], value->i8[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_UINT8:
                for (i = 0; i < attr.count; i++) {
                    if (value->ui8[i] < attr.data.ui8[0] ||
                            value->ui8[i] > attr.data.ui8[1]) {
                        printf(
                                "**** object: '%s', outside range in attribute '%s'[%d]\n", attr.obj_nm, attr.access_nm, i);
                        printf("valid range: %d - %d\nread value: %d\n",
                                attr.data.ui8[0], attr.data.ui8[1], value->ui8[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_INT16:
                for (i = 0; i < attr.count; i++) {
                    if (value->i16[i] < attr.data.i16[0] ||
                            value->i16[i] > attr.data.i16[1]) {
                        printf(
                                "**** object: '%s', outside range in attribute '%s'[%d]\n",
                                attr.obj_nm, attr.access_nm, i);
                        printf("valid range: %d - %d\nread value: %d\n",
                                attr.data.i16[0], attr.data.i16[1], value->i16[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            case DFNT_INT32:
                for (i = 0; i < attr.count; i++) {
                    if (value->i32[i] < attr.data.i32[0] ||
                            value->i32[i] > attr.data.i32[1]) {
                        printf(
                                "**** object: '%s', outside range in attribute '%s'[%d]\n",
                                attr.obj_nm, attr.access_nm, i);
                        printf("valid range: %d - %d\nread value: %d\n",
                                attr.data.i32[0], attr.data.i32[1], value->i32[i]);
                        fmt_status = fmt_status | 2;
                    }
                }
                break;
            default:
                printf("************* DEFAULT CASE OF VALUE CHECKING\n");
                printf("Program error\n");
                fmt_status = fmt_status | 1;
                break;
            }
        }
    }
    return 0;
}
