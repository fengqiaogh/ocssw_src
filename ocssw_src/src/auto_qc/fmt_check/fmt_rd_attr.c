#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmt_check.h"
extern int fmt_status; /* format check status, see fmt_check */

static char *att_rd_str[] = {"NOREAD", "READ_NOCK", "READ_ONE_VAL",
    "READ_INCLUSIVE"};

int fmt_rd_attr(char *file, FILE *fid, fmt_str *fmt)
/*******************************************************************

   fmt_rd_attr

   purpose: read in the attribute portion of the format description for 
            a dataset

   Returns type: int - 0 if all went well, 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      file name containing the format
                                                table description
      FILE *            fid              I      format table file handle
      fmt_str *         fmt             I/O     format structure

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 29 Mar 2005     Original development
                                        (a portion of the old fmt_read)

 *******************************************************************/
 {
    char line[500], line_sav[500], str_typ[200], *str;
    int i, j, k, ifound, ndata, ierr;

    /*
     *  allocate space for the attribute descriptions
     */
    if ((fmt->att = malloc(fmt->n_attr * sizeof ( attr_str))) == NULL) {
        printf("**************Program error\n");
        printf("in: %s, attrib allocation failed\n", __FILE__);
        fmt_status = fmt_status | 1;
        return -1;
    }
    /*
     *  read in each line of attribute description to the proper place
     */
    for (i = 0; i < fmt->n_attr; i++) {
        if (get_line(line, 500, fid, '#') == NULL) {
            printf("**************Program error\n");
            printf("unable to read attribute description # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line);
            fmt_status = fmt_status | 1;
            return -1;
        }

        strcpy(line_sav, line); /* save an in-tact copy of line */
        fmt->att[i].dim_index = -1;

        /*
         *  use the s_parse so that any quoted strings get taken together
         */
        if ((str = s_parse(line, '\"')) == NULL) {
            printf("**************Program error\n");
            printf(
                    "unable to read attr ctl value 1 of description # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        strcpy(fmt->att[i].obj_nm, str);

        for (j = 0; j < 5; j++) {
            if ((str = s_parse(NULL, '\"')) == NULL) {
                printf("**************Program error\n");
                printf(
                        "unable to read attr ctl value %d of description # %d of file:\n'%s'\n",
                        j + 2, i, file);
                printf("last line read:\n'%s'\n", line_sav);
                fmt_status = fmt_status | 1;
                return -1;
            }
            switch (j) {
            case 0:
            {
                strcpy(fmt->att[i].access_nm, str);
                break;
            }
            case 1:
            {
                strcpy(fmt->att[i].int_nm, str);
                break;
            }
            case 2:
            {
                strcpy(str_typ, str);
                break;
            }
            case 3:
            {
                if (sscanf(str, "%d", &(fmt->att[i].count)) == EOF) {
                    printf("**************Program error\n");
                    printf(
                            "unable to decode att ctl value count in desc # %d of file:\n'%s'\n",
                            i, file);
                    printf("last line read:\n'%s'\n", line_sav);
                    fmt_status = fmt_status | 1;
                    return -1;
                }
                break;
            }
            case 4:
            {
                ifound = 0;
                for (k = 0; k < 4; k++) {
                    if (strcmp(str, att_rd_str[k]) == 0) {
                        fmt->att[i].read = k;
                        ifound = 1;
                        break;
                    }
                }
                if (ifound == 0) {
                    printf("**************Program error\n");
                    printf("Improper read code for the attribute\n");
                    printf("desc # %d of file:\n'%s'\n", i, file);
                    printf("last line read:\n'%s'\n", line_sav);
                    fmt_status = fmt_status | 1;
                    return -1;
                }
                break;
            }
            }
        }
        /*
         *  we are still left with the data to be read.  As it can have
         *  multiple formats and variable # values
         *
         *  set up the # of attribute values that need to be read.  This is
         *  the count for read flag!= 2 (a set value for each attribute) and 
         *  2 for read flag = 2 (min and max)
         */
        ndata = fmt->att[i].count;
        if (fmt->att[i].read == ATT_RD_READ_INCLUSIVE) ndata = 2;

        if (strcmp(str_typ, "DFNT_CHAR") == 0) {
            fmt->att[i].type = DFNT_CHAR;
            ndata = 1;
        } else if (strcmp(str_typ, "DFNT_INT8") == 0) {
            fmt->att[i].type = DFNT_INT8;
        } else if (strcmp(str_typ, "DFNT_UINT8") == 0) {
            fmt->att[i].type = DFNT_UINT8;
        } else if (strcmp(str_typ, "DFNT_INT16") == 0) {
            fmt->att[i].type = DFNT_INT16;
        } else if (strcmp(str_typ, "DFNT_INT32") == 0) {
            fmt->att[i].type = DFNT_INT32;
        } else if (strcmp(str_typ, "DFNT_FLOAT32") == 0) {
            fmt->att[i].type = DFNT_FLOAT32;
        } else if (strcmp(str_typ, "DFNT_FLOAT64") == 0) {
            fmt->att[i].type = DFNT_FLOAT64;
        } else {
            fmt->att[i].type = -9999; /* unhandled, checked later on in var_decode */
        }
        /*
         *  in the case where no read or no check is done, no more arguments
         *  reqquired on a line (ndata can be 0 for these cases)
         */
        if ((fmt->att[i].read != ATT_RD_NOREAD) &
                (fmt->att[i].read != ATT_RD_READ_NOCK)) {
            for (j = 0; j < ndata; j++) {
                if ((str = s_parse(NULL, '\"')) == NULL) {
                    printf("**************Program error\n");
                    printf("unable to find data #%d in desc # %d of file:\n'%s'\n",
                            j, i, file);
                    printf("last line read:\n'%s'\n", line_sav);
                    fmt_status = fmt_status | 1;
                    return -1;
                }
                if ((ierr = var_decode(str, fmt->att[i].type,
                        (void *) fmt->att[i].data.i32, j, 0)) != 0) {
                    if (ierr == -1) {
                        printf("**************Program error\n");
                        printf(
                                "unable to decode #%d in attribute description # %d of file:\n'%s'\n",
                                j, i, file);
                        printf("last line read:\n'%s'\n", line_sav);
                        fmt_status = fmt_status | 1;
                        return -1;
                    }
                    if (ierr == -2) {
                        printf("**************Program error\n");
                        printf(
                                "Value # %d is out of range in attr descrip # %d of file:\n'%s'\n",
                                j, i, file);
                        printf("last line read:\n'%s'\n", line_sav);
                        fmt_status = fmt_status | 1;
                        return -1;
                    }
                    if (ierr == -3) {
                        printf("**************Program error\n");
                        printf(
                                "data type of: %d is not handled by fmt_read currently\n",
                                fmt->att[i].type);
                        printf(
                                "In attr descrip # %d of file:\n'%s'\n", i, file);
                        fmt_status = fmt_status | 1;
                        return -1;
                    }
                }
            }
        }
    }
    return 0;
}
