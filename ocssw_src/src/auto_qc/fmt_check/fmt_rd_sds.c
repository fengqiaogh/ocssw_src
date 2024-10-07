#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmt_check.h"
extern int fmt_status; /* format check status, see fmt_check */

int fmt_rd_sds(char *file, FILE *fid, fmt_str *fmt)
/*******************************************************************

   fmt_rd_sds

   purpose: read in the SDS description section of the format description for 
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
    char line[500], line_sav[500], str_typ[200], *str, *str2;
    int i, j, k, idim, ifound, ierr;

    /*
     *  allocate space for the attribute descriptions
     */
    if ((fmt->sds_info = malloc(fmt->n_sds * sizeof ( sds_info_str)))
            == NULL) {
        printf("**************Program error\n");
        printf("in: %s, SDS description storage allocation failed\n", __FILE__);
        fmt_status = fmt_status | 1;
        return -1;
    }
    /*
     *  read in each line of SDS description to the proper place
     */
    for (i = 0; i < fmt->n_sds; i++) {
        if (get_line(line, 500, fid, '#') == NULL) {
            printf("**************Program error\n");
            printf("unable to read SDS description # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line);
            fmt_status = fmt_status | 1;
            return -1;
        }
        strcpy(line_sav, line); /* save original line */

        /*
         *  SDS name field
         */
        if ((str = s_parse(line, '\"')) == NULL) {
            printf("**************Program error\n");
            printf("unable to read SDS ctl line %d, entry #1\n", i);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        strcpy(fmt->sds_info[i].name, str);

        /*
         *  The next 3 fields: data type, # attributes and rank
         */
        for (j = 0; j < 3; j++) {
            if ((str = s_parse(NULL, '\"')) == NULL) {
                printf("**************Program error\n");
                printf(
                        "unable to read SDS ctl value %d of line # %d of file:\n'%s'\n",
                        j + 2, i, file);
                printf("last line read:\n'%s'\n", line_sav);
                fmt_status = fmt_status | 1;
                return -1;
            }
            switch (j) {
            case 0: /* data type  */
                strcpy(str_typ, str);
                break;

            case 1: /* # attributes  */
                if (sscanf(str, "%d", &(fmt->sds_info[i].n_attr)) == EOF) {
                    printf("**************Program error\n");
                    printf(
                            "Cannot decode SDS n_attr ctl, line %d, file:\n'%s'\n",
                            i, file);
                    printf("last line read:\n'%s'\n", line_sav);
                    fmt_status = fmt_status | 1;
                    return -1;
                }
                break;

            case 2: /*  rank (# dimensions)  */
                if (sscanf(str, "%d", &(fmt->sds_info[i].rank)) == EOF) {
                    printf("**************Program error\n");
                    printf("Cannot decode SDS rank ctl, line %d, file:\n'%s'\n",
                            i, file);
                    printf("last line read:\n'%s'\n", line_sav);
                    fmt_status = fmt_status | 1;
                    return -1;
                }
                break;
            }
        }
        /*
         *  depending on the rank, there are up to 3 pairs of dimension
         *  type / size fields
         */
        for (idim = 0; idim < fmt->sds_info[i].rank; idim++) {
            if ((str = s_parse(NULL, '\"')) == NULL) {
                printf("**************Program error\n");
                printf(
                        "unable to read SDS field %d of line # %d of file:\n'%s'\n",
                        idim * 2 + 5, i, file);
                printf("last line read:\n'%s'\n", line_sav);
                fmt_status = fmt_status | 1;
                return -1;
            }
            /*
             *  read the next value and have it ready
             */
            if ((str2 = s_parse(NULL, '\"')) == NULL) {
                printf("**************Program error\n");
                printf(
                        "unable to read SDS field %d of line # %d of file:\n'%s'\n",
                        idim * 2 + 6, i, file);
                printf("last line read:\n'%s'\n", line_sav);
                fmt_status = fmt_status | 1;
                return -1;
            }
            /*
             *  check the dimension type and deal with size field appropriately
             */
            if (strcmp(str, "DSIZE_NOCK") == 0)
                fmt->sds_info[i].e_ranges[idim] = 0;

            else if (strcmp(str, "DSIZE_CK") == 0) {
                if (sscanf(str2, "%d", &(fmt->sds_info[i].e_ranges[idim]))
                        == EOF) {
                    printf("**************Program error\n");
                    printf("unable to get dimension size from format table\n");
                    printf("for dimension %d, SDS field %d, of line %d of file:"
                            "\n'%s'\n", idim, idim * 2 + 6, i, file);
                    printf("last line read:\n'%s'\n", line_sav);
                    fmt_status = fmt_status | 1;
                    return -1;
                }
            } else if (strcmp(str, "DSIZE_DIM_NM") == 0) {
                /*
                 *  the index to the dimension structure array (+1 and negated)
                 *  is placed in the e_range
                 */
                ifound = 0;
                for (k = 0; k < fmt->n_dim_defs; k++) {
                    if (strcmp(str2, fmt->dim_id[k].att_short) == 0) {
                        ifound = 1;
                        fmt->sds_info[i].e_ranges[idim] = -(k + 1);
                        break;
                    }
                }
                if (ifound == 0) {
                    printf("**************Program error\n");
                    printf(
                            "No definition of %s was found in the array definition section\n",
                            str2);
                    printf("For dimension %d, SDS field %d, of line %d of file:"
                            "\n'%s'\n", idim, idim * 2 + 6, i, file);
                    printf("last line read:\n'%s'\n", line_sav);
                    fmt_status = fmt_status | 1;
                    return -1;
                }
            }
        }
        /*
         *  get the type put in
         */
        if (strcmp(str_typ, "DFNT_CHAR") == 0) {
            fmt->sds_info[i].type = DFNT_CHAR;
            fmt->sds_info[i].byt_per_val = 1;
        } else if (strcmp(str_typ, "DFNT_INT16") == 0) {
            fmt->sds_info[i].type = DFNT_INT16;
            fmt->sds_info[i].byt_per_val = 2;
        } else if (strcmp(str_typ, "DFNT_INT32") == 0) {
            fmt->sds_info[i].type = DFNT_INT32;
            fmt->sds_info[i].byt_per_val = 4;
        } else if (strcmp(str_typ, "DFNT_FLOAT32") == 0) {
            fmt->sds_info[i].type = DFNT_FLOAT32;
            fmt->sds_info[i].byt_per_val = 4;
        } else if (strcmp(str_typ, "DFNT_FLOAT64") == 0) {
            fmt->sds_info[i].type = DFNT_FLOAT64;
            fmt->sds_info[i].byt_per_val = 4;
        } else if (strcmp(str_typ, "DFNT_INT8") == 0) {
            fmt->sds_info[i].type = DFNT_INT8;
            fmt->sds_info[i].byt_per_val = 1;
        } else if (strcmp(str_typ, "DFNT_UINT8") == 0) {
            fmt->sds_info[i].type = DFNT_UINT8;
            fmt->sds_info[i].byt_per_val = 1;
        } else {
            printf("**************Program error\n");
            printf(
                    "Unplanned SDS format descript encountered for line %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        /*
         *  for optional range checking, get the range checking keyword
         *  to decide if any needs to be done
         */
        if ((str = s_parse(NULL, '\"')) == NULL) {
            printf("**************Program error\n");
            printf("range checking keyword must be provided\n");
            printf("at field %d, line %d of file:\n'%s'\n",
                    7 + 2 * fmt->sds_info[i].rank, i, file);
            fmt_status = fmt_status | 1;
            return -1;
        } else {
            if (strcmp(str, "NO_MN_MX") == 0) {
                fmt->sds_info[i].flg_rng = 0;
            } else if (strcmp(str, "MN_MX_CK") == 0) {
                fmt->sds_info[i].flg_rng = 1;

                for (j = 0; j < 2; j++) {
                    if ((str = s_parse(NULL, '\"')) == NULL) {
                        printf("**************Program error\n");
                        printf(
                                "unable to find SDS range value #%d in desc # %d of file:\n'%s'\n",
                                j, i, file);
                        printf("last line read:\n'%s'\n", line_sav);
                        fmt_status = fmt_status | 1;
                        return -1;
                    }
                    if ((ierr = var_decode(str, fmt->sds_info[i].type,
                            (void *) fmt->sds_info[i].sds_rng.i32, j, 1)) != 0) {
                        if (ierr == -1) {
                            printf("**************Program error\n");
                            printf("unable to decode SDS range value #%d"
                                    " in description # %d of file: n'%s'\n", j, i, file);
                            printf("last line read:\n'%s'\n", line_sav);
                            fmt_status = fmt_status | 1;
                            return -1;
                        }
                        if (ierr == -2) {
                            printf("**************Program error\n");
                            printf("SDS range value # %d is out of range in descrip #"
                                    " %d of file:\n'%s'\n", j, i, file);
                            printf("last line read:\n'%s'\n", line_sav);
                            fmt_status = fmt_status | 1;
                            return -1;
                        }
                        if (ierr == -3) {
                            printf("**************Program error\n");
                            printf("data type of: %d is not handled by fmt_read currently\n",
                                    fmt->att[i].type);
                            printf("In descrip # %d of file:\n'%s'\n", i, file);
                            fmt_status = fmt_status | 1;
                            return -1;
                        }
                    }
                }
            } else {
                printf("**************Program error\n");
                printf(
                        "Range checking keyword must be either 'NO_MN_MX' or 'MN_MX_CK'\n");
                printf("at field %d, line %d of file:\n'%s'\n",
                        7 + 2 * fmt->sds_info[i].rank, i, file);
                fmt_status = fmt_status | 1;
                return -1;
            }
        }
    }
    return 0;
}
