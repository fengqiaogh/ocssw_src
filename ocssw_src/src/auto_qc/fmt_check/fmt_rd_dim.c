#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmt_check.h"
extern int fmt_status; /* format check status, see fmt_check */

int fmt_rd_dim(char *file, FILE *fid, fmt_str *fmt)
/*******************************************************************

   fmt_rd_dim

   purpose: read in the array dimension portion of the format description for 
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

 *******************************************************************/
 {
    char line[500], line_sav[500], *str;
    int i, j, ifound;

    /*
     *  allocate space for the dimension descriptions
     */
    if ((fmt->dim_id = malloc(fmt->n_dim_defs * sizeof ( attr_str))) == NULL) {
        printf("**************Program error\n");
        printf("in: %s, array definition structure allocation failed\n",
                __FILE__);
        fmt_status = fmt_status | 1;
        return -1;
    }
    /*
     *  read in each line of attribute description to the proper place
     */
    for (i = 0; i < fmt->n_dim_defs; i++) {
        if (get_line(line, 500, fid, '#') == NULL) {
            printf("**************Program error\n");
            printf(
                    "unable to read array dimension description # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line);
            fmt_status = fmt_status | 1;
            return -1;
        }

        strcpy(line_sav, line); /* keep an in-tact copy of the line */

        /*
         *  use the s_parse so that any quoted strings get taken together
         *  1st, the attrib name, then, the short name
         */
        if ((str = s_parse(line, '\"')) == NULL) {
            printf("**************Program error\n");
            printf("Unable to get <attribute name> (col 1) portion of dimension "
                    "description\n");
            printf("line %d, fmt table: %s\n", i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        strcpy(fmt->dim_id[i].att_nm, str);

        if ((str = s_parse(NULL, '\"')) == NULL) {
            printf("**************Program error\n");
            printf("Unable to get <short name> (col 2) portion of dimension "
                    "description\n");
            printf("line %d, fmt table: %s\n", i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        strcpy(fmt->dim_id[i].att_short, str);

        /*
         *  make sure the attribute name is one of the global attributes
         */
        ifound = 0;
        for (j = 0; j < fmt->n_attr; j++) {
            if (strcmp(fmt->att[j].obj_nm, "gbl") == 0) {
                if (strcmp(fmt->att[j].int_nm, fmt->dim_id[i].att_nm) == 0) {
                    ifound = 1;
                    /*
                     *  set the dim_index in the sttr struct to the location in here
                     */
                    fmt->att[j].dim_index = i;
                    /*
                     *  also make sure the attrib is not char, read in and a single
                     *  value
                     */
                    if (fmt->att[j].type == DFNT_CHAR) {
                        printf("**************Program error\n");
                        printf("dimension descriptor long name: '%s'\n",
                                fmt->dim_id[i].att_nm);
                        printf("cannot be a char type\n");
                        printf("line %d, fmt table: %s\n", i, file);
                        printf("last line read:\n'%s'\n", line_sav);
                        fmt_status = fmt_status | 1;
                        return -1;
                    }
                    if (fmt->att[j].read == ATT_RD_NOREAD) {
                        printf("**************Program error\n");
                        printf("dimension descriptor long name: '%s'\n",
                                fmt->dim_id[i].att_nm);
                        printf("must be read in in attribute section\n");
                        printf("(cannot have read keyword of NOREAD)\n");
                        printf("line %d, fmt table: %s\n", i, file);
                        printf("last line read:\n'%s'\n", line_sav);
                        fmt_status = fmt_status | 1;
                        return -1;
                    }
                    if (fmt->att[j].count != 1) {
                        printf("**************Program error\n");
                        printf("dimension descriptor long name: '%s'\n",
                                fmt->dim_id[i].att_nm);
                        printf("must be for attribute with a single value\n");
                        printf("(number of values must be 1)\n");
                        printf("line %d, fmt table: %s\n", i, file);
                        printf("last line read:\n'%s'\n", line_sav);
                        fmt_status = fmt_status | 1;
                        return -1;
                    }
                    break;
                }
            }
        }
        if (ifound != 1) {
            printf("**************Program error\n");
            printf("No attribute named: '%s' was found\n", fmt->dim_id[i].att_nm);
            printf("in the global attributes, from line %d, file %s\n", i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        /*
         *  fill in the dimension length with -1 to show it is unfilled as of yet
         */
        fmt->dim_id[i].dim_size = -1;
    }
    return 0;
}
