#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmt_check.h"
extern int fmt_status; /* format check status, see fmt_check */

int fmt_rd_ras(char *file, FILE *fid, fmt_str *fmt)
/*******************************************************************

   fmt_rd_ras

   purpose: read in the raster definition section of the format 
            description for a dataset

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
    int i, k, ifound;

    /*
     *  allocate space for theraster definitions 
     */
    if ((fmt->ras = malloc(fmt->n_raster * sizeof ( ras_str)))
            == NULL) {
        printf("**************Program error\n");
        printf("in: %s, raster definition storage allocation failed\n", __FILE__);
        fmt_status = fmt_status | 1;
        return -1;
    }
    /*
     *  read in each line of raster definition to the proper place
     */
    for (i = 0; i < fmt->n_raster; i++) {
        if (get_line(line, 500, fid, '#') == NULL) {
            printf("**************Program error\n");
            printf("unable to read raster definition # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line);
            fmt_status = fmt_status | 1;
            return -1;
        }
        strcpy(line_sav, line); /* save original line */

        /*
         *  raster label name
         */
        if ((str = s_parse(line, '\"')) == NULL) {
            printf("**************Program error\n");
            printf("unable to read raster label name, line %d, entry #1\n", i);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        strcpy(fmt->ras[i].lbl_ras, str);

        /*
         *  The next field is the palette label
         */
        if ((str = s_parse(NULL, '\"')) == NULL) {
            printf("**************Program error\n");
            printf(
                    "unable to read raster palette label, line # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        strcpy(fmt->ras[i].lbl_pal, str);

        /*
         *  the next 2 fields are the dimensions of the image, they are linked
         *  to the attributes by the short name
         */
        if ((str = s_parse(NULL, '\"')) == NULL) {
            printf("**************Program error\n");
            printf(
                    "unable to read raster # pixels short name, line # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        ifound = 0;
        for (k = 0; k < fmt->n_dim_defs; k++) {
            if (strcmp(str, fmt->dim_id[k].att_short) == 0) {
                ifound = 1;
                fmt->ras[i].npix_indx = k;
                break;
            }
        }
        if (ifound == 0) {
            printf("**************Program error\n");
            printf(
                    "No definition of %s was found in the array definition section\n",
                    str);
            printf("For raster pixel dimension of line %d of file:"
                    "\n'%s'\n", i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }

        if ((str = s_parse(NULL, '\"')) == NULL) {
            printf("**************Program error\n");
            printf(
                    "unable to read raster # lines short name, line # %d of file:\n'%s'\n", i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
        ifound = 0;
        for (k = 0; k < fmt->n_dim_defs; k++) {
            if (strcmp(str, fmt->dim_id[k].att_short) == 0) {
                ifound = 1;
                fmt->ras[i].nlin_indx = k;
                break;
            }
        }
        if (ifound == 0) {
            printf("**************Program error\n");
            printf(
                    "No definition of %s was found in the array definition section\n",
                    str);
            printf("For raster line dimension of line %d of file:"
                    "\n'%s'\n", i, file);
            printf("last line read:\n'%s'\n", line_sav);
            fmt_status = fmt_status | 1;
            return -1;
        }
    }
    return 0;
}
