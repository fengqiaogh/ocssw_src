#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmt_check.h"
extern int fmt_status; /* format check status, see fmt_check */

int fmt_rd_l3vg(char *file, FILE *fid, fmt_str *fmt)
/*******************************************************************

   fmt_rd_l3vg

   purpose: read in the level 3 Vgroup description section of the format 
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
    int i, j, k, ifound, n_tot_vg, nrows;
    int valid_resolve[5] = {1, 2, 4, 9, 36};
    int numrows[5] = {17280, 8640, 4320, 2160, 540};

    /*
     *  allocate space for the Vgroup descriptions
     */
    n_tot_vg = fmt->n_vgroup + 3;
    if ((fmt->vg_info = malloc(n_tot_vg * sizeof ( vg_info_str)))
            == NULL) {
        printf("**************Program error\n");
        printf("in: %s, Vgroup description storage allocation failed\n",
                __FILE__);
        fmt_status = fmt_status | 1;
        return -1;
    }
    /*
     *  read the # of bins short name
     */
    i = 1;
    if (get_line(line, 500, fid, '#') == NULL) {
        printf("**************Program error\n");
        printf("unable to read Vgroup description # %d of file:\n'%s'\n",
                i, file);
        printf("last line read:\n'%s'\n", line);
        fmt_status = fmt_status | 1;
        return -1;
    }
    strcpy(line_sav, line);

    if ((str = s_parse(line, '\"')) == NULL) {
        printf("**************Program error\n");
        printf(
                "unable to read Vgroup description, line # %d of file:\n'%s'\n",
                i, file);
        printf("last line read:\n'%s'\n", line);
        fmt_status = fmt_status | 1;
        return -1;
    }

    ifound = 0;
    for (k = 0; k < fmt->n_dim_defs; k++) {
        if (strcmp(str, fmt->dim_id[k].att_short) == 0) {
            ifound = 1;
            fmt->vg_nbin_indx = k;
            break;
        }
    }
    if (ifound == 0) {
        printf("**************Program error\n");
        printf(
                "No definition of %s was found in the array definition section\n",
                str);
        printf("For Vgroup 'nbins' definition, line %d of file:"
                "\n'%s'\n", i, file);
        printf("last line read:\n'%s'\n", line);
        fmt_status = fmt_status | 1;
        return -1;
    }

    /*
     *  get the resolve and check it is a valid value
     */
    if ((str = s_parse(NULL, '\"')) == NULL) {
        printf("**************Program error\n");
        printf(
                "unable to read Vgroup resolve value, line # %d of file:\n'%s'\n",
                i, file);
        printf("last line read:\n'%s'\n", line_sav);
        fmt_status = fmt_status | 1;
        return -1;
    }

    if (sscanf(str, "%d", &(fmt->resolve)) == EOF) {
        printf("**************Program error\n");
        printf(
                "Cannot decode Vgroup resolve value, line # %d of file:\n'%s'\n",
                i, file);
        printf("last line read:\n'%s'\n", line_sav);
        fmt_status = fmt_status | 1;
        return -1;
    }

    ifound = 0;
    for (i = 0; i < 5; i++) {
        if (fmt->resolve == valid_resolve[i]) {
            ifound = 1;
            nrows = numrows[i];
            break;
        }
    }
    if (ifound != 1) {
        printf("**************Program error\n");
        printf(
                "Vgroup resolve is not a valid value.  It must be\n"
                "1, 2, 4, 9, or 36\n");
        printf("last line read:\n'%s'\n", line_sav);
        fmt_status = fmt_status | 1;
        return -1;
    }

    /*
     *  set up the first 3 Vdatas with standard information
     */
    fmt->vg_info[0].n_fields = 7;
    fmt->vg_info[0].num_typs[0] = DFNT_INT32;
    fmt->vg_info[0].num_typs[1] = DFNT_INT32;
    fmt->vg_info[0].num_typs[2] = DFNT_INT32;
    fmt->vg_info[0].num_typs[3] = DFNT_FLOAT64;
    fmt->vg_info[0].num_typs[4] = DFNT_FLOAT64;
    fmt->vg_info[0].num_typs[5] = DFNT_FLOAT64;
    fmt->vg_info[0].num_typs[6] = DFNT_FLOAT64;
    fmt->vg_info[0].nrec = 1;
    strcpy(fmt->vg_info[0].name, "SEAGrid");
    strcpy(fmt->vg_info[0].class, "Geometry");
    strcpy(fmt->vg_info[0].field_list,
            "registration,straddle,bins,radius,max_north,max_south,seam_lon");

    fmt->vg_info[1].n_fields = 7;
    fmt->vg_info[1].num_typs[0] = DFNT_INT32;
    fmt->vg_info[1].num_typs[1] = DFNT_FLOAT64;
    fmt->vg_info[1].num_typs[2] = DFNT_FLOAT64;
    fmt->vg_info[1].num_typs[3] = DFNT_INT32;
    fmt->vg_info[1].num_typs[4] = DFNT_INT32;
    fmt->vg_info[1].num_typs[5] = DFNT_INT32;
    fmt->vg_info[1].num_typs[6] = DFNT_INT32;
    fmt->vg_info[1].nrec = nrows;
    strcpy(fmt->vg_info[1].name, "BinIndex");
    strcpy(fmt->vg_info[1].class, "Index");
    strcpy(fmt->vg_info[1].field_list,
            "row_num,vsize,hsize,start_num,begin,extent,max");

    fmt->vg_info[2].n_fields = 7;
    fmt->vg_info[2].num_typs[0] = DFNT_INT32;
    fmt->vg_info[2].num_typs[1] = DFNT_INT16;
    fmt->vg_info[2].num_typs[2] = DFNT_INT16;
    fmt->vg_info[2].num_typs[3] = DFNT_INT16;
    fmt->vg_info[2].num_typs[4] = DFNT_FLOAT32;
    fmt->vg_info[2].num_typs[5] = DFNT_INT8;
    fmt->vg_info[2].num_typs[6] = DFNT_INT32;
    fmt->vg_info[2].nrec = -1;
    strcpy(fmt->vg_info[2].name, "BinList");
    strcpy(fmt->vg_info[2].class, "DataMain");
    strcpy(fmt->vg_info[2].field_list,
            "bin_num,nobs,nscenes,time_rec,weights,sel_cat,flags_set");
    /*
     *  The following lines are the names of the products in the L3
     *  set up the structures for each of those products
     */
    for (i = 0; i < fmt->n_vgroup; i++) {
        j = i + 3;
        if (get_line(line, 500, fid, '#') == NULL) {
            printf("**************Program error\n");
            printf("unable to read Vgroup product description # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line);
            fmt_status = fmt_status | 1;
            return -1;
        }

        if ((str = s_parse(line, '\"')) == NULL) {
            printf("**************Program error\n");
            printf(
                    "unable to read Vgroup product description, line # %d of file:\n'%s'\n",
                    i, file);
            printf("last line read:\n'%s'\n", line);
            fmt_status = fmt_status | 1;
            return -1;
        }
        /*
         *  OK, set up the general and product-specific structure info
         */
        fmt->vg_info[j].n_fields = 2;
        fmt->vg_info[j].num_typs[0] = DFNT_FLOAT32;
        fmt->vg_info[j].num_typs[1] = DFNT_FLOAT32;
        fmt->vg_info[j].nrec = -1;
        strcpy(fmt->vg_info[j].name, str);
        strcpy(fmt->vg_info[j].class, "DataSubordinate");
        sprintf(fmt->vg_info[j].field_list, "%s_sum,%s_sum_sq", str, str);
    }
    return 0;
}
