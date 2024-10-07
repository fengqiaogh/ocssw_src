#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmt_check.h"
extern int fmt_status; /* format check status, see fmt_check */

int fmt_read(char *file, fmt_str *fmt)
/*******************************************************************

   fmt_read

   purpose: read in the format description for a dataset
            from a file

   Returns type: int - 0 if all went well, 

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      file name containing the format
                                                table description
      fmt_str *         fmt              O      format structure

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       21-Feb-1995     Original development
      W. Robinson       13-Mar-1995     adapt to use the s_parse routine
      W. Robinson        8-May-1995     take fmt files from $SWFTBL/fmt_check
      W. Robinson        7-Jun-1995     add code to handle UINT8 - unsigned 
                                        8-bit bytes
      W. Robinson       30-Jun-1995     set fmt_status if problems
      W. Robinson       26-Sep-1996     add logic for new attribute read 
                                        flag of 2: read in the minimum and 
                                        maximum.  Also, read in min, max
                                        range for sds range checking
      W. Robinson       31-Oct-1996     upgrade to handle float64 data
      W. Robinson       18-Feb-1997     changed call to var_decode
      W. Robinson, SAIC 29 Mar 2005     update to fully table driven model
                                        and call I/O for sections from here

 *******************************************************************/
 {
    FILE *fid;
    char line[500], *str;
    int ival, isec;
    char *sec_ids[] = {"**ATTR", "**ARR_DIM", "**SDS", "**RASTER", "**VGROUP"};

    enum sec_code {
        SEC_CODE_ATTR, SEC_CODE_ARR_DIM, SEC_CODE_SDS,
        SEC_CODE_RASTER, SEC_CODE_VGROUP
    };
    int nsec = 5; /* # of sections in fmt table and identifiers (above) */

    /*
     *  Open the file
     */
    if ((fid = fopen(file, "r")) == NULL) {
        printf("**************Program error\n");
        printf("unable to open format file:\n'%s'\n", file);
        fmt_status = fmt_status | 1;
        return -1;
    }

    /*
     *  process each section description, each must be there in order
     */
    for (isec = 0; isec < nsec; isec++) {
        if (get_line(line, 500, fid, '#') == (char *) NULL) {
            printf("**************Format table read error\n");
            printf("table: %s\n", file);
            printf("for section # %d, name: %s\n", isec, sec_ids[isec]);
            printf("unable to read section description line\n");
            fmt_status = fmt_status | 1;
            return -1;
        }
        /*
         * extract the section keyword and # following entries
         */
        if ((str = s_parse(line, '\"')) == NULL) {
            printf("**************Format table read error\n");
            printf("table: %s\n", file);
            printf("for section # %d, name: %s\n", isec, sec_ids[isec]);
            printf("unable to find section keyword\n");
            fmt_status = fmt_status | 1;
            return -1;
        }
        if (strcmp(str, sec_ids[isec]) != 0) {
            printf("**************Format table read error\n");
            printf("table: %s\n", file);
            printf("for section # %d, name: %s\n", isec, sec_ids[isec]);
            printf("improper section name: %s, found\n", str);
            fmt_status = fmt_status | 1;
            return -1;
        }
        if ((str = s_parse(NULL, '\"')) == NULL) {
            printf("**************Format table read error\n");
            printf("table: %s\n", file);
            printf("for section # %d, name: %s\n", isec, sec_ids[isec]);
            printf("no # of section descriptions found\n");
            fmt_status = fmt_status | 1;
            return -1;
        }
        if (sscanf(str, "%d", &ival) == EOF) {
            printf("**************Format table read error\n");
            printf("table: %s\n", file);
            printf("for section # %d, name: %s\n", isec, sec_ids[isec]);
            printf("could not read # description lines\n");
            fmt_status = fmt_status | 1;
            return -1;
        }
        /*
         *  for each section, call a reader
         */
        switch (isec) {
        case SEC_CODE_ATTR:
            printf("%s - fmt_rd call for %s\n", __FILE__, sec_ids[isec]);
            fmt->n_attr = ival;
            if (fmt_rd_attr(file, fid, fmt) != 0) return -1;
            break;

        case SEC_CODE_ARR_DIM:
            printf("%s - fmt_rd call for %s\n", __FILE__, sec_ids[isec]);
            fmt->n_dim_defs = ival;
            if (fmt_rd_dim(file, fid, fmt) != 0) return -1;
            break;

        case SEC_CODE_SDS:
            printf("%s - fmt_rd call for %s\n", __FILE__, sec_ids[isec]);
            fmt->n_sds = ival;
            if (fmt_rd_sds(file, fid, fmt) != 0) return -1;
            break;

        case SEC_CODE_RASTER:
            printf("%s - fmt_rd call for %s\n", __FILE__, sec_ids[isec]);
            fmt->n_raster = ival;
            if (fmt_rd_ras(file, fid, fmt) != 0) return -1;
            break;

        case SEC_CODE_VGROUP:
            printf("%s - fmt_rd call for %s\n", __FILE__, sec_ids[isec]);
            fmt->n_vgroup = ival;
            if (ival > 0)
                if (fmt_rd_l3vg(file, fid, fmt) != 0) return -1;
            break;

        default:
            printf("************** Program error\n");
            printf("%s - This shouldn't happen, isec = %d\n", __FILE__, isec);
            fmt_status = fmt_status | 1;
            return -1;
            break;
        }
    }
    /*
     *  some verification: array dimension section info
     */
    /*
      printf( "list of the array dimension descriptors, # = %d\n", 
        fmt->n_dim_defs );
      printf( "'attribute name'  'short name'\n" );
      for( i = 0; i < fmt->n_dim_defs; i++ )
        {
        printf( "%s    %s\n", fmt->dim_id[i].att_nm, fmt->dim_id[i].att_short );
        }
     */
    /*
     *  SDS description section info
     */
    /*
      printf( "\n\nList of SDS section info, # entries: %d\n", fmt->n_sds );
      printf( "SDS name  type  # attr  rank  dim 1 code, size   dim 2  dim 3  flg_rng\n" );
      for( i = 0; i < fmt->n_sds; i++ )
        {
        printf( "%s  %d  %d  %d    %d    %d    %d    %d\n",
          fmt->sds_info[i].name, fmt->sds_info[i].type, fmt->sds_info[i].n_attr,
          fmt->sds_info[i].rank, fmt->sds_info[i].e_ranges[0], 
          fmt->sds_info[i].e_ranges[1], fmt->sds_info[i].e_ranges[2], 
          fmt->sds_info[i].flg_rng );
        }
     */

    fclose(fid);
    return 0;
}
