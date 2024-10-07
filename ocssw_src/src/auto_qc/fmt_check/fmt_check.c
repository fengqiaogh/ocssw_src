#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mfhdf.h"
#include "fmt_check.h"
int verbose; /*  switch to write the global attrib info: 0 no, 1 yes */
int fmt_status; /* status of format checking: 0 all good, 1 program problem,
                   2 format problem, 3 both problems  */

int main(int argc, char **argv)
/*******************************************************************

   fmt_check

   purpose: open the SeaWiFS dataset and check the format

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               argc            I       count of input args
      char *            argv[]          I       args: [1] hdf file to check
                                                [2] table to check against
                                                [3] verbose switch - 0 or
                                                blank, none, 1 - be more 
                                                verbose

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       16-Feb-1995     Original development
      W. Robinson       9-Mar-1995      update to v2.6
                                        1 - Title is now a constant
                                        2 - new fields added: Scene Center 
                                            Solar Zenith, csol_z, entry_year, 
                                            entry_day
      W. Robinson       16-Mar-1995     run the checking of all levels
                                        and all versions
      W. Robinson       18-Apr-1995     add the level 2 and 3
      W. Robinson       30-Jun-1995     add verbose switch to call args
      W. Robinson       29-Aug-1995     update for added v2.7 formats
      W. Robinson       20-Nov-1995     update for added v2.8 formats
      L. Kumar		19-Dec-1995     add nrt met and ozone
      W. Robinson       21-Oct-1996     update for v3.0 formats
      W. Robinson        3-Mar-1997     update for v3.1 met format
      W. Robinson        9-Dec-1999     update for new repro 3 formats
      W. Robinson       15-Mar-2001     update argv for linux use
      W. Robinson, SAIC 29 Mar 2005     update to fully table driven model
                                        and call I/O for sections from here

 *******************************************************************/
 {
    int i, j, ierr, num, a_count, index;
    l1info_struct l1info;
    void fmt_exit(int);
    char file[400], table[400];
    fmt_str fmt; /* master format structure */
    u_data value;

    fmt_status = 0;
    /*
     *  Start check the input line arguments
     */
    if (argc != 3 && argc != 4) {
        printf("*** format: fmt_check <file> <table> <verbose>\n");
        printf("    <file> is the hdf file to check\n");
        printf("    <table> is the table to use in the format checking\n");
        printf("    <verbose> is blank or 0 to list attributes 1 to omit\n");
        fmt_status = 1;
        fmt_exit(fmt_status);
    }
    /*
     *  set global to activate attr_disp if 4th arg is 1 (verbose output)
     */
    verbose = (argc == 4 && argv[3][0] == '1') ? 0 : 1;
    strcpy(file, argv[1]);
    strcpy(table, argv[2]);

    printf("Format Check of Dataset: '%s'\n", file);
    printf("with control table: \n  '%s'\n\n", table);
    /*
     *  hopen will open the file and read data descriptor blocks
     *  to memory
     */
    if ((l1info.fid = Hopen(file, DFACC_RDONLY, 0)) < 0) {
        printf("fmt_check: Failed on the Hopen of \n'%s'\n", file);
        fmt_status = 1;
        fmt_exit(fmt_status);
    }

    /*
     *  SDstart opens the hdf interface and inits the SD interface, wierd.
     *  apparently, if both SD and HDF are used, both these need to
     *  be called.
     */
    if ((l1info.sdfid = SDstart(file, DFACC_RDONLY)) < 0) {
        printf("fmt_check: Failure at SDstart of \n'%s'\n", file);
        l1io_close(l1info);
        fmt_status = 2;
        fmt_exit(fmt_status);
    }

    /*
     *  ok, it's open. Get the format descriptions fromthe table
     */
    if (fmt_read(table, &fmt) != 0) {
        printf("**** Exiting due to failure in control file read\n");
        printf("     control file: '%s'\n", table);
        fmt_status = fmt_status | 1;
        fmt_exit(fmt_status);
    }
    printf("\n\nChecking Global Attributes...\n\n");

    for (i = 0; i < fmt.n_attr; i++) {
        if (strcmp(fmt.att[i].obj_nm, "gbl") == 0) {
            ierr = get_attr(l1info.sdfid, fmt.att[i], &value, &a_count);
            if (ierr == 0) {
                /*
                 *  display the attribute
                 */
                attr_disp(fmt.att[i], value, a_count);

                /*
                 *  also, save the dimension sizes desired
                 */
                if ((index = fmt.att[i].dim_index) >= 0) {
                    switch (fmt.att[i].type) {
                    case DFNT_FLOAT32:
                        fmt.dim_id[ index ].dim_size = (int) value.f32[0];
                        break;
                    case DFNT_FLOAT64:
                        fmt.dim_id[ index ].dim_size = (int) value.f64[0];
                        break;
                    case DFNT_INT8:
                        fmt.dim_id[ index ].dim_size = (int) value.i8[0];
                        break;
                    case DFNT_UINT8:
                        fmt.dim_id[ index ].dim_size = (int) value.ui8[0];
                        break;
                    case DFNT_INT16:
                        fmt.dim_id[ index ].dim_size = (int) value.i16[0];
                        break;
                    case DFNT_INT32:
                        fmt.dim_id[ index ].dim_size = (int) value.i32[0];
                        break;
                    default:
                        printf(
                                "****Program error, dimension getting reached impossible condition,\n");
                        printf("routine %s\n", __FILE__);
                        fmt_status = fmt_status | 1;
                        fmt_exit(fmt_status);
                        break;
                    }
                }
            }
        }
    }
    /* temporary, list the dim_id stuff */
    printf("\n\ndimension structure info:\n");
    printf("attrib name, short name, dim size, # entries: %d\n", fmt.n_dim_defs);
    for (i = 0; i < fmt.n_dim_defs; i++) {
        printf("%s  %s   %d\n", fmt.dim_id[i].att_nm, fmt.dim_id[i].att_short, fmt.dim_id[i].dim_size);
    }
    /*
     *  next, proceed to check the science datasets
     */
    printf("\n\nChecking Science datasets...\n\n");

    /*
     *  loop through the sds's
     */
    for (i = 0; i < fmt.n_sds; i++) {
        /*
         *  Find the indixes of the attributes for this sds
         */
        num = 0;
        for (j = 0; j < fmt.n_attr; j++) {
            if (strcmp(fmt.att[j].obj_nm, fmt.sds_info[i].name) == 0) {
                fmt.sds_info[i].attr_ind[ num++ ] = j;
            }
        }
        if (num != fmt.sds_info[i].n_attr) {
            printf("********************\n");
            printf("Problem with fmt_check, for SDS: '%s',\n",
                    fmt.sds_info[i].name);
            printf(" mismatch in # attributes in list and expected #\n");
            fmt_status = fmt_status | 2;
        }
        /*
         *  Do most of the checking of the sds size and its attributes
         */
        if ((ierr = chk_sds(l1info.sdfid, &fmt, i)) == 0) {
            /*
             *  code for checking the ranges that could not be done above
             */
        }
    }

    /*
     *  check the raster data if requested
     */
    if (fmt.n_raster > 0) {
        printf("\n\nChecking Raster data...\n\n");
        for (i = 0; i < fmt.n_raster; i++) {
            hdf_ras_chk(file, fmt.ras[i].lbl_ras, fmt.ras[i].lbl_pal,
                    fmt.dim_id[ fmt.ras[i].npix_indx ].dim_size,
                    fmt.dim_id[ fmt.ras[i].nlin_indx ].dim_size);
        }
    }

    /*
     *  check the L3 Vgroups if needed
     */
    if (fmt.n_vgroup > 0) {
        ck_v_l3(l1info.fid, &fmt);
    }

    l1io_close(l1info);
    fmt_exit(fmt_status);
    return 0;
}

void fmt_exit(int status)
/*******************************************************************

   fmt_exit

   purpose: provide a common exit from the fmt_check program and
            error summary info.

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               status          I       format check status

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       30-jun-1995     Original development

 *******************************************************************/ {
    if (status == 0)
        printf("\n\nSuccessfull format check, no format errors\n");
    else {
        if (status & 1)
            printf("\n\nFailure of format check due to program error\n");
        if (status & 2)
            printf("\n\nFailure of format check due to dataset format error\n");
    }
    exit(status);
}
