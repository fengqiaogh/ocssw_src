#include "viirs_sim_sdr.h"
#include <stdio.h>
#include <string.h>

int wr_attr_seq(h5io_str *fid, int n_attr, h5attr_struc *attrs)
/*-----------------------------------------------------------------------------
    Program:   mk_geo_sdr.c

    Description:  create the geolocation file

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        h5io_str *  fid         I    dataset or group id to place attributes in
        int         n_attr      I    number of attributes to wtite
        h5attr_struc *  attrs   I    Atribute description structure

    Modification history:

    W. Robinson, SAIC  20 Oct 2008  Original development

----------------------------------------------------------------------------*/ {
    int i;
    /*
     *  Loop through the attribs and output
     */
    for (i = 0; i < n_attr; i++) {
        /* only output if the express value is on */
        if (attrs[i].express == 1) {
            if (attrs[i].type == 0) {
                if (h5io_wr_attr(fid, attrs[i].name, attrs[i].typ,
                        attrs[i].ndim, attrs[i].dim_siz, (void *) attrs[i].data)
                        != 0) {
                    printf("%s %d: Trouble writing non-string attr # %d\n",
                            __FILE__, __LINE__, i);
                    return 1;
                }
            } else {
                if (h5io_wr_attr_str(fid, attrs[i].name, attrs[i].ndim,
                        attrs[i].dim_siz, attrs[i].str_len, (void *) attrs[i].data)
                        != 0) {
                    printf("%s %d: Trouble writing string attr # %d\n",
                            __FILE__, __LINE__, i);
                    return 1;
                }
            }
        }
    }
    return 0;
}
