#include "l1io.h"
#include <mfhdf.h>

extern int fmt_status; /*  format checker status, see fmt_check */

int hdf_ras_chk(char *file, char *lbl_ras, char *lbl_pal, int npix, int nlin)
/*******************************************************************

   hdf_ras_chk

   purpose: check the raster and associated palette for a file

   Returns type: int for possible future use

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      hdf file to use
      char *            lbl_ras          I      expected raster image
                                                label, if it says '*NONE*',
                                                there is no raster
      char *            lbl_pal          I      expected palette label
      int               npix             I      expected # of pixels in
                                                image
      int               nlin             I      expected # of lines in image

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       15-Mar-1995     Original development
      W. Robinson       2-May-1995      adapt for no raster but a palettel
      W. Robinson       30-Jun-1995     set fmt_status if problems
      W. Robinson       13-jun-1996     DFANlablist expects char* even tho
                                        it should be char** so cast it so
      W. Robinson       4-aug-1997      disable the # raster images check
                                        (hdfpack makes this 2 for some reason)

 *******************************************************************/
 {
    intn nimg, label_len = 300, pal_flg;
    int nlbl, list_len = 300, ipix, ilin;
    uint16 ref_list[300];
    char label_list[300][300];
    /*
     *  get the labels for the raster and palette and check them
     *  we only expect 1 image and palette
     */
    if (strcmp(lbl_ras, "*NONE*") != 0) {
        nlbl = DFANlablist(file, DFTAG_RIG, ref_list, (char *) label_list,
                list_len, label_len, 1);
        if (nlbl != 1) {
            printf("****dataset has improper # of raster images = %d\n",
                    nlbl);
            fmt_status = fmt_status | 2;
        }

        /*
         *  check label
         */
        if (strcmp(label_list[0], lbl_ras) != 0) {
            printf(
                    "****Raster image has improper label\nread    : '%s'\nexpected: '%s'\n",
                    label_list[0], lbl_ras);
            fmt_status = fmt_status | 2;
        }
    }

    nlbl = DFANlablist(file, DFTAG_LUT, ref_list, (char *) label_list,
            list_len, label_len, 1);
    if (nlbl != 1) {
        printf("****dataset has improper # of palettes = %d\n",
                nlbl);
        fmt_status = fmt_status | 2;
    }

    /*
     *  check label
     */
    if (strcmp(label_list[0], lbl_pal) != 0) {
        printf(
                "****palette has improper label\nread    : '%s'\nexpected: '%s'\n",
                label_list[0], lbl_pal);
        printf("However, this is alright (There should be no label)\n");
    }
    /*
     *  Ok, use DFR8nimages to see if any images are there
     *      use DFR8getdims to get the image info
     *      use DFR8getimage to finally get the image
     */
    if (strcmp(lbl_ras, "*NONE*") != 0) {
        nimg = DFR8nimages(file);
        if (nimg != 1) {
            printf("****dataset has improper # of raster images = %d\n",
                    nimg);
            /*
             *  WDR 4aug97 the hdf pack makes this 2 so disable it for now
             fmt_status = fmt_status | 2;
             return -1;
             */
            printf("However, this will be ignored while hdfpack is checked\n");
        }

        if (DFR8getdims(file, (int32 *) & ipix, (int32 *) & ilin, &pal_flg) == -1) {
            printf("****For raster image, DFR8getdims got an error\n");
            fmt_status = fmt_status | 2;
            return -1;
        }
        if (!pal_flg) {
            printf("****DFR8getdims found no palette\n");
            fmt_status = fmt_status | 2;
        }

        if (ipix != npix) {
            printf("****raster image has %d pixels, expected %d pixels\n",
                    ipix, npix);
            fmt_status = fmt_status | 2;
        }

        if (ilin != nlin) {
            printf("****raster image has %d lines, expected %d lines\n",
                    ilin, nlin);
            fmt_status = fmt_status | 2;
        }
    }
    return 0;
}
