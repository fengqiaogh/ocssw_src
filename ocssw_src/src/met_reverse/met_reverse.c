#include "hio.h"

int main(int argc, char *argv[])
/*******************************************************************

   met_reverse

   purpose: reverse the lines of a met file in the latitude direction - to fix
     a problem in ancnrt that happened starting Jul 2009

   Returns type: status of 0 is all good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               argc             I      count of arguments
      char *[]          argv             I      [1] is the met file to reverse 
                                                the fields of

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       18-Sep-2009     Original development, 

 *******************************************************************/ {
    hio_struct hinfo;
    float *data, *xfr_arr;
    int nlin = 181, npix = 360, ilin, iret = 0, irc, ifld, nfld = 5;
    char *b_data, *b_xfr;
    char *met_fields[] = {"z_wind", "m_wind", "press", "p_water", "rel_hum"};
    char *qc_fields[] = {"z_wind_QC", "m_wind_QC", "press_QC",
        "p_water_QC", "rel_hum_QC"};
    /*
     *  open the file
     */
    if (hio_open(argv[1], DFACC_RDWR, &hinfo) != 0) {
        printf("%s: hio_open failure\n", __FILE__);
        iret = 1;
        exit(iret);
    }
    /*
     *  do all fields, parm and QC
     */
    irc = 0;
    if ((xfr_arr = (float *) malloc(npix * sizeof ( float))) == NULL)
        irc++;
    if ((b_xfr = (char *) malloc(npix * sizeof ( char))) == NULL)
        irc++;
    if ((data = (float *) malloc(npix * nlin * sizeof ( float))) == NULL)
        irc++;
    if ((b_data = (char *) malloc(npix * nlin * sizeof ( char))) == NULL)
        irc++;
    if (irc != 0) {
        printf("%s: Transfer space allocation failed\n", __FILE__);
        iret = 1;
    } else {
        for (ifld = 0; ifld < nfld; ifld++) {
            /*
             *  read a SDS
             */
            if (hio_r_sds(hinfo, met_fields[ifld],
                    DFNT_FLOAT32, (void *) data) != 0) {
                printf("%s: hio_r_sds failure\n", __FILE__);
                iret = 1;
                exit(iret);
            }
            /*
             *  reverse the data
             */
            for (ilin = 0; ilin < nlin / 2; ilin++) {
                memcpy((void *) xfr_arr, (void *) (data + ilin * npix),
                        npix * sizeof ( float));
                memcpy((void *) (data + ilin * npix),
                        (void *) (data + (nlin - 1 - ilin) * npix),
                        npix * sizeof ( float));
                memcpy((void *) (data + (nlin - 1 - ilin) * npix),
                        (void *) xfr_arr, npix * sizeof ( float));
            }
            /*
             *  write it out again
             */
            if (hio_i_sds(hinfo, met_fields[ifld],
                    DFNT_FLOAT32, (void *) data) != 0) {
                printf("%s: hio_i_sds failure\n", __FILE__);
                iret = 1;
                exit(iret);
            }
            /*
             *  do the same for the QC byte arrays
             */
            if (hio_r_sds(hinfo, qc_fields[ifld],
                    DFNT_INT8, (void *) b_data) != 0) {
                printf("%s: hio_r_sds failure\n", __FILE__);
                iret = 1;
                exit(iret);
            }

            for (ilin = 0; ilin < nlin / 2; ilin++) {
                memcpy((void *) b_xfr, (void *) (b_data + ilin * npix),
                        npix * sizeof ( char));
                memcpy((void *) (b_data + ilin * npix),
                        (void *) (b_data + (nlin - 1 - ilin) * npix),
                        npix * sizeof ( char));
                memcpy((void *) (b_data + (nlin - 1 - ilin) * npix),
                        (void *) b_xfr, npix * sizeof ( char));
            }

            if (hio_i_sds(hinfo, qc_fields[ifld],
                    DFNT_INT8, (void *) b_data) != 0) {
                printf("%s: hio_i_sds failure\n", __FILE__);
                iret = 1;
                exit(iret);
            }
        }
    }
    /*
     *  close the file
     */
    if (hio_close(hinfo) != 0) {
        printf("%s: hio_close failure\n", __FILE__);
        iret = 1;
    }
    free(xfr_arr);
    free(b_xfr);
    free(data);
    free(b_data);

    exit(iret);
}
