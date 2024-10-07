#include "hio.h"
#include <mfhdf.h>

int32_t hio_open(char *fname, int32_t mode, hio_struct *info)
/*******************************************************************

   hio_open

   purpose: open a general hdf file

   Returns type: int32 - return status: 0 is good
                        -1 if any other error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            fname            I      name of file to open
      int32             mode             I      open mode: DFACC_RDONLY,
                                                  DFACC_CREATE, DFACC_RDWR
      struct hio_struct * info           O      information  struct
                                                about the opened file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       5-Feb-1997      Original development

 *******************************************************************/
 {

    if (((mode == DFACC_RDONLY) || (mode == DFACC_RDWR)) &&
            !Hishdf(fname)) {
        printf("hio_open: Failure at Hishdf of \n%s\n", fname);
        return -1;
    }
    /*
     *  SDstart opens the hdf interface and inits the SD interface, wierd.
     *  apparently, if both SD and HDF are used, both these need to
     *  be called.
     */
    if ((info->sdfid = SDstart(fname, mode)) < 0) {
        printf("hio_open: Failure at SDstart of \n%s\n", fname);
        hio_close(*info);
        return -1;
    }

    /*
     */
    if ((info->fid = Hopen(fname, DFACC_RDONLY, 0)) < 0) {
        printf("hio_open: Failed on the Hopen of \n%s\n", fname);
        return -1;
    }
    /*
     *  ok, it's open.
     */
    return 0;

}
