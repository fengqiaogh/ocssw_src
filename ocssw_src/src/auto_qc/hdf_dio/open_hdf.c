#include "l1io.h"

#include <mfhdf.h>

int32 open_hdf(char *fname, l1info_struct *l1info)
/*******************************************************************

   open_hdf

   purpose: open a general hdf file

   Returns type: int32 - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            fname            I      name of file to open
      struct l1info_struct * l1info      O      information  struct
                                                about the opened file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       5-Feb-1997      Original development

 *******************************************************************/
 {

    /*
     *  hopen will open the file and read data descriptor blocks
     *  to memory
     */
    if ((l1info->fid = Hopen(fname, DFACC_RDONLY, 0)) < 0) {
        printf("open_hdf: Failed on the Hopen of \n%s\n", fname);
        return -1;
    }

    /*
     *  SDstart opens the hdf interface and inits the SD interface, wierd.
     *  apparently, if both SD and HDF are used, both these need to
     *  be called.
     */
    if ((l1info->sdfid = SDstart(fname, DFACC_RDONLY)) < 0) {
        printf("open_hdf: Failure at SDstart of \n%s\n", fname);
        l1io_close(*l1info);
        return -1;
    }

    /*
     *  ok, it's open.
     */
    return 0;

}
