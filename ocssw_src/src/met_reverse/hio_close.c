#include "hio.h"

int32_t hio_close(hio_struct info)
/*******************************************************************

   hio_close

   purpose: close the hdf dataset 

   Returns type: int 0 if OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      struct hio_struct  info            I      information  struct
                                                about the opened file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       30-Sep-2005     from l1io_close

 *******************************************************************/
 {
    int iret;
    /*
     *  just close the file up
     */
    iret = Hclose(info.fid);
    if (iret == 0)
        iret = SDend(info.sdfid);
    return iret;
}
