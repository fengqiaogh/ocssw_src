#include "l1io.h"
#include <mfhdf.h>

void l1io_close(l1info_struct l1info)
/*******************************************************************

   get_l1a_close

   purpose: close the level 1 dataset 

   Returns type: void - nothing

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      struct l1info_struct   l1info      O      information  struct
                                                about the opened file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       16-Aug-1994     Original development

 *******************************************************************/
 {

    /*
     *  just close the file up
     */
    Hclose(l1info.fid);
    SDend(l1info.sdfid);
}
