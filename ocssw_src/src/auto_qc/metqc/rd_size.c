#include <mfhdf.h>
#include "l1io.h"

int rd_size(char *file, int *npix, int *nlin)
/*******************************************************************

   rd_size

   purpose: open a hdf file, read the size from the attributes
            and close it  (assume size is stored in the 'Number of Rows'
            and 'Number of Columns' attributes)

   Returns type: int - 0 if all went well, 
                       -1 can't find the attribute

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      file to work with
      int *             npix             O      # pixels (columns)
      int *             nlin             O      # lines (rows)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       5-Feb-1997      Original development

 *******************************************************************/
 {
    l1info_struct hdf_info;
    int32 n_type, count;

    n_type = DFNT_INT32;
    count = 1;
    /*
     *  open the file
     */
    if (open_hdf(file, &hdf_info) != 0) {
        return -1;
    }
    /*
     *  read the 2 attributes
     */
    if (read_g_attr(hdf_info, "Number of Rows", &n_type, &count, nlin) != 0) {
        return -1;
    }
    if (read_g_attr(hdf_info, "Number of Columns", &n_type, &count, npix) != 0) {
        return -1;
    }
    /*
     *  close the file
     */
    l1io_close(hdf_info);

    return 0;
}
