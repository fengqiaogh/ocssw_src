#include <string.h>
//#include "netcdf.h"
#include "mfhdf.h"
#include "l1io.h"

int32 nav_read(l1info_struct l1info, int irec,
        navblockType *nav)
/*******************************************************************

   nav_read

   purpose: read just a record of nav info from a hdf dataset 

   Returns type: int32 - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      struct l1info_struct   l1info      I      information  struct
                                                about the opened file
      int               irec             I      line, record # to read
      struct navblockType *   nav        O      navigation info struct

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       10-Feb-1995     Original development

 *******************************************************************/
 {
    int i, j;
    int32 index, rank, sdid, numbertype, nattrs, data_dims[5];
    int32 start[5], edge[5];
    void *ptr; /*  generic data pointer  */
    char name[H4_MAX_NC_NAME]; /*  MAX_NC_NAME is max # chars for this */

    /*   This is the list of names of the data elements we will be getting  */
    static char *namelist[4] ={"orb_vec", "sun_ref", "sen_mat", "scan_ell"};

    /*
     *  now, SDnametoindex is used to zero in on the data we want,
     *  ie. get the index of the data
     */

    for (i = 0; i < 4; i++) {
        if ((index = SDnametoindex(l1info.sdfid, namelist[i])) < 0) {
            printf("l1io_read: couldn't find %s in the qc dataset\n",
                    namelist[i]);
            return -1;
        }

        /*
         *  next, get the SD ID for the data
         */
        if ((sdid = SDselect(l1info.sdfid, index)) < 0) {
            printf("l1io_read: Failed in SDselect for item %s\n",
                    namelist[i]);
            return -1;
        }

        /*
         *  now, get the size information of the data item
         */
        if (SDgetinfo(sdid, name, &rank, data_dims, &numbertype,
                &nattrs) < 0) {
            printf("l1io_read: Failed in SDgetinfo for item %s\n",
                    namelist[i]);
            return -1;
        }

        /*
         *  compute the read controls
         */
        for (j = 0; j < rank; j++) {
            start[j] = 0;
            edge[j] = 1; /*  test WDR */
            edge[j] = data_dims[j];
        }
        start[0] = irec; /* the first dim is the row dim so set for the rec */
        edge[0] = 1;

        /*
         *  set the pointer to read into for the 8 sds arrays being read
         */

        switch (i) {
        case 0: /*  "orb_vec"  */
        {
            ptr = (void *) (nav->orb_vec);
            break;
        }
        case 1: /*  "sun_ref"  */
        {
            ptr = (void *) (nav->sun_ref);
            break;
        }
        case 2: /*  "sen_mat"  */
        {
            ptr = (void *) (nav->sen_mat);
            break;
        }
        case 3: /*  "scan_ell"  */
        {
            ptr = (void *) (nav->scan_ell);
            break;
        }
        }
        /*
         *  read the data
         */
        if (SDreaddata(sdid, start, NULL, edge, ptr) < 0) {
            printf("l2io_read: failure to read data for item %s\n", namelist[i]);
            return -1;
        }
    }
    return 0;
}
