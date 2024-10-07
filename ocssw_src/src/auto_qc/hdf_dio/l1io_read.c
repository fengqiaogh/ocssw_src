#include <string.h>
#include "l1io.h"
#include <mfhdf.h>

int32_t l1io_read(l1info_struct l1info, int irec,
        int16_t *data, navblockType *nav)
/*******************************************************************

   get_l1a_read

   purpose: read a record of some information from the level 1 dataset 

   Returns type: int32 - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      struct l1info_struct   l1info      I      information  struct
                                                about the opened file
      int               irec             I      line, record # to read
      int16 *           data             O      array of data: 1 line
                                                with 8 channels
      struct navblockType *   nav        O      navigation info struct

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       16-Aug-1994     Original development

 *******************************************************************/
 {
    int i, j;
    int32 index, rank, sdid, numbertype, nattrs, data_dims[5];
    int32 start[5], edge[5];
    void *ptr; /*  generic data pointer  */
    char name[H4_MAX_NC_NAME]; /*  MAX_NC_NAME is max # chars for this */

    /*   for use after 15 Aug 94 I/O spec v3.0  */
    static char *namelist[8] ={"l1a_data", "orb_vec", "l_vert", "sun_ref",
        "att_ang", "sen_mat", "scan_ell", "nflag"};

    /*  for use before 15 Aug 94 I/O spec v3.0
       static char *namelist[8] =
          { "scan_line", "orb_vec", "l_vert", "sun_ref",
            "att_ang", "sen_mat", "scan_ell", "nflag" };
     */

    /*
     *  now, SDnametoindex is used to zero in on the data we want,
     *  ie. get the index of the data
     */

    for (i = 0; i < 8; i++) {
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
        case 0: /*  "l1a_data"   */
        {
            ptr = (void *) data;
            break;
        }
        case 1: /*  "orb_vec"  */
        {
            ptr = (void *) (nav->orb_vec);
            break;
        }
        case 2: /*  "l_vert"  */
        {
            ptr = (void *) (nav->l_vert);
            break;
        }
        case 3: /*  "sun_ref"  */
        {
            ptr = (void *) (nav->sun_ref);
            break;
        }
        case 4: /*  "att_ang"  */
        {
            ptr = (void *) (nav->att_ang);
            break;
        }
        case 5: /*  "sen_mat"  */
        {
            ptr = (void *) (nav->sen_mat);
            break;
        }
        case 6: /*  "scan_ell"  */
        {
            ptr = (void *) (nav->scan_ell);
            break;
        }
        case 7: /*  "nflag"  */
        {
            ptr = (void *) (nav->nflag);
            break;
        }
        }
        /*
         *  read the data
         */
        if (SDreaddata(sdid, start, NULL, edge, ptr) < 0)
            /*      if( SDreaddata( sdid, start, stride, edge, ptr ) < 0 )  */ {
            printf("l1io_read: failure to read data for item %s\n", namelist[i]);
            return -1;
        }
    }
    return 0;
}
