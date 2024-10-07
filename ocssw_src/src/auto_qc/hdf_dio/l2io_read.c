#include <string.h>
#include "l1io.h"
#include <mfhdf.h>

int32 l2io_read(l1info_struct l1info, int irec,
        float *data, int16 *l2_flags, navblockType *nav)
/*******************************************************************

   l2io_read

   purpose: read a record of some information from the level 2 dataset 

   Returns type: int32 - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      struct l1info_struct   l1info      I      information  struct
                                                about the opened file
      int               irec             I      line, record # to read
      float *           data             O      array of data: 1 line
                                                with 12 channels, all vals
                                                for 1 pixel are together
      int16 *           flags            O      level 2 flag values
      struct navblockType *   nav        O      navigation info struct

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       25-Jan-1995     Original development

 *******************************************************************/
 {
    int i, j, ipix;
    int32 index, rank, sdid, numbertype, nattrs, data_dims[5];
    int32 start[5], edge[5];
    void *ptr; /*  generic data pointer  */
    char name[H4_MAX_NC_NAME]; /*  MAX_NC_NAME is max # chars for this */
    /*  temporary storage for the unscaled data  */
    int16 *i2temp;
    uint8 *i1temp;
    /*  the bias and slope are stored here  */
    static float l2_bias[12] ={0., 0., 0., 0., 0., 0., 0., 32., 32., 0., 0., 0.};
    static float l2_slope[12] ={.001, .001, .001, .001, .001, .002,
        .001, .001, .001, .0002, .01, .005};

    /*   This is the list of names of the data elements we will be getting  */
    static char *namelist[20] ={"nLw_412", "nLw_443", "nLw_490", "nLw_510",
        "nLw_555", "La_670", "La_865", "CZCS_pigment",
        "chlor_a", "K_490", "eps_78", "tau_865",
        "l2_flags", "orb_vec", "l_vert", "sun_ref",
        "att_ang", "sen_mat", "scan_ell", "nflag"};

    /*
     *  before the festivities start, allocate space for the temp storage
     */
    i2temp = malloc(248 * 2);
    i1temp = malloc(248);

    /*
     *  now, SDnametoindex is used to zero in on the data we want,
     *  ie. get the index of the data
     */

    for (i = 0; i < 20; i++) {
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
        case 0: /*  This is the 10 channels of 2-byte values  */
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        {
            ptr = (void *) i2temp;
            break;
        }
        case 10: /*  This is the 2 channels of byte data  */
        case 11:
        {
            ptr = (void *) i1temp;
            break;
        }
        case 12:
        {
            ptr = (void *) l2_flags;
            break;
        }
        case 13: /*  "orb_vec"  */
        {
            ptr = (void *) (nav->orb_vec);
            break;
        }
        case 14: /*  "l_vert"  */
        {
            ptr = (void *) (nav->l_vert);
            break;
        }
        case 15: /*  "sun_ref"  */
        {
            ptr = (void *) (nav->sun_ref);
            break;
        }
        case 16: /*  "att_ang"  */
        {
            ptr = (void *) (nav->att_ang);
            break;
        }
        case 17: /*  "sen_mat"  */
        {
            ptr = (void *) (nav->sen_mat);
            break;
        }
        case 18: /*  "scan_ell"  */
        {
            ptr = (void *) (nav->scan_ell);
            break;
        }
        case 19: /*  "nflag"  */
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
            printf("l2io_read: failure to read data for item %s\n", namelist[i]);
            return -1;
        }
        /*
         *  for the data that needs to be scaled to floats, do the work here
         */
        if (i < 12) {
            if (i < 10) /*  do the 2 byte values here */ {
                for (ipix = 0; ipix < 248; ipix++) {
                    *(data + i + 12 * ipix) =
                            l2_bias[i] + l2_slope[i] * i2temp[ipix];
                }
            } else /*  do the 1 byte values here */ {
                for (ipix = 0; ipix < 248; ipix++) {
                    *(data + i + 12 * ipix) =
                            l2_bias[i] + l2_slope[i] * i1temp[ipix];
                }
            }
        }
    }
    /*
     *  free that temp space
     */
    free(i2temp);
    free(i1temp);
    return 0;
}
