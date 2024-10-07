#include <string.h>
#include "hio.h"

int32 hio_r_sds(hio_struct info, char *arr_name, int32 exp_ntyp,
        void *array)
/*******************************************************************

   hio_r_sds

   purpose: read a full science dataset to an array

   Returns type: int32 - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      struct hio_struct info             I      information  struct
                                                about the opened file
      char *            arr_name         I      name of science dataset 
                                                to read
      int32             exp_ntyp         I      expected number type
      void*             array            O      array of data from the sds
                                                already allocated outside

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       16-Sep-2005      Original development

 *******************************************************************/
 {
    int j;
    int32 index, rank, sdid, numbertype, nattrs, data_dims[5];
    int32 start[5], edge[5];
    char name[H4_MAX_NC_NAME]; /*  MAX_NC_NAME is max # chars for this */

    /*
     *  SDnametoindex is used to zero in on the data we want,
     *  ie. get the index of the data
     */
    if ((index = SDnametoindex(info.sdfid, arr_name)) < 0) {
        printf("%s: couldn't find %s in the dataset\n", __FILE__,
                arr_name);
        return -1;
    }

    /*
     *  next, get the SD ID for the data
     */
    if ((sdid = SDselect(info.sdfid, index)) < 0) {
        printf("%s: Failed in SDselect for item %s\n", __FILE__,
                arr_name);
        return -1;
    }

    /*
     *  now, get the size information of the data item
     */
    if (SDgetinfo(sdid, name, &rank, data_dims, &numbertype,
            &nattrs) < 0) {
        printf("%s: Failed in SDgetinfo for item %s\n", __FILE__,
                arr_name);
        return -1;
    }

    if (numbertype != exp_ntyp) {
        printf(
                "%s: Type of data to be read is not type expected for item %s\n",
                __FILE__, arr_name);
        printf("actual type: %d, expected type: %d\n", numbertype, exp_ntyp);
        return -1;
    }

    /*
     *  compute the read controls for whole array
     */

    for (j = 0; j < rank; j++) {
        start[j] = 0;
        edge[j] = data_dims[j];
    }

    /*
     *  read the data
     */
    if (SDreaddata(sdid, start, NULL, edge, array) < 0) {
        printf("%s: failure to read data for item %s\n", __FILE__, arr_name);
        return -1;
    }
    return 0;
}
