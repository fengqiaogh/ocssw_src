/*******************************************************************

   hdfio.c - general routines for hdf file input, output

   Contents:
     hdfio_open - open the file
     hdfio_rd_gattr - read a global attribute
     hdfio_rd_sd - read a science dataset
     hdfio_close - close a hdf file

 *******************************************************************/

#include <string.h>
#include "hdfio.h"

int hdfio_open(char *fname, hdfio_struc *hdfinfo)
/*******************************************************************

   hdfio_open

   purpose: open a general hdf file

   Returns type: int - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            fname            I      name of file to open
      hdfio_struc *     hdfinfo          O      information  struct
                                                about the opened file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 5 Aug 2004      Original development

 *******************************************************************/
 {

    if ((hdfinfo->fid = Hopen(fname, DFACC_RDONLY, 0)) < 0) {
        printf("hdfio_open: Failed on the Hopen of \n%s\n", fname);
        return -1;
    }

    /*
     *  SDstart opens the hdf interface and inits the SD interface, wierd.
     *  apparently, if both SD and HDF are used, both these need to
     *  be called.
     */
    if ((hdfinfo->sdfid = SDstart(fname, DFACC_RDONLY)) < 0) {
        printf("hdfio_open: Failure at SDstart of \n%s\n", fname);
        hdfio_close(*hdfinfo);
        return -1;
    }

    /*
     *  ok, it's open.
     */
    return 0;
}

int hdfio_rd_gattr(hdfio_struc hdfinfo, char *name, int32 n_type,
        int32 count, void *data)
/*******************************************************************

   hdfio_rd_gattr

   purpose: read a global attribute to a waiting array (if size correct
              and type is correct)

   Returns type: int - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      hdfio_struc       hdfinfo          I      information  struct
                                                about the opened file
      char *            name             I      name of glabal attr to read
      int32             n_type           I      number type expected using
                                                HDF nomenclature
      int32             count            I      count of values
      void*             data             O      array of data from the
                                                global attribute 

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 5 Aug 2004      Original development

 *******************************************************************/
 {
    int32 attr_index;
    char a_name[200];
    int32 a_n_type, a_count;

    /*
     *  get the index for that attribute name
     */
    if ((attr_index = SDfindattr(hdfinfo.sdfid, name))
            == -1) {
        printf("hdfio_rd_gattr: Error in SDfindattr for attribute: '%s'\n",
                name);
        return -1;
    }
    /*
     *  find out about the attribute and be sure it is the type and size expected
     */
    if (SDattrinfo(hdfinfo.sdfid, attr_index, a_name, &a_n_type, &a_count)
            == -1) {
        printf("hdfio_rd_gattr: Error in SDattrinfo for attribute: '%s'\n",
                name);
        return -1;
    }

    if (a_n_type != n_type) {
        printf(
                "hdfio_rd_gattr: attribute: '%s', difference in number type found\n",
                name);
        printf("          expected: %d, read: %d\n", (int) n_type, (int) a_n_type);
        return -1;
    }

    if (((n_type != DFNT_CHAR) && (a_count != count)) ||
            ((n_type == DFNT_CHAR) && (count < a_count))) {
        printf("hdfio_rd_gattr: attribute: '%s', difference in count found\n",
                name);
        printf("             expected: %d, read: %d\n", (int) count, (int) a_count);
        return -1;
    }

    /*
     *  read in the stuff
     */
    if (SDreadattr(hdfinfo.sdfid, attr_index, data) == -1) {
        printf("hdfio_rd_gattr: attribute: '%s', Error reading data\n",
                name);
        return -1;
    }
    return 0;
}

int hdfio_rd_sd(hdfio_struc hdfinfo, char *arr_name, void *array)
/*******************************************************************

   hdfio_rd_sd

   purpose: read a full science dataset to a float array

   Returns type: int - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      hdfio_struc       hdfinfo          I      hdf information  struct
                                                about the opened file
      char *            arr_name         I      name of science dataset 
                                                to read
      void *            array            O      array of data from the sds
                                                already allocated outside

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       9-Feb-1995      Original development
      W. Robinson, SAIC  5 Aug 2004     update this version to do more general 
                                        void array type

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
    if ((index = SDnametoindex(hdfinfo.sdfid, arr_name)) < 0) {
        printf("hdfio_rd_sd: couldn't find %s in the dataset\n",
                arr_name);
        return -1;
    }

    /*
     *  next, get the SD ID for the data
     */
    if ((sdid = SDselect(hdfinfo.sdfid, index)) < 0) {
        printf("hdfio_rd_sd: Failed in SDselect for item %s\n",
                arr_name);
        return -1;
    }

    /*
     *  now, get the size information of the data item
     */
    if (SDgetinfo(sdid, name, &rank, data_dims, &numbertype,
            &nattrs) < 0) {
        printf("hdfio_rd_sd: Failed in SDgetinfo for item %s\n",
                arr_name);
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
        printf("hdfio_rd_sd: failure to read data for item %s\n", arr_name);
        return -1;
    }
    return 0;
}

void hdfio_close(hdfio_struc hdfinfo)
/*******************************************************************

   hdfio_close

   purpose: close the hdf dataset 

   Returns type: void - nothing

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      hdfio_struc       hdfinfo          O      information  struct
                                                about the opened file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC, 5 Aug 2004     Original development

 *******************************************************************/
 {

    /*
     *  just close the file up
     */
    Hclose(hdfinfo.fid);
    SDend(hdfinfo.sdfid);
}
