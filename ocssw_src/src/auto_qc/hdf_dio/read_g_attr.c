#include <string.h>
//#include "netcdf.h"
#include "mfhdf.h"
#include "l1io.h"

int32 read_g_attr(l1info_struct l1info, char *name, int32 *n_type,
        int32 *count, void *data)
/*******************************************************************

   read_g_attr

   purpose: read a global attribute to a waiting array (if size correct
              and type is correct)

   Returns type: int32 - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      struct l1info_struct   l1info      I      information  struct
                                                about the opened file
      char *            name             I      name of glabal attr to read
      int32             n_type           I      number type expected using
                                                HDF nomenclature
      int32 *           count            I      count of values
      void*             data             O      array of data from the
                                                global attribute 

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       5-Feb-1997      Original development

 *******************************************************************/
 {
    int32 attr_index;
    char a_name[200];
    int32 a_n_type, a_count;

    /*
     *  get the index for that attribute name
     */
    if ((attr_index = SDfindattr(l1info.sdfid, name))
            == -1) {
        printf("read_g_attr: Error in SDfindattr for attribute: '%s'\n",
                name);
        return -1;
    }
    /*
     *  find out about the attribute and be sure it is the type and size expected
     */
    if (SDattrinfo(l1info.sdfid, attr_index, a_name, &a_n_type, &a_count)
            == -1) {
        printf("read_g_attr: Error in SDattrinfo for attribute: '%s'\n",
                name);
        return -1;
    }

    if (a_n_type != *n_type) {
        printf(
                "read_g_attr: attribute: '%s', difference in number type found\n",
                name);
        printf("          expected: %d, read: %d\n", *n_type, a_n_type);
        return -1;
    }

    if (a_count != *count) {
        printf("read_g_attr: attribute: '%s', difference in count found\n",
                name);
        printf("             expected: %d, read: %d\n", *count, a_count);
    }

    /*
     *  read in the stuff
     */
    if (SDreadattr(l1info.sdfid, attr_index, data) == -1) {
        printf("read_g_attr: attribute: '%s', Error reading data\n",
                name);
        return -1;
    }
    return 0;
}
