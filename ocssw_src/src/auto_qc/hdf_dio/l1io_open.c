#include "l1io.h"
#include <mfhdf.h>

int32_t l1io_open(char *fname, l1info_struct *l1info,
        int32_t *npix, int32_t *nlin)
/*******************************************************************

   get_l1a_open

   purpose: open the level 1 dataset and send back the # lines, pixels

   Returns type: int32 - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            fname            I      name of file to open
      struct l1info_struct * l1info      O      information  struct
                                                about the opened file
      int32 *           npix             O      # of pixels in data
      int32 *           nlin             O      # of lines in data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       16-Aug-1994     Original development

 *******************************************************************/
 {

    int32 attr_index;
    /*  hopen will open the file and read data descriptor blocks
        to memory
     */
    if ((l1info->fid = Hopen(fname, DFACC_RDONLY, 0)) < 0) {
        printf("l1io_open: Failed on the Hopen of \n%s\n", fname);
        return -1;
    }

    /*
     *  SDstart opens the hdf interface and inits the SD interface, wierd.
     *  apparently, if both SD and HDF are used, both these need to
     *  be called.
     */
    if ((l1info->sdfid = SDstart(fname, DFACC_RDONLY)) < 0) {
        printf("l1io_open: Failure at SDstart of \n%s\n", fname);
        l1io_close(*l1info);
        return -1;
    }

    /*
     *  ok, it's open.
     *  the # samples and # lines are global attributes, get here
     *  first get the index for the name
     */
    /*  for use before 15 Aug 94 I/O v3.0
       if( ( attr_index = SDfindattr( l1info->sdfid, "Pixels per scan line" ) )
           == -1 )
     */

    /*   for use after 15 Aug 94 I/O v3.0  */
    if ((attr_index = SDfindattr(l1info->sdfid, "Pixels per Scan line"))
            == -1)
        /*  
         *  Well, there is another possibility, and I don't wish to exclude 
         *  the older name yet, so check for this possibility
         */ {
        if ((attr_index = SDfindattr(l1info->sdfid, "Pixels per Scan Line"))
                == -1) {
            printf("l1io_open: Error trying to get global attr Pixels per Scan line (or ...Line)\n");
            l1io_close(*l1info);
            return -1;
        }
    }
    /*
     *  for an unknown attribute, you would call SDattrinfo so you
     *  could allocate the right type and space for the data but
     *  we know this already (sic) so
     *
     *  read in the attribute
     */
    if (SDreadattr(l1info->sdfid, attr_index, npix) == -1) {
        printf("l1io_open: Error trying to get data for: Pixels per Scan line\n");
        l1io_close(*l1info);
        return -1;
    }

    /*
     *  get the # lines also the same way
     */
    /*   for use after 15 Aug 94 I/O v3.0  */
    if ((attr_index = SDfindattr(l1info->sdfid, "Number of Scan Lines"))
            == -1)
        /*  for use before 15 Aug 94 I/O v3.0
           if( ( attr_index = SDfindattr( l1info->sdfid, "Number of scan lines" ) )
               == -1 )
         */ {
        printf("l1io_open: Error trying to get global attr: Number of Scan Lines\n");
        l1io_close(*l1info);
        return -1;
    }

    if (SDreadattr(l1info->sdfid, attr_index, nlin) == -1) {
        printf("l1io_open: Error trying to get data for Number of Scan Lines\n");
        l1io_close(*l1info);
        return -1;
    }
    /*
     *  and return successfully
     */
    return 0;

}
