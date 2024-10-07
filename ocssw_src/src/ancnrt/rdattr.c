/***************************************************************** 
 * File:          rdattr  
 *
 * Purpose:       read global attributes from an HDF file
 *
 * Description:   string arrays of label/value pairs read from the 
 *                HDF file.
 *
 * Input parms: 
 * char *filename HDF file name basename (extension '.hdf' added) 
 *
 * Output parms:
 * GLOBAL struct char *annot - array of header labels and values read
 *
 * Returns:       -1 on failure, 0 on success
 *
 * Subs called:   HDF library routines
 *
 * History:       based on DAAC "Metadata Submission Guide" 2/93
 *
 * Author:        Brian D. Schieber, GSC, 2/93
 *
 * Notes:         Uses header file (*.h) descriptors for maximum label,  
 *                value and filename lengths.
 *
 * Mod history:   
 *   BDS, 9/1/93	- modified to support new HDF formats.
 *
 *****************************************************************/

#include <mfhdf.h>
#include "ancil.h"

int rdattr(char *filename) {

    /*  
     * local variables 
     */

    int32 dfile;
    int status;
    int attnum;
    char name[50];
    int32 nt;
    int32 count = 0;

    /**** open header file ***/

    dfile = SDstart(filename, DFACC_RDONLY);
    if (dfile == FAIL) return -1;

    /****  read the first label and description  ******/

    attnum = 1;
    status = SDattrinfo(dfile, attnum, name, &nt, &count);
    if (status < 0) printf("< 0  status from SDattrinfo\n");
    if (status < 0) return (-1);

    printf("SDattrinfo returned:\n");
    printf("attnum: %d\n", attnum);
    printf("name:   [%s]\n", name);
    printf("nt:     %d\n", nt);
    printf("count:  %d\n", count);

#if 0
    /******  read the rest of the labels and descriptions *******/

    while (llength >= 0) {
        llength = DFANgetfidlen(dfile, NOTFIRST);
        if (llength > 0) {
            count++;
            lreturn = DFANgetfid(dfile, annot[count].label, MAXLABLEN, NOTFIRST);

            dlength = DFANgetfdslen(dfile, NOTFIRST);

            dreturn = DFANgetfds(dfile, annot[count].descr, MAXDESCLEN, NOTFIRST);
        }
    }

    SDend(dfile);
    if (result < 0) return (-1);
#endif

    return 0;

} /* rdattr */

