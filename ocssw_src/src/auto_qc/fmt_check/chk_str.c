#include "fmt_check.h"
extern int fmt_status; /* status of format check, see fmt_check */

void chk_str(attr_str attr, char *str, int32 len)
/*******************************************************************

   chk_str

   purpose: make sure that a string read from a hdf dataset has 
     only one null at the end

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      attr_str          attr             I      attribute description
                                                structure.
      char *            str              I      string to check
      int32             len              I      length reported by HDF

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       20-Nov-1995     Original development
      L. Kumar		19-Dec-1995     Modified an output (printf) stmt
      W. Robinson       17-Sep-1996     Fix strings so they print if
                                        null is imbedded somewhere (noted
                                        here but done in attr_disp because 
                                        string can be used and should not 
                                        be altered)

 *******************************************************************/
 {
    int32 i;

    /*
     *  only do this if it is char string read in
     */
    if ((attr.type == DFNT_CHAR) && (attr.read != -1)) {
        /*
         *  loop through the bytes and look for the null at the right place
         */
        for (i = 0; i < len; i++) {
            if ((str[i] == 0) && (i != len - 1)) {
                printf("**** object: '%s', attribute: '%s', String has \n",
                        attr.obj_nm, attr.access_nm);
                printf("     a null at a location ( %d ) other than the end\n",
                        i);
                printf("     For printing, replaced with '?'\n");
                printf("     (NOTE, no error condition set currently\n");
                /*   WDR temporarily disable for processing
                            fmt_status = fmt_status | 2;
                 */
            } else if ((str[i] != 0) && (i == len - 1)) {
                printf("**** object: '%s', attribute: '%s', String lacks \n",
                        attr.obj_nm, attr.access_nm);
                /*
                            printf( "     a null at the end location ( %d )\n", i ); 
                 */
                printf("     a null at the end location ( %d )\n", i + 1); /*LK*/
                fmt_status = fmt_status | 2;
            }
        }
    }
}
