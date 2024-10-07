#include <mfhdf.h>
#include "ancil.h"
#include "ancnrt_proto.h"

/***************************************************************** 
 * Function: 			rdsdsid
 * Purpose:				return HDF Scientific Data Set (SDS) id value
 * Description:		read HDF Scientific Data Set (SDS) and return 
 *                	Reads one SDS at a time.
 * Input parameters: 		
 * Output parameters: 		
 * Returns:
 * Local variables:	
 * Subs called:		HDF library routines
 * History:     		none
 * Note:					similar to rdsds routine but only returns 
 *							refid for a simple (single Vgroup level) SDS.
 *							Primarily used as a workaround to IDL not yet
 *							supporting the HDF structures.  Perhaps this 
 *							can be removed in the future.
 *
 * Author:				Brian D. Schieber, GSC, 2/94
 * Modification history:	
 *****************************************************************/

int rdsdsid(fid, sdfid, vgname, sdslabel)
int32 fid, sdfid;
char *vgname;
char *sdslabel;
{
    int j;
    int32 vid, vg;
    int32 sdsindex, sdsid;
    int32 rank, numbertype;
    int32 tagarray1[255], refarray1[255], nrefs1, nattrs1;
    int32 dimsizes[2];
    char lname[MAXNAMELNG];

    /* Open the Vgroup by name and fileid, attach to it and get the number 
     * and tag type (SDG, VH, etc) in it, loop over each item and find 
     * the ones that match the desired SDS and return SDSid.
     */

    if ((vid = Vfind(fid, vgname)) <= 0) return (ERROR);
    vg = Vattach(fid, vid, "r");
    nrefs1 = Vntagrefs(vg);
    nattrs1 = nrefs1;
    Vgettagrefs(vg, tagarray1, refarray1, nattrs1);

    for (j = 0; j < nrefs1; j++) {
        if ((tagarray1[j] == DFTAG_NDG) ||
                (tagarray1[j] == DFTAG_SD)) { /* SDS found */
            sdsindex = SDreftoindex(sdfid, refarray1[j]);
            sdsid = SDselect(sdfid, sdsindex);
            SDgetinfo(sdsid, lname, &rank, dimsizes, &numbertype,
                    &nattrs1);
            if (!strcmp(sdslabel, lname)) {
                Vdetach(vg);
                return (sdsid); /* SUCCESS */
            }
        } /* tagarray1[j] */
    } /* for j */
    Vdetach(vg);
    pwarning("SDS data id not found");
    return (-1);
} /* rdsdsid */

