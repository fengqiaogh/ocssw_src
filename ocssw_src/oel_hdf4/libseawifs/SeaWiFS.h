/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.4/L012_Util/SeaWiFS.h,v 4.12 1995/12/07 18:50:40 seawifsd Exp seawifsd $
$Log: SeaWiFS.h,v $
Revision 4.12  1995/12/07 18:50:40  seawifsd
increased max number of objects in HDF file from 120 to 130
due to the added objects in PROD_SPEC_V2.8. Added macro definition
NAME_DELIMETER.

Revision 4.11  1995/02/24 14:55:54  seawifsd
added support of new object DFANObj.

Revision 4.10  1995/01/17 20:59:37  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:13:21  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:49:18  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:03:15  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.1  1994/10/03 18:36:26  frank
increased the MaxObjtotal from 100 to 120 because the changes in spec.

Revision 1.2  1994/05/10 18:33:47  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:19:48  seawifst
Initial revision

Included objmax.h to get Max*total definitions.
Norman Kuring		7-Nov-1996

 */


#ifndef SEAWIFS_H_
#define SEAWIFS_H_

#include "objmax.h"

typedef struct SeaWiFS_HDFStruct {
    int n; /* total records		*/
    int fid; /* HDF file id for SD routines	*/
    int fid2; /* HDF file id for all others	*/
    char *fname; /* HDF data file name		*/
    int ngattr; /* number of global attributes	*/
    int nvgrp; /* number of VGRP object	*/
    int nmfsd; /* number of MFSD object	*/
    int nsfsd; /* number of SFSD object	*/
    int nvset; /* number of VSET object	*/
    int ndfan; /* number of DFAN object	*/
    ATTRObj attr[MaxATTRtotal]; /* buffer for all ATTR objects	*/
    VGRPObj vgrp[MaxVGRPtotal]; /* buffer for all VGRP objects	*/
    MFSDObj mfsd[MaxMFSDtotal]; /* buffer for all MFSD objects	*/
    SFSDObj sfsd[MaxSFSDtotal]; /* buffer for all SFSD objects	*/
    VSETObj vset[MaxVSETtotal]; /* buffer for all VSET objects	*/
    DFANObj dfan[MaxDFANtotal]; /* buffer for all DFAN objects	*/
    Objidx oidx[MaxObjtotal]; /* index reference for all objects*/
} SeaWiFS_HDFType;


/* this is used to indicate the next record will be read		*/
#define GET_NEXT_RECORD -1

/* define the delimeter character for name list				*/
#define NAME_DELIMETER ","

#endif /* SEAWIFS_H_ */
