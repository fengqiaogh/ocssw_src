/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/hdf/WIFSHDF.h,v 4.10 1995/01/17 19:58:14 seawifsd Exp seawifsd $
$Log: WIFSHDF.h,v $
Revision 4.10  1995/01/17 19:58:14  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:14:30  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:46:33  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:31  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.3  1994/10/04 15:48:40  frank
added the proper header files to be compiled under ANSI.

Revision 1.1.1.2  1994/05/23 15:27:51  frank
included prototype declaration header file 'WIFSHDF_proto.h' so that
compiling with ANSI will work correctly.

Revision 1.2  1994/05/10 18:48:55  seawifst
May 6, 1994 version 1.2
        Summary 1.1.1.1 to 1.1.1.1(total 1 revisions):
 * added BY_LABEL definition

Revision 1.1.1.1  1994/04/20 19:10:39  frank
added BY_LABEL definition

Revision 1.1  1994/04/19 13:32:33  seawifst
Initial revision


 */


#ifndef WIFSHDF_H_
#define WIFSHDF_H_

#define ACCESS_BY_NAME 1
#define ACCESS_BY_INDEX 2
#define ACCESS_BY_ID 3
#define ACCESS_BY_REF 4
#define ACCESS_BY_LABEL 5

#define BY_NAME  ACCESS_BY_NAME
#define BY_INDEX ACCESS_BY_INDEX
#define BY_ID  ACCESS_BY_ID
#define BY_REF  ACCESS_BY_REF
#define BY_LABEL ACCESS_BY_LABEL

#define BY_NAME_STR "BY_NAME"
#define BY_INDEX_STR "BY_INDEX"
#define BY_ID_STR "BY_ID"
#define BY_REF_STR "BY_REF"
#define BY_LABEL_STR "BY_LABEL"

//#include "usrhdr.h"
#include "usrmac.h"
//#include "hdfhdr.h"
//#include "hdfmac.h"
#include "SeaWiFS.h"
#include "WIFSHDF_proto.h"
#endif /* WIFSHDF_H_ */
