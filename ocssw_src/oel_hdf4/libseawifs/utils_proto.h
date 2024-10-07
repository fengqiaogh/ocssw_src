/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.3/L012_Util/utils_proto.h,v 4.12 1995/05/04 17:56:22 seawifsd Exp seawifsd $
$Log: utils_proto.h,v $
Revision 4.12  1995/05/04 17:56:22  seawifsd
replaced all references of 'exist' to 'generic'.

Revision 4.11  1995/01/25 18:03:59  seawifsd
added include eng_qual_proto.h

Revision 4.10  1995/01/17 20:59:40  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:13:23  seawifsd
Jan. 9, 1994, 4.0

Revision 3.4  1994/12/23 15:41:04  seawifsd
commented out unused include file cdl_attr_proto.h

Revision 3.3  1994/11/08 18:49:22  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:03:18  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.2  1994/10/03 18:41:43  frank
included exist_proto.h and get_norad_orbnum_proto.h.

Revision 1.1.1.1  1994/05/19 15:40:07  frank
took out newtime_proto.h include statement. This is a duplicate of
utiltime_proto.h.

Revision 1.2  1994/05/10 18:54:16  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:52:38  seawifst
Initial revision


 */


#include "generic_proto.h"
#include "filesize_proto.h"
#include "utiltime_proto.h"
#include "MFSD_proto.h"
#include "WIFSHDF_proto.h"
#include "cdl_object_proto.h"
/*
#include	"cdl_attr_proto.h"
 */
#include "hdf_object_proto.h"
#include "getinfo_proto.h"
#include "soh_proto.h"
#include "tlm_proto.h"
#include "eng_qual_proto.h"
#include "mopsaccess_proto.h"
#include "datalevel_proto.h"
#include "dataname_proto.h"
#include "datatime_proto.h"
#include "datatype_proto.h"
#include "get_norad_orbnum_proto.h"
#include "ffmaccess_proto.h"
#include "gain_tdi_proto.h"
#include "sc_id_proto.h"
#include "tilt_proto.h"
#ifdef NEED_DATABASE
#include "dbaccess_proto.h"
#endif /* NEED_DATABASE */
