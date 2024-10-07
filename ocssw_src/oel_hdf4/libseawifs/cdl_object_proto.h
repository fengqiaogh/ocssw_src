#ifndef CDL_OBJECT_PROTO_H_
#define CDL_OBJECT_PROTO_H_

extern int wifstype_atoi
PROTO((char *wifstypestr));

extern char *wifstype_itoa
PROTO((int wifstype));

extern int list_attrobj
PROTO((ATTRObjType *attr));

extern int list_vgrpobj
PROTO((VGRPObjType *vgrp));

extern int list_mfsdobj
PROTO((MFSDObjType *mfsd));

extern int list_sfsdobj
PROTO((SFSDObjType *sfsd));

extern int list_vsetobj
PROTO((VSETObjType *vset));

extern int list_dfanobj
PROTO((DFANObjType *dfan));

extern int list_CDL
PROTO((uint32 MFSaddr));


#endif /* CDL_OBJECT_PROTO_H_ */
