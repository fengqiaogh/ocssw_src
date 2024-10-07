#ifndef HDF_OBJECT_PROTO_H_
#define HDF_OBJECT_PROTO_H_

extern int32 create_ATTRObj
PROTO((int32 oid, ATTRObjType *attr));

extern int32 get_ATTRObj
PROTO((int32 oid, ATTRObj *attr));

extern int32 create_VGRPObj
PROTO((int32 fid, VGRPObjType *vgrp));

extern int32 create_MFSDObj
PROTO((int32 fid, MFSDObjType *mfsd));

extern uaddr find_VGRP
PROTO((char *group, char *class, uaddr MFSaddr));

extern int32 attach_MFSD_to_VGRP
PROTO((int32 fid, MFSDObjType *mfsd, VGRPObjType *vgrp));

extern int open_WIFSfile
PROTO((char *fname, int mode, int32 *fid, int32 *fid2));

extern int close_WIFSfile
PROTO((int32 fid, int32 fid2));

extern int32 create_WIFSfile
PROTO((uaddr MFSaddr));


#endif /* HDF_OBJECT_PROTO_H_ */
