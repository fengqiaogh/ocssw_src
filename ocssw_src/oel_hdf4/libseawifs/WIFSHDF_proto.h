#ifndef WIFSHDF_PROTO_H_
#define WIFSHDF_PROTO_H_

extern char *access_str
PROTO((int queryaccess));

extern int VGnumber
PROTO((int fid2));

extern int VGnametoid
PROTO((int fid2, char *name, char *class));

extern int VGinqtagref
PROTO((int fid2, int gid, int32 tag, int32 ref));

extern int VGidtoname
PROTO((int fid2, int gid, char *name, char *class));

extern int VGfileinfo
PROTO((int fid2, int *nvgrp, int *gidarray));

extern int get_VGRPObjInfo
PROTO((int fid2, VGRPObj *VGRP, int queryaccess));

extern int free_VGRPObj
PROTO((VGRPObj *VGRP));

extern int get_ATTRObjInfo
PROTO((int32 oid, ATTRObj *ATTR, int queryaccess));

extern int free_ATTRObj
PROTO((ATTRObj *ATTR));

extern int get_MFSDObjInfo
PROTO((int32 fid, MFSDObj *MFSD, int queryaccess));

extern int free_MFSDObj
PROTO((MFSDObj *MFSD));

extern int is_valid_gid
PROTO((int fid2, int gid));

extern int is_valid_VGRPObj
PROTO((int fid2, VGRPObj *VGRP));

extern int is_valid_ATTRObj
PROTO((int32 oid, ATTRObj *ATTR, int nattr));

extern int is_valid_MFSDObj
PROTO((int fid, MFSDObj *MFSD));

extern int get_ATTRObjDatasize
PROTO((ATTRObj *ATTR));

extern int get_ATTRObjData
PROTO((int oid, ATTRObj *ATTR));

extern int get_MFSDObjDataitem
PROTO((int fid, MFSDObj *MFSD));

extern int get_MFSDObjDatasize
PROTO((int fid, MFSDObj *MFSD));

extern int get_MFSDObjDataRecitem
PROTO((int fid, MFSDObj *MFSD));

extern int get_MFSDObjDataRecsize
PROTO((int fid, MFSDObj *MFSD));

extern int allocMFSDObjData
PROTO((int fid, MFSDObj *MFSD));

extern int allocMFSDObjDataRec
PROTO((int fid, MFSDObj *MFSD, int nrec));

extern int get_MFSDObjData
PROTO((int32 fid, MFSDObj *MFSD));

extern int get_MFSDObjDataRec
PROTO((int fid, MFSDObj *MFSD, int roff, int rcount, void *data));

extern int get_MFSDObjDataSlice
PROTO((int32 fid, MFSDObj *MFSD,
        int32 start[], int32 stride[], int32 count[], void *data));

extern int build_objindex
PROTO((SeaWiFS_HDFType *WIFS, char *labels[], char *objtype[]));

extern int get_WIFSinfo
PROTO((SeaWiFS_HDFType *WIFS, char *labels[], char *objtype[]));

extern int build_hdf
PROTO((char *fname, SeaWiFS_HDFType *WIFS));


#endif /* WIFSHDF_PROTO_H_ */
