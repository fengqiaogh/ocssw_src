/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/hdf/cdl_object.h,v 4.12 1995/02/24 15:34:52 seawifsd Exp seawifsd $
$Log: cdl_object.h,v $
Revision 4.12  1995/02/24 15:34:52  seawifsd
added support of new object DFANObj.

Revision 4.11  1995/02/21 16:46:08  seawifsd
implemented CAL_OFFSET feature.

Revision 4.10  1995/01/17 19:58:15  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:14:31  seawifsd
Jan. 9, 1994, 4.0

Revision 3.4  1994/12/23 19:00:17  seawifsd
commented out un-used header file cdl_attr.h.
deleted commented out #ifdef/#endif that was calling compiling errors in
some unknown situations when -cckr was used. No problem for -ansi.

Revision 3.3  1994/11/08 18:46:36  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:34  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.2  1994/10/04 15:51:20  frank
added code to support MFSD_MINMAX feature for future use.
added code to support REC_FLAG and NSAMP_DIMS features in the future so
that generic read routine can find out whether the SDS is a 'record'
variable or contains 'nsamp' dimension.

Revision 1.1.1.1  1994/05/23 14:52:04  frank
changed 'int' to 'int32' on MFSDObjStruct so that compiling with ANSI
will be correct.

Revision 1.2  1994/05/10 18:49:16  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:33:00  seawifst
Initial revision


 */


#ifndef _CDL_OBJECT_H_
#define _CDL_OBJECT_H_
/*
#include	"cdl_attr.h"
 */
#ifndef MAXDIMS
#define MAXDIMS 3
#endif

#define WIFSATTRTYPE 1
#define WIFSVGRPTYPE 2
#define WIFSMFSDTYPE 3
#define WIFSSFSDTYPE 4
#define WIFSVSETTYPE 5
#define WIFSDFANTYPE 6
#define WIFSATTRTYPESTR "Attribute"
#define WIFSVGRPTYPESTR "Vgroup"
#define WIFSMFSDTYPESTR "MF_SDS"
#define WIFSSFSDTYPESTR "SF_SDS"
#define WIFSVSETTYPESTR "Vset"
#define WIFSDFANTYPESTR "DFAN"

/* this structure must match the beginning part of the SeaWiFS data	*/

/* file structure. For example: level_1a_gac_datStruct			*/
typedef struct cdl_hdrStruct {
    int n; /* total number of records		*/
    int fid; /* HDF file id for SD routines		*/
    /* return by SDstart()			*/
    int fid2; /* HDF file id for all other routines	*/
    /* return by Hopen()			*/
    char *fname; /* HDF data file name			*/
    int ngattr; /* number of global attributes		*/
    int nvgrp; /* number of VGRP object		*/
    int nmfsd; /* number of MFSD object		*/
    int nsfsd; /* number of SFSD object		*/
    int nvset; /* number of VSET object		*/
    int ndfan; /* number of DFAN object		*/
} cdl_hdrType;

/* CDL variable attribute Object_type = "Attribute"			*/
/* this type will be implemented as HDF multifile SDS global attribute	*/
/* Also, the attributes for MFSD variables will also be implemented in	*/

/* this structure.							*/
typedef struct ATTRObjStruct {
    int idx; /* HDF Attribute index				*/
    char *name; /* CDL variable name, not in HDF data		*/
    int DFNT; /* HDF global attribute data type		*/
    int count; /* dimension, number of data item		*/
    char *label; /* name of HDF global attribute, CDL attribute	*/
    /* 'long_name'					*/
    void *data; /* value of HDF global attribute		*/
} ATTRObjType;

typedef struct SFSDObjStruct {
    int ref;
} SFSDObjType;

typedef struct MFSDObjStruct {
    int32 sid; /* SDS id			*/
    int ref; /* SDS reference number		*/
    int idx; /* HDF SDS index		*/
    char *name; /* name of the SDS		*/
    int DFNT; /* data type			*/
    int32 rank; /* rank				*/
    int32 dims[MAXDIMS]; /* dimension			*/
    char *dimname[MAXDIMS]; /* name of each dimension	*/
    int32 dimid[MAXDIMS]; /* dimension id			*/
    void *data; /* data pointer			*/
    char *label; /* label(long_name) of SDS	*/
    /* redundant of attr 'long_name'*/
    unsigned char rec_flag; /* whether it's record var T/F	*/
    int32 nsamp_dims[MAXDIMS]; /* link-list pointed to the	*/
    /* dimension that is 'nsamp'	*/
    void *valid_max;
    void *valid_min;
    void *max;
    void *min;

#ifndef NO_CAL_OFFSET
    /* the relationship between the actual value 'y' and the value	*/
    /* 'x' that is stored in a data set is defined as:		*/
    /*								*/
    /* y = cal*(x - cal_offset)					*/
    /*								*/
    /* which is not the same as slope/intercept defined in SeaWiFS:	*/
    /*								*/
    /* y = slope * x + intercept					*/
    /*								*/
    /* Thus								*/
    /*								*/
    /*	cal = slope ; cal_offset = -(intercept/slope)		*/
    /*								*/
    /* if uncal_num_type is defined as DFNT_NONE, it means that	*/
    /* there is no slope/intercept(cal/cal_offset) been defined	*/
    /* for this SDS							*/
    int32 uncal_num_type; /* number type of uncalibrated data */
    float64 cal; /* calibration factor		*/
    float64 cal_err; /* calibration error		*/
    float64 cal_offset; /* uncalibrated offset		*/
    float64 cal_offset_err; /* uncalibrated offset error	*/
#endif /* !NO_CAL_OFFSET */

    char *group; /* group name for this SDS	*/
    /* value from corresponding VGRP*/
    char *class; /* class name for this SDS	*/
    /* value from corresponding VGRP*/
    int nattr; /* total number of attributes	*/
    ATTRObjType *attrs; /* list of attribute pointers	*/
} MFSDObjType;


/* CDL variable attribute Object_type = "Vgroup"			*/

/* this will be implemented in HDF Vgroup				*/
typedef struct VGRPObjStruct {
    int gid; /* Vgroup ID					*/
    char *name; /* variable name, not in HDF			*/
    char *group; /* Vgroup name					*/
    char *class; /* Vgroup class					*/
} VGRPObjType;

typedef struct VSETObjStruct {
    int ref;
} VSETObjType;

typedef struct DFANObjStruct {
    int dfantype; /* DFTAG_FID,DFTAG_FD,DFTAG_DIL,DFTAG_DIA*/
    int32 tag; /* target tag(DFTAG_DIL,DFTAG_DIA)	*/
    int32 ref; /* target ref(DFTAG_DIL,DFTAG_DIA)	*/
    int seq; /* sequence number(DFTAG_FID,DFTAG_FD)	*/
    char *data;
    int len; /* data size(DFTAG_FD,DFTAG_DIA)	*/
} DFANObjType;

typedef MFSDObjType MFSDObj;
typedef ATTRObjType ATTRObj;
typedef VGRPObjType VGRPObj;
typedef SFSDObjType SFSDObj;
typedef VSETObjType VSETObj;
typedef DFANObjType DFANObj;

typedef struct ObjIndexStruct {
    int objtype; /* HDF Object type,WIFSATTRTYPE,...	*/
    int objindex; /* 0,1,2,...				*/
} ObjIndexType;

typedef ObjIndexType Objidx;


//#include "usrhdr.h"
//#include "usrmac.h"
#include "cdl_object_proto.h"
#endif /* _CDL_OBJECT_H_ */
