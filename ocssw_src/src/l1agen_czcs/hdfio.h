/*******************************************************************

   hdfio.h

   purpose: include file for the use of hdf I/O routines

   Parameters: none

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       15-Feb-1993     Original development
      W. Robinson, SAIC, 28 Jan 2004    add include of mfhdf.h
      W. Robinson, SAIC, 5 Aug 2004     updated for generic use

 *******************************************************************/

/*
 *  Note that hdf.h is needed for this and navigation defs
 */
#include "hdf.h"
#include "mfhdf.h"

/*
 *  the l1info_struct structure has file ids...
 */
struct hdfio_struc_d {
    int32 fid;
    int32 sdfid;
};

typedef struct hdfio_struc_d hdfio_struc;

/*
 *  prototypes
 */
int hdfio_rd_sd(hdfio_struc, char *, void *);
int hdfio_open(char *, hdfio_struc *);
int hdfio_rd_gattr(hdfio_struc, char *, int32, int32, void *);
void hdfio_close(hdfio_struc);
