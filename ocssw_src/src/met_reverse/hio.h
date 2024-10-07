/*******************************************************************

   hio.h

   purpose: include file for the use of the hdf I/O routine aids

   Parameters: 
      Type              Name            I/O     Description
      ----              ----            ---     -----------

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 30-Sep-2005     rom l1io.h

 *******************************************************************/

/*
 *  Note that hdf.h is needed for this and navigation defs
 */
#include "hdf.h"
#include "mfhdf.h"

/*
 *  the hio_struct structure has file ids...
 */
typedef struct hio_struct_d {
    int32 fid;
    int32 sdfid;
} hio_struct;

/*
 *  prototypes
 */
int32 hio_open(char *, int32, hio_struct*);
int32 hio_w_sds(hio_struct, int, int32 *, int32, char *, void *);
int32 hio_close(hio_struct);
int32 hio_r_sds(hio_struct, char *, int32, void *);
int32 hio_i_sds(hio_struct, char *, int32, void *);
int32 hio_r_sds_sl(hio_struct, char *, int32, int32 *, int32 *, void *);
int32 hio_r_sds_sls(hio_struct, char *, int32, int32 *, int32 *, int32 *,
        void *);
int32 hio_rg_attr(hio_struct, char *, int32, int32, void *);
