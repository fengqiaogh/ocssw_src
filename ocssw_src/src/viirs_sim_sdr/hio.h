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

#include <stdint.h>

/*
 *  the hio_struct structure has file ids...
 */
typedef struct hio_struct_d {
    int32_t fid;
    int32_t sdfid;
} hio_struct;

/*
 *  prototypes
 */
int32_t hio_open(char *, int32_t, hio_struct*);
//int32_t hio_w_sds(hio_struct, int, int32_t *, int32_t, char *, void *);
int32_t hio_close(hio_struct);
int32_t hio_r_sds(hio_struct, char *, int32_t, void *);
//int32_t hio_r_sds_sl(hio_struct, char *, int32_t, int32_t *, int32_t *, void *);
//int32_t hio_r_sds_sls(hio_struct, char *, int32_t, int32_t *, int32_t *, int32_t *,
//        void *);
//int32_t hio_rg_attr(hio_struct, char *, int32_t, int32_t, void *);
