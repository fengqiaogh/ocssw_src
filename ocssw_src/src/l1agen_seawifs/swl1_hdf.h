#ifndef _SWL1_HDF_H_
#define _SWL1_HDF_H_

#include <stdio.h>
#include "swl0_struc.h"
#include "swl1_struc.h"
#include "mfhdf.h"
#include "passthebuck.h"


#define MALLOC(ptr,typ,num) {    \
  (ptr) = (typ *)malloc((num) * sizeof(typ));  \
  if((ptr) == NULL){     \
    fprintf(stderr,     \
    "-E- %s line %d: Memory allocation failure.\n", \
    __FILE__,__LINE__);     \
    return(MEMORY_ALLOCATION_ERROR);   \
  }       \
}

#define CALLOC(ptr,typ,num) {                                           \
  (ptr) = (typ *)calloc((num) , sizeof(typ));                           \
  if((ptr) == NULL){                                                    \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n",      \
    __FILE__,__LINE__);                                                 \
    return(MEMORY_ALLOCATION_ERROR);                                   \
  }                                                                     \
}

char * L1aFilename(swl0ctl *l0ctl, double, unsigned char);
char * DataTypeString(swl0ctl *l0ctl, unsigned char);
char * DTypeString(unsigned char);
int CreateL1aFile(char *, swl0scene *, char *, char *, swl0ctl *);
int CloseL1aFile(l1met *);
int CreateScanData(int32 ns, int32 np);
int WriteScanData(int32, swl1rec *);
void DecomposeTime(double, int16 *, int16 *, int32 *);
int AddCalData(void);
int AddTiltData(int32, int16 f[20], int16 r[20][2],
        float32 lat[20][2][2], float32 lon[20][2][2]);
int MakeVgroups(void);

#endif /* _SWL1_HDF_H_ */
