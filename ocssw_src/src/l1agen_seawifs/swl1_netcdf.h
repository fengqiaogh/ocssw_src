#ifndef _SWL1_NETCDF_H_
#define _SWL1_NETCDF_H_

#include "swl0_struc.h"
#include "swl1_struc.h"
#include "mfhdf.h"
#include "passthebuck.h"

#ifdef __cplusplus
extern "C" {
#endif


char * L1aFilename_netcdf(swl0ctl *l0ctl, double, unsigned char);
char * DataTypeString(swl0ctl *l0ctl, unsigned char);
char * DTypeString(unsigned char);
int CreateL1aFile_netcdf(char *, swl0scene *, char *, char *, swl0ctl *);
int CloseL1aFile_netcdf(l1met *);
int CreateScanData(int32 ns, int32 np);
int WriteScanData_netcdf(int32, swl1rec *);
void DecomposeTime(double, int16 *, int16 *, int32 *);
int AddCalData(void);
int AddTiltData(int32, int16 f[20], int16 r[20][2],
        float32 lat[20][2][2], float32 lon[20][2][2]);
int MakeVgroups(void);



#ifdef __cplusplus
}
#endif


#endif /* _SWL1_NETCDF_H_ */
