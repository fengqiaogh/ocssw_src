//#include <VL1_Params.h>
//#include <VcstViirsCal.h>

#ifdef __cplusplus
extern "C" {
#endif

void VcstViirsCal_initialize(char *l1aFilename, char *l1aCalParFilename);
void VcstViirsCal_calibrateMOD(int iScan, int nbands, float **l1bptrs);

#ifdef __cplusplus
}
#endif
