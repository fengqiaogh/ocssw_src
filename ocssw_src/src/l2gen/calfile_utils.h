
#ifndef CALFILE_UTILS_H
#define CALFILE_UTILS_H
extern float ***aots;
extern float ***nlws;

#ifdef __cplusplus
extern "C" {
#endif

extern int16_t fileID; // updated by file appender function...
//arrays to contain the values needed to avg over for inversion


#include "l1.h"
#include "aer_struc.h"
#include "target_struc.h"
#include <dfutils.h>

typedef enum {
    DET2DET, CROSSCAL, BINMATCH
} caltype;

typedef struct cal_struct {
    int32_t sensorID;
    int32_t year;
    int32_t day;
    int32_t msec;
    int16_t iscan;
    int32_t pixnum;
    uint8_t detnum;
    /* move to vars1D
    float    lon;
    float    lat;
     */
    uint8_t mside;
    float *vars1D; //vars1D[nvars1d]
    float **Lt; // Lt[nband][ndets]
    float **vLt; // vLt[nband][ndets]
    float **data; //data[nproducts][ndets]
} calstr;

idDS calfile_open(char *ofile, int sensorID, int ydim, int xdim, int nprods, int nvars1d, char l2prods[L1_MAXPROD][32],
        char vars1Dnames[L1_MAXPROD][32], long* numExistingRuns, caltype ctype);
int calfile_close(idDS ds_id);
int calfile_create(char *ofile, idDS *ds_id, int sensorID, int ydim, int xdim, int nprods, int nvars1d,
        char l2prods[L1_MAXPROD][32], char vars1Dnames[L1_MAXPROD][32], caltype ctype);
int calfile_write(idDS ds_id, calstr *calrec, int recnum, int ydim, int xdim, int nprods, int nbands, int nvars1d,
        char l2prods[L1_MAXPROD][32], char vars1Dnames[L1_MAXPROD][32], caltype ctype);
void inversion_init(long xdim, long iscan, int nbands, long ipix, aestr *aerec, tgstr *tgrec);

calstr* alloc_calrec(int xdim, int nbands, int nprods, int nvar1d);
void free_calrec(calstr *calrec, int nbands, int nprods);

#ifdef __cplusplus
}
#endif
#endif
