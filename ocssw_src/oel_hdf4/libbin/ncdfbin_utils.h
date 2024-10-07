#ifndef  NETCDFBINUTILS_H
#define  NETCDFBINUTILS_H

#ifdef __cplusplus
extern "C" {
#endif
#pragma GCC diagnostic ignored "-Wpadded"

typedef struct {
    uint32_t binnum;
    int16_t nobs;
    int16_t nscenes;
    float weights;
    float time_rec;
} binListStruct_nc;

typedef struct {
    uint64_t binnum;
    int16_t nobs;
    int16_t nscenes;
    float weights;
    float time_rec;
} binListStruct64_nc;

typedef struct {
    uint32_t start_num;
    uint32_t begin;
    uint32_t extent;
    uint32_t max;
} binIndexStruct_nc;

typedef struct {
    uint64_t start_num;
    uint64_t begin;
    uint32_t extent;
    uint64_t max;
} binIndexStruct64_nc;


int defineBinList_nc(int32_t deflate, int32_t grpid);
int defineBinList64_nc(int32_t deflate, int32_t grpid);
int writeBinList_nc(int32_t grpid, int32_t nbins_to_write, const void *data);
int defineBinData_nc(int32_t deflate, int32_t grpid, int32_t nprods, char** prodnames);
int writeBinData_nc(int32_t grpid, int32_t nbins_to_write, int32_t iprod, const void *data);
int defineBinIndex_nc(int32_t deflate, int32_t grpid);
int defineBinIndex64_nc(int32_t deflate, int32_t grpid);
int writeBinIndex_nc(int32_t grpid, int32_t n_write, const void *data);
int defineQuality_nc(int32_t deflate, int32_t grpid);
int writeQuality_nc(int32_t grpid, int32_t nbins_to_write, const void *data);

#ifdef __cplusplus
}
#endif

#endif
