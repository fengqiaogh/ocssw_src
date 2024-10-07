/* 
 * File:   dfutils.h
 * Author: dshea
 *
 * Created on May 11, 2015, 10:44 AM
 */

#ifndef DFUTILS_H
#define DFUTILS_H

#include <hdf4utils.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    DS_HDF, // use a HDF file
    DS_NCDF // use a netCDF file
} ds_format_t;

typedef enum {
    DS_READ, // open file for reading
    DS_WRITE // open/create file for writing
} ds_access_t;

typedef struct {
    int32_t fid;
    int32_t sid;
    ds_format_t fftype;
    int32_t deflate;
} idDS;


int s2u(const char *in, char *out);
int getProdlist(const char *fname, char **prodlist, int32_t *l2_flags_type);

idDS startDS(const char *filename, ds_format_t format,
        ds_access_t accessmode, int32_t deflate);
idDS openDS(const char *filename);
int endDS(idDS ds_id);

#if 0
int createDS(
        idDS ds_id,
        const char *sname, /* short name */
        const char *lname, /* long name */
        const char *standard_name, /* NetCDF standard name (not set if passed NULL or "") */
        char *reference;
        const char *units, /* units (not set if passed NULL or "") */
        double low, /* low end of valid range */
        double high, /* high end of range (no range set if low >= high) */
        float slope, /* scale factor (not set if 0)  */
        float offset, /* scaling offset (not set if 0)  */
        int32 nt, /* HDF number type */
        int32 rank, /* number of dimensions (must be <= 3) */
        int32 d0, /* size of 1st dimension */
        int32 d1, /* size of 2nd dimension */
        int32 d2, /* size of 3rd dimension (1 if rank < 2) */
        const char *dn0, /* name of 1st dimension */
        const char *dn1, /* name of 2nd dimension (NULL if rank < 2) */
        const char *dn2 /* name of 3rd dimension (NULL if rank < 3) */
        );
#endif

int createDS(
        idDS ds_id, int sensorId,
        const char *sname, /* short name */
        int32_t dm[3], /* dimension sizes */
        const char dm_name[3][80] /* dimension names */
        );

/** create dataset using the productInfo structure */
int createDS2(idDS ds_id, /**< DS structure for netCDF/HDF file */
        const char *sname, /**< short name */
        productInfo_t* pinfo, /**< productInfo structure defining the product */
        int32_t dm[3], /**< dimension sizes */
        const char dm_name[3][80] /**< dimension names */
        );

int32_t selectDS(idDS ds_id, const char *l2_prod_names);
int32_t checkDS(idDS ds_id, const char *l2_prod_name);
int readDS(idDS ds_id, const char *name, int32_t *start, int32_t *stride,
        int32_t *count, void* data);
int writeDS(idDS ds_id, const char *name, const void* data, int32_t s0, int32_t s1,
        int32_t s2, int32_t e0, int32_t e1, int32_t e2);
int endaccessDS(idDS ds_id);

int fileInfo(idDS ds_id, int32_t *n_datasets, int32_t *n_globalattr);
int getDimsDS(idDS ds_id, const char sdsname[], int32_t dims[]);
int getTypeDS(idDS ds_id, const char sdsname[HDF4_UTILS_MAX_NAME], int32_t *dtype);

int setAttr(idDS ds_id, const char *nam, int32_t typ, int32_t cnt, const void* data);
int8_t findAttr(idDS ds_id, const char *nam);
int readAttr(idDS ds_id, const char *nam, void* data);
char* readAttrStr(idDS ds_id, const char *name);
int infoAttr(idDS ds_id, const char *nam, int32_t *dtype, int32_t *count);

int SetChrGA(idDS ds_id, const char *name, const char *value);
int SetF32GA(idDS ds_id, const char *name, float value);
int SetF64GA(idDS ds_id, const char *name, double value);
int SetI8GA(idDS ds_id, const char *name, uint8_t value);
int SetI16GA(idDS ds_id, const char *name, int16_t value);
int SetI32GA(idDS ds_id, const char *name, int32_t value);

int16_t getDataTypeInt(productInfo_t *p_info);
int16_t *float2int16(float fbuf[], int32_t spix, int32_t npix, int32_t incr,
        float slope, float offset);
int8_t *float2int8(float fbuf[], int32_t spix, int32_t npix, int32_t incr,
        float slope, float offset);
uint16_t *float2uint16(float fbuf[], int32_t spix, int32_t npix, int32_t incr,
        float slope, float offset);
uint8_t *float2uint8(float fbuf[], int32_t spix, int32_t npix, int32_t incr,
        float slope, float offset);
void *scale_sds(float *data, productInfo_t *p, int32_t npix);
float *unscale_sds(void *data, productInfo_t *p, int32_t spix, int32_t npix,
        int incr);

// netcdf specific routines
int CreateNCDF(
        idDS ds_id,
        const char *sname, /* short name */
        const char *lname, /* long name */
        const char *standard_name, /* NetCDF standard name (not set if passed NULL or "") */
        const char *reference,
        const char *comment,
        const char *units, /* units (not set if passed NULL or "") */
        double low, /* low end of valid range */
        double high, /* high end of range (no range set if low >= high) */
        float scale_factor, /* scale factor (not set if 0)  */
        float add_offset, /* scaling offset (not set if 0)  */
        int32_t fillValue, /* fill value */
        int32_t nt, /* NCDF number type */
        int32_t rank, /* number of dimensions (must be <= 3) */
        int32_t dimids[3]/* dimension ids */
        );
void nc_init_compress(int32_t nc_id, int32_t var_id, int32_t *dimids, int32_t rank,
        size_t *chunksize, int deflate_level);

#ifdef __cplusplus
}
#endif

#endif /* DFUTILS_H */

