#ifndef SRC_GET_FORMAT_GET_FORMAT_H_
#define SRC_GET_FORMAT_GET_FORMAT_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Use enum to give each input file format a unique number */
typedef enum {
    FT_INVALID = -1,
    FT_UNKNOWN,
    FT_AVIRIS,
    FT_CLASSAVHRR,
    FT_CZCSL1A,
    FT_GOCIL1B,
    FT_HICOL1B,
    FT_L1BNCDF,
    FT_L1HDF,
    FT_L1XCAL,
    FT_L2HDF,
    FT_L2NCDF,
    FT_L3BIN,
    FT_L3MAP,
    FT_MERISCC,
    FT_MERISL1B,
    FT_MERISL2,
    FT_MERISL1BSAFE,
    FT_MODISGEO, // MODIS Geolocation (hdf4)
    FT_MODISL1B, // MODIS L1B (hdf4, all bands)
    FT_MOSL1B,
    FT_OCM2L1B,
    FT_OCML1B,
    FT_OCML1BDB,
    FT_OCTSL1A,
    FT_OCTSL1ANC,
    FT_OCTSL1B,
    FT_OLCI,
    FT_OLCIGEO,
    FT_OLIL1B,
    FT_OCIA,
    FT_OCIL1B,
    FT_OCIS, // OCI simulated data, L1B
    FT_OSMIL1A,
    FT_PRISM,
    FT_SEAWIFSL1A,
    FT_SEAWIFSL1ANC, // Seawifs L1A (Netcdf)
    FT_VIIRSGEO, // VIIRS Geolocation (hdf5)
    FT_VIIRSGEONC, // VIIRS Geolocation (NetCDF4)
    FT_VIIRSL1A, // VIIRS Level-1A (NetCDF4)
    FT_VIIRSL1B, // VIIRS M-band (hdf5)
    FT_VIIRSL1BNC, // VIIRS M-band (NetCDF4)
    FT_SGLI, // SGLI files (hdf5)
    FT_L5TML1B,
    FT_L7ETML1B,
    FT_MSIL1C, /* Sentinel MSI L1C */
    FT_HAWKEYEL1A,
    FT_MISR,
    FT_SEABASSRRS,
    FT_SPEXONE,
    FT_HARP2,
    FT_HARP,
    FT_HKT,
    FT_L1C,
    FT_L1CANC
} file_type;

typedef struct file_format {
    file_type type;
    int32_t sensor_id;
    int32_t subsensor_id;
} file_format;

file_format getFormat(char *filename);
file_type   getFormatType(char *filename);

file_format chk_oli_geo(char *filename);

#ifdef __cplusplus
}
#endif // C++

#endif /* SRC_GET_FORMAT_GET_FORMAT_H_ */
