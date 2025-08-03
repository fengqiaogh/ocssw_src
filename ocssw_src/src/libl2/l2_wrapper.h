#ifndef READL2SCAN_H
#define READL2SCAN_H

#include <dfutils.h>
#include <stdint.h>
#include <stdint.h>
#define MAXNUMBERPRODUCTS 2000
#define MAXNFILES 512
#ifdef __cplusplus
extern "C" {
#endif

typedef struct l2prod_struct {
    int32_t fileindex;          // Index of the current file being processed

    char filename[256];         // Name/path of the L2 file
    char oformat[32];          // Output format specification

    int32_t nrec;              // Number of scan records/lines
    int32_t nsamp;             // Number of samples per scan

    int16_t syear;             // Start year of data
    int16_t sday;              // Start day of year
    int32_t smsec;             // Start milliseconds of day
    int16_t eyear;             // End year of data
    int16_t eday;              // End day of year
    int32_t emsec;             // End milliseconds of day
    int32_t orbit;             // Orbit number
    char dtype[8];             // Data type string

    int32_t ntilts;            // Number of sensor tilt angles
    int16_t tilt_flags[20];    // Flags for each tilt angle
    int16_t tilt_ranges[2][20];// Min/max ranges for each tilt

    char *flagnames;           // Names of quality flags

    int32_t year;              // Current scan year
    int32_t day;               // Current scan day of year
    int32_t msec;              // Current scan milliseconds

    float *longitude;          // Array of longitude values
    float *latitude;           // Array of latitude values

    int32_t nprod;             // Number of data products
    char *prodname[MAXNUMBERPRODUCTS];  // Names of data products
    float **l2_data;           // 2D array of product data values
    int32_t *l2_flags;         // Quality control flags

    uint8_t *mside;            // Mirror side identifier
    uint8_t *detnum;           // Detector number
    int32_t *pixnum;           // Pixel number within scan

} l2_prod;

typedef struct meta_l2Struct {
    char *product_name; /* ATTR Product name(file name)         */
    char *title;
    char *mission;     /* ATTR mission                         */
    char *sensor_name; /* ATTR sensor name                     */
    char *sensor;      /* ATTR sensor                          */
    float northlat;    /* ATTR Northernmost latitude           */
    float southlat;    /* ATTR Southernmost latitude           */
    float westlon;     /* ATTR Westernmost longitude           */
    float eastlon;     /* ATTR Easternmost longitude           */
} meta_l2Type;

/**
 * @brief Converts between HDF and NetCDF data types
 * @param dtype Input data type (HDF format)
 * @param fileformat File format (DS_NCDF or other)
 * @return Converted data type for NetCDF
 */
int32_t get_dtype(int32_t dtype, ds_format_t fileformat);

/**
 * @brief Opens and initializes an L2 file for reading
 * @param fname Path to L2 file
 * @param plist Product list string (comma-separated list or "ALL")
 * @param l2_str Pointer to l2_prod structure to populate
 * @return 0 on success
 */
int32_t openL2(const char *fname, const char *plist, l2_prod *l2_str);

/**
 * @brief Reopens an L2 file that was previously opened
 * @param fileindex Index of file to reopen
 * @param l2_str Pointer to l2_prod structure
 * @return 0 on success
 */
int32_t reopenL2(int32_t fileindex, l2_prod *l2_str);

/**
 * @brief Reads data from an L2 file for a specific scan line
 * @param l2_str Pointer to l2_prod structure
 * @param ifile File index
 * @param recnum Record/scan number to read
 * @param iprod Product index (-1 for all products)
 * @param scan_in_rowgroup Scan row group flag
 */
int32_t readL2(l2_prod *l2_str, int32_t ifile, int32_t recnum, int32_t iprod,
               unsigned char *scan_in_rowgroup);


/**
 * @brief Reads latitude and longitude to l2_prod structure
 * @param l2_str Pointer to l2_prod structure
 * @param start start index of slice
 * @param edges size of the slise
 * @param ifile ifile number
 * @param scan_in_rowgroup Scan row group flag
 */               
int32_t readlonlat(l2_prod *l2_str, int32_t ifile, int32_t *start, int32_t *edges,
                   unsigned char *scan_in_rowgroup);

/**
 * @brief Closes an open L2 file
 * @param l2_str Pointer to l2_prod structure
 * @param ifile File index
 * @return 0 on success
 */
int32_t closeL2(l2_prod *l2_str, int32_t ifile);

/**
 * @brief Finds the index of a product in the L2 file
 * @param l2_str Pointer to l2_prod structure
 * @param prodname Name of product to find
 * @return Product index if found, -1 if not found
 */
int32_t findprod(l2_prod *l2_str, char *prodname);

/**
 * @brief Reads metadata from an L2 file
 * @param meta_l2 Pointer to meta_l2Type structure to populate
 * @param ifile File index
 * @return 0 on success
 */
int32_t readL2meta(meta_l2Type *meta_l2, int32_t ifile);

/**
 * @brief Frees memory allocated for L2 metadata structures
 * @param meta_l2 Pointer to meta_l2Type structure
 * @return 0 on success
 */
int32_t freeL2meta(meta_l2Type *meta_l2);

/**
 * @brief Frees memory allocated for L2_prod structures
 * @param l2_str Pointer to L2_prod structure
 * @return 0 on success
 */
int32_t freeL2(l2_prod *l2_str);

/**
 * @brief Gets the units for a specified L3B product
 * @param l2_str Pointer to l2_prod structure
 * @param ifile File index
 * @param l3b_prodname Name of L3B product
 * @param units String to store units
 * @return 0 on success, -1 if product not found
 */
int32_t getL3units(l2_prod *l2_str, int32_t ifile, char *l3b_prodname, char *units);

#ifdef __cplusplus
}
#endif
#endif
