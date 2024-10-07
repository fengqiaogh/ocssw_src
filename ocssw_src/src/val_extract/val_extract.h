/** @file val_extract.h
    @brief Process a small section of a Level-2 NetCDF file.
*/

#ifndef ___VAL_EXTRACT_H_
#define ___VAL_EXTRACT_H_

#include "argpar.h"
#include "olog.h"
#include "shash.h"

#include <netcdf.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#ifndef VALEXTRACT_API_VERSION
/** @brief Current API version.  If already defined, use API defined first. */
#define VALEXTRACT_API_VERSION 2007005
/** @brief Current API version as a string. */
#define VALEXTRACT_API_VERSION_STR "2.7.5"
#endif

/** @brief Returned on successful processing. */
#define VALEXTRACT_ERR_NONE EXIT_SUCCESS
/** @brief Returned when the desired point is not in the file boundaries. */
#define VALEXTRACT_ERR_POINT_NOT_FOUND 99
/** @brief Returned for unexpected errors like malloc failures or, possibly, permissions
    problems and the like. */
#define VALEXTRACT_ERR_UNKNOWN 100
/** @brief Returned when the NetCDF file can't be opened (due to errors or corruption). */
#define VALEXTRACT_ERR_NCFILE_ERR 101
/** @brief Returned when the NetCDF file isn't in the format expected (not an L2, etc). */
#define VALEXTRACT_ERR_NCFILE_INVALID 102
/** @brief Returned when the something goes wrong processing l2_flags. */
#define VALEXTRACT_ERR_FLAG 103
/** @brief Returned when the something goes wrong processing a product or when finding a
    product specified on the command line. */
#define VALEXTRACT_ERR_VARIABLE 104
/** @brief Returned when given bad or no arguments. */
#define VALEXTRACT_ERR_INPUT 105
/** @brief Not a real error, but returned when the L2QC step says it's a poor quality file. */
#define VALEXTRACT_ERR_L2QC 106

/** @brief argpar structure used for making programs that inherit options from this library. */
extern argpar val_extract_argpar;

/** @brief Passed to parser when beginning the processing. */
#define VALEXTRACT_KEY_INIT    1
/** @brief Passed to parser when first opening the NetCDF file.  The parser is passed a pointer
    to an nc_region structure as input. */
#define VALEXTRACT_KEY_FILE    2
/** @brief Passed to parser when a variable is processed.  The parser is passed a pointer
    to an nc_var structure as input. */
#define VALEXTRACT_KEY_VAR     3
/** @brief Passed to parser when the processing is finished and no errors were encountered. */
#define VALEXTRACT_KEY_SUCCESS 4
/** @brief Passed to parser when the processing is finished and an error was encountered. */
#define VALEXTRACT_KEY_ERROR   5
/** @brief Passed to parser when ending the processing. */
#define VALEXTRACT_KEY_FINI    6

/** @brief Returned from a parser to stop processing. */
#define VALEXTRACT_CMD_STOP    2

/** @brief Initial value for location arguments (lat/lon/line/pixl) before argument parsing. */
#define VALEXTRACT_UNSET -255

/** @brief Pointer to a callback function to call for each argument parsed.

    @param[in] key One of the VALEXTRACT_KEY_ macros.
    @param[in] nc_input Either NULL, an nc_region or nc_var pointer, depending on they key. All
    pointers should persist until after the VALEXTRACT_KEY_FINI call.
    @param[in] user_input Input pointer from the val_extract_arguments structure.

    @return 0 on success, 1 on error.
*/
typedef int (*val_extract_parser)(int key, void *nc_input, void *user_input);

/** @typedef typedef struct val_extract_valid_range
    @brief For the lazy. */
/** @struct val_extract_valid_range
    @brief Given to val_extract_arguments, these contain the valid ranges, inclusive, of variables.
    Any values outside the range(s) will not be used for statistics.
*/
typedef struct val_extract_valid_range {
    char *prefix;
    double min, max;
} val_extract_valid_range;

/** @typedef typedef struct val_extract_arguments
    @brief For the lazy. */
/** @struct val_extract_arguments
    @brief Passed into the library function to control the processing.  Many of the fields will
    be unspecified.

    For the box definition, it accepts either 1) box_size (box_size OR box_size_km) and a center
    point of start_line and start_pixel or start_lat and start_lon or 2) a region defined by
    start_ and end_ line and pixel or start_ and end_ lat and lon or 3) nothing, and the whole
    file is processed.
*/
typedef struct val_extract_arguments {
    /** @brief Path to Level 2 NetCDF file. */
    const char *ifile;

    /** @brief Box width (and height) for extract, in pixels. */
    int box_size;

    /** @brief Box radius for extract, in kilometers. */
    double box_size_km;

    /** @brief Whether to use lines and pixels instead of latitudes and longitudes. */
    bool line_and_pixel;
    /** @brief Either center line for a box extract or the beginning line for regions. */
    int start_line;
    /** @brief End line for regions. */
    int end_line;
    /** @brief Either center pixel for a box extract or the beginning pixel for regions. */
    int start_pixel;
    /** @brief End pixel for regions. */
    int end_pixel;

    /** @brief Whether to use latitudes and longitudes instead of lines and pixels. */
    bool lat_and_lon;
    /** @brief Either center latitude for a box extract or the beginning latitude for regions. */
    double start_lat;
    /** @brief End latitude for regions. */
    double end_lat;
    /** @brief Either center longitude for a box extract or the beginning longitude for regions. */
    double start_lon;
    /** @brief End longitude for regions. */
    double end_lon;

    /** @brief Only process geospatial variables, those with lines and pixels as the only dimensions. */
    bool geospatial_only;
    /** @brief Read all geospatial variables as doubles, regardless of original type, useful for
        avoiding copy-and-pasting the processing code. */
    bool geospatial_to_double;

    /** @brief If given, only process the NetCDF variables listed. */
    char **products;
    /** @brief How many NetCDF variables should be processed. */
    int product_count;

    /** @brief If given, ignore all pixels with ANY of these flag bits sets. */
    uint32_t ignore_flags_mask;
    /** @brief If given, ignore all pixels with ANY of these flags.  This is added to ignore_flags_mask. */
    char **ignore_flags;
    /** @brief How many ignore_flags are present. */
    int ignore_flag_count;

    /** @brief If given, count all pixels with ANY of these flag bits sets. */
    uint32_t count_flags_mask;
    /** @brief If given, count all pixels with ANY of these flags.  This is added to count_flags_mask. */
    char **count_flags;
    /** @brief How many count_flags are present. */
    int count_flag_count;

    /** @brief If given, use these flags for the L2QC process, paired with the l2qc_thresholds input. */
    char **l2qc_flags;
    /** @brief If given, use these thresholds for the L2QC process, aborting processing if any are
        thresholds are surpassed. */
    double *l2qc_thresholds;
    /** @brief How many l2qc_thresholds are present. */
    int l2qc_threshold_count;

    /** @brief If given, use these to filter variables' data used for statistics. */
    val_extract_valid_range *valid_ranges;
    /** @brief How many ignore_flags are present. */
    int valid_range_count;

    /** @brief If given, override default list of products where statistics are not meaningful (such as
        l2_flags). */
    char **skip_stats;
    /** @brief How many skip_stats products are present. */
    int skip_stats_count;

    /** @brief Create .global file containing all global attributes. */
    bool global_atts;

    /** @brief Include all variables attributes in output files. */
    bool variable_atts;

    /** @brief olog object to use for messages, defaults to global logger. */
    olog *log;

    /** @brief If given, OPTICS algorithm threshold for preventing reaching over data gaps. Lower = stricter. */
    double optics_threshold;

    /** @brief Required, this is the function that is given the nc_region and nc_var structures after
        they are processed. */
    val_extract_parser val_extract_parser;
    /** @brief If given, this is passed to the val_extract_parser for arbitrary use by the user. */
    void *user_input;
} val_extract_arguments;

/** @brief Process a small section of a Level-2 NetCDF file.

    @param[in] arguments val_extract_argument structure.

    @return 0 on success, one of VALEXTRACT_ERR_ on error.
*/
int val_extract(val_extract_arguments *arguments);

/** @brief Clean up stuff malloc'd by the argpar callback.  Should be called at the end of processing
    if val_extract_argpar is used. */
int val_extract_clean(val_extract_arguments *arguments);

/** @brief Clear an argument struct to make it suitable for using the library without using argpar.

    This is sometimes necessary because a freshly malloc'd struct can have valid values for arguments
    that a user may not want; for example, if start and end coords aren't explicitly unset, omitting
    them from the arguments won't grab the entire file, as would be expected.

    @param[in] arguments val_extract_argument structure.
*/
void val_extract_clear_args(val_extract_arguments *arguments);

/** @brief Given dimension lengths from an nc_file and a one dimensional index, find the corresponding
    n-dimension index. */
void unflatten_index(int index, int ndims, const int *dims, int *result);

/** @brief Returns a string representation of the library's version number, without a label. */
const char *val_extract_version();

/** @brief Returns a string representation of the library's implemented API version number, without a
    label. */
const char *val_extract_api_version();

/* The following should go into an nc_read library */

typedef struct nc_point {
    size_t line, pixel;
} nc_point;

typedef struct nc_box {
    size_t start[2], count[2];
} nc_box;

typedef struct nc_file {
    const char *file_path;
    int ndims, *dim_lengths;
    int ngrps, *gids;
    int line_dimid, pixel_dimid, pixel_control_dimid, control_point_size;
    int flag_count;
    uint32_t *flag_masks;
    char *flag_string;
    char **flag_meanings;
    shash *attributes;
} nc_file;

typedef struct nc_region {
    nc_file *file;
    nc_point center;
    int box_count;
    nc_box *boxes;
    int pixel_count, unflagged_pixel_count, flagged_pixel_count;
    int *flag_counts;
    uint32_t *pixel_flags;
    double lat, lon;
    time_t utime;
    char ascii_time[24];
} nc_region;

typedef struct nc_var_stats {
    bool initialized;
    int count;
    double min, max, mean, median, stddev, rms;
    double center_value;
} nc_var_stats;

typedef struct nc_var {
    nc_file *file;
    nc_region *region;
    char *name, *group_name;
    int gid, varid;
    bool is_geospatial, uses_control_points;
    bool has_scale, has_offset, has_fill;
    double scale, offset, fill;
    int ndims, *dim_ids;
    shash *attributes;
    nc_var_stats stats, filtered_stats, iqr_stats;
    nc_type data_type;
    size_t data_sizeof;
    double *valid_data;
    int valid_data_count;
    void *data;
} nc_var;

size_t nc_get_region_size(nc_region *region);
size_t nc_get_region_dim_size(nc_region *region, int dim_index);

#define declare_nc_get_varr(type, suffix) \
    int nc_get_varr ## suffix(int ncid, int varid, nc_region *region, type *p); \
    int nc_get_varrd ## suffix(int ncid, int varid, nc_region *region, int dim_index, type *p);

declare_nc_get_varr(char, _text)
declare_nc_get_varr(unsigned char, _uchar)
declare_nc_get_varr(signed char, _schar)
declare_nc_get_varr(short, _short)
declare_nc_get_varr(int, _int)
declare_nc_get_varr(long, _long)
declare_nc_get_varr(float, _float)
declare_nc_get_varr(double, _double)
declare_nc_get_varr(unsigned short, _ushort)
declare_nc_get_varr(unsigned int, _uint)
declare_nc_get_varr(long long, _longlong)
declare_nc_get_varr(unsigned long long, _ulonglong)
declare_nc_get_varr(char *, _string)
declare_nc_get_varr(void,)

/* The following should go into a math library */

#define declare_compare(type, suffix) int compare ## suffix(const void *a, const void *b)
declare_compare(int8_t, _int8);
declare_compare(uint8_t, _uint8);
declare_compare(int16_t, _int16);
declare_compare(uint16_t, _uint16);
declare_compare(int32_t, _int32);
declare_compare(uint32_t, _uint32);
declare_compare(int64_t, _int64);
declare_compare(uint64_t, _uint64);
declare_compare(float, _float);
declare_compare(double, _double);

#endif /* ___VAL_EXTRACT_H_ */
