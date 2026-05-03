#ifndef bin_util_h
#define bin_util_h

#pragma GCC diagnostic ignored "-Wpadded"
#ifndef PI
#define PI  3.141592653589793
#endif
#include <limits>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include "hdf5.h"
#include "L3Shape.h"
#include "genutils.h"

#define MAXNPROD 1024
#define MAXNVDATA MAXNPROD+3
#define MAXNUMATTR 16
#define MAXATTRLEN 32
#define MAXATTRNAMELEN 32

#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )

enum AveragingScheme { ARITHMETIC_MEAN, GEOMETRIC_MEAN, HARMONIC_MEAN, QUADRATIC_MEAN, UNDEFINED };
// product attributes read from l3bin file
class ProductL3Attributes {
   public:
    std::string averaging_scheme_name{};
    AveragingScheme averaging_scheme{UNDEFINED};
    float valid_max = std::numeric_limits<float>::max();
    float valid_min = std::numeric_limits<float>::min();
    float wavelength{BAD_FLT};
};
/**
 * @brief Get the averaging scheme from a string representation
 * @param averaging_scheme The string representation of the averaging scheme (e.g., "arithmetic", "geometric", "harmonic", "quadratic")
 * @return The corresponding AveragingScheme enum value
 */
AveragingScheme get_averaging_scheme_from_string(const std::string& averaging_scheme);

/**
 *  @brief Applies the specified averaging scheme to a value, reverse
 * @param val The value to which the averaging scheme is applied
 * @param method The averaging scheme to apply (arithmetic, geometric, harmonic and quadratic)
 * @return
 */
template <typename T>
T apply_reverse_averaging_scheme(T val, AveragingScheme method) {
    if (val == BAD_FLT)
        return val;
    switch (method) {
        case ARITHMETIC_MEAN:
            return val;
        case GEOMETRIC_MEAN:
            return std::exp(val);
        case HARMONIC_MEAN:
            if (val == 0)
                return BAD_FLT;
            else
                return 1.0 / val;
        case QUADRATIC_MEAN:
            if (val < 0)
                return BAD_FLT;
            return std::sqrt(val);
        default:
            std::cerr << "-E-: unidentified averaging scheme. See " << __FILE__ << ":" << __LINE__ << std::endl;
            exit(EXIT_FAILURE);
    }
}

/**
 * @brief Applies the specified averaging scheme to a value
 * @param val The value to which the averaging scheme is applied
 * @param method The averaging scheme to apply (arithmetic, geometric, harmonic and quadratic)
 */
float apply_averaging_scheme(float val, AveragingScheme method);

/**
 * @brief Finds the mean method for each product based on the input specification. If the input is empty, all products will use arithmetic mean.
 * @param averaging_scheme_specification The input specification for mean methods, in the format "prod1:method1,prod2:method2,...". Method can be arithmetic, geometric, harmonic or quadratic. If the input is empty, all products will use arithmetic mean.
 * @return A map of product names to their corresponding mean methods. If a product is not specified in the input, it will not be included in the map (the caller should apply the default method in that case).
 */
std::unordered_map<std::string, std::pair<AveragingScheme, std::string>> get_averaging_scheme_per_product(const std::string & averaging_scheme_specification);

namespace Hdf {

struct binIndexStruct {
    int32_t row_num; // 0-based
    double vsize;
    double hsize;
    int32_t basebin;
    int32_t beg;
    int32_t ext;
    int32_t numbin;
};

struct binIndexStruct_cdf4 {
    int32_t basebin;
    int32_t beg;
    int32_t ext;
    int32_t numbin;
};

struct binIndexStruct64_cdf4 {
    int64_t basebin;
    int64_t beg;
    int32_t ext;
    int64_t numbin;
};

struct binListStruct {
    int32_t bin_num;
    int16_t nobs;
    int16_t nscenes;
    int16_t time_rec;
    float weights;
    uint8_t sel_cat;
    int32_t flags_set;
    float lat;
    float lon;
};

struct binListStruct_hdf5 {
    int32_t bin_num;
    short nobs;
    short nscenes;
    float weights;
    int64_t flags_set;
};

struct binListStruct_cdf4 {
    uint32_t bin_num;
    short nobs;
    short nscenes;
    float weights;
    float time_rec;
};

struct binListStruct64_cdf4 {
    uint64_t bin_num;
    short nobs;
    short nscenes;
    float weights;
    float time_rec;
};

int create_vdata(int32_t file_id, int32_t vg_id, int32_t *vdata_id,
        const char *vdata_name, const char *class_name,
        int32_t n_flds, char const * const fldname[], int32_t type[],
        int32_t noext, int32_t *aid);

int32_t write_vdata(int vdata_id, int32_t n_recs_to_write, void *data);

int read_binList(int n_elem, int32_t vdata_id_binlist, binListStruct* binList, l3::L3Shape *shape);

int write_binList(int n_elem, int32_t vdata_id_binlist, binListStruct* binList);

int write_prodData(int n_elem, int32_t vdata_id_proddata, float *data,
        binListStruct* binList);
int copy_prodData(int n_elem, int32_t *binsToCopy,
        char const * const fldname3[],
        int32_t in_vdata_id_proddata,
        int32_t out_vdata_id_proddata);
int create_compound(hid_t group_id, const char *dataset_name, hid_t *dataset_id,
        hid_t *type_id, size_t typesize, int32_t n_flds,
        char const * const fldname[], size_t offset[], hid_t type[],
        hid_t *filespace, hid_t dataspace);
}

#endif
