#ifndef PRODUCT_INFO_H
#define PRODUCT_INFO_H

#include <stdint.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif

/* default values for the product info structure */

#define PRODUCT_DEFAULT_description "no description"
#define PRODUCT_DEFAULT_units ""
#define PRODUCT_DEFAULT_palette "default"
#define PRODUCT_DEFAULT_paramDesignator "none"
#define PRODUCT_DEFAULT_paramWaveMin -1
#define PRODUCT_DEFAULT_paramWaveMax -1
#define PRODUCT_DEFAULT_standardName NULL
#define PRODUCT_DEFAULT_category "Miscellaneous"
#define PRODUCT_DEFAULT_dataType "float"
#define PRODUCT_DEFAULT_prefix NULL
#define PRODUCT_DEFAULT_suffix NULL
#define PRODUCT_DEFAULT_algorithmName ""
#define PRODUCT_DEFAULT_productName ""
#define PRODUCT_DEFAULT_cat_ix -1
#define PRODUCT_DEFAULT_prod_ix -1
#define PRODUCT_DEFAULT_rank 2
#define PRODUCT_DEFAULT_fillValue -32767
#define PRODUCT_DEFAULT_fillValue_byte -1
#define PRODUCT_DEFAULT_fillValue_ubyte 0xFF
#define PRODUCT_DEFAULT_fillValue_ushort 0xFFFF
#define PRODUCT_DEFAULT_fillValue_uint 0xFFFFFFFF
#define PRODUCT_DEFAULT_validMin (-FLT_MAX)
#define PRODUCT_DEFAULT_validMax FLT_MAX
#define PRODUCT_DEFAULT_displayScale "linear"
#define PRODUCT_DEFAULT_displayMin 0
#define PRODUCT_DEFAULT_displayMax 0
#define PRODUCT_DEFAULT_scaleFactor 1
#define PRODUCT_DEFAULT_addOffset 0
#define PRODUCT_DEFAULT_reference NULL
#define PRODUCT_DEFAULT_comment NULL
#define PRODUCT_DEFAULT_titleFormat "no format string"

/**
   Product data structure for product XML file.
 */
typedef struct productInfo_str {
    char *description; /**< description of the product  */
    char *units; /**< product units */
    char *palette; /**< palette file name (without .pal extention) */
    char *paramDesignator; /**< wavelengths (none,wave,band,int) */
    int paramWaveMin; /**< parameter minimum wavelength */
    int paramWaveMax; /**< parameter maximum wavelength */
    char *standardName; /**< CF compliant name for product */
    char *category; /**< category for product grouping */
    char *dataType; /**< data type (byte,ubyte,short,ushort,int,uint,float,double) */
    char *prefix; /**< first part of full product name (before wavelength) */
    char *suffix; /**< last part of product name (after wavelength) */
    char *algorithmName; /**< algorithm name, usually the default suffix */
    char *productName; /**< product name, usually the default prefix */
    int cat_ix; /**< catalog index for l2gen */
    int prod_ix; /**< param in prod name (none,wavelength,band#,int) */
    int rank; /**< how many dimensions in the product array */
    double fillValue; /**< missing data value */
    double validMin; /**< valid data range min */
    double validMax; /**< valid data range max */
    char *displayScale; /**< scaling for display (linear,log,arctan) */
    double displayMin; /**< suggested display min */
    double displayMax; /**< suggested display max */
    double scaleFactor; /**< physicalVal = fileVal * scaleFactor + addOffset */
    double addOffset; /**< physicalVal = fileVal * scaleFactor + addOffset */
    char *reference; /**< reference to a paper describing the algorithm */
    char *comment; /**< comments to provide additional useful information */
    char *titleFormat; /**< format string used to make the description */
} productInfo_t;


/*
 * function prototypes
 */

void clearProductInfo(productInfo_t* info);
void initProductInfo(productInfo_t* info);
productInfo_t* allocateProductInfo();
void freeProductInfo(productInfo_t* info);
void copyProductInfo(productInfo_t* dest, const productInfo_t* src);
void getFirstProductInfo(productInfo_t* info);
int getNextProductInfo(productInfo_t* info);
int findProductInfo(const char* name, int sensorId, productInfo_t* info);
char* getProductNameFull(productInfo_t* info);
void printProductInfo(const char* productFullName, const productInfo_t* info);
const char* getGCMDKeywords(const char* suite);

#ifdef __cplusplus
}

#include <string>

struct ProductAlias {
    std::string name;
    std::string prefix;
    std::string suffix;
};

bool findProductAlias(std::string productName, int sensorId, ProductAlias &productAlias);

#endif


#endif
