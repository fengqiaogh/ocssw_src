#include <productInfo.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <pugixml.hpp>
#include <genutils.h>
#include <sensorInfo.h>

#include <string>
#include <map>
#include <fstream>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
 
#define XML_STRING_SIZE 512

using namespace std;
using namespace pugi;
using namespace rapidjson;

/*
 * global variables
 */
static bool productInitialized = false;
static char productXMLFileName[FILENAME_MAX];
static xml_document rootNode;
static xml_node productsNode;
static xml_node productNode;
static xml_node algorithmNode;

static productInfo_t* productInfo = NULL;
static productInfo_t* algorithmInfo = NULL;

static string aliasFileName = "product_alias.json";

static map<string, ProductAlias> aliasMap;
static int aliasSensorId = -1;


/*
 * dup string if src is not NULL.
 */
static char* duplicateString(const char* src) {
    if (src == NULL)
        return NULL;
    else
        return strdup(src);
}

/**
 * set the product structure to defaults freeing memory from the current values
 *
 * @param info destination product structure
 */
extern "C" void clearProductInfo(productInfo_t* info) {
    if (info->description)
        free(info->description);
    info->description = duplicateString(PRODUCT_DEFAULT_description);
    if (info->units)
        free(info->units);
    info->units = duplicateString(PRODUCT_DEFAULT_units);
    if (info->palette)
        free(info->palette);
    info->palette = duplicateString(PRODUCT_DEFAULT_palette);
    if (info->paramDesignator)
        free(info->paramDesignator);
    info->paramDesignator = duplicateString(PRODUCT_DEFAULT_paramDesignator);
    info->paramWaveMin = PRODUCT_DEFAULT_paramWaveMin;
    info->paramWaveMax = PRODUCT_DEFAULT_paramWaveMax;
    if (info->standardName)
        free(info->standardName);
    info->standardName = duplicateString(PRODUCT_DEFAULT_standardName);
    if (info->category)
        free(info->category);
    info->category = duplicateString(PRODUCT_DEFAULT_category);
    if (info->dataType)
        free(info->dataType);
    info->dataType = duplicateString(PRODUCT_DEFAULT_dataType);
    if (info->prefix)
        free(info->prefix);
    info->prefix = duplicateString(PRODUCT_DEFAULT_prefix);
    if (info->suffix)
        free(info->suffix);
    info->suffix = duplicateString(PRODUCT_DEFAULT_suffix);
    if (info->algorithmName)
        free(info->algorithmName);
    info->algorithmName = duplicateString(PRODUCT_DEFAULT_algorithmName);
    if (info->productName)
        free(info->productName);
    info->productName = duplicateString(PRODUCT_DEFAULT_productName);
    info->cat_ix = PRODUCT_DEFAULT_cat_ix;
    info->prod_ix = PRODUCT_DEFAULT_prod_ix;
    info->rank = PRODUCT_DEFAULT_rank;
    info->fillValue = PRODUCT_DEFAULT_fillValue;
    info->validMin = PRODUCT_DEFAULT_validMin;
    info->validMax = PRODUCT_DEFAULT_validMax;
    if (info->displayScale)
        free(info->displayScale);
    info->displayScale = duplicateString(PRODUCT_DEFAULT_displayScale);
    info->displayMin = PRODUCT_DEFAULT_displayMin;
    info->displayMax = PRODUCT_DEFAULT_displayMax;
    info->scaleFactor = PRODUCT_DEFAULT_scaleFactor;
    info->addOffset = PRODUCT_DEFAULT_addOffset;
    if (info->reference)
        free(info->reference);
    info->reference = duplicateString(PRODUCT_DEFAULT_reference);
    if (info->comment)
        free(info->comment);
    info->comment = duplicateString(PRODUCT_DEFAULT_comment);
    if (info->titleFormat)
        free(info->titleFormat);
    info->titleFormat = duplicateString(PRODUCT_DEFAULT_titleFormat);
}

/**
 * set the product structure to defaults ignoring the current values
 *
 * @param info destination product structure
 */
extern "C" void initProductInfo(productInfo_t* info) {
    bzero(info, sizeof (productInfo_t));
    clearProductInfo(info);
}

/**
 * allocate memory for the product into structure and init to defaults.
 * structure should be freed using freeProductInfo()
 *
 * @return pointer to newly allocated structure
 */
extern "C" productInfo_t* allocateProductInfo() {
    productInfo_t* info = (productInfo_t*) allocateMemory(sizeof (productInfo_t),
            "productInfo");
    clearProductInfo(info);
    return info;
}

/**
 * free all the internal memory and the productInfo structure memory.
 *
 * @param info pointer to the product structure
 */
extern "C" void freeProductInfo(productInfo_t* info) {
    if (info->description)
        free(info->description);
    if (info->units)
        free(info->units);
    if (info->palette)
        free(info->palette);
    if (info->paramDesignator)
        free(info->paramDesignator);
    if (info->standardName)
        free(info->standardName);
    if (info->category)
        free(info->category);
    if (info->dataType)
        free(info->dataType);
    if (info->prefix)
        free(info->prefix);
    if (info->suffix)
        free(info->suffix);
    if (info->algorithmName)
        free(info->algorithmName);
    if (info->productName)
        free(info->productName);
    if (info->displayScale)
        free(info->displayScale);
    if (info->reference)
        free(info->reference);
    if (info->comment)
        free(info->comment);
    if (info->titleFormat)
        free(info->titleFormat);
    free(info);
}

/**
 * copy product info structure, just header data
 *
 * @param dest destination product structure
 * @param src source product structure
 */
void copyProductInfoHeader(productInfo_t* dest, const productInfo_t* src) {
    if (dest->productName)
        free(dest->productName);
    dest->productName = duplicateString(src->productName);
    if (dest->paramDesignator)
        free(dest->paramDesignator);
    dest->paramDesignator = duplicateString(src->paramDesignator);
    dest->paramWaveMin = src->paramWaveMin;
    dest->paramWaveMax = src->paramWaveMax;
}

/**
 * copy product info structure
 *
 * @param dest destination product structure
 * @param src source product structure
 */
extern "C" void copyProductInfo(productInfo_t* dest, const productInfo_t* src) {
    if (dest->description && strlen(dest->description) > 0)
        free(dest->description);
    dest->description = duplicateString(src->description);
    if (dest->units)
        free(dest->units);
    dest->units = duplicateString(src->units);
    if (dest->palette)
        free(dest->palette);
    dest->palette = duplicateString(src->palette);
    if (dest->paramDesignator)
        free(dest->paramDesignator);
    dest->paramDesignator = duplicateString(src->paramDesignator);
    dest->paramWaveMin = src->paramWaveMin;
    dest->paramWaveMax = src->paramWaveMax;
    if (dest->standardName)
        free(dest->standardName);
    dest->standardName = duplicateString(src->standardName);
    if (dest->category)
        free(dest->category);
    dest->category = duplicateString(src->category);
    if (dest->dataType)
        free(dest->dataType);
    dest->dataType = duplicateString(src->dataType);
    if (dest->prefix)
        free(dest->prefix);
    dest->prefix = duplicateString(src->prefix);
    if (dest->suffix)
        free(dest->suffix);
    dest->suffix = duplicateString(src->suffix);
    if (dest->algorithmName)
        free(dest->algorithmName);
    dest->algorithmName = duplicateString(src->algorithmName);
    if (dest->productName)
        free(dest->productName);
    dest->productName = duplicateString(src->productName);
    dest->cat_ix = src->cat_ix;
    dest->prod_ix = src->prod_ix;
    dest->rank = src->rank;
    dest->fillValue = src->fillValue;
    dest->validMin = src->validMin;
    dest->validMax = src->validMax;
    if (dest->displayScale)
        free(dest->displayScale);
    dest->displayScale = duplicateString(src->displayScale);
    dest->displayMin = src->displayMin;
    dest->displayMax = src->displayMax;
    dest->scaleFactor = src->scaleFactor;
    dest->addOffset = src->addOffset;
    if (dest->reference)
        free(dest->reference);
    dest->reference = duplicateString(src->reference);
    if (dest->comment)
        free(dest->comment);
    dest->comment = duplicateString(src->comment);
    if (dest->titleFormat)
        free(dest->titleFormat);
    dest->titleFormat = duplicateString(src->titleFormat);
}

/**
 * load XML file into local XML structure.
 */
void initXmlFile() {

    // bail if root node is already set
    if (productInitialized)
        return;
    productInitialized = true;

    char *dataRoot;
    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        exit(1);
    }
    strcpy(productXMLFileName, dataRoot);
    strcat(productXMLFileName, "/common/product.xml");

    xml_parse_result result = rootNode.load_file(productXMLFileName);
    if (!result) {
        printf("%s Line %d: could not open product XML file = %s\n",
                __FILE__, __LINE__, productXMLFileName);
        printf("    %s", result.description());
        exit(1);
    }
    productsNode = rootNode.child("products");
    if (!productsNode) {
        printf("%s Line %d: could not find products tag in XML file = %s\n",
                __FILE__, __LINE__, productXMLFileName);
        exit(1);
    }

    productInfo = allocateProductInfo();
    algorithmInfo = allocateProductInfo();

}

/**
 * check the product node to make sure the name == "product"
 */
void checkProductNode() {
    const char *tmpStr;

    tmpStr = productNode.name();
    if (strcmp(tmpStr, "product") != 0) {
        printf("%s Line %d: \"product\" node expected, found \"%s\" in file %s\n",
                __FILE__, __LINE__, tmpStr, productXMLFileName);
        exit(1);
    }
}

/**
 * check the algorithm node to make sure the name == "algorithm"
 */
void checkAlgorithmNode() {
    const char *tmpStr;

    tmpStr = algorithmNode.name();
    if (strcmp(tmpStr, "algorithm") != 0) {
        printf("%s Line %d: \"algorithm\" node expected, found \"%s\" in file %s\n",
                __FILE__, __LINE__, tmpStr, productXMLFileName);
        exit(1);
    }
}

/**
 * read the paramDesignator information out of the XML file
 */
void readParamDesignator(productInfo_t* info, xml_node node) {
    char *tmpStr;
    const char *tmpStr1;
    xml_node tmpNode;

    tmpNode = node.child("paramDesignator");

    // reset the param designator to none if there is one
    if (tmpNode) {
        free(info->paramDesignator);
        info->paramDesignator = strdup(PRODUCT_DEFAULT_paramDesignator);
        info->paramWaveMin = PRODUCT_DEFAULT_paramWaveMin;
        info->paramWaveMax = PRODUCT_DEFAULT_paramWaveMax;
    }

    while (tmpNode) {

        // grab the text

        tmpStr1 = tmpNode.child_value();
        tmpStr = trimBlanksDup(tmpStr1);

        if (strcmp(tmpStr, "none") == 0 ||
                strcmp(tmpStr, "band") == 0 ||
                strcmp(tmpStr, "int") == 0) {
            free(info->paramDesignator);
            info->paramDesignator = tmpStr;
            info->paramWaveMin = PRODUCT_DEFAULT_paramWaveMin;
            info->paramWaveMax = PRODUCT_DEFAULT_paramWaveMax;
            return;
        } else if (strcmp(tmpStr, "uv") == 0) {
            if (strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if (info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 100)
                info->paramWaveMin = 100;
            if (info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 400)
                info->paramWaveMax = 400;
        } else if (strcmp(tmpStr, "visible") == 0) {
            if (strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if (info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 400)
                info->paramWaveMin = 400;
            if (info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 725)
                info->paramWaveMax = 725;
        } else if (strcmp(tmpStr, "nir") == 0) {
            if (strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if (info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 725)
                info->paramWaveMin = 725;
            if (info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 1400)
                info->paramWaveMax = 1400;
        } else if (strcmp(tmpStr, "swir") == 0) {
            if (strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if (info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 1400)
                info->paramWaveMin = 1400;
            if (info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 3000)
                info->paramWaveMax = 3000;
        } else if (strcmp(tmpStr, "emissive") == 0) {
            if (strcmp(info->paramDesignator, "wave") != 0) {
                free(info->paramDesignator);
                info->paramDesignator = strdup("wave");
            }
            if (info->paramWaveMin == PRODUCT_DEFAULT_paramWaveMin ||
                    info->paramWaveMin > 3000)
                info->paramWaveMin = 3000;
            if (info->paramWaveMax == PRODUCT_DEFAULT_paramWaveMax ||
                    info->paramWaveMax < 15000)
                info->paramWaveMax = 15000;
        } else {
            printf("%s Line %d: \"paramDesignator\" node has illegal value \"%s\" in file %s\n",
                    __FILE__, __LINE__, tmpStr, productXMLFileName);
            exit(1);
        }

        free(tmpStr);

        // grab next node
        tmpNode = tmpNode.next_sibling("paramDesignator");
    }
}

/**
 * read the range node and write the data into the productInfo structure
 *
 * @param productInfo product structure to write the range data
 * @param rangeNode XML node to read the info from
 */
void readSingleRange(productInfo_t* productInfo, xml_node rangeNode) {
    xml_node node;
    if ((node = rangeNode.child("validMin")))
        productInfo->validMin = atof(node.child_value());
    if ((node = rangeNode.child("validMax")))
        productInfo->validMax = atof(node.child_value());
    if ((node = rangeNode.child("displayMin")))
        productInfo->displayMin = atof(node.child_value());
    if ((node = rangeNode.child("displayMax")))
        productInfo->displayMax = atof(node.child_value());
    if ((node = rangeNode.child("scaleFactor")))
        productInfo->scaleFactor = atof(node.child_value());
    if ((node = rangeNode.child("addOffset")))
        productInfo->addOffset = atof(node.child_value());
}

/**
 * read the range structures given the paramVal
 *
 * @param productInfo product structure to write the range data
 * @param productNode XML node to read the info from
 * @param paramVal integer to match to the range's min and max
 */
void readRange(productInfo_t* productInfo, xml_node productNode, int paramVal) {
    int min, max;
    xml_node rangeNode;
    xml_attribute attr;

    // set the default
    rangeNode = productNode.child("range");
    do {

        // make sure min and max are not defined
        if (rangeNode.attribute("min"))
            continue;
        if (rangeNode.attribute("max"))
            continue;

        readSingleRange(productInfo, rangeNode);

    } while ((rangeNode = rangeNode.next_sibling("range")));

    // now search for the matching range
    rangeNode = productNode.child("range");
    do {
        // make sure min and max are defined
        if ((attr = rangeNode.attribute("min"))) {
            min = attr.as_int();
            if ((attr = rangeNode.attribute("max"))) {
                max = attr.as_int();
                if (min <= paramVal && paramVal <= max) {
                    readSingleRange(productInfo, rangeNode);
                    break;
                }
            }
        }
    } while ((rangeNode = rangeNode.next_sibling("range")));

}

/**
 * read product header from the XML productNode into productInfo
 */
void readProductHeader() {

    // clear and set productName
    if (productInfo->productName)
        free(productInfo->productName);
    productInfo->productName = trimBlanksDup(productNode.attribute("name").value());

    // clear and set paramDesignator
    if (productInfo->paramDesignator)
        free(productInfo->paramDesignator);
    productInfo->paramDesignator = duplicateString(PRODUCT_DEFAULT_paramDesignator);
    readParamDesignator(productInfo, productNode);
}

/**
 * read product data from the XML productNode into productInfo
 */
void readProduct(int paramVal) {
    xml_node node;

    // clear product structure to defaults
    clearProductInfo(productInfo);

    if (productInfo->productName)
        free(productInfo->productName);
    productInfo->productName = trimBlanksDup(productNode.attribute("name").value());

    if ((node = productNode.child("standardName"))) {
        if (productInfo->standardName)
            free(productInfo->standardName);
        productInfo->standardName = trimBlanksDup(node.child_value());
    }

    if ((node = productNode.child("units"))) {
        if (productInfo->units)
            free(productInfo->units);
        productInfo->units = trimBlanksDup(node.child_value());
    }


    if ((node = productNode.child("palette"))) {
        if (productInfo->palette)
            free(productInfo->palette);
        productInfo->palette = trimBlanksDup(node.child_value());
    }

    node = productNode.child("category");
    if (!node) {
        printf("-E- Need to define \"category\" in product \"%s\" in product.xml\n",
                productInfo->productName);
        exit(EXIT_FAILURE);
    }
    if (productInfo->category)
        free(productInfo->category);
    productInfo->category = trimBlanksDup(node.child_value());

    if ((node = productNode.child("displayScale"))) {
        if (productInfo->displayScale)
            free(productInfo->displayScale);
        productInfo->displayScale = trimBlanksDup(node.child_value());
    }

    if ((node = productNode.child("type"))) {
        if (productInfo->dataType)
            free(productInfo->dataType);
        productInfo->dataType = trimBlanksDup(node.child_value());
        if(!strcmp(productInfo->dataType, "byte")) {
            if(productInfo->fillValue == PRODUCT_DEFAULT_fillValue) {
                productInfo->fillValue = PRODUCT_DEFAULT_fillValue_byte;
            }
        } else if(!strcmp(productInfo->dataType, "ubyte")) {
            if(productInfo->fillValue == PRODUCT_DEFAULT_fillValue) {
                productInfo->fillValue = PRODUCT_DEFAULT_fillValue_ubyte;
            }
        } else if(!strcmp(productInfo->dataType, "ushort")) {
            if(productInfo->fillValue == PRODUCT_DEFAULT_fillValue) {
                productInfo->fillValue = PRODUCT_DEFAULT_fillValue_ushort;
            }
        } else if(!strcmp(productInfo->dataType, "uint")) {
            if(productInfo->fillValue == PRODUCT_DEFAULT_fillValue) {
                productInfo->fillValue = PRODUCT_DEFAULT_fillValue_uint;
            }
        }
    }

    if ((node = productNode.child("reference"))) {
        if (productInfo->reference)
            free(productInfo->reference);
        productInfo->reference = trimBlanksDup(node.child_value());
    }

    if ((node = productNode.child("comment"))) {
        if (productInfo->comment)
            free(productInfo->comment);
        productInfo->comment = trimBlanksDup(node.child_value());
    }

    readParamDesignator(productInfo, productNode);
    if (strcmp(productInfo->paramDesignator, "none") != 0) {
        productInfo->prod_ix = paramVal;
    }

    readRange(productInfo, productNode, paramVal);
}

/**
 * read algorithm header from the XML algorithmNode and productNode
 * into algorithmInfo
 */
void readAlgorithmHeader() {

    xml_attribute attr;
    xml_node node;

    // copy header info from the productInfo
    copyProductInfoHeader(algorithmInfo, productInfo);

    // clear algorithm header
    if (algorithmInfo->algorithmName) {
        free(algorithmInfo->algorithmName);
        algorithmInfo->algorithmName = NULL;
    }
    if (algorithmInfo->prefix) {
        free(algorithmInfo->prefix);
        algorithmInfo->prefix = NULL;
    }
    if (algorithmInfo->suffix) {
        free(algorithmInfo->suffix);
        algorithmInfo->suffix = NULL;
    }

    // populate the fields
    if ((attr = algorithmNode.attribute("name"))) {
        algorithmInfo->algorithmName = trimBlanksDup(attr.value());
    } else {
        algorithmInfo->algorithmName = strdup("");
    }

    readParamDesignator(algorithmInfo, algorithmNode);

    char defaultPrefix[XML_STRING_SIZE];
    if (strcmp(algorithmInfo->paramDesignator, "none") == 0) {
        strcpy(defaultPrefix, algorithmInfo->productName);
    } else {
        strcpy(defaultPrefix, algorithmInfo->productName);
        strcat(defaultPrefix, "_");
    }

    if ((node = algorithmNode.child("prefix"))) {
        algorithmInfo->prefix = trimBlanksDup(node.child_value());
    } else {
        algorithmInfo->prefix = strdup(defaultPrefix);
    }
    if ((node = algorithmNode.child("suffix"))) {
        algorithmInfo->suffix = trimBlanksDup(node.child_value());
    } else {
        if (strlen(algorithmInfo->algorithmName) > 0) {
            algorithmInfo->suffix = (char*) malloc(strlen(algorithmInfo->algorithmName) + 2);
            strcpy(algorithmInfo->suffix, "_");
            strcat(algorithmInfo->suffix, algorithmInfo->algorithmName);
        } else {
            algorithmInfo->suffix = strdup("");
        }
    }
}

/**
 * read algorithm data from the XML algorithmNode into algorithmInfo
 */
void readAlgorithm(int paramVal) {
    char tmpStr[XML_STRING_SIZE];
    xml_node node;

    // start with the productInfo as default values
    copyProductInfo(algorithmInfo, productInfo);
    readAlgorithmHeader();

    algorithmInfo->cat_ix = atoi(algorithmNode.child_value("cat_ix"));
    if ((node = algorithmNode.child("rank")))
        algorithmInfo->rank = atoi(node.child_value());


    if ((node = algorithmNode.child("units"))) {
        if (algorithmInfo->units)
            free(algorithmInfo->units);
        algorithmInfo->units = trimBlanksDup(node.child_value());
    }

    if ((node = algorithmNode.child("fillValue")))
        algorithmInfo->fillValue = atof(node.child_value());

    node = algorithmNode.child("description");
    if (!node) {
        printf("-E- Need to define \"description\" in product \"%s_%s\" in product.xml\n",
                productInfo->productName, productInfo->algorithmName);
        exit(EXIT_FAILURE);

    }
    strcpy(tmpStr, node.child_value());
    trimBlanks(tmpStr);

    if (algorithmInfo->titleFormat)
        free(algorithmInfo->titleFormat);
    algorithmInfo->titleFormat = strdup(tmpStr);
    if (algorithmInfo->description)
        free(algorithmInfo->description);
    if (strcmp(algorithmInfo->paramDesignator, "none") == 0) {
        if ((node = algorithmNode.child("prod_ix")))
            algorithmInfo->prod_ix = atoi(node.child_value());
        algorithmInfo->description = trimBlanksDup(tmpStr);
    } else {
        algorithmInfo->description = (char*) malloc(strlen(tmpStr) + 64);
        algorithmInfo->prod_ix = paramVal;
        sprintf(algorithmInfo->description, tmpStr, algorithmInfo->prod_ix);
    }

    if ((node = algorithmNode.child("reference"))) {
        if (algorithmInfo->reference)
            free(algorithmInfo->reference);
        algorithmInfo->reference = trimBlanksDup(node.child_value());
    }

    if ((node = algorithmNode.child("comment"))) {
        if (algorithmInfo->comment)
            free(algorithmInfo->comment);
        algorithmInfo->comment = trimBlanksDup(node.child_value());
    }

    if ((node = algorithmNode.child("type"))) {
        if (algorithmInfo->dataType)
            free(algorithmInfo->dataType);
        algorithmInfo->dataType = trimBlanksDup(node.child_value());
        if(!strcmp(algorithmInfo->dataType, "byte")) {
            if(algorithmInfo->fillValue == PRODUCT_DEFAULT_fillValue) {
                algorithmInfo->fillValue = PRODUCT_DEFAULT_fillValue_byte;
            }
        } else if(!strcmp(algorithmInfo->dataType, "ubyte")) {
            if(algorithmInfo->fillValue == PRODUCT_DEFAULT_fillValue) {
                algorithmInfo->fillValue = PRODUCT_DEFAULT_fillValue_ubyte;
            }
        } else if(!strcmp(algorithmInfo->dataType, "ushort")) {
            if(algorithmInfo->fillValue == PRODUCT_DEFAULT_fillValue) {
                algorithmInfo->fillValue = PRODUCT_DEFAULT_fillValue_ushort;
            }
        } else if(!strcmp(algorithmInfo->dataType, "uint")) {
            if(algorithmInfo->fillValue == PRODUCT_DEFAULT_fillValue) {
                algorithmInfo->fillValue = PRODUCT_DEFAULT_fillValue_uint;
            }
        }
    }

    readRange(algorithmInfo, algorithmNode, paramVal);
}

/**
 * find first algorithm node header
 */
void findFirstAlgorithmHeader() {
    productNode = productsNode.child("product");
    if (!productNode) {
        printf("%s Line %d: could not find first product tag in XML file = %s\n",
                __FILE__, __LINE__, productXMLFileName);
        exit(1);
    }
    checkProductNode();
    readProductHeader();

    algorithmNode = productNode.child("algorithm");
    if (!algorithmNode) {
        printf("%s Line %d: could not find first algorithm tag for product=%s in XML file = %s\n",
                __FILE__, __LINE__, productInfo->productName, productXMLFileName);
        exit(1);
    }
    checkAlgorithmNode();
    readAlgorithmHeader();
}

/**
 * find next algorithm node header
 *
 * @return 1 if algorithm found, 0 if reached end of file
 */
int findNextAlgorithmHeader() {
    algorithmNode = algorithmNode.next_sibling("algorithm");
    if (!algorithmNode) {

        // try next product
        productNode = productNode.next_sibling("product");
        if (!productNode)
            return 0;
        checkProductNode();
        readProductHeader();

        algorithmNode = productNode.child("algorithm");
        if (!algorithmNode) {
            printf("%s Line %d: could not find first algorithm tag for product=%s in XML file = %s\n",
                    __FILE__, __LINE__, productInfo->productName, productXMLFileName);
            exit(1);
        }
    }
    checkAlgorithmNode();
    readAlgorithmHeader();
    return 1;
}

/**
 * compare algorithmInfo header to product string
 *
 * @param productFullName full name of the product to compare
 * @param sensorId sensor ID to use for wavelength comparison
 * @param paramVal pointer to an int where parameter value will be written
 * @return 1 if matches, 0 if not
 */
int compareAlgorithmHeader(const char* productFullName, int sensorId, int* paramVal) {

    // see if prefix matches
    if (strncmp(productFullName, algorithmInfo->prefix, strlen(algorithmInfo->prefix)))
        return 0;

    char tmpStr[XML_STRING_SIZE];
    if (strcmp(algorithmInfo->paramDesignator, "none") == 0) {
        strcpy(tmpStr, algorithmInfo->prefix);
        strcat(tmpStr, algorithmInfo->suffix);
        if (strcmp(tmpStr, productFullName) == 0) {
            *paramVal = PRODUCT_DEFAULT_prod_ix;
            return 1;
        } else {
            return 0;
        }
    } else {
        int nameLen = strlen(productFullName);
        int prefixLen = strlen(algorithmInfo->prefix);
        int suffixLen = strlen(algorithmInfo->suffix);
        int paramLen = nameLen - prefixLen - suffixLen;

        // need at least one char to be a valid param
        if (paramLen < 1)
            return 0;

        // check the suffix
        if (strcmp(algorithmInfo->suffix, productFullName + nameLen - suffixLen) != 0)
            return 0;

        // get the param value
        char paramStr[XML_STRING_SIZE];
        if (paramLen > XML_STRING_SIZE - 1) // just in case
            paramLen = XML_STRING_SIZE - 1;
        paramStr[0] = 0;
        strncat(paramStr, productFullName + prefixLen, paramLen);

        if (!isValidInt(paramStr))
            return 0;

        int i = atoi(paramStr);
        if (strcmp(algorithmInfo->paramDesignator, "wave") == 0) {
            if (algorithmInfo->paramWaveMin <= i && i <= algorithmInfo->paramWaveMax) {
                int32_t numBands;
                int32_t *iwave;
                int waveIndex;
                int old_verbose = want_verbose;
                want_verbose = 0;
                numBands = rdsensorinfo(sensorId, 0, "iwave", (void**) &iwave);
                numBands += rdsensorinfo(sensorId, 0, "NbandsIR", NULL);
                want_verbose = old_verbose;
                if (numBands != -1) {
                    for (waveIndex = 0; waveIndex < numBands; waveIndex++) {
                        if (i == iwave[waveIndex]) {
                            *paramVal = i;
                            return 1;
                        }
                    }
                }

            }
        } else {
            *paramVal = i;
            return 1;
        }
    }
    return 0;
}

/**
 * find first product and fill in the product structure.
 *
 * @param info pre allocated product structure to fill in
 */
extern "C" void getFirstProductInfo(productInfo_t* info) {
    initXmlFile();
    findFirstAlgorithmHeader();
    readProduct(-1);
    readAlgorithm(-1);
    copyProductInfo(info, algorithmInfo);
}

/**
 * find next product and fill in the product structure.
 *
 * @param info pre allocated product structure to fill in
 * @return 1 if product found, 0 if no more products
 */
extern "C" int getNextProductInfo(productInfo_t* info) {
    if (findNextAlgorithmHeader()) {
        readProduct(-1);
        readAlgorithm(-1);
        copyProductInfo(info, algorithmInfo);
        return 1;
    }
    return 0;
}

/* json alias file
{
    "alias" : {
	    "name" : "realProduct",
	    "prefix" : "prefixReplacement",
	    "suffix" : "suffixReplacement"
    }
}
*/

void readProductAliasFile(string fileName, map<string, ProductAlias> &aliasMap) {
    char fileNameFixed[FILENAME_MAX];
    parse_file_name(fileName.c_str(), fileNameFixed);

    // does the file exist
    if (access(fileNameFixed, R_OK) != -1) {
        ifstream ifs(fileNameFixed);
        IStreamWrapper isw(ifs);
        Document document;
        document.ParseStream(isw);

        for (Value::ConstMemberIterator itr = document.MemberBegin(); itr != document.MemberEnd(); ++itr) {
            ProductAlias alias;
            alias.name = itr->value["name"].GetString();
            alias.prefix = itr->value["prefix"].GetString();
            alias.suffix = itr->value["suffix"].GetString();

            aliasMap[itr->name.GetString()] = alias;
        }
    }
}

void initProductAliasMap(int sensorId) {
    if(sensorId != aliasSensorId) {
        aliasSensorId = sensorId;
        aliasMap.clear();
        
        // read aliases file from common
        string fileName = (string)"$OCDATAROOT/common/" + aliasFileName;
        readProductAliasFile(fileName, aliasMap);

        // read aliases file from sensor
        fileName = (string)"$OCDATAROOT/" + sensorId2SensorDir(sensorId) + "/" + aliasFileName;
        readProductAliasFile(fileName, aliasMap);

        // read aliases file from sub-sensor
        int subsensorId = sensorId2SubsensorId(sensorId);
        if(subsensorId != -1) {
            fileName = (string)"$OCDATAROOT/" + sensorId2SensorDir(sensorId) + "/" + subsensorId2SubsensorDir(subsensorId) + "/" + aliasFileName;
            readProductAliasFile(fileName, aliasMap);
        }
    }
}

bool findProductAlias(string productName, int sensorId, ProductAlias &productAlias) {
    initProductAliasMap(sensorId);

    // search for alias
    try {
        productAlias = aliasMap.at(productName);
        return true;
    } catch (const std::out_of_range& oor) {
        //std::cerr << "Out of Range error: " << oor.what() << '\n';
    }
    return false;
}

/**
 * find product and fill in the product structure.
 *
 * @param productName product to find
 * @param sensorId sensor ID to use when looking up the product
 * @param info pre allocated product structure to fill (allocateProductInfo)
 * @return 1 if product found, 0 if not found
 */
extern "C" int findProductInfo(const char* productName, int sensorId, productInfo_t* info) {
    int paramVal;
    const char* tmpProductName = productName;

    // search for alias
    ProductAlias productAlias;
    bool isAlias = findProductAlias(productName, sensorId, productAlias);
    if(isAlias) {
        tmpProductName = productAlias.name.c_str();
    }

    initXmlFile();
    findFirstAlgorithmHeader();
    do {
        if (compareAlgorithmHeader(tmpProductName, sensorId, &paramVal)) {
            readProduct(paramVal);
            readAlgorithm(paramVal);
            copyProductInfo(info, algorithmInfo);

            // fix up the record for aliases
            if(isAlias) {
                free(info->prefix);
                info->prefix = strdup(productAlias.prefix.c_str());
                free(info->suffix);
                info->suffix = strdup(productAlias.suffix.c_str());
            }

            return 1;
        }
    } while (findNextAlgorithmHeader());

    return 0;
}

/**
 *
 * @param info product structure to read
 * @return full product name.  Pointer to internal memory.
 */
extern "C" char* getProductNameFull(productInfo_t* info) {
    static char name[XML_STRING_SIZE];

    if (strcmp(info->paramDesignator, "none") == 0) {
        strcpy(name, info->prefix);
        strcat(name, info->suffix);
    } else {
        sprintf(name, "%s%d%s", info->prefix, info->prod_ix, info->suffix);
    }
    return name;
}

/**
 * print product structure.
 *
 * @param productFullName product name
 * @param info structure to print
 */
extern "C" void printProductInfo(const char* productFullName, const productInfo_t* info) {
    printf("fullProductName=%s\n", productFullName);
    printf("description=%s\n", info->description);
    printf("units=%s\n", info->units);
    printf("palette=%s\n", info->palette);
    printf("paramDesignator=%s\n", info->paramDesignator);
    printf("paramWaveMin=%d\n", info->paramWaveMin);
    printf("paramWaveMax=%d\n", info->paramWaveMax);
    printf("standardName=%s\n", info->standardName);
    printf("category=%s\n", info->category);
    printf("dataType=%s\n", info->dataType);
    printf("prefix=%s\n", info->prefix);
    printf("suffix=%s\n", info->suffix);
    printf("algorithmName=%s\n", info->algorithmName);
    printf("productName=%s\n", info->productName);
    printf("cat_ix=%d\n", info->cat_ix);
    printf("prod_ix=%d\n", info->prod_ix);
    printf("rank=%d\n", info->rank);
    printf("fillValue=%g\n", info->fillValue);
    printf("validMin=%g\n", info->validMin);
    printf("validMax=%g\n", info->validMax);
    printf("displayScale=%s\n", info->displayScale);
    printf("displayMin=%g\n", info->displayMin);
    printf("displayMax=%g\n", info->displayMax);
    printf("scaleFactor=%g\n", info->scaleFactor);
    printf("addOffset=%g\n", info->addOffset);
    printf("reference=%s\n", info->reference);
    printf("comment=%s\n", info->comment);
    printf("titleFormat=%s\n", info->titleFormat);
}
