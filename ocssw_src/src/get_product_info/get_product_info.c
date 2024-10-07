#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <genutils.h>
#include <clo.h>
#include <sensorInfo.h>
#include <productInfo.h>

#define VERSION "3.0"
#define PROG_NAME "get_product_info"

void printProductList() {
    productInfo_t* productInfo = allocateProductInfo();

    //printf("chlor_a\n");
    getFirstProductInfo(productInfo);
    do {
        if (strcmp(productInfo->paramDesignator, "none") == 0) {
            printf("%s%s\n", productInfo->prefix, productInfo->suffix);
        } else {
            printf("%snnn%s\n", productInfo->prefix, productInfo->suffix);
        }
    } while (getNextProductInfo(productInfo));
    freeProductInfo(productInfo);
}

void printProductListSensor(int sensor) {

    // get sensor wavelengths
    int32_t* iwave;

    // stop verbose output
    int old_verbose = want_verbose;
    want_verbose = 0;
    int numWavelengths = rdsensorinfo(sensor, 0, "iwave", (void**) &iwave);
    want_verbose = old_verbose;

    if (numWavelengths == -1) {
        printf("-E- Could not lookup sensor %d wavelengths\n", sensor);
        exit(1);
    }

    //printf("chlor_a\n");
    productInfo_t* productInfo = allocateProductInfo();
    getFirstProductInfo(productInfo);
    do {
        if (strcmp(productInfo->paramDesignator, "none") == 0) {
            printf("%s%s\n", productInfo->prefix, productInfo->suffix);
        } else if (strcmp(productInfo->paramDesignator, "wave") == 0) {
            int i;
            for (i = 0; i < numWavelengths; i++) {
                if((iwave[i] > productInfo->paramWaveMin) && (iwave[i] < productInfo->paramWaveMax)) {
                    printf("%s%d%s\n", productInfo->prefix, iwave[i], productInfo->suffix);
                }
            }
        } else {
            printf("%snnn%s\n", productInfo->prefix, productInfo->suffix);
        }
    } while (getNextProductInfo(productInfo));
    freeProductInfo(productInfo);
}

void printProductCSV(productInfo_t* productInfo) {
    // name
    if (strcmp(productInfo->paramDesignator, "none") == 0) {
        printf("%s%s,", productInfo->prefix, productInfo->suffix);
    } else {
        printf("%snnn%s,", productInfo->prefix, productInfo->suffix);
    }

    printf("%s,", productInfo->paramDesignator);
    printf("%d,", productInfo->paramWaveMin);
    printf("%d,", productInfo->paramWaveMax);

    printf("%s,", productInfo->dataType);
    printf("%g,", productInfo->scaleFactor);
    printf("%g,", productInfo->addOffset);
    printf("%g,", productInfo->fillValue);
    printf("%g,", productInfo->validMin);
    printf("%g,", productInfo->validMax);
  
    printf("%s,", productInfo->displayScale);
    printf("%g,", productInfo->displayMin);
    printf("%g,", productInfo->displayMax);

    printf("\"%s\",", productInfo->units);
    printf("%s,", productInfo->standardName);
    printf("%s,", productInfo->palette);
    printf("%s,", productInfo->category);
    printf("\"%s\",", productInfo->reference);
    printf("\"%s\",", productInfo->description);
    printf("\"%s\",", productInfo->comment);
    
    printf("\n");
}

void printProductListCSV() {
    productInfo_t* productInfo = allocateProductInfo();

    // print headings
    printf("name,");
    printf("paramDesignator,");
    printf("paramWaveMin,");
    printf("paramWaveMax,");

    printf("dataType,");
    printf("scaleFactor,");
    printf("addOffset,");
    printf("fillValue,");
    printf("validMin,");
    printf("validMax,");

    printf("displayScale,");
    printf("displayMin,");
    printf("displayMax,");

    printf("units,");
    printf("standardName,");
    printf("palette,");
    printf("category,");
    printf("reference,");
    printf("description,");
    printf("comment");
    
    printf("\n");
    
    getFirstProductInfo(productInfo);
    do {
        printProductCSV(productInfo);
    } while (getNextProductInfo(productInfo));
    freeProductInfo(productInfo);
}

int main(int argc, char *argv[]) {
    // setup clo
    clo_setEnablePositionOptions(1);
    clo_setEnableParOption(0);
    clo_setVersion2(PROG_NAME, VERSION);
    char helpStr[2048];
    strcpy(helpStr, "Program to list products or detailed information about a single product\n");
    strcat(helpStr, "Usage: ");
    strcat(helpStr, PROG_NAME);
    strcat(helpStr, " [option=val] <productName>\nOptions:");
    clo_setHelpStr(helpStr);
    clo_optionList_t* list = clo_createList();

    clo_addOption(list, "l", CLO_TYPE_BOOL, "false", "list all of the products");
    clo_addOption(list, "r", CLO_TYPE_BOOL, "false", "list all of the products recursing through\n        the wavelengths of the sensor specified");
    clo_addOption(list, "csv", CLO_TYPE_BOOL, "false", "list all of the products, output in CSV format");
    clo_addOption(list, "sensor", CLO_TYPE_STRING, "MODISA", "sensor name (or ID) to use for wavelength\n        expansion or product lookup");
    clo_addOption(list, "<product>", CLO_TYPE_HELP, NULL, "<productName> product to print detailed information about");

    clo_readArgs(list, argc, argv);
    int numPositionOptions = clo_getPositionNumOptions(list);
    if (numPositionOptions > 1) {
        printf("-E- Too many parameters on the command line\n");
        exit(1);
    }

    // print usage if no args given
    if (argc == 1) {
        clo_printUsage(list);
        exit(0);
    }

    // list all products
    if (clo_getBool(list, "l")) {
        printProductList();
        exit(0);
    }

    // get the sensor ID
    const char* sensorName = clo_getString(list, "sensor");
    int sensorId = sensorName2SensorId(sensorName);
    if (sensorId == -1) {
        printf("-E- Could not find sensor \"%s\"\n", sensorName);
        exit(1);
    }
    // just in case a sensor ID was passed in
    sensorName = sensorId2SensorName(sensorId);

    // list all products recursively
    if (clo_getBool(list, "r")) {
        printProductListSensor(sensorId);
        exit(0);
    }

    // list all products recursively
    if (clo_getBool(list, "csv")) {
        printProductListCSV();
        exit(0);
    }

    // list detailed info for a single product
    if (numPositionOptions == 0) {
        printf("-E- A product name needs to be given on the command line\n");
        exit(1);
    }

    char* productName = clo_getPositionString(list, 0);

    productInfo_t* productInfo = allocateProductInfo();
    if (findProductInfo(productName, sensorId, productInfo)) {
        printf("sensorName=%s\n", sensorName);
        printProductInfo(productName, productInfo);
        return 0;
    }

    printf("-E- Product \"%s\" is not a valid product\n", productName);
    return 1;
}
