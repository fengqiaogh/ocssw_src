/*
    Modification history:
    Programmer       Organization      Date      Description of change
    --------------   ------------    --------    ---------------------
    Joel Gales       Futuretech      08/10/03    Original Development

 */

/*
  Revision 2.0.2 09/13/22
  Add support for CLO
  B. Yang

  Revision 2.0.2 12/06/17
  Add support for prodlist for netcdf4 files
  J. Gales

  Revision 1.000 01/08/13
  Don't read control points if lonlatinterp == 0
  J. Gales

  Revision 0.993 05/21/10
  Fix east/westmost lon for extracts crossing dateline
  J. Gales

  Revision 0.992 09/26/08
  Only print prodlist (arg 9) if it exists
  J. Gales

  Revision 0.991 06/18/08
  Don't extract on "pxls" dimention for OCTS
  J. Gales

  Revision 0.990 12/05/06
  Fix slat,clat,...,slon,clon,...,& lon/lat boundary metadata
  J. Gales

  Revision 0.983 08/08/06
  Make cntl_pt_cols 1-based
  J. Gales

  Revision 0.982 06/15/06
  Set octs only if Level 1
  J. Gales

  Revision 0.981 06/01/06
  Make sure sscan is odd and escan is even if OCTS
  J. Gales

  Revision 0.98 05/30/06
  Fix problem with OCTS L1A extraction
  J. Gales
 */


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <genutils.h>

#include "l2extract.h"
#include <timeutils.h>

clo_optionList_t *optionList = NULL;

int main(int argc, char *argv[]) {
    int32_t sScan, eScan;
    int32_t sPixel, ePixel;
    static char productList[2048];
    static char waveList[2048];
    int32_t numPositionOptions;

    char *inFile, *outFile;
    char *productListStr = "";
    char *waveStr = "";

    int32_t pixelSub, scanSub;

    char *version = "3.2";

    printf("This is version %s of %s (compiled on %s %s)\n", version, "l2extract", __DATE__, __TIME__);

    /*** load input parameters into local variables */

    clo_setEnablePositionOptions(1);

    optionList = clo_createList();

    l2extractInitOptions(optionList, version);
    if (argc == 1) {
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }

    l2extractReadOptions(optionList, argc, argv);
    numPositionOptions = clo_getPositionNumOptions(optionList);

    if (numPositionOptions == 0) {
        inFile = clo_getString(optionList, "ifile");
        outFile = clo_getString(optionList, "ofile");
        sPixel = clo_getInt(optionList, "spix");
        ePixel = clo_getInt(optionList, "epix");
        sScan = clo_getInt(optionList, "sline");
        eScan = clo_getInt(optionList, "eline");
        pixelSub = 1;
        scanSub = 1;

        if (clo_isSet(optionList, "product")) {
            productListStr = clo_getRawString(optionList, "product");
        }
        if (clo_isSet(optionList, "wavelist")) {
            waveStr = clo_getRawString(optionList, "wavelist");
        }
    } else if (numPositionOptions == 8) {
        inFile = clo_getPositionString(optionList, 0);
        sPixel = atoi(clo_getPositionString(optionList, 1));
        ePixel = atoi(clo_getPositionString(optionList, 2));
        sScan = atoi(clo_getPositionString(optionList, 3));
        eScan = atoi(clo_getPositionString(optionList, 4));
        pixelSub = atoi(clo_getPositionString(optionList, 5));
        scanSub = atoi(clo_getPositionString(optionList, 6));
        outFile = clo_getPositionString(optionList, 7);
    } else if (numPositionOptions == 9) {
        inFile = clo_getPositionString(optionList, 0);
        sPixel = atoi(clo_getPositionString(optionList, 1));
        ePixel = atoi(clo_getPositionString(optionList, 2));
        sScan = atoi(clo_getPositionString(optionList, 3));
        eScan = atoi(clo_getPositionString(optionList, 4));
        pixelSub = atoi(clo_getPositionString(optionList, 5));
        scanSub = atoi(clo_getPositionString(optionList, 6));
        outFile = clo_getPositionString(optionList, 7);
        productListStr = clo_getPositionString(optionList, 8);
    } else {
        printf("\nNote: Number of Positional Arguments should be 0, 8 or 9!\n");
        exit(EXIT_FAILURE);
    }

    if (pixelSub != 1 || scanSub != 1) {
        printf("Subsampling not yet implemented.\n");
        exit(EXIT_FAILURE);
    }

    productList[0] = 0;
    if (productListStr[0] != 0) {
        strcpy(productList, ",l2_flags,");
        strcat(productList, productListStr);
        strcat(productList, ",");
        printf("product: %s\n", productListStr);
    }

    if (clo_getBool(optionList, "verbose")) {
        want_verbose = 1;
    }
    if (waveStr[0] != 0) {
        strcat(waveList, waveStr);
        if (want_verbose)
            printf("Wavelength requested by the user: %s\n", waveStr);
    }

    return extractNetCDF(inFile, outFile, sPixel, ePixel, sScan, eScan, productList, waveList);
}
