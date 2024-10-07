/******************************************************************************
 * 
 * NAME: VcstGeoRctnglStruct.h
 *
 * DESCRIPTION: Contains interpolation rectangle structure and decimated
 *              geolocation data structure as well as constants used in these.
 *
 *
 ******************************************************************************/
#ifndef VCST_GEO_RCTNGL_STRUCT_H
#define VCST_GEO_RCTNGL_STRUCT_H

#include <VcstGeoDataStructs.h>
#include <VcstGeoParameters.h>
#include <vector>


const int REC_ROWS = 3;
const int REC_COLS = REC_ROWS;

// Geolocation data for the decimated interpolation rectangles.
// NOTE:  The GEO field height is not in the Decimated Geo structure because
// it is set from the Terrain Info in the storeGranule() method for IMG data
// and copied directly into the DMS buffer in calcModFromImg() for MOD data.

typedef struct {
    double lat[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    double lon[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    double satazm[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    double satzen[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    double sunazm[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    double sunzen[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    double moonazm[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    double moonzen[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    double range[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    //Add pixel phase and fraction.
    //0.0 to 100, percentage of Moon illuminated
    double moonFrac[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    //In radians, 0.0=full to +pi=new
    double moonPhase[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    unsigned char qFlags[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL][REC_ROWS][REC_COLS];
    unsigned char rectQFlags[VIIRS_SCANS][VIIRSI_MAX_NRCTNGL];
} ViirsGeoDecimType;


#endif

