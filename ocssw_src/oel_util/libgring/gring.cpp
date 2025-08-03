// Includes from main_l1info.c
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdbool.h>
#include <fcntl.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include <unistd.h>
#include <sstream>
#include <algorithm>
#include <utility>

// Additional includes
#include "gring.h"
#include "allocate2d.h"
#include "genutils.h"
#include "piutils.h"

void Gring::initialize() {
    centerLats = std::vector<float>(MAX_GEOBOX_CNT);

    // State vars
    foundGoodScan = false;
    prevTrackDirection = 0;
    centerLatitude = -999;
    geoboxCount = 0;
    directionCheckPrevious = -1;
    latitudeMin = 90;   // All valid values will be lower than this
    latitudeMax = -90;  // All valid values will be higher than this
    longitudeMin = 359.999;
    longitudeMax = 0;
    previousLon = 0;
    invalidPixelMask = 0;
}

/**
 * @brief Initializes a Gring with the latitude delta set to 5, the scan direction set to right-to-left,
 * and the direction check interval set to 1
 */
Gring::Gring() {
    initialize();
    deltaLatDegrees = 20;
    scanDirection = -1;
    directionCheckInterval = 1;
}

Gring::Gring(float latPointDelta, int scanDirection, size_t directionCheckInterval) {
    initialize();
    this->deltaLatDegrees = latPointDelta;
    this->scanDirection = scanDirection;
    this->directionCheckInterval = directionCheckInterval;
}

/**
 * @brief Instantiate a Gring given a sensor ID
 * @param sensorId The identification number of the sensor
 * @param deltaLatDegrees The chosen interval at which this Gring will add scans to itself
 */
Gring::Gring(int sensorId, float deltaLatDegrees) {
    initialize();
    bool sensorScansLeftToRight =
        sensorId == CZCS || sensorId == OCI || sensorId == HARP2 || sensorId == SPEXONE;

    if (sensorScansLeftToRight)
        scanDirection = 1;
    else
        scanDirection = -1;

    std::string instrumentName = sensorId2InstrumentName(sensorId);
    if (instrumentName == "MODIS")
        directionCheckInterval = 10;
    else if (instrumentName == "VIIRS")
        directionCheckInterval = 16;
    else
        directionCheckInterval = 5;

    this->deltaLatDegrees = deltaLatDegrees;
}

Gring::~Gring() {
}

void Gring::setDeltaLat(float newDeltaLat) {
    deltaLatDegrees = newDeltaLat;
}

void Gring::setDirectionCheckInterval(size_t newInterval) {
    directionCheckInterval = newInterval;
}

void Gring::setScanDirection(int newDirection) {
    scanDirection = newDirection;
}

void Gring::setInvalidPixelMask(int32_t newMask) {
    invalidPixelMask = newMask;
}

float Gring::getDeltaLat() {
    return deltaLatDegrees;
}

size_t Gring::getDirectionCheckInterval() {
    return directionCheckInterval;
}

int Gring::getScanDirection() {
    return scanDirection;
}

/**
 * @brief Return all of the geospatial extremes variables
 * @return An array of latitude min, latitude max, longitude min, longitude max
 */
std::array<float, 4> Gring::getGeospatialExtremes() {
    return {latitudeMin, latitudeMax, longitudeMin, longitudeMax};
}

float Gring::getGeospatialLatitudeMin() {
    return latitudeMin;
}
float Gring::getGeospatialLatitudeMax() {
    return latitudeMax;
}
float Gring::getGeospatialLongitudeMin() {
    return longitudeMin;
}
float Gring::getGeospatialLongitudeMax() {
    return longitudeMax;
}

bool Gring::isValidLat(float lat) {
    return (lat != BAD_FLT && -90.0 <= lat && lat <= 90.0);
}

bool Gring::isValidLon(float lon) {
    return (lon != BAD_FLT && -180.0 <= lon && lon <= 180.0);
}

// Expects a lon in 180 degree space and returns that same value mapped to 0-360
float normalizeLonTo360Space(float lon) {
    return fmod(lon + 360, 360);
}

/**
 * @brief Adds a scan's lon/lat values as long as the maximum number of gring points is not exceeded
 * @param l1rec A struct containing a pointer to lat, lon, flags, and number of pixels values
 * @param recnum The current scan index
 * @param escan The index of the last scan
 */
void Gring::processScan(float *lat, float *lon, int32_t npix, size_t recnum, size_t escan, int32_t *flags) {
    bool scanShouldBeIncluded;
    static ScanBounds scanBounds;

    scanShouldBeIncluded = shouldIncludeScan(lat, lon, npix, recnum, escan, flags, scanBounds);
    if (scanShouldBeIncluded) {
        if (geoboxCount >= (MAX_GEOBOX_CNT - 1)) {
            printf("Error - Max number of gring points exceeded.\n");
            exit(EXIT_FAILURE);
        }

        //  This scan should be included in the gring, so caller stuffs scanBoundsOut into a geoboxes array.
        //  Note that a geoboxes (borrowed from l2_generic.c) is just the endpoints (slon,slat,elon,elat) of a
        //  scan line.
        geoboxes[0][geoboxCount] = scanBounds.start.lon;  // Starting lon of this scan
        geoboxes[1][geoboxCount] = scanBounds.start.lat;  // Starting lat of this scan
        geoboxes[2][geoboxCount] = scanBounds.end.lon;    // Ending lon of this scan
        geoboxes[3][geoboxCount] = scanBounds.end.lat;    // Ending lat of this scan
        centerLats[geoboxCount] = scanBounds.center.lat;  // save the center lat also
        geoboxCount++;

        updateExtremes(scanBounds);
    }
}

void Gring::updateExtremes(const ScanBounds &scanBounds) {
    const LatLonPoint &start = scanBounds.start;
    const LatLonPoint &end = scanBounds.end;

    if (start.lat > latitudeMax)
        latitudeMax = start.lat;
    if (start.lat < latitudeMin)
        latitudeMin = start.lat;
    if (end.lat > latitudeMax)
        latitudeMax = end.lat;
    if (end.lat < latitudeMin)
        latitudeMin = end.lat;

    float lonMin360 = normalizeLonTo360Space(longitudeMin);
    float lonMax360 = normalizeLonTo360Space(longitudeMax);
    float startLon360 = normalizeLonTo360Space(start.lon);
    float endLon360 = normalizeLonTo360Space(end.lon);

    // Update longitude min/max as well as local copies
    if (startLon360 < lonMin360) {
        longitudeMin = start.lon;
        lonMin360 = startLon360;
    }
    if (startLon360 > lonMax360) {
        longitudeMax = start.lon;
        lonMax360 = startLon360;
    }
    if (endLon360 < lonMin360) {
        longitudeMin = end.lon;
        lonMin360 = endLon360;
    }
    if (endLon360 > lonMax360) {
        longitudeMax = end.lon;
        lonMax360 = endLon360;
    }
    
    float lonSpan360 = fmod((lonMax360 - lonMin360) + 360.0, 360.0);
    if (lonSpan360 >= 180.0) { // Antimeridian check
        float temp = longitudeMax;
        longitudeMax = longitudeMin;
        longitudeMin = temp;
    }

    previousLon = end.lon;  // TODO: This assumes a scan direction
}

/**
 * This function is called per scan. It returns true if this (recnum) scan should be included in a gring.
 * It borrows code from main_l1info.c (e.g. for finding good spix,endPixel,cpix) and from l2_generic.c (e.g.
 * for building the gring from a geoboxes).
 *
 * NB: Allocation (and freeing) of memory for the output struct scanBounds is outside the scope of this
 * function. If the return value is false, the output scanBounds should be ignored.
 *
 * @param lats Latitude values
 * @param lons Longitude values
 * @param numPixels
 * @param currScan Index of this scan
 * @param lastScan Index of last scan
 * @param scanBounds Output, a structure containing the good endpoints--if any--and center points of the scan
 * line
 * @param flags Indicators of whether a given pixel is invalid. Should be checked against a predefined mask,
 * and expected to be numPixels long
 * @return Whether the given scan should be included in this Gring
 */
bool Gring::shouldIncludeScan(float *lats, float *lons, size_t numPixels, size_t currScan, size_t lastScan,
                              int32_t *flags, ScanBounds &scanBounds) {
    int goodStartPixel, goodEndPixel, goodCenterPixel;
    int thisScanDirection = 0;

    // Initialize
    size_t endPixel = numPixels - 1;
    size_t centerPixel = numPixels / 2;
    goodStartPixel = -1;
    goodEndPixel = -1;
    goodCenterPixel = -1;

    // Find the good start, center, and end pixels for this scan
    for (size_t ip = 0; ip <= endPixel; ip++) {
        bool pixelIsValid = false;

        if (flags != nullptr) {
            pixelIsValid =
                (flags[ip] & invalidPixelMask) == 0 && isValidLat(lats[ip]) && isValidLon(lons[ip]);
        } else {  // Flags came in unset
            pixelIsValid = isValidLat(lats[ip]) && isValidLon(lons[ip]);
        }

        if (pixelIsValid) {
            if (goodStartPixel == -1) {
                goodStartPixel = ip;
            }
            goodEndPixel = ip;  // Eventually this will truly be the last pixel
            // Store this good pixel in centerPixel. Eventually (somewhere in the middle), it will be last
            // good one at or before centerPixel.
            if ((ip <= centerPixel) || (goodCenterPixel == -1)) {
                goodCenterPixel = ip;
            }
        }
    }

    // If this is a good scan, set the scan bounds (i.e. endpoints), and test to see if it should be included
    // in the gring
    bool scanIsValid = (goodStartPixel != -1) && (goodEndPixel != -1) && (goodCenterPixel != -1);
    if (scanIsValid) {
        scanBounds.start.lat = lats[goodStartPixel];
        scanBounds.start.lon = lons[goodStartPixel];
        scanBounds.center.lat = lats[goodCenterPixel];
        scanBounds.center.lon = lons[goodCenterPixel];
        scanBounds.end.lat = lats[goodEndPixel];
        scanBounds.end.lon = lons[goodEndPixel];

        // CONDITION 1: If this is the FIRST good scan, return scanBounds, after resetting "state
        // information".
        if (!foundGoodScan) {
            foundGoodScan = true;
            centerLatitude = lats[goodCenterPixel];
            directionCheckPrevious = currScan;
            return true;
        } else {
            // This is a good scan, but not the first good one. So, we can start tracking direction; i.e. no
            // sense in tracking direction prior to the first good scan. wait for a few lines before checking
            // direction
            bool needToCheckDirection = (currScan - directionCheckPrevious) > directionCheckInterval;
            if (needToCheckDirection) {
                directionCheckPrevious = currScan;
                if ((lats[goodCenterPixel] - centerLatitude) > 0) {
                    thisScanDirection = 1;
                } else {
                    thisScanDirection = -1;
                }
                if (prevTrackDirection == 0)
                    prevTrackDirection = thisScanDirection;
            }
        }
        // CONDITION 2: If the latitude of this scan is more than deltaLatDegrees from the centerLatitude,
        // return scanBounds, after resetting "state information". Also CONDITION 1 must have already been
        // satisfied; i.e. centerLatitude has a non-default value.
        if (foundGoodScan && (fabs(centerLatitude - lats[goodCenterPixel]) > deltaLatDegrees)) {
            prevTrackDirection = thisScanDirection;
            centerLatitude = lats[goodCenterPixel];
            return true;
        }
        // CONDITION 3: If the satellite passes over a pole, detect the change from e.g. ascending to
        // descending (or vice versa) and return scanBounds after resetting "state information". Also
        // CONDITION 1 must have already been satisfied; i.e. centerLatitude has a non-default value.  Only
        // look in high latitudes.
        if (foundGoodScan && (prevTrackDirection == -thisScanDirection) && (thisScanDirection != 0) &&
            (fabs(lats[goodCenterPixel]) > 75)) {
            prevTrackDirection = thisScanDirection;
            centerLatitude = lats[goodCenterPixel];
            return true;
        }
    }

    // If we get this far, this is not a scan to include in the gring. So, just return 0. However, if this is
    // the last scan, we want to use the last GOOD scan.
    if (currScan == lastScan) {
        // If we got this far (without returning) AND this is the last scan, use the last GOOD scan (still
        // stored in scanBounds).
        return true;
    } else {
        return false;
    }
}

void Gring::deleteExtraLine() {
    if (geoboxCount > 2) {
        if (fabs(centerLats[geoboxCount - 1] - centerLats[geoboxCount - 2]) < deltaLatDegrees / 2) {
            geoboxCount--;
            geoboxes[0][geoboxCount - 1] = geoboxes[0][geoboxCount];
            geoboxes[1][geoboxCount - 1] = geoboxes[1][geoboxCount];
            geoboxes[2][geoboxCount - 1] = geoboxes[2][geoboxCount];
            geoboxes[3][geoboxCount - 1] = geoboxes[3][geoboxCount];
            centerLats[geoboxCount - 1] = centerLats[geoboxCount];
        }
    }
}

/**
 * After accumulating the desired number of scanBounds into a geoboxes array (of dimension 4 x N), the caller
 * (of shouldIncludeScan) can then call creatGring  to return a gring_t struct. This code is (mostly) borrowed
 * from l2_generic.c.
 *
 * To create the clockwise gring, it goes forward through the ending lats/lons and backwards through the
 * starting lats/lons. (This approach works if the geoboxes is organized such that the scans are always in the
 * -x direction, where the path is in the +y direction; i.e. the instrument scans "to the left". This function
 * calls normalize to guarantee that this is true for the geoboxesOrdered.)
 *
 * @param geoboxes                geoboxes[4][MAX_GEOBOX_CNT] is a 2D array (borrowed from l2_generic.c),
 * where geoboxes[0][i] is the starting longitude (slon) of the ith scan line to be included in the gring
 * geoboxes[1][i] is the starting latitude  (slat) of the ith scan line to be included in the gring
 * geoboxes[2][i] is the ending   longitude (elon) of the ith scan line to be included in the gring
 *                                  geoboxes[3][i] is the ending   latitude  (elat) of the ith scan line to be
 * included in the gring
 * @param scanDirection        Indicates whether the instrument scans to the left (-1) or to the right (+1).
 * @param gringOut              Output, a struct having arrays gring_lat,gring_lon,gring_seq.
 */
void Gring::createGring(gring_t *gringOut) {
    deleteExtraLine();

    // Declarations for packing gring struct
    std::array<std::array<float, MAX_GEOBOX_CNT>, 4> geoboxesOrdered;
    size_t i, j;

    // Order the geoboxes
    normalize(geoboxesOrdered);

    // Create the gring from the ordered geoboxes
    j = 1;
    // Start at first slon
    gringOut->lons[0] = geoboxesOrdered[0][0];
    // Get elons, going forward
    for (i = 0; i < geoboxCount; i++) {
        gringOut->lons[j++] = geoboxesOrdered[2][i];
    }
    // Get slons, going backwards
    for (i = 0; i < geoboxCount - 1; i++) {
        gringOut->lons[j++] = geoboxesOrdered[0][geoboxCount - 1 - i];
    }
    // gring latitides and sequence numbers
    j = 1;
    // Start at first slat
    gringOut->lats[0] = geoboxesOrdered[1][0];
    gringOut->pointSequence[0] = j;
    // Get elats, going forward
    for (i = 0; i < geoboxCount; i++) {
        gringOut->pointSequence[j] = j + 1;
        gringOut->lats[j++] = geoboxesOrdered[3][i];
    }
    // Get slats, going backwards
    for (i = 0; i < geoboxCount - 1; i++) {
        gringOut->pointSequence[j] = j + 1;
        gringOut->lats[j++] = geoboxesOrdered[1][geoboxCount - 1 - i];
    }
    // Verify
    for (i = 0; i < j; i++) {
        // printf("INFO: createGring: lat=%f lon=%f\n",gringOut->gring_lat[i],gringOut->gring_lon[i]);
    }
    // num_gring_pts
    gringOut->numPoints = j;
}

/**
 * This function creates an ordered geobox (geoboxesOrdered) from the possibly unordered geobox.
 * By ordered, we mean that the endpoints of each scan are ordered in a way that will make it
 * easy for createGring to create a gring in a clockwise order; i.e. geoboxesOrdered is organized such that
 * the scans are always in the -x direction, where the path is in the +y direction; i.e. the instrument scans
 * "to the left".
 *
 * @param normalizedGeoboxes A copy of geoboxes whose directions are all right-to-left
 */
void Gring::normalize(std::array<std::array<float, MAX_GEOBOX_CNT>, 4> &normalizedGeoboxes) {
    // Initialize normalizedGeoboxes = geoboxes
    for (size_t scanNum = 0; scanNum < geoboxCount; scanNum++) {
        for (size_t j = 0; j < 4; j++) {
            normalizedGeoboxes[j][scanNum] = geoboxes[j][scanNum];
        }
    }

    // For certain missions, the instrument scans to the right, so flip the endpoints
    if (scanDirection > 0) {
        for (size_t scanNum = 0; scanNum < geoboxCount; scanNum++) {
            //  Lons
            normalizedGeoboxes[0][scanNum] = geoboxes[2][scanNum];
            normalizedGeoboxes[2][scanNum] = geoboxes[0][scanNum];
            // Lats
            normalizedGeoboxes[1][scanNum] = geoboxes[3][scanNum];
            normalizedGeoboxes[3][scanNum] = geoboxes[1][scanNum];
        }
    }
}

int Gring::getGeoboxCount() {
    deleteExtraLine();
    return (int)geoboxCount;
}

/**
 *
 * @param longitudes: csv of gring longitudes, formatted
 * @param latitudes: csv of gring latitudes
 * @param sequence: csv of gring sequence numbers
 * @param max_string_length: Caller specifies max length for any of the above strings. So, caller must know
 * the formats used in convertFloatArrayToCsv, convertIntArrayToCsv
 * @return 0 on success, -1 on failure
 */
int Gring::getGringStrings(std::string &longitudes, std::string &latitudes, std::string &sequence) {
    int numPoints, currGringPoint;

    if (geoboxCount > 1) {
        // Enough scans for a gring, so keep going
        gring_t *gring = new gring_t;
        createGring(gring);
        // Declare arrays
        float *localLats, *localLons;
        size_t *localPointSequence;
        // Allocate memory for arrays
        numPoints = gring->numPoints;
        localLats = (float *)calloc(numPoints, sizeof(float));
        localLons = (float *)calloc(numPoints, sizeof(float));
        localPointSequence = (size_t *)calloc(numPoints, sizeof(size_t));

        // Unpack the gring_t struct
        for (currGringPoint = 0; currGringPoint < numPoints; currGringPoint++) {
            localLons[currGringPoint] = gring->lons[currGringPoint];
            localLats[currGringPoint] = gring->lats[currGringPoint];
            localPointSequence[currGringPoint] = gring->pointSequence[currGringPoint];
            // printf("INFO: get_gring_strings: lat=%f
            // lon=%f\n",localLats[currGringPoint],localLons[currGringPoint]);
        }

        // Convert arrays into strings and print
        convertFloatArrayToCsv(longitudes, localLons, numPoints);
        convertFloatArrayToCsv(latitudes, localLats, numPoints);
        convertIntArrayToCsv(sequence, localPointSequence, numPoints);

        // Clean up
        delete (gring);
        free(localLats);
        free(localLons);
        free(localPointSequence);
        // return success
        return (EXIT_SUCCESS);
    } else {
        // NOT enough scans for a gring, so return -1
        return (-1);
    }
}

std::vector<std::string> tokenizeCsv(std::string csvList) {
    std::vector<std::string> tokenizedList;
    std::string token;
    std::istringstream inputStream(csvList);

    while (getline(inputStream, token, ',')) {
        tokenizedList.push_back(token);
    }

    return tokenizedList;
}

/**
 * @brief Produce a polygon string that lists the given latitudes and longitudes as lat/lon points
 * @param latitudeCsv A comma-separated list of latitude values
 * @param longitudeCsv A comma-separated list of longitude values
 * @param sequence The ordering of the points
 * @return "POLYGON ((<lat, lon> ...))"
 */
std::string Gring::getWktString(std::string &latitudeCsv, std::string &longitudeCsv, std::string &sequence) {
    std::vector<std::string> latitudeTokens = tokenizeCsv(latitudeCsv);
    std::vector<std::string> longitudeTokens = tokenizeCsv(longitudeCsv);

    if (latitudeTokens.size() != longitudeTokens.size()) {
        std::invalid_argument up("Latitude and longitude CSV lists must be of equal length");
        throw up;
    }

    std::string result = "POLYGON ((";
    for (size_t i = 0; i < latitudeTokens.size(); i++) {
        if (i > 0)
            result += ", ";

        result += latitudeTokens[i] + " " + longitudeTokens[i];
    }
    result += " ))";

    return result;
}

/**
 * @brief Produce a polygon string that lists the given latitudes and longitudes as lat/lon points
 * @return "POLYGON ((<lat, lon> ...))"
 */
std::string Gring::getWktString() {
    std::string longitudesCsv, latitudesCsv, sequence;
    int returnStatus = Gring::getGringStrings(longitudesCsv, latitudesCsv, sequence);

    if (returnStatus != EXIT_SUCCESS) {
        std::runtime_error up("Unable to get Gring strings");
        throw up;
    }

    return getWktString(latitudesCsv, longitudesCsv, sequence);
}

/**
 * To help display gring, convert an array of floats into a csv string, without whitespace, formatted like
 * %6.2f
 * @param myString
 * @param myArray
 * @param array_length
 */
void Gring::convertFloatArrayToCsv(std::string &myString, float *myArray, int array_length) {
    int j;
    char tmp[64];  // One value of myArray, formatted like %6.2f

    myString.clear();
    // Iterate through all but the last, adding a delimiter at the end
    for (j = 0; j < array_length; j++) {
        if (j != 0)
            myString.append(",");
        sprintf(tmp, "%.5f", myArray[j]);
        myString.append(tmp);
    }
}

/**
 * To help display gring, convert an array of ints into a csv string, formatted like %d
 * @param myString
 * @param myArray
 * @param array_length
 */
void Gring::convertIntArrayToCsv(std::string &myString, size_t *myArray, int array_length) {
    int j;
    char tmp[64];

    myString.clear();
    // Iterate through all
    for (j = 0; j < array_length; j++) {
        if (j != 0)
            myString.append(",");
        sprintf(tmp, "%d", (int)myArray[j]);
        myString.append(tmp);
    }
}
