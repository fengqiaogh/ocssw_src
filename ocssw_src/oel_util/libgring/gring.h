#ifndef _GRING_HELPER_H
#define _GRING_HELPER_H

#ifndef MAX_GEOBOX_CNT
#define MAX_GEOBOX_CNT 100
#endif

#include <string>
#include <vector>
#include <array>

//  The member functions of this class can be used by the CALLER to create a gring
//  Call sequence:
//      CALLER (instantiation)  > Gring
//      CALLER (per scan)       > processScan      > shouldIncludeScan    > setScanBounds
//      CALLER                  > getGeoboxCount
//      CALLER                  > getGringStrings > createGring           > normalize
//                                                  > convertFloatArrayToCsv
//                                                  > convertIntArrayToCsv

class Gring {
   public:
    Gring();
    Gring(float latPointDelta, int scanDirection, size_t directionCheckInterval);
    Gring(int sensorID, float deltaLatDegrees);
    virtual ~Gring();  // Destructor
    void processScan(float *lat, float *lon, int32_t npix, size_t recnum, size_t escan,
                     int32_t *flags = nullptr);
    void processScan(std::vector<float> &lat, std::vector<float> &lon, std::vector<int32_t> &flags,
                     int32_t npix, size_t recnum, size_t escan);
    int getGeoboxCount();
    int getGringStrings(std::string &lons, std::string &lats, std::string &sequence);
    std::string getWktString(std::string &latitudeCsv, std::string &longitudeCsv, std::string &sequence);
    std::string getWktString();

    void setDeltaLat(float newDeltaLat);
    void setDirectionCheckInterval(size_t);
    void setScanDirection(int newScanDirection);
    void setInvalidPixelMask(int32_t newMask);

    float getDeltaLat();
    size_t getDirectionCheckInterval();
    int getScanDirection();

    std::array<float, 4> getGeospatialExtremes();
    float getGeospatialLatitudeMin();
    float getGeospatialLatitudeMax();
    float getGeospatialLongitudeMin();
    float getGeospatialLongitudeMax();

   private:
    struct LatLonPoint {
        float lat;
        float lon;
    };
    // The scanBounds structure is used to mark the end-points (and mid-point) of a scan that can contribute
    // to a gring
    struct ScanBounds {
        LatLonPoint start;
        LatLonPoint center;
        LatLonPoint end;
    };

    // gring_t struct
    struct gring_t {
        float lats[2 * MAX_GEOBOX_CNT];
        float lons[2 * MAX_GEOBOX_CNT];
        size_t pointSequence[2 * MAX_GEOBOX_CNT];
        int numPoints;
    };

    // State information
    bool foundGoodScan;
    int prevTrackDirection;
    float centerLatitude;  // Declared at class scope to keep the number of geoboxes down
    std::array<std::array<float, MAX_GEOBOX_CNT>, 4> geoboxes;  // Array 4 x MAX_GEOBOX_CNT
    std::vector<float> centerLats;
    size_t geoboxCount;
    size_t directionCheckPrevious;  // last line where we checked the direction
    size_t directionCheckInterval;  // number of line to wait before checking the direction
    float latitudeMin;
    float latitudeMax;
    float longitudeMin; // Easternmost longitude, unless this crosses 180
    float longitudeMax; // Westernmost longitude, unless this crosses 180
    float previousLon; // Helps with determining whether this crossed 180

    // Config
    float deltaLatDegrees;
    int scanDirection;
    void initialize();
    int32_t invalidPixelMask;

    bool isValidLat(float lat);
    bool isValidLon(float lon);

    bool shouldIncludeScan(float *lat, float *lon, size_t numPixels, size_t currScan, size_t lastScan,
                           int32_t *flags, ScanBounds &scanBounds);

    void deleteExtraLine();

    void createGring(gring_t *this_gring);

    void normalize(std::array<std::array<float, MAX_GEOBOX_CNT>, 4> &geoboxOrdered);

    void convertFloatArrayToCsv(std::string &myString, float *myArray, int array_length);

    void convertIntArrayToCsv(std::string &myString, size_t *myArray, int array_length);

    void updateExtremes(const ScanBounds &scanBounds);
};

#endif
