#include <stddef.h>
#include <string>
#include <array>
#include <vector>
#include <boost/geometry.hpp>
#include "boost/geometry/geometries/point.hpp"

#ifndef __GRING_H__
#define __GRING_H__

class Gring {
   public:
    /**
     * @brief Instantiates a default Gring object tocompute the geospatial extrema and bounding polygon of a
     * given set of data
     *
     * @param scanInclusionRate The change in latitude across the data necessary for a scan to be included
     * @param invalidPixelMask A bitwise mask to filter pixels out of the extrema and bounding polygon.
     */
    Gring(double scanInclusionRate = 20.0, size_t flowCheckInterval = 20,
          int32_t invalidPixelMask = 0);

    /**
     * @brief Attempts to include a scan in the current G-Ring boundary and updates the lat/lon extrema
     *
     * This method examines the given latitude and longitude data. First it updates the granule-wide
     * lat/lon extrema of all non-flagged pixels. Then, it determines whether this scan should contribute any
     * vertices to the bounding polygon (G-Ring) of the granule.
     *
     * @tparam real A floating-point type representing the precision of the latitude
     *              and longitude arrays (i.e., float or double).
     * @param lats Pointer to an array of latitude values for the current scan line.
     * @param lons Pointer to an array of longitude values for the current scan line.
     * @param numPixels The number of valid pixels in the current scan.
     * @param currScan The zero-based index of the current scan being processed.
     * @param flags Optional pointer to an array of per-pixel flags, which may be used
     *              to skip or mask invalid data. If omitted, all pixels will be considered.
     *
     * @note
     * The granule-wide lat/lon extrema (computed from non-flagged pixels) may lie inside the
     * granule and therefore do not represent its outline. In contrast, the vertices of the G-Ring (bounding
     * polygon) originate from the geospatial edges (typically the first and last valid pixels) of included
     * scans, so the polygon describes the external perimeter instead of the mere min/max.
     */
    template <typename real>
    void tryIncludeScan(real* lats, real* lons, size_t numPixels, size_t currScan, int32_t* flags = nullptr);

    double getGeospatialLongitudeMin() const;
    double getGeospatialLongitudeMax() const;
    double getGeospatialLatitudeMin() const;
    double getGeospatialLatitudeMax() const;
    /** @return latitude min, latitude max, longitude min, longitude max */
    std::array<double, 4> getGeospatialExtremes() const;

    /**
     * @brief Returns the inferred satellite orbital direction based on processed scan lines.
     *
     * The direction is computed incrementally by counting ascending and descending scans over the granule (or
     * a configurable subset of scans). A positive value indicates a net ascending (northbound) motion, a
     * negative value indicates a net descending (southbound) motion, and zero indicates an equal number of
     * ascending and descending scans—typically meaning the granule perfectly spans a pole or the direction
     * cannot be determined yet.
     * @note The direction check interval is configurable; by default, direction is
     *       evaluated every 20 scans.
     * @return A signed integer whose sign indicates the dominant orbital direction:
     *         > 0 for ascending, < 0 for descending, 0 if evenly divided or undetermined.
     */
    int getSatelliteDirection() const;
    int32_t getInvalidPixelMask() const;

    std::string getGeospatialBounds();
    int getGringStrings(std::string& longitudes, std::string& latitudes, std::string& sequence);

    void setDirectionCheckInterval(size_t newInterval);
    void setDeltaLat(double newRate);
    void setInvalidPixelMask(int32_t mask);

   private:
    double scanInclusionInterval;  // Based on latitude
    int32_t invalidPixelMask;
    int satelliteFlow;         // Positive for ascending, negative for descending
    size_t flowCheckInterval;  // Scans between flow checks. Used to include scans when going over poles
    size_t previousFlowCheck;  // The last time the satellite's flow was checked
    size_t previousScanIncluded;
    double prevIncludedLat;  // Set on the first call to processScan and every time the
                             // difference between this and any lat in a give scan meets or
                             // exceeds scanInclusionRate

    // The start, middle, and end of the valid portion of this scan
    std::array<int, 3> validPortionOfScan; // Class level for more intelligent inclusion of scan bounds

    double longitudeMin;
    double longitudeMax;
    double latitudeMin;
    double latitudeMax;

    template <typename real>
    bool shouldIncludeScan(const real* lats, const real* lons, size_t numPixels, size_t currScan,
                           int32_t* flags = nullptr);

    template <typename real>
    void updateExtremes(const real* lats, const real* lons, size_t numPixels);
    template <typename real>
    void updateExtremes(const real* lats, const real* lons, const int32_t* flags, size_t numPixels);

    void addBounds();

    using point_t =
        boost::geometry::model::point<double, 2, boost::geometry::cs::geographic<boost::geometry::degree>>;
    // Polygon of point_t's in a CCW fashion
    using wktPolygon_t = boost::geometry::model::polygon<point_t, false>;
    // Unclosed polygon of point_t's in a clockwise fashion
    using gringPath_t = boost::geometry::model::polygon<point_t, true, false>;

    // lons, lats of the start of each scan that is included
    std::vector<point_t> scanStartPoints;
    // lons, lats of the end of each scan that is included
    std::vector<point_t> scanEndPoints;
    point_t currScanStart;  // lon, lat
    point_t currScanEnd;    // lon, lat
};

#endif