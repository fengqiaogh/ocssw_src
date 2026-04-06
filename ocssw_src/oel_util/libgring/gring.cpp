
#include "gring.h"
#include <algorithm>
#include <cfloat>
#include <sstream>
#include <genutils.h>

using namespace boost::geometry;

Gring::Gring(double scanInclusionRate, size_t flowCheckInterval, int32_t invalidPixelMask) {
    this->scanInclusionInterval = scanInclusionRate;
    this->flowCheckInterval = flowCheckInterval;
    this->invalidPixelMask = invalidPixelMask;
    satelliteFlow = 0;

    latitudeMin = std::numeric_limits<float>::max();
    latitudeMax = std::numeric_limits<float>::lowest();
    longitudeMin = std::numeric_limits<float>::max();
    longitudeMax = std::numeric_limits<float>::lowest();
}

template <typename real>
std::array<int, 3> findValidScanBounds(const real* lats, const real* lons, size_t numPixels,
                                       const int32_t* flags, int32_t invalidPixelMask) {
    int firstGoodPixel = -1;
    int lastGoodPixel = -1;
    int middleGoodPixel = -1;
    size_t endPixel = numPixels - 1;
    size_t centerPixel = numPixels / 2;  // Imprecision with integer division doesn't matter here

    auto isValidLat = [&](const real& lat) { return lat != BAD_FLT && (-90.0 <= lat && lat <= 90.0); };
    auto isValidLon = [&](const real& lon) { return lon != BAD_FLT && (-180.0 <= lon && lon <= 180.0); };

    for (size_t i = 0; i <= endPixel; i++) {
        bool pixelIsValid = false;

        if (flags != nullptr) {
            pixelIsValid = (flags[i] & invalidPixelMask) == 0 && isValidLat(lats[i]) && isValidLon(lons[i]);
        } else {
            pixelIsValid = isValidLat(lats[i]) && isValidLon(lons[i]);
        }

        if (pixelIsValid) {
            if (firstGoodPixel == -1) {
                firstGoodPixel = i;
            }

            lastGoodPixel = i;  // Eventually this will truly be the last good pixel

            if ((i <= centerPixel) || middleGoodPixel == -1) {
                middleGoodPixel = i;
            }
        }
    }

    return {firstGoodPixel, middleGoodPixel, lastGoodPixel};
}


template <typename real>
void Gring::tryIncludeScan(real* lats, real* lons, size_t numPixels, size_t currScan, int32_t* flags) {
    
    // Used in extrema update and in scan inclusion heuristics
    validPortionOfScan = findValidScanBounds(lats, lons, numPixels, flags, invalidPixelMask);
    if (validPortionOfScan[0] == -1 && validPortionOfScan[1] == -1 && validPortionOfScan[2] == -1) {
        return; // Nothing can be done with this scan
    }

    if (flags) {
        updateExtremes(lats, lons, flags, numPixels);
    } else {
        updateExtremes(lats, lons, numPixels);
    }

    bool shouldIncludeScan = Gring::shouldIncludeScan(lats, lons, numPixels, currScan, flags);

    // Storing then conditionally including allows us to catch the final scan without taking in a parameter for
    // the final scan.
    currScanStart = {lons[validPortionOfScan.front()], lats[validPortionOfScan.front()]};
    currScanEnd = {lons[validPortionOfScan.back()], lats[validPortionOfScan.back()]};

    if (shouldIncludeScan) {
        previousScanIncluded = currScan;
        addBounds();
    }
}
template void Gring::tryIncludeScan(float* lats, float* lons, size_t numPixels, size_t currScan,
                                    int32_t* flags);
template void Gring::tryIncludeScan(double* lats, double* lons, size_t numPixels, size_t currScan,
                                    int32_t* flags);

/**
 * @brief Returns true if the given array is non-strictly monotonic (entirely non-decreasing or
 * non-increasing)
 */
template <typename real>
bool arrayIsMonotonic(const real* array, size_t numPixels) {
    if (numPixels == 0 || numPixels == 1) {
        return true;
    }

    float epsilon = 1e-5;
    bool increasing = true;
    bool decreasing = true;
    for (size_t i = 0; i < numPixels - 1; i++) {
        if (array[i] > array[i + 1] + epsilon)
            increasing = false;

        if (array[i] < array[i + 1] + epsilon)
            decreasing = false;

        if (!increasing && !decreasing)
            return false;
    }
    return increasing || decreasing;
}
template bool arrayIsMonotonic<float>(const float*, size_t);
template bool arrayIsMonotonic<double>(const double*, size_t);

template <typename real>
auto assignExtrema = [](double& currentMin, double& currentMax, real newMin, real newMax) {
    currentMax = std::max(currentMax, static_cast<double>(newMax));
    currentMin = std::min(currentMin, static_cast<double>(newMin));
};

template <typename real>
void Gring::updateExtremes(const real* lats, const real* lons, size_t numPixels) {
    auto nonFillCompare = [&](const auto& a, const auto& b) {
        if (a == BAD_FLT)
            return false;
        if (b == BAD_FLT)
            return true;
        return a < b;
    };

    // We don't mind the weirdness that comes at the poles
    double scanMaxLat = *std::max_element(lats, lats + numPixels, nonFillCompare);
    double scanMinLat = *std::min_element(lats, lats + numPixels, nonFillCompare);
    assignExtrema<real>(latitudeMin, latitudeMax, scanMinLat, scanMaxLat);

    bool lonsAreMonotonic = arrayIsMonotonic(lons, numPixels);
    if (lonsAreMonotonic) {  // Does not cross 180
        float scanMaxLon = *std::max_element(lons, lons + numPixels, nonFillCompare);
        float scanMinLon = *std::min_element(lons, lons + numPixels, nonFillCompare);
        assignExtrema<real>(longitudeMin, longitudeMax, scanMinLon, scanMaxLon);
    } else {  // crosses 180
        if (lons[0] > lons[numPixels - 1]) {
            // Scan goes west-to-east, and longitude min should be the 'leftmost' lon with North being up
            assignExtrema<real>(longitudeMin, longitudeMax, lons[validPortionOfScan.front()],
                                lons[validPortionOfScan.back()]);
        } else {  // lons[0] <= lons[numPixels - 1]
            // Scan goes east-to-west, and longitude min should be the 'leftmost' lon with North being up
            assignExtrema<real>(longitudeMin, longitudeMax, lons[validPortionOfScan.back()],
                                lons[validPortionOfScan.front()]);
        }
    }
}
template void Gring::updateExtremes<float>(const float*, const float*, size_t);
template void Gring::updateExtremes<double>(const double*, const double*, size_t);

template <typename real>
void Gring::updateExtremes(const real* lats, const real* lons, const int32_t* flags, size_t numPixels) {
    for (size_t i = 0; i < numPixels; i++) {
        if (lats[i] < latitudeMin && (flags[i] & invalidPixelMask) == 0) {
            latitudeMin = lats[i];
        }
        if (lats[i] > latitudeMax && (flags[i] & invalidPixelMask) == 0) {
            latitudeMax = lats[i];
        }
    }

    if (arrayIsMonotonic(lons, numPixels)) {
        for (size_t i = 0; i < numPixels; i++) {
            if (lons[i] < longitudeMin && (flags[i] & invalidPixelMask) == 0) {
                longitudeMin = lons[i];
            }
            if (lons[i] > longitudeMax && (flags[i] & invalidPixelMask) == 0) {
                longitudeMax = lons[i];
            }
        }
    } else {
        real firstValidLon = lons[validPortionOfScan.front()];
        real lastValidLon = lons[validPortionOfScan.back()];

        if (firstValidLon > lastValidLon) {
            // Scan goes west-to-east, and longitude min should be the 'leftmost' lon with North being up
            assignExtrema<real>(longitudeMin, longitudeMax, firstValidLon, lastValidLon);
        } else {
            // Scan goes east-to-west, and longitude min should be the 'leftmost' lon with North being up
            assignExtrema<real>(longitudeMin, longitudeMax, lastValidLon, firstValidLon);
        }
    }
}
template void Gring::updateExtremes<float>(const float*, const float*, const int32_t*, size_t);
template void Gring::updateExtremes<double>(const double*, const double*, const int32_t*, size_t);

template <typename real>
bool Gring::shouldIncludeScan(const real* lats, const real* lons, size_t numPixels, size_t currScan,
                              int32_t* flags) {
    bool noScansIncludedYet = scanStartPoints.empty() || scanEndPoints.empty();

    bool foundValidPortionOfScan =
        validPortionOfScan[0] != -1 && validPortionOfScan[1] != -1 && validPortionOfScan[2] != -1;

    if (!foundValidPortionOfScan) {
        return false;
    }

    // From here on out, foundValidPortionOfScan is always true

    double middleLatitude = static_cast<double>(lats[validPortionOfScan[1]]);

    if (noScansIncludedYet) {
        prevIncludedLat = middleLatitude;
        previousFlowCheck = currScan;
        return true;
    } else {  // At least one scan included
        if ((currScan - previousFlowCheck) > flowCheckInterval) {
            previousFlowCheck = currScan;

            if ((middleLatitude - prevIncludedLat) > 0) {
                satelliteFlow++;  // Scan ascends
            } else {
                satelliteFlow--;  // Scan descends
            }
        }
    }

    // Either south or north
    bool crossedLatInterval = std::abs(middleLatitude - prevIncludedLat) >= scanInclusionInterval;

    if (crossedLatInterval) {
        prevIncludedLat = middleLatitude;
        return true;
    }

    return false;  // By default, don't include a scan
}
template bool Gring::shouldIncludeScan(const double* lats, const double* lons, size_t numPixels,
                                       size_t currScan, int32_t* flags);
template bool Gring::shouldIncludeScan(const float* lats, const float* lons, size_t numPixels,
                                       size_t currScan, int32_t* flags);

void Gring::addBounds() {
    bool someBoundsExist = (!scanStartPoints.empty() || !scanEndPoints.empty());

    auto arraysAreEqual = [](point_t compareTo, point_t compared) {
        bool index0Same = std::abs(compareTo.get<0>() - compared.get<0>()) <= DBL_EPSILON;
        bool index1Same = std::abs(compareTo.get<1>() - compared.get<1>()) <= DBL_EPSILON;
        return index0Same && index1Same;
    };

    if (someBoundsExist && (arraysAreEqual(scanStartPoints.back(), currScanStart) &&
                            arraysAreEqual(scanEndPoints.back(), currScanEnd))) {
        // Duplicates should not be included
        return;
    }

    scanStartPoints.push_back(currScanStart);
    scanEndPoints.push_back(currScanEnd);
}

std::string Gring::getGeospatialBounds() {
    addBounds();  // Make sure the final scan is included

    if (scanStartPoints.size() < 2 || scanEndPoints.size() < 2) {
        std::cerr << "-W- Unable to make WKT string: Not enough included scans" << std::endl;
        return "";
    }

    wktPolygon_t polygon;

    // The separation of the following two loops ensure the polygon doesn't intersect itself
    for (point_t point : scanStartPoints) {
        append(polygon, point_t(point.get<0>(), point.get<1>()));
    }
    for (auto it = scanEndPoints.rbegin(); it != scanEndPoints.rend(); it++) {
        point_t point = *it;
        append(polygon, point_t(point.get<0>(), point.get<1>()));
    }

    correct(polygon);

    std::string reason;
    if (!is_valid(polygon, reason)) {
        std::cerr << "-W- Unable to make WKT string: " << reason << std::endl;
        return "";
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(5) << wkt(polygon);

    return oss.str();
}

int Gring::getGringStrings(std::string& longitudes, std::string& latitudes, std::string& sequence) {
    addBounds();  // Make sure the final scan is included

    if (scanStartPoints.size() < 2 || scanEndPoints.size() < 2) {
        std::cerr << "-W- Unable to make Gring strings: Not enough included scans" << std::endl;

        longitudes = "";
        latitudes = "";
        sequence = "";
        return -1;
    }

    gringPath_t polygon;

    // The separation of the following two loops ensure the polygon doesn't intersect itself
    for (point_t point : scanStartPoints) {
        append(polygon, point_t(point.get<0>(), point.get<1>()));
    }
    for (auto it = scanEndPoints.rbegin(); it != scanEndPoints.rend(); it++) {
        point_t point = *it;
        append(polygon, point_t(point.get<0>(), point.get<1>()));
    }

    correct(polygon);

    std::string reason;
    if (!is_valid(polygon, reason)) {
        std::cerr << "-W- Unable to make Gring points: " << reason << std::endl;
        return -1;
    }

    std::ostringstream longitudeOss;
    longitudeOss << std::fixed << std::setprecision(5);
    std::ostringstream latitudeOss;
    latitudeOss << std::fixed << std::setprecision(5);
    for (size_t i = 0; i < polygon.outer().size(); i++) {
        point_t point = polygon.outer()[i];

        longitudeOss << get<0>(point);
        latitudeOss << get<1>(point);
        sequence += std::to_string(i + 1);

        if (i < polygon.outer().size() - 1) {  // Don't want a comma at the end
            longitudeOss << ',';
            latitudeOss << ',';
            sequence += ',';
        }
    }

    longitudes = longitudeOss.str();
    latitudes = latitudeOss.str();

    return 0;
}

double Gring::getGeospatialLongitudeMin() const {
    return longitudeMin;
}
double Gring::getGeospatialLongitudeMax() const {
    return longitudeMax;
}
double Gring::getGeospatialLatitudeMin() const {
    return latitudeMin;
}
double Gring::getGeospatialLatitudeMax() const {
    return latitudeMax;
}

std::array<double, 4> Gring::getGeospatialExtremes() const {
    return {latitudeMin, latitudeMax, longitudeMin, longitudeMax};
}

int Gring::getSatelliteDirection() const {
    return satelliteFlow;
}

int32_t Gring::getInvalidPixelMask() const {
    return invalidPixelMask;
}

void Gring::setDirectionCheckInterval(size_t newInterval) {
    flowCheckInterval = newInterval;
}

void Gring::setDeltaLat(double newRate) {
    scanInclusionInterval = newRate;
}

void Gring::setInvalidPixelMask(int32_t mask) {
    invalidPixelMask = mask;
}
