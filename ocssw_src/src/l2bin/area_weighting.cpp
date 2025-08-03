#include "area_weighting.h"

AreaWeighting::AreaWeighting(int area_mode) : area_mode(area_mode) {
    switch (area_mode) {
        case 0:
            mode = L2PixelOff;
            break;
        case 1:
            mode = L2PixelDelta;
            break;
        default:
            mode = L2PixelCorner;
            break;
    }
}

int AreaWeighting::get_mode() const {
    return area_mode;
}

void AreaWeighting::set_scan(size_t iscan, std::vector<uint8_t> &mask) {
    if (iscan > 0) {
        mask[iscan - 1] = scan_within_the_group;
        std::vector<size_t> start = {iscan - 1, 0};
        std::vector<size_t> count = {1, nsamp};
        l2file->readLatitude(lastLat.data(), start, count, mask);
        l2file->readLongitude(lastLon.data(), start, count, mask);
    } else {
        size_t forward_line = 1;  // 1
        mask[forward_line] = scan_within_the_group;
        std::vector<size_t> start = {forward_line, 0};
        std::vector<size_t> count = {1, nsamp};
        l2file->readLatitude(lastLat.data(), start, count, mask);
        l2file->readLongitude(lastLon.data(), start, count, mask);
    }
    lastLine = iscan;
    l2file->readLatitudeScan(&latitude, iscan, mask);
    l2file->readLongitudeScan(&longitude, iscan, mask);
    if (mode == L2PixelCorner) {
        // Find center sample index
        int i_center = nsamp / 2;
        
        // Set initial scale factor for latitude differences
        float scale_factor = scale_factor_min_2;
        
        // Use larger scale factor for high latitudes
        if (std::abs(latitude[i_center]) > high_lat_value)
            scale_factor = scale_factor_max_2;
            
        bool use_last_delta_lat = false;
        
        // Check if latitude differences are too large compared to previous deltas
        for (size_t i = 0; i < nsamp; i++) {
            if (std::abs(dlat[i]) * scale_factor < std::abs(latitude[i] - lastLat[i]) && dlat[i] != 0) {
                lat2Valid = false;
                use_last_delta_lat = true;
                break;
            }
        }
        
        if (use_last_delta_lat) {
            // If differences too large, calculate new last lat/lon using previous deltas
            for (size_t i = 0; i < nsamp; i++) {
                lastLat[i] = latitude[i] - dlat[i];
                lastLon[i] = longitude[i] - dlon[i];
            }
        } else {
            // Otherwise calculate new deltas from current differences
            for (size_t i = 0; i < nsamp; i++) {
                dlat[i] = latitude[i] - lastLat[i];
                dlon[i] = longitude[i] - lastLon[i];
            }
        }
    }
}

void AreaWeighting::set_corners(size_t iscan, std::vector<uint8_t> &mask) {
    // Check if L2 file has been set
    if (l2file == nullptr) {
        EXIT_LOG(std::cerr << "-Error-: L2 file has not been set")
    }
    // Return early if pixel mode is off
    if (mode == L2PixelOff)
        return;
    // Initialize arrays if this is first scan
    if (lastLat.empty()) {
        size_t nrec;
        l2file->getDimensions(nrec, nsamp);
        lastLat.resize(nsamp, BAD_FLT);
        lastLon.resize(nsamp, BAD_FLT);

        if (mode == L2PixelCorner) {
            // For corner mode, allocate corner arrays and initialize deltas
            upperCornerLat.resize(nsamp + 1);
            upperCornerLon.resize(nsamp + 1);
            lowerCornerLat.resize(nsamp + 1);
            lowerCornerLon.resize(nsamp + 1);
            dlat.resize(nsamp,median_start_value);
            dlon.resize(nsamp,median_start_value);
        } else if (mode == L2PixelDelta) {
            // For delta mode, allocate delta arrays and initialize center value
            delta_lat.resize(nsamp + 1);
            delta_lon.resize(nsamp + 1);
            delta_lat[nsamp / 2] = median_start_value;
        }
    }
    // Set scan data based on whether this is first scan
    if (!lastLine.has_value()) {
        set_scan(iscan, mask);
    } else {
        // Check if scan lines are sequential
        if (iscan - lastLine.value() != 1)
            lat2Valid = false;
        set_scan(iscan, mask);
    }
    if (mode == L2PixelCorner) {
        // Handle first scan differently since no previous scan exists
        if (iscan == 0) {
            // For first scan, extrapolate upper corners and interpolate lower corners
            extrapolatePixelCorners(lastLat.data(), lastLon.data(), latitude, longitude, upperCornerLat.data(),
                                    upperCornerLon.data(), nsamp);
            interpolatePixelCorners(latitude, longitude, lastLat.data(), lastLon.data(), lowerCornerLat.data(),
                                    lowerCornerLon.data(), nsamp);
        } else {
            // For subsequent scans, reuse previous corners if valid
            if (lat2Valid) {
                // Swap previous lower corners to upper corners
                upperCornerLat.swap(lowerCornerLat);
                upperCornerLon.swap(lowerCornerLon);
            } else {
                // Otherwise interpolate new upper corners
                interpolatePixelCorners(lastLat.data(), lastLon.data(), latitude, longitude, upperCornerLat.data(),
                                        upperCornerLon.data(), nsamp);
            }
            // Extrapolate new lower corners
            extrapolatePixelCorners(lastLat.data(), lastLon.data(), latitude, longitude, lowerCornerLat.data(),
                                    lowerCornerLon.data(), nsamp);
        }
    } else {
        // For delta mode, calculate pixel deltas between scans
        calculatePixelDeltas(latitude, longitude, lastLat.data(), lastLon.data(), delta_lat.data(),
                             delta_lon.data(), nsamp);
    }

    lat2Valid = true;
}

bool AreaWeighting::valid_geolocation(size_t ipixel) const {
    if (mode == L2PixelOff)
        return true;
    if (mode == L2PixelDelta) {
        if (!valid_lat(lastLat[ipixel]))
            return false;
        if (ipixel < nsamp - 1 && !valid_lon(longitude[ipixel + 1]))
            return false;
    }
    if (mode == L2PixelCorner) {
        if (!valid_lat(upperCornerLat[ipixel]))
            return false;
        if (!valid_lon(upperCornerLon[ipixel]))
            return false;
        if (!valid_lat(lowerCornerLat[ipixel]))
            return false;
        if (!valid_lon(lowerCornerLon[ipixel]))
            return false;
        if (!valid_lat(upperCornerLat[ipixel + 1]))
            return false;
        if (!valid_lon(upperCornerLon[ipixel + 1]))
            return false;
        if (!valid_lat(lowerCornerLat[ipixel + 1]))
            return false;
        if (!valid_lon(lowerCornerLon[ipixel + 1]))
            return false;
    }
    return true;
}
void clampDeltaLon(float *deltaLon) {
    if (*deltaLon > 90) {
        if (*deltaLon > 270)
            *deltaLon -= 360.0;
        else
            *deltaLon -= 180.0;
    } else if (*deltaLon < -90) {
        if (*deltaLon < -270)
            *deltaLon += 360.0;
        else
            *deltaLon += 180.0;
    }
}


void AreaWeighting::interpolatePixelCorners(float *lat0, float *lon0, float *lat1, float *lon1, float *latOut,
                                            float *lonOut, int32_t numPoints) {
    // calc the delta for first point
    float dLat = (lat1[0] - lat0[0] + lat0[1] - lat0[0]) / 2.0;
    float dLon = (lon1[0] - lon0[0] + lon0[1] - lon0[0]) / 2.0;
    clampDeltaLon(&dLon);

    // first calc 0 point
    latOut[0] = lat1[0] - dLat;
    lonOut[0] = lon1[0] - dLon;

    // calc the rest of the points
    for (int i = 0; i < numPoints - 1; i++) {
        dLat = (lat1[i] - lat0[i] + lat0[i + 1] - lat0[i]) / 2.0;
        dLon = (lon1[i] - lon0[i] + lon0[i + 1] - lon0[i]) / 2.0;
        clampDeltaLon(&dLon);

        latOut[i + 1] = lat0[i] + dLat;
        lonOut[i + 1] = lon0[i] + dLon;
    }

    // calc last point
    latOut[numPoints] = lat0[numPoints - 1] + dLat;
    lonOut[numPoints] = lon0[numPoints - 1] + dLon;
}


void AreaWeighting::extrapolatePixelCorners(float *lat0, float *lon0, float *lat1, float *lon1, float *latOut,
                                            float *lonOut, int32_t numPoints) {
    // calc the delta for first point
    float dLat = (lat1[1] - lat1[0] + lat0[0] - lat1[0]) / 2.0;
    float dLon = (lon1[1] - lon1[0] + lon0[0] - lon1[0]) / 2.0;
    clampDeltaLon(&dLon);

    // first calc 0 point
    latOut[0] = lat1[0] - dLat;
    lonOut[0] = lon1[0] - dLon;

    // calc the rest of the points
    for (int i = 0; i < numPoints - 1; i++) {
        dLat = (lat1[i] - lat0[i] + lat0[i + 1] - lat0[i]) / 2.0;
        dLon = (lon1[i] - lon0[i] + lon0[i + 1] - lon0[i]) / 2.0;
        clampDeltaLon(&dLon);

        latOut[i + 1] = lat1[i] + dLat;
        lonOut[i + 1] = lon1[i] + dLon;
    }

    // calc last point
    latOut[numPoints] = lat1[numPoints - 1] + dLat;
    lonOut[numPoints] = lon1[numPoints - 1] + dLon;
}


void AreaWeighting::calculatePixelDeltas(float *lat0, float *lon0, float *lat1, float *lon1, float *latOut,
                                         float *lonOut, int32_t numPoints) {
    // bail if the new delta is more than 2x bigger than the previous line
    int i_center = numPoints / 2;
    float scale_factor = scale_factor_min;
    if (fabs(lat0[i_center]) > high_lat_value)
        scale_factor = scale_factor_max;
    if (fabsf(lat0[i_center] - lat1[i_center]) > latOut[i_center] * scale_factor && latOut[i_center] != 0)
        return;

    // calc all points except the last one
    for (int i = 0; i < numPoints - 1; i++) {
        latOut[i] = fabsf(lat0[i] - lat1[i]) / 2.0;
        lonOut[i] = fabsf(lon0[i] - lon0[i + 1]) / 2.0;
        if (lonOut[i] > 90) {
            lonOut[i] = 180.0 - lonOut[i];
        }
    }

    // calc last point
    latOut[numPoints - 1] = latOut[numPoints - 2];
    lonOut[numPoints - 1] = lonOut[numPoints - 2];
}

float AreaWeighting::get_latitude(size_t ipixel) const {
    return latitude[ipixel];
};
float AreaWeighting::get_longitude(size_t ipixel) const {
    return longitude[ipixel];
};
float AreaWeighting::get_upperCornerLat(size_t ipixel) const {
    return upperCornerLat[ipixel];
};
float AreaWeighting::get_upperCornerLon(size_t ipixel) const {
    return upperCornerLon[ipixel];
};
float AreaWeighting::get_lowerCornerLat(size_t ipixel) const {
    return lowerCornerLat[ipixel];
};
float AreaWeighting::get_lowerCornerLon(size_t ipixel) const {
    return lowerCornerLon[ipixel];
};

void AreaWeighting::reset() {
    latitude = nullptr;
    longitude = nullptr;
    free_vector(upperCornerLat);
    free_vector(upperCornerLon);
    free_vector(lowerCornerLat);
    free_vector(lowerCornerLon);
    free_vector(lastLat);
    free_vector(lastLon);
    lastLine.reset();
    free_vector(delta_lat);
    free_vector(delta_lon);
    lat2Valid = false;
}
float AreaWeighting::get_delta_lat(size_t ipixel) const {
    return delta_lat[ipixel];
}
float AreaWeighting::get_delta_lon(size_t ipixel) const {
    return delta_lon[ipixel];
}
