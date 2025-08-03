#ifndef AREA_WEIGHTING_H
#define AREA_WEIGHTING_H
#include <vector>
#include <optional>
#include <l2_reader.hpp>
class AreaWeighting {
    // area weighting variables
    enum L2PixelMode_t { L2PixelOff, L2PixelCorner, L2PixelDelta };
    L2PixelMode_t mode{L2PixelOff};
    std::vector<float> delta_lat;           // delta lat for corner of pixel
    std::vector<float> delta_lon;           // delta lon for corner of pixel
    std::vector<float> upperCornerLat;           // lat1 for corner of pixel
    std::vector<float> upperCornerLon;           // lon1 for corner of pixel
    std::vector<float> lowerCornerLat;           // lat2 for corner of pixel
    std::vector<float> lowerCornerLon;           // lon2 for corner of pixel
    bool lat2Valid{false};             // does lat2 and lon2 contain valid coordinates
    std::optional<size_t> lastLine;  // if no value the last line does not have good lat/lon
    size_t nsamp{0};
    std::vector<float> lastLat;
    std::vector<float> lastLon;
    L2_Reader *l2file{nullptr};
    float median_start_value{1000000};
    std::vector<float> dlat;
    std::vector<float> dlon;
    float high_lat_value{73.0f};
    float scale_factor_min{4.0f};
    float scale_factor_min_2{16.0};
    float scale_factor_max{30.0f};
    float scale_factor_max_2{60.0f};
    float *latitude{nullptr}, *longitude{nullptr};
    int area_mode{0};
    /**
     * @brief Sets scan line data by reading latitude and longitude values
     * @param iscan Index of the scan line to read
     * @param mask Vector of mask values for the scan line
     */
    void set_scan(size_t iscan, std::vector<uint8_t> &mask);

    /**
     * @brief Calculates pixel delta values for latitude and longitude
     * @param lat0 First row latitude values
     * @param lon0 First row longitude values  
     * @param lat1 Second row latitude values
     * @param lon1 Second row longitude values
     * @param latOut Output latitude delta values
     * @param lonOut Output longitude delta values
     * @param numPoints Number of points in each row
     */
    void calculatePixelDeltas(float *lat0, float *lon0, float *lat1, float *lon1, float *latOut,
                              float *lonOut, int32_t numPoints);

    /**
     * @brief Extrapolates pixel corner coordinates
     * @param lat0 First row latitude values
     * @param lon0 First row longitude values
     * @param lat1 Second row latitude values
     * @param lon1 Second row longitude values
     * @param latOut Output latitude corner values (has one more element than input)
     * @param lonOut Output longitude corner values (has one more element than input)
     * @param numPoints Number of points in each input row
     */
    void extrapolatePixelCorners(float *lat0, float *lon0, float *lat1, float *lon1, float *latOut,
                                 float *lonOut, int32_t numPoints);

    /**
     * @brief Interpolates pixel corner coordinates
     * @param lat0 First row latitude values
     * @param lon0 First row longitude values
     * @param lat1 Second row latitude values
     * @param lon1 Second row longitude values
     * @param latOut Output latitude corner values (has one more element than input)
     * @param lonOut Output longitude corner values (has one more element than input)
     * @param numPoints Number of points in each input row
     */
    void interpolatePixelCorners(float *lat0, float *lon0, float *lat1, float *lon1, float *latOut,
                                 float *lonOut, int32_t numPoints);

   public:
    AreaWeighting() = default;
    AreaWeighting(int area_mode);
    
    /**
     * @brief Sets the l2 reader pointer
     * @param reference to an L2 reader
     */
    void set_l2_reader(L2_Reader &l2file) {
        this->l2file = &l2file;
    }
    /**
     * @brief Get the current pixel mode
     * @return The pixel mode as an integer (0=Off, 1=Corner, 2=Delta)
     */
    int get_mode() const;

    /**
     * @brief Check if geolocation data is valid for a given pixel
     * @param ipixel Index of the pixel to check
     * @return true if geolocation is valid, false otherwise
     */
    bool valid_geolocation(size_t ipixel) const;

    /**
     * @brief Calculate and set the corner coordinates for pixels in a scan line
     * @param iscan Index of the scan line
     * @param mask Vector of mask values for the scan line
     */
    void set_corners(size_t iscan, std::vector<uint8_t> &mask);

    /**
     * @brief Get the latitude value for a pixel
     * @param ipixel Index of the pixel
     * @return The latitude value
     */
    float get_latitude(size_t ipixel) const;

    /**
     * @brief Get the longitude value for a pixel
     * @param ipixel Index of the pixel
     * @return The longitude value
     */
    float get_longitude(size_t ipixel) const;

    /**
     * @brief Get the first latitude corner value for a pixel
     * @param ipixel Index of the pixel
     * @return The lat1 corner value
     */
    float get_upperCornerLat(size_t ipixel) const;

    /**
     * @brief Get the first longitude corner value for a pixel
     * @param ipixel Index of the pixel
     * @return The lon1 corner value
     */
    float get_upperCornerLon(size_t ipixel) const;

    /**
     * @brief Get the second latitude corner value for a pixel
     * @param ipixel Index of the pixel
     * @return The lat2 corner value
     */
    float get_lowerCornerLat(size_t ipixel) const;

    /**
     * @brief Get the second longitude corner value for a pixel
     * @param ipixel Index of the pixel
     * @return The lon2 corner value
     */
    float get_lowerCornerLon(size_t ipixel) const;

    /**
     * @brief Get the latitude delta value for a pixel
     * @param ipixel Index of the pixel
     * @return The latitude delta value
     */
    float get_delta_lat(size_t ipixel) const;

    /**
     * @brief Get the longitude delta value for a pixel
     * @param ipixel Index of the pixel
     * @return The longitude delta value
     */
    float get_delta_lon(size_t ipixel) const;

    /**
     * @brief Reset all internal buffers and state
     */
    void reset();
};
#endif  // AREA_WEIGHTING_H
