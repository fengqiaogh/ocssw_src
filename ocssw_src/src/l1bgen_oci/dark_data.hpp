#ifndef __L1B_DARK__
#define __L1B_DARK__

struct DarkData {
    uint16_t numScansAvg = 1;  // Number of dark scans to average
    uint16_t numPixSkip = 0;   // Number of dark pixels to skip at the beginning of the file
    uint16_t numPix;
    int16_t darkZone = -1;  // Start of dark collect in terms of spatial zones
    int16_t spatialAgg;
    uint32_t ***data;                 // The actual dark data
    std::vector<double> corrections;  // Should be sized in accordance with the number of instrument bands
};

/**
 * @brief Filters dark noise from the given dark values.
 *
 * This function processes the dark values to remove noise, improving the quality of the dark data.
 *
 * @param numScans The number of scans in the dark data.
 * @param numBands The number of spectral bands in the dark data.
 * @param numPixels The number of pixels per scan line.
 * @param fillValue The value used to represent missing or invalid data.
 * @param darkValues A 3D array containing the dark values to be filtered.
 *                   The dimensions are [numScans][numBands][numPixels].
 */
void filterDarkNoise(const size_t numScans, const size_t numBands, const size_t numPixels,
                     const uint32_t fillValue, uint32_t ***darkValues);

#endif