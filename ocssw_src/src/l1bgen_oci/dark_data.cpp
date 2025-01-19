#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include "dark_data.hpp"

using namespace std;

void filterDarkNoise(const size_t numScans, const size_t numBands, const size_t numPixels,
                     const uint32_t fillValue, uint32_t ***darkValues) {
    const float REJECTION_FACTOR = 4.5;  // Equivalent to 6 sigma

    for (size_t i = 0; i < numBands; ++i) {
        // Flatten the 2D array (dark[*,i,*]) into a 1D vector to make use of some standard functions
        vector<uint32_t> bandData;
        for (size_t j = 0; j < numScans; ++j) {
            for (size_t k = 0; k < numPixels; ++k) {
                bandData.push_back(darkValues[j][i][k]);
            }
        }

        // The data over which we are finding the median are integers, so
        // we don't need to handle the size % 2 == 0 case.
        size_t medianIndex = bandData.size() / 2;
        sort(bandData.begin(), bandData.end());
        uint32_t darkMedian = bandData[medianIndex];

        int32_t size = bandData.size();
        int32_t q1Index = size / 4;
        int32_t q3Index = 3 * size / 4;
        int32_t q1 = bandData[q1Index];
        int32_t q3 = bandData[q3Index];
        int32_t interQuartileRange = q3 - q1;

        float outlierThreshold = REJECTION_FACTOR * interQuartileRange;

        // Filter out noise
        for (size_t j = 0; j < numScans; ++j) {
            for (size_t k = 0; k < numPixels; ++k) {
                // Floating point to avoid mixing integer and floating point arithmetic
                uint32_t medianDeviation = abs((int32_t)darkValues[j][i][k] - (int32_t)darkMedian);

                if (medianDeviation > outlierThreshold) {
                    darkValues[j][i][k] = fillValue;
                }
            }
        }
    }
}
