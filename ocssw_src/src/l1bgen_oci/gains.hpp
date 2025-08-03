#ifndef __L1B_GAINS__
#define __L1B_GAINS__
#include <stdint.h>
#include <array>
#include "types.hpp"

enum GainDims {
    T,  // Time
    TEMP,
    TEMP_CORR,      // Temperature correction
    RVS,            // Response Versus Scan
    NL,             // Nonlinearity
    MS,             // Mirror side. Which of either side of the half angle mirror
    GAIN_DIMS_SIZE  // Sentinel value
};

typedef struct Gains {
    std::array<uint16_t, GAIN_DIMS_SIZE> dimensions = {0, 0, 0, 0, 0, 0};  // The shape of this LUT
    vec2D<float> k1K2;
    vec3D<float> k3Coefs;
    vec3D<float> k4Coefs;
    vec2D<double> k5Coefs;
    std::vector<uint32_t> saturationThresholds;
} Gains;

#endif