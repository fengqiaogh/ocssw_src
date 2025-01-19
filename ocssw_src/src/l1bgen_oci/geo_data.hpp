#ifndef __GEO_DATA_H__
#define __GEO_DATA_H__

#include <stdlib.h>
#include <vector>

// Mechanism Control Electronics Telemetry
struct MceTlm {
    double comRotRate;
    double lineRate;  // Lines per second
    double **hamEncData;  // Half angle mirror assembly data
    double **rtaEncData;  // Rotating telescope assembly data
    int pprOffset;
    short mceBoardId;
    size_t numEncoderChannels;
    size_t numMceScans;
    std::vector<int32_t> mceSpinIds;
    std::vector<uint8_t> encoderSampleCounts;
};

// The side effects of geolocation, used in many places
struct GeoData {
    bool isDark;
    double auCorrection;   // The Earth's current distance from the Sun as a multiplier of 1 AU
    double earthViewTimeOffset; // Time of first pixel that sees the Earth (?)
    double unixTimeEnd;    // of this L1B file, derived from input L1A
    double unixTimeStart;  // of this L1B file, derived from input L1A
    size_t numCcdPix;  // Aka numHyperSciPix
    size_t numGoodScans;
    size_t numSwirPix;

    MceTlm mceTelem;
    std::vector<double> earthViewTimes;
    std::vector<double> scanStartTimes;
    std::vector<double> sciPixOffset;
    std::vector<double> swirPixOffset;
    std::vector<double> ccdScanAngles;
    std::vector<double> swirScanAngles;
    std::vector<float> pixelLatitudes;
    std::vector<float> pixelLongtitudes;
    std::vector<float> scienceLines;
    std::vector<float> swirLines;
    std::vector<int32_t> spinIds;
    std::vector<short> height; // ASL per pixel
    std::vector<short> sensorAzimuths;
    std::vector<short> sensorZeniths;
    std::vector<short> solarAzimuths;
    std::vector<short> solarZeniths;
    std::vector<uint8_t> hamSides;
    std::vector<uint8_t> qualityFlag;
};

#endif