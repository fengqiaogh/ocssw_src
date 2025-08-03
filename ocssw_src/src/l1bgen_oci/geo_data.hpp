#ifndef __GEO_DATA_H__
#define __GEO_DATA_H__

#include <stdlib.h>
#include <vector>
#include <cstdint>
#include <netcdf>

// Mechanism Control Electronics Telemetry
struct MceTlm {
    double comRotRate;
    double lineRate;      // Lines per second
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
    struct GeoLut {
        double masterClock;  // OCI clock
        double mceClock;

        double craftToTilt[3][3];  // Spacecraft to tilt base transformation
        double tiltAxis[3];        // In tilt reference frame
        double tiltAngles[2];      // Tilt angles at either the aft or forward positions
        double tiltHome;
        double tiltToOciMech[3][3];    // Tilt platform to OCI mechanical transformation
        double ociMechToOciOpt[3][3];  // OCI mechanical to optical transformation
        double rtaAxis[3];  // Rotating telescope assembly (RTA) rotation axis in OCI optical reference frame
        double hamAxis[3];  // HAM rotation axis in OCI optical reference frame
        double hamAlongTrackAngles[2];
        double hamCrossTrackAngles[2];
        double rtaEncoderScale;  // RTA encoder conversion to arcseconds
        double hamEncoderScale;  // HAM encoder conversion to arcseconds

        int32_t rtaNadir[2];  // Pulse per revolution (PPR) offset from RTA nadir angle in encoder counts

        double alongTrackPlanarity[5];   // The plane that intersects `acrossTrackPlanarity` at the spacecraft
        double acrossTrackPlanarity[5];  // The plane that intersects `alongTrackPlanarity` at the spacecraft

        netCDF::NcFile *file;

    };  // GeoLut

    bool isDark;
    double auCorrection;         // The Earth's current distance from the Sun as a multiplier of 1 AU
    double earthViewTimeOffset;  // Time of first pixel that sees the Earth (?)
    double unixTimeEnd;          // of this L1B file, derived from input L1A
    double unixTimeStart;        // of this L1B file, derived from input L1A
    size_t numCcdPix;            // Aka numHyperSciPix
    size_t numGoodScans;
    size_t numSwirPix;

    MceTlm mceTelem;
    GeoLut geoLut;

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
    std::vector<short> height;  // ASL per pixel
    std::vector<float> sensorAzimuths;
    std::vector<float> sensorZeniths;
    std::vector<float> solarAzimuths;
    std::vector<float> solarZeniths;
    std::vector<uint8_t> hamSides;
    std::vector<uint8_t> qualityFlag;
};

#endif