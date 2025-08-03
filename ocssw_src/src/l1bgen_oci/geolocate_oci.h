/**
 * @file geolocate_oci.h
 * @brief Header file containing functions for OCI geolocation
 *
 * This file contains function declarations for various geolocation operations
 * specific to the Ocean Color Instrument (OCI). It includes functions for
 * coordinate transformations, interpolations, sun vector calculations,
 * and pixel geolocation.
 *
 * @authors Joel Gales (SAIC), Jakob C Lindo (SSAI)
 * @date Aug 2024
 */

/**
 * @file geolocate_oci.h
 * @brief Header file containing functions for OCI geolocation
 *
 * This file contains function declarations for various geolocation operations
 * specific to the Ocean Color Instrument (OCI). It includes functions for
 * coordinate transformations, interpolations, sun vector calculations,
 * and pixel geolocation.
 *
 * @authors Joel Gales (SAIC), Jakob C Lindo (SSAI)
 * @date Aug 2024
 */

#ifndef _GEOLOCATE_OCI_H_
#define _GEOLOCATE_OCI_H_

#include "l1b_file.hpp"
#include "geo_data.hpp"
#include <genutils.h>
#include "types.hpp"
#include <array>
#include "l1b_options.hpp"
#include "processing_tracker.hpp"

double constexpr RADIANS_TO_ARCSECONDS = OEL_RADEG * 3600;
const int constexpr MAX_ENC_COUNT = 0x20000;

enum DataType {
    NO_DATA,
    EARTH,
    DARK_CALIBRATION,
    SOLAR_DAILY_CALIBRATION,
    SOLAR_MONTHLY_CALIBRATION,
    RESPONSE_CURVE,
    LUNAR,
    DIAGNOSTIC,
    EARTH_SPECTRAL,
    NO_PROCESSING,
    EXTERNAL_SNAPSHOP_TRIGGER,  // sic
    INTERNAL_SNAPSHOP_TRIGGER,  // sic
    SPCA,
    LUNAR_STARE,
    NON_BASELINE_SPECTRAL_AGGREGATION
};

using namespace netCDF;
using namespace netCDF::exceptions;

typedef float quaternion[4];
typedef std::array<float, 3> orbArray;

/**
 * @brief Expands environment variables in a given string.
 *
 * This function searches for environment variables in the input string,
 * denoted by a '$' prefix and delimited by '/'. If found, it replaces
 * the variable with its corresponding value from the environment.
 *
 * @param sValue Pointer to the string to be processed.
 * @return int Returns 0 on success. Exits the program if fails
 * @note If the environment variable is not defined, the function prints an error message and exits.
 */
inline int expandEnvVar(std::string *sValue) {
    if ((*sValue).find_first_of("$") == std::string::npos)
        return 0;
    std::string::size_type posEndIdx = (*sValue).find_first_of("/");
    if (posEndIdx == std::string::npos)
        return 0;
    const std::string envVar = sValue->substr(1, posEndIdx - 1);
    char *envVar_str = getenv(envVar.c_str());
    if (envVar_str == 0x0) {
        printf("Environment variable: %s not defined.\n", sValue->c_str());
        exit(1);
    }
    *sValue = envVar_str + (*sValue).substr(posEndIdx);

    return 0;
}

/**
 * @brief Geolocate an L1A file, producing an L1B file.
 * @param l1aFile The L1A file to be geolocated
 * @param l1bFile The L1B file object to which geolocation data will be written
 * @param locatingContext Indicator as to which kind of location will happen
 * @param processingTracker A ProcessingTracker
 * @return GeoData struct
 */
GeoData locateOci(netCDF::NcFile *l1aFile, Level1bFile &l1bFile, const oel::L1bOptions &options,
                  LocatingContext locatingContext, oel::ProcessingTracker &processingTracker);

/**
 * @brief Read mechanism control electronics (MCE) Telemetry and compute related variables
 *
 * This function reads MCE telemetry data from the L1A file and computes various parameters
 * necessary for geolocation. It processes data related to the rotating telescope assembly (RTA)
 * and half-angle mirror (HAM).
 *
 * @param l1aFile Pointer to the L1A netCDF file
 * @param geoLut Reference to the GeoLut structure containing geolocation lookup table data
 * @param egrData NcGroup object representing the engineering_data group in PACE OCI L1A files
 *
 * @return MceTlm structure
 *
 * @note This function populates the GeoLut structure with relevant data read from the L1A file.
 */
MceTlm readMceTelemetry(netCDF::NcFile *l1aFile, GeoData::GeoLut &geoLut, netCDF::NcGroup egrData);
/**
 * @brief Determine Earth view pixels and time offsets from Pulse per revolution (PPR) using scan mode table
 * fields to determine pixel offset
 * @param comRotRate Commanded rotation rate of the telescope in revolutions per second
 * @param dataTypes The types of
 * @param spatialZoneLines
 * @param spatialAggregation
 * @param numHyperSciPix Number of hyperspectral science pixels
 * @param numSwirPixels Number of short-wave infrared science pixels
 * @param earthViewTimeOffset Time offset from Pulse per revolution (PPR) to center of Earth View
 * @param scienceLines Line numbers corresponding to science pixel centers
 * @param swirLines Line numbers corresponding to short-wave infrared pixel centers
 * @param sciencePixelOffset Time offset from Earth view mid time to each science pixel
 * @param swirPixelOffset Time offset from Earth View mid time to each short-wave infrared pixel
 * @param isDark Indicates whether this collection is dark
 * @param returnStatus 0 if normal, 1 if no Earth view
 * @return 0 if normal, 1 if no Earth view
 */
int getEarthView(double comRotRate, int16_t *dataTypes, int16_t *spatialZoneLines,
                 int16_t *spatialAggregation, uint16_t &numHyperSciPix, uint16_t &numSwirPixels,
                 double &earthViewTimeOffset, float *scienceLines, float *swirLines,
                 double *sciencePixelOffset, double *swirPixelOffset, bool isDark);

/**
 * @brief Generates OCI Earth view vectors for one spin using MCE telemetry and encoder data.
 *
 * This function calculates the Earth view vectors for a single spin of the Ocean Color Instrument (OCI)
 * using Mechanism Control Electronics (MCE) telemetry and encoder data. It applies corrections based on
 * the Half-Angle Mirror (HAM) and Rotating Telescope Assembly (RTA) encoder data to determine accurate
 * pixel line-of-sight vectors.
 *
 * @param geoData Structure containing geolocation data
 * @param geoLut Structure containing geolocation lookup table data
 * @param numPixels Number of pixels in Earth view
 * @param earthViewTimeOffset Time offset from Pulse per revolution (PPR) to center of Earth View
 * @param delt Vector of time offsets in seconds from Earth View time to pixel times
 * @param scanNum The current scan number being processed
 * @param vectors [out] View vectors in the instrument frame (pre-allocated 2D array)
 * @param scanAngles [out] Scan angles for science pixels
 *
 * @return int Returns EXIT_SUCCESS upon success, EXIT_FAILURE if no MCE data for the spin
 *
 * @note This function implements the algorithm described in "Use of OCI Telemetry to Determine Pixel
 * Line-of-Sight" by F. Patt, 2020-05-18
 */
int getEarthViewVectors(GeoData &geoData, const GeoData::GeoLut &geoLut, const int device,
                        const double earthViewTimeOffset, const std::vector<double> &delt,
                        const size_t scanNum, std::vector<std::array<float, 3>> &vectors);

/**
 * @brief Interpolate encoder data to pixel times and add to scan angles
 * @param numPix Number of pixels in the scan
 * @param geoData Structure containing geolocation data
 * @param geoLut Structure containing geolocation lookup table data
 * @param timeOffsets Vector of time offsets for each pixel
 * @param mceSpinId MCE spin ID for the current scan
 * @param currScan Current scan number
 * @param thetaCorrections [out] Vector to store the calculated theta corrections
 */
void getThetaCorrections(const size_t numPix, const GeoData &geoData, const GeoData::GeoLut &geoLut,
                         const std::vector<double> &timeOffsets, int mceSpinId, size_t currScan,
                         std::vector<double> &thetaCorrections);

/**
 * @brief Check for any missing scan start times and interpolate as necessary.
 * @param scanStartTimes The times that will be checked and modified
 * @param sfl flag array which also might get modified
 */
void interpolateMissingScanTimes(std::vector<double> &scanStartTimes, std::vector<short> &sfl);

/**
 * @brief Get J2000 to ECEF transformation matrix
 * @param year (I) describing 4-digit year
 * @param dayOfYear (I) describing the day of the year
 * @param secondsOfDay (I) describing seconds since the start of the day
 * @param transformationMatrix (O), a 3x3 transformation matrix between J2000 and ECEF reference frames
 * @return 0
 */
int j2000ToEcr(int32_t year, int32_t dayOfYear, double secondsOfDay, double transformationMatrix[3][3]);

/**
 * @brief Get J2000 to MOD (precession) transformation
 * @param iyr (I) describing 4-digit year
 * @param idy (I) describing the day of the year
 * @param sec (I) describing seconds since the start of the day
 * @param ecmat (O), a 3x3 transformation matrix between j2000 and MOD reference frames
 * @return 0
 */
int j2000ToMod(int32_t iyr, int32_t idy, double sec, double j2mod[3][3]);

/**
 * @brief Get the nutation matrix for a given date
 * @param year The year (4-digit)
 * @param dayOfYear The day of the year (1-366)
 * @param nutationMatrix [out] 3x3 nutation matrix
 * @return int Status code (0 for success, non-zero for error)
 */
int getNutation(int32_t year, int32_t dayOfYear, double nutationMatrix[3][3]);

/**
 * @brief Get the UT1-UTC difference for a given date
 * @param year The year (4-digit)
 * @param dayOfYear The day of the year (1-366)
 * @param ut1ToUtc [out] The UT1-UTC difference in seconds
 * @return int Status code (0 for success, non-zero for error)
 */
int getUt1ToUtc(int32_t year, int32_t dayOfYear, double &ut1ToUtc);

/**
 * @brief Convert a 3x3 rotation matrix to a quaternion
 * @param rm Input 3x3 rotation matrix
 * @param q [out] Output quaternion (4-element array)
 * @return int Status code (0 for success, non-zero for error)
 */
int matrixToQuaternion(double rm[3][3], double q[4]);

/**
 * @brief Multiply two quaternions (double precision * float precision)
 * @param q1 First quaternion (double precision)
 * @param q2 Second quaternion (float precision)
 * @param result [out] Result quaternion (double precision)
 * @return 0 always
 */
int multiplyQuaternions(double q1[4], float q2[4], double result[4]);

/**
 * @brief Multiply two quaternions (float precision * float precision)
 * @param q1 First quaternion (float precision)
 * @param q2 Second quaternion (float precision)
 * @param result [out] Result quaternion (float precision)
 * @return 0 always
 */
int multiplyQuaternions(float q1[4], float q2[4], float result[4]);

/**
 * @brief Interpolate orbit position and velocity vectors to scan line times using cubic interpolation.
 * @param numOrbRec Number of orbit records. Used to find vectors bracketing each scan
 * @param sdim Number of scans
 * @param orbitTimes Array of orbit vector times in seconds of day
 * @param positions Array of position vectors for each time in orbitTimes
 * @param velocities Array of velocity vectors for each time in orbitTimes
 * @param scanTimes Array of times in seconds of day for each scan line
 * @param posI Array of interpolated positions
 * @param velI Array of interpolated velocities
 */
void interpolateOrbitVectors(size_t numOrbRec, size_t numScans, const std::vector<double> &orbitTimes,
                             const std::vector<orbArray> &positions, const std::vector<orbArray> &velocities,
                             const std::vector<double> &scanTimes, std::vector<orbArray> &posi,
                             std::vector<orbArray> &veli);

/**
 * @brief Interpolate quaternions to scan line times
 * @param numAttRec Number of attitude records
 * @param numScans Number of scans
 * @param quatTimes Array of times for each quaternion
 * @param quats Quaternions to be interpolated
 * @param scanTimes Array of scan start times
 * @param quatsInterpolated Interpolated quaternions
 */
int interpolateAttitudeForScanTimes(size_t numAttRec, size_t numScans, double *quatTimes, quaternion *quats,
                                    double *scanTimes, quaternion *quatsInterpolated);

/**
 * @brief Interpolate tilt angles to scan line times
 * @param numTiltRecords Number of tilt records
 * @param numScans Number of scans
 * @param tiltTimes Array of times for each tilt measurement
 * @param tiltAngles Array of input tilt angles
 * @param scanTimes Array of scan start times
 * @param interpolatedTilts Array of interpolated tilt angles
 */
int interpolateTilt(size_t numTiltRecords, size_t numScans, double *tiltTimes, float *tiltAngles,
                    double *scanTimes, float *interpolatedTilts);

/**
 * @brief Computes the Sun vector in Earth-Centered Rotating (ECR) coordinates
 * @param numScans Number of scans
 * @param year Year of the observation
 * @param dayOfYear Day of the year (1-366)
 * @param secondsOfDay Array of seconds of the day for each scan
 * @param sunVector Array of unit Sun vectors in ECR coordinates
 * @param sunEarthDistance Array of distances between Sun and Earth in astronomical units
 */
int getEcrSunVector(size_t numScans, int32_t year, int32_t dayOfYear, double *secondsOfDay,
                    orbArray *sunVector, double *sunEarthDistance);

/**
 * @brief Computes the Sun vector in geocentric inertial (equatorial) coordinates.
 *
 * This function uses the model referenced in The Astronomical Almanac for 1984,
 * Section S (Supplement) and documented for the SeaWiFS Project in
 * "Constants and Parameters for SeaWiFS Mission Operations", in TBD.
 * The accuracy of the Sun vector is approximately 0.1 arcminute.
 *
 * @param numScans Number of scans
 * @param year Year of the observation
 * @param dayOfYear Day of the year (1-366)
 * @param secondsOfDay Seconds of the day for each scan
 * @param[out] sunVector Array of unit Sun vectors in geocentric inertial coordinates
 * @param[out] sunEarthDistance Array of distances between Sun and Earth in astronomical units
 * @return int Status code (0 for success, non-zero for error)
 */
int getInertialSunVector(size_t numScans, int32_t year, int32_t dayOfYear, double *secondsOfDay,
                         orbArray *sunVector, double *sunEarthDistance);

/**
 * @brief Converts a quaternion to a rotation matrix
 *
 * This function takes a quaternion and converts it to a 3x3 rotation matrix.
 *
 * @param quaternion Input quaternion (float array of size 4)
 * @param rotationMatrix Output rotation matrix (double 3x3 array)
 * @return 0 always
 */
int quaternionToMatrix(float quaternion[4], double rotationMatrix[3][3]);

/**
 * @brief Calculate the coefficients which represent the Earth scan track in the sensor frame. Uses an
 * equatorial radius of 6378.137 km and a flattening factor of 1/298.257 (WGS 1984).
 *
 * @param pos (I) ECR orbit position vector (km)
 * @param sOMat (I) describing sensor orientation matrix
 * @param locatingContext the locating context
 * @param spCoef (O) describing scan path coefficients
 */
int getEllipsoidScanTrackCoefs(orbArray pos, double sOMat[3][3], LocatingContext locatingContext,
                               double spCoefs[10]);

/** @brief This subroutine performs navigation of a scanning sensor on the
 * surface of an ellipsoid based on an input orbit position vector and
 * spacecraft orientation matrix.  It uses a closed-form algorithm for
 * determining points on the ellipsoidal surface which involves
 * determining the intersection of the scan plan with the ellipsoid.
 * The sensor view vectors in the sensor frame are passed in as a 3xN array.
 * The reference ellipsoid is set according to the scan
 * intersection coefficients in the calling sequence an equatorial
 * radius of 6378.137 km. and a flattening factor of 1/298.257 are
 * used by both the Geodetic Reference System (GRS) 1980 and the
 * World Geodetic System (WGS) 1984.
 * It then computes geometric parameters using the pixel locations on
 * the Earth, the spacecraft position vector and the unit Sun vector in
 * the geocentric rotating reference frame.  The outputs are arrays of
 * geodetic latitude and longitude, solar zenith and azimuth and sensor
 * zenith and azimuth.  The azimuth angles are measured from local
 * North toward East.  Flag values of 999. are returned for any pixels
 * whose scan angle is past the Earth's horizon.
 *
 * @param[in] demFile Digital Elevation Map file path
 * @param[in] position ECR orbit position vector in km at scan mid-time
 * @param[in] velocity ECR orbit velocity vector in km/sec
 * @param[in] sensorOrientationMatrix Sensor Orientation Matrix
 * @param[in] scanPathCoefs Scan path coefficients
 * @param[in] locatingContext Indicates which celestial body is being observed (GEO, HELIO, SELENO, or NONE)
 * @param[in] sunUnitVector Sun unit vector in ECR
 * @param[in] checkForEclipse Boolean flag to check for eclipse conditions
 * @param[in] earthToMoon Vector from Earth's center to Moon's center
 * @param[in] earthSunDistance Distance between Earth and Sun in astronomical units
 * @param[in] sensorViewVectors Sensor view vectors
 * @param[in] numPixels Number of pixels to geolocate
 * @param[in] deltaT Change in time
 * @param[out] qualityFlags Quality flags for each pixel
 * @param[in] currScan Current scan number
 * @param[out] latitudes Geodetic latitude of each pixel
 * @param[out] longitudes Geodetic longitude of each pixel
 * @param[out] solarZeniths Solar zeniths for each pixel
 * @param[out] solarAzimuths Solar azimuths for each pixel
 * @param[out] sensorZeniths Sensor zeniths for each pixel
 * @param[out] sensorAzimuths Sensor azimuths for each pixel
 * @param[out] terrainHeights Terrain height for each pixel
 * @param[out] geoBox Geospatial lat/lon mins/maxes
 * @return int Status code
 */
int locatePixelsOci(const char *demFile, orbArray position, orbArray velocity,
                    double sensorOrientationMatrix[3][3], double scanPathCoefs[10],
                    LocatingContext locatingContext, orbArray sunUnitVector, bool checkForEclipse,
                    orbArray earthToMoon, double earthSunDistance,
                    std::vector<std::array<float, 3>> &sensorViewVectors, size_t numPixels, double *deltaT,
                    std::vector<uint8_t> &qualityFlags, size_t currScan, float *latitudes, float *longitudes,
                    float *solarZeniths, float *solarAzimuths, float *sensorZeniths, float *sensorAzimuths,
                    short *terrainHeights);

/**
 * @brief Calculate the attitude angles given the ECEF orbit vector and attiude matrix. Rotation order is
 * (1, 2, 3). The reference ellipsoid uses an equatorial radius of 6378.137 km and a flattening factor of
 * 1/298.257 (WGS 1984)
 *
 * @param pos (I) Orbit position vector (ECEF)
 * @param vel (I) Orbit vel vector (ECEF)
 * @param smat (I) Sensor attitude matrix (ECEF to sensor)
 * @param rpy (O) Attitude angles (roll, pitch, yaw)
 *
 */
int getAttitudeAngles(orbArray pos, orbArray vel, double smat[3][3], orbArray &rpy);

#endif  // _GEOLOCATE_OCI_H_
