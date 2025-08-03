/**
 * @file geolocate_oci.cpp
 * @brief Geolocation module for OCI instrument data
 *
 * This file contains functions for geolocation of OCI (Ocean Color Instrument)
 * data, including orbit interpolation, attitude interpolation, pixel geolocation,
 * and various coordinate transformations.
 *
 * @authors Joel Gales (SAIC), Jakob C Lindo (SSAI)
 * @date Aug 2024
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/tensor.hpp>

#include <netcdf>

#include "geolocate_oci.h"
#include "nc4utils.h"
#include "timeutils.h"
#include "global_attrs.h"
#include <terrain.h>
#include <libnav.h>
#include <sun2000.h>
#include <genutils.h>
#include "aggregate.hpp"
#include <allocate2d.h>
#include <array>
#include <gring.h>

// These two J2000 constants are defined here because the compiler doesn't like them in a higher scope
constexpr double J_2000 = 2451545.0;           // A reference epoch for astronomical calculations
constexpr double J_2000_MODIFIED = 2451545.5;  // Some calculations prefer starting from midnight

#define ERROR_EXIT(dim)                                                                                    \
    {                                                                                                      \
        if (dim == 0) {                                                                                    \
            printf("--Error--: the dimension %s is zero. Exiting. See line %d in file %s", #dim, __LINE__, \
                   __FILE__);                                                                              \
            exit(EXIT_FAILURE);                                                                            \
        }                                                                                                  \
    }

// Earth ellipsoid parameters
constexpr double EARTH_RADIUS = 6378.137;         // Kilometers, at the equator
constexpr double EARTH_FLATTENING = 1 / 298.257;  // Describes how 'squished' the Earth is
constexpr double EARTH_OMF2 = (1 - EARTH_FLATTENING) * (1 - EARTH_FLATTENING);

constexpr double MOON_RADIUS = 1738.1;  // Kilometers
constexpr double MOON_FLATTENING = 1 / 827.667;
constexpr double MOON_OMF2 = (1 - MOON_FLATTENING) * (1 - MOON_FLATTENING);

constexpr double omegaE = 7.29211585494e-5;

// TODO: Define J2000 = 2451545 somwehere both calibrations and geolocation can see it

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     SAIC           09/13/21 0.01  Original development
//  Joel Gales     SAIC           06/14/22 0.03  Add support for tilt change
//  Gwyn Fireman   SAIC           07/10/23 0.07  Read global metadata from json file
//  Joel Gales     SAIC           08/08/23 0.08  Check scan times
//  Joel Gales     SAIC           08/18/23 0.80  Add earth_sun distance
//                                               correction attribute, support
//                                               for reflectance output
//  Joel Gales     SAIC           08/29/23 0.81  Change fill values
//  Joel Gales     SAIC           10/03/23 0.85  Change tilt output name to
//                                               tilt_angle
//                                               Change to callable routine
//                                               rather than standalone program

using namespace std;
namespace ublas = boost::numeric::ublas;

void setGslVector(gsl_vector *vector, orbArray orbVector) {
    gsl_vector_set(vector, 0, static_cast<double>(orbVector[0]));
    gsl_vector_set(vector, 1, static_cast<double>(orbVector[1]));
    gsl_vector_set(vector, 2, static_cast<double>(orbVector[2]));
}

void setGslVector(gsl_vector *vector, orbArray orbVector, double scalar) {
    gsl_vector_set(vector, 0, static_cast<double>(orbVector[0]) * scalar);
    gsl_vector_set(vector, 1, static_cast<double>(orbVector[1]) * scalar);
    gsl_vector_set(vector, 2, static_cast<double>(orbVector[2]) * scalar);
}

gsl_vector *getGslVector(orbArray vector) {
    gsl_vector *result = gsl_vector_alloc(3);
    setGslVector(result, vector);
    return result;
}

gsl_vector *getGslVector(orbArray vector, double scalar) {
    gsl_vector *result = gsl_vector_alloc(3);
    setGslVector(result, vector, scalar);
    return result;
}

/**
 * @brief Scale all elements of an array of double of size 3 in place.
 * @note orbArrays are typedef'd pointers, so there's no need to make the orbArray parameter a reference or a
 * pointer
 * @note Such things should be implemented by orbArray, but that will result in several changes, and we plan
 * to move to an ppropriate data structure soon
 *
 * @param orbArray The vector to be scaled
 * @param scalar The quantity by which to scale the whole orbArray
 */
void scaleElements(double vec[3], const double scalar) {
    double temp[3];
    temp[0] = vec[0] * scalar;
    temp[1] = vec[1] * scalar;
    temp[2] = vec[2] * scalar;

    vec[0] = temp[0];
    vec[1] = temp[1];
    vec[2] = temp[2];
}

/**
 * @brief Computes coordinate transformation matrix corresponding to Euler sequence; assumes order of
 * rotations is 3, 1, 3
 * @param angles The Euler angles
 * @return The rotation matrix
 */
gsl_matrix *createRotationMatrix(const std::array<double, 3> &angles) {
    double cosX = cos(angles[0] / OEL_RADEG);
    double sinX = sin(angles[0] / OEL_RADEG);
    double cosY = cos(angles[1] / OEL_RADEG);
    double sinY = sin(angles[1] / OEL_RADEG);
    double cosZ = cos(angles[2] / OEL_RADEG);
    double sinZ = sin(angles[2] / OEL_RADEG);

    // Create individual rotation matrices
    gsl_matrix *rotX = gsl_matrix_alloc(3, 3);
    gsl_matrix *rotY = gsl_matrix_alloc(3, 3);
    gsl_matrix *rotZ = gsl_matrix_alloc(3, 3);

    gsl_matrix_set(rotX, 0, 0, cosX);
    gsl_matrix_set(rotX, 0, 1, sinX);
    gsl_matrix_set(rotX, 1, 0, -sinX);
    gsl_matrix_set(rotX, 1, 1, cosX);
    gsl_matrix_set(rotX, 2, 2, 1.0);

    gsl_matrix_set(rotY, 0, 0, 1.0);
    gsl_matrix_set(rotY, 1, 1, cosY);
    gsl_matrix_set(rotY, 1, 2, sinY);
    gsl_matrix_set(rotY, 2, 1, -sinY);
    gsl_matrix_set(rotY, 2, 2, cosY);

    gsl_matrix_set(rotZ, 0, 0, cosZ);
    gsl_matrix_set(rotZ, 0, 1, sinZ);
    gsl_matrix_set(rotZ, 1, 0, -sinZ);
    gsl_matrix_set(rotZ, 1, 1, cosZ);
    gsl_matrix_set(rotZ, 2, 2, 1.0);

    // Compute total rotation as rotZ * rotX * rotZ
    gsl_matrix *xzRotation = gsl_matrix_alloc(3, 3);
    gsl_matrix *rotationMatrix = gsl_matrix_alloc(3, 3);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rotX, rotZ, 0.0, xzRotation);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rotZ, xzRotation, 0.0, rotationMatrix);

    // Clean up temporary matrices
    gsl_matrix_free(rotX);
    gsl_matrix_free(rotY);
    gsl_matrix_free(rotZ);
    gsl_matrix_free(xzRotation);

    return rotationMatrix;
}

/**
 * @brief Get ECR to selenographic transformation matrix
 * @param year Year, four digits
 * @param dayOfYear Day of year
 * @param secondsOfDay Seconds of day
 * @param selenocentricMatrix ECR to selenographic transformation matrix
 */
void ecrToSeleno(int year, int dayOfYear, double secondsOfDay, gsl_matrix *selenocentricMatrix) {
    using namespace boost;

    const float MOON_INCLINATION = -1.535;  // Angle off the Earth's ecliptic. AKA ximeq
    double utcTime;

    (void)getUt1ToUtc(year, dayOfYear, utcTime);

    double greenwichHourAngle;
    double dayOfYearInSeconds = dayOfYear + (secondsOfDay + utcTime) / SECONDS_IN_DAY;
    gha2000(year, dayOfYearInSeconds, &greenwichHourAngle);

    // The IDL version of gha2000 includes calls to ephparms() and nutate() and outputs their results
    double meanLongitudeMoon;    // AKA xlm
    double meanLongitudeSun;     // AKA xls
    double meanAnomalySun;       // AKA gs
    double ascNodeMoon;          // Ascending node of the Moon's mean orbit. AKA omega
    double eclipticObliquity;    // True obliquity of the Earth's ecliptic. AKA eps
    double nutationInLongitude;  // AKA dpsi
    double fractionOfDayOfYear = dayOfYearInSeconds - dayOfYear;
    double julianDay = jday(year, 1, (int16_t)dayOfYear);
    double daysSinceJ2000 = julianDay - J_2000_MODIFIED + fractionOfDayOfYear;
    ephparms(daysSinceJ2000, &meanLongitudeSun, &meanAnomalySun, &meanLongitudeMoon, &ascNodeMoon);
    nutate(daysSinceJ2000, meanLongitudeSun, meanAnomalySun, meanLongitudeMoon, ascNodeMoon,
           &nutationInLongitude, &eclipticObliquity);

    double cosGhaDegrees = cos(greenwichHourAngle / OEL_RADEG);
    double sinGhaDegrees = sin(greenwichHourAngle / OEL_RADEG);
    gsl_matrix *ghaRotation = gsl_matrix_alloc(3, 3);
    gsl_matrix_set(ghaRotation, 0, 0, cosGhaDegrees);
    gsl_matrix_set(ghaRotation, 1, 1, cosGhaDegrees);
    gsl_matrix_set(ghaRotation, 2, 2, 1.0);
    gsl_matrix_set(ghaRotation, 0, 1, sinGhaDegrees);
    gsl_matrix_set(ghaRotation, 1, 0, -sinGhaDegrees);

    double cosTodDegrees = cos(eclipticObliquity / OEL_RADEG);
    double sinTodDegrees = sin(eclipticObliquity / OEL_RADEG);
    gsl_matrix *eclipticRotation = gsl_matrix_alloc(3, 3);
    gsl_matrix_set(eclipticRotation, 0, 0, 1.0);
    gsl_matrix_set(eclipticRotation, 1, 1, cosTodDegrees);
    gsl_matrix_set(eclipticRotation, 2, 2, cosTodDegrees);
    gsl_matrix_set(eclipticRotation, 1, 2, sinTodDegrees);
    gsl_matrix_set(eclipticRotation, 2, 1, -sinTodDegrees);

    double fm = meanLongitudeMoon - ascNodeMoon - 180;
    // ecl2sel
    gsl_matrix *moonRotation = createRotationMatrix({ascNodeMoon, MOON_INCLINATION, fm});

    gsl_matrix *eclipticToMoon = gsl_matrix_alloc(3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, moonRotation, eclipticRotation, 0.0, eclipticToMoon);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, eclipticToMoon, ghaRotation, 0.0, selenocentricMatrix);

    gsl_matrix_free(ghaRotation);
    gsl_matrix_free(eclipticRotation);
    gsl_matrix_free(moonRotation);
    gsl_matrix_free(eclipticToMoon);
}

/**
 * @brief Calculates the Moon vector in geocentric *inertial* (equatorial) coordinates
 * @param[in] year The 4-digit year
 * @param[in] day Between 1 and 366
 * @param[in] secondsOfDay Between 0 and 86400
 * @param[out] moonVecInertial  Moon vector in geocentric rotating coordinates (3 elements). Accurate to
 * approximately 1 arcminute
 * @param[out] lunarDistance Earth-Moon distance in kilometers. Expect ~360,000 to ~400,000
 *
 * @ref "Low-precision Formula for Planetary Positions", T.C. Van Flandern and K.F. Pulkkinen,
 *       Ap. J. Supplement Series, 41:391-411, 1979.
 */
void moon2000(int year, int dayOfYear, double secondsOfDay, orbArray &moonVecInertial,
              double &lunarDistance) {
    int month = 1;                                 // Passing a month of 0 to jday() will cause issues
    int julianDay = jday(year, month, dayOfYear);  // jday
    constexpr double SECONDS_IN_HALF_DAY = SECONDS_IN_DAY / 2;  // Removes a magic number
    // Days since Jan 2 1980 with a fractional component in terms of seconds. t
    double fractionalJulianDay = julianDay - J_2000 + (secondsOfDay - SECONDS_IN_HALF_DAY) / SECONDS_IN_DAY;

    double meanLongitudeSun;   // xls
    double meanAnomalySun;     // In degrees. gs
    double meanLongitudeMoon;  // xlm
    double ascendingNodeMoon;  // In degrees. omega
    ephparms(fractionalJulianDay, &meanLongitudeSun, &meanAnomalySun, &meanLongitudeMoon, &ascendingNodeMoon);

    double nutation;           // Actual nutation value in degrees. dpsi
    double eclipticObliquity;  // In terms of degrees. Includes nutation. Axial tilt of the Earth. eps
    nutate(fractionalJulianDay, meanLongitudeSun, meanAnomalySun, meanLongitudeMoon, ascendingNodeMoon,
           &nutation, &eclipticObliquity);

    double meanAnomalyMoon = fmod(134.96292 + 13.06499295 * fractionalJulianDay, 360);  // gm
    double argumentOfLatitudeMoon = fmod(meanLongitudeMoon - ascendingNodeMoon, 360);   // fm
    double meanElongationMoon = fmod(meanLongitudeMoon - meanLongitudeSun, 360);        // dm

    // rm = ...
    lunarDistance = 60.36298 - 3.27746 * cos(meanAnomalyMoon / OEL_RADEG) -
                    .57994 * cos((meanAnomalyMoon - 2.0 * meanElongationMoon) / OEL_RADEG) -
                    0.46357 * cos(2.0 * meanElongationMoon / OEL_RADEG) -
                    0.08904 * cos(2.0 * meanAnomalyMoon / OEL_RADEG) +
                    0.03865 * cos((2.0 * meanAnomalyMoon - 2.0 * meanElongationMoon) / OEL_RADEG) -
                    0.03237 * cos((2.0 * meanElongationMoon - meanAnomalySun) / OEL_RADEG) -
                    0.02688 * cos((meanAnomalyMoon + 2.0 * meanElongationMoon) / OEL_RADEG) -
                    0.02358 * cos((meanAnomalyMoon - 2.0 * meanElongationMoon + meanAnomalySun) / OEL_RADEG) -
                    0.02030 * cos((meanAnomalyMoon - meanAnomalySun) / OEL_RADEG) +
                    0.01719 * cos(meanElongationMoon / OEL_RADEG) +
                    0.01671 * cos((meanAnomalyMoon + meanAnomalySun) / OEL_RADEG) +
                    0.01247 * cos((meanAnomalyMoon - 2.0 * argumentOfLatitudeMoon) / OEL_RADEG);

    // dlm = ...
    double lunarLongitudeCorrection =
        22640.0 * sin(meanAnomalyMoon / OEL_RADEG) -
        4586.0 * sin((meanAnomalyMoon - 2.0 * meanElongationMoon) / OEL_RADEG) +
        2370.0 * sin(2.0 * meanElongationMoon / OEL_RADEG) + 769.0 * sin(2.0 * meanAnomalyMoon / OEL_RADEG) -
        668.0 * sin(meanAnomalySun / OEL_RADEG) - 412.0 * sin(2.0 * argumentOfLatitudeMoon / OEL_RADEG) -
        212.0 * sin((2.0 * meanAnomalyMoon - 2.0 * meanElongationMoon) / OEL_RADEG) -
        206.0 * sin((meanAnomalyMoon - 2.0 * meanElongationMoon + meanAnomalySun) / OEL_RADEG) +
        192.0 * sin((meanAnomalyMoon + 2.0 * meanElongationMoon) / OEL_RADEG) +
        165.0 * sin((2.0 * meanElongationMoon - meanAnomalySun) / OEL_RADEG) +
        148.0 * sin((meanAnomalyMoon - meanAnomalySun) / OEL_RADEG) -
        125.0 * sin(meanElongationMoon / OEL_RADEG) -
        110.0 * sin((meanAnomalyMoon + meanAnomalySun) / OEL_RADEG) -
        55.0 * sin((2.0 * argumentOfLatitudeMoon - 2.0 * meanElongationMoon) / OEL_RADEG) -
        45.0 * sin((meanAnomalyMoon + 2.0 * argumentOfLatitudeMoon) / OEL_RADEG) +
        40.0 * sin((meanAnomalyMoon - 2.0 * argumentOfLatitudeMoon) / OEL_RADEG) -
        38.0 * sin((meanAnomalyMoon - 4.0 * meanElongationMoon) / OEL_RADEG) +
        36.0 * sin(3.0 * meanAnomalyMoon / OEL_RADEG) -
        31.0 * sin((2.0 * meanAnomalyMoon - 4.0 * meanElongationMoon) / OEL_RADEG) +
        28.0 * sin((meanAnomalyMoon - 2.0 * meanElongationMoon - meanAnomalySun) / OEL_RADEG) -
        24.0 * sin((2.0 * meanElongationMoon + meanAnomalySun) / OEL_RADEG) +
        19.0 * sin((meanAnomalyMoon - meanElongationMoon) / OEL_RADEG) +
        18.0 * sin((meanElongationMoon + meanAnomalySun) / OEL_RADEG) +
        15.0 * sin((meanAnomalyMoon + 2.0 * meanElongationMoon - meanAnomalySun) / OEL_RADEG) +
        14.0 * sin((2.0 * meanAnomalyMoon + 2.0 * meanElongationMoon) / OEL_RADEG) +
        14.0 * sin((4.0 * meanElongationMoon) / OEL_RADEG) -
        13.0 * sin((3.0 * meanAnomalyMoon - 2.0 * meanElongationMoon) / OEL_RADEG);

    // xlma = xlm + dlm / 3600.d0 + dpsi
    double lunarLongitude = meanLongitudeMoon + lunarLongitudeCorrection / 3600.0 + nutation;

    // beta = ...
    double lunarEclipticLatitude =
        18461.0 * sin(argumentOfLatitudeMoon / OEL_RADEG) +
        1010.0 * sin((meanAnomalyMoon + argumentOfLatitudeMoon) / OEL_RADEG) +
        1000.0 * sin((meanAnomalyMoon - argumentOfLatitudeMoon) / OEL_RADEG) -
        624.0 * sin((argumentOfLatitudeMoon - 2.0 * meanElongationMoon) / OEL_RADEG) -
        199.0 * sin((meanAnomalyMoon - argumentOfLatitudeMoon - 2.0 * meanElongationMoon) / OEL_RADEG) -
        167.0 * sin((meanAnomalyMoon + argumentOfLatitudeMoon - 2.0 * meanElongationMoon) / OEL_RADEG) +
        117.0 * sin((argumentOfLatitudeMoon + 2.0 * meanElongationMoon) / OEL_RADEG) +
        62.0 * sin((2.0 * meanAnomalyMoon + argumentOfLatitudeMoon) / OEL_RADEG) +
        33.0 * sin((meanAnomalyMoon - argumentOfLatitudeMoon + 2.0 * meanElongationMoon) / OEL_RADEG) +
        32.0 * sin((2.0 * meanAnomalyMoon - argumentOfLatitudeMoon) / OEL_RADEG) -
        30.0 * sin((argumentOfLatitudeMoon - 2.0 * meanElongationMoon + meanAnomalySun) / OEL_RADEG) -
        16.0 * sin((2.0 * meanAnomalyMoon + argumentOfLatitudeMoon - 2.0 * meanElongationMoon) / OEL_RADEG) +
        15.0 * sin((meanAnomalyMoon + argumentOfLatitudeMoon + 2.0 * meanElongationMoon) / OEL_RADEG) +
        12.0 * sin((argumentOfLatitudeMoon - 2.0 * meanElongationMoon - meanAnomalySun) / OEL_RADEG);
    lunarEclipticLatitude /= 3600.0;

    // Compute unit Moon vector in ecliptic plane
    double lunarEclipticLatitudeDegrees = lunarEclipticLatitude / OEL_RADEG;
    double cosLunarEclipticLatitude = cos(lunarEclipticLatitudeDegrees);
    // xme
    double moonVecEcliptic[3] = {cos(lunarLongitude / OEL_RADEG) * cosLunarEclipticLatitude,
                                 sin(lunarLongitude / OEL_RADEG) * cosLunarEclipticLatitude,
                                 sin(lunarEclipticLatitudeDegrees)};

    // Rotate vector to equatorial plane
    double eclipObliqDegrees = eclipticObliquity / OEL_RADEG;
    // xmn = xme ...
    double cosEclipObliq = cos(eclipObliqDegrees);
    double sinEclipObliq = sin(eclipObliqDegrees);
    moonVecInertial[0] = moonVecEcliptic[0];
    moonVecInertial[1] = moonVecEcliptic[1] * cosEclipObliq - moonVecEcliptic[2] * sinEclipObliq;
    moonVecInertial[2] = moonVecEcliptic[1] * sinEclipObliq + moonVecEcliptic[2] * cosEclipObliq;

    lunarDistance *= EARTH_RADIUS;
}

/**
 * @brief Computes the unit Moon vector in geocentric *rotating* coordinates, using subprograms to compute
 * inertial Moon vector and Greenwhich hour angle
 * @note Referred to in reference IDL as l_moon
 * @param[in] year The 4-digit year
 * @param[in] dayOfYear Between 1 and 366
 * @param[in] secondsOfDay Between 0 and 86400
 * @param[out] moonVecRotating Unit Moon vector in geocentric rotating coordinates (3 elements). Magnitude
 * should be close to the Earth-Moon distance at the given time
 * @param[out] lunarDistance Earth-Moon distance in kilometers. Expect ~360,000 to ~400,000
 */
void getMoonVector(const int year, const int dayOfYear, const double secondsOfDay, orbArray &moonVecRotating,
                   double &lunarDistance) {
    orbArray moonVecInertial;  // AKA xm
    moon2000(year, dayOfYear, secondsOfDay, moonVecInertial, lunarDistance);

    // Get Greenwich mean sidereal angle
    double fractionalDay = dayOfYear + secondsOfDay / SECONDS_IN_DAY;  // Turns day-of-year XXX into XXX.YYY
    double greenwichHourAngleDegrees;
    gha2000(year, fractionalDay, &greenwichHourAngleDegrees);
    double greenwichHourAngleRadians = greenwichHourAngleDegrees / OEL_RADEG;

    // Transform Moon vector into geocentric rotating frame
    double cosGhar = cos(greenwichHourAngleRadians);
    double sinGhar = sin(greenwichHourAngleRadians);
    moonVecRotating[0] = moonVecInertial[0] * cosGhar + moonVecInertial[1] * sinGhar;
    moonVecRotating[1] = moonVecInertial[1] * cosGhar - moonVecInertial[0] * sinGhar;
    moonVecRotating[2] = moonVecInertial[2];
}

vector<string> tokenize(istringstream stringStream) {
    vector<string> tokens;
    string token;
    while (stringStream >> token) {
        tokens.push_back(token);
    }

    return tokens;
}

/**
 * @brief Read a predicted PACE ephemeris file and convert vectors to ECR.
 * @param[in] ephFile Name of the ephemeris file. Expected to be non-null and not empty.
 * @param[in] recordTimes Out-parameter indicating time of each record.
 * @param[out] positions In/out-parameter indicating position of each record.
 * @param[out] velocities In/out-parameter indicating velocity of each record.
 * @param[in] l1aEpoch Unix time of the start of the L1A file's first day
 */
void readDefinitiveEphemeris(string ephFile, vector<double> &recordTimes, vector<vector<double>> &positions,
                             vector<vector<double>> &velocities, double l1aEpoch) {
    ifstream inStream(ephFile, ios::in);

    if (!inStream.is_open()) {
        cerr << "-E- Unable to open definitive ephemeris file\n";
        exit(EXIT_FAILURE);
    }

    vector<double> tempPositions(3);
    vector<double> tempVelocities(3);
    string ephStr;
    vector<string> tokens;  // Datetime, 3d position, 3d velocity.
    double prevTime = 0;

    while (inStream.peek() != EOF) {
        getline(inStream, ephStr, '\n');

        istringstream lineStream(ephStr);
        string token;  // Each individual token in this line

        tokens.clear();
        while (lineStream >> token) {
            tokens.push_back(token);
        }

        if (tokens.size() != 7 || tokens[0].find("T") == string::npos) {
            continue;
        }

        // skip duplicate records
        double theTime = isodate2unix(tokens[0].c_str());
        if (prevTime == theTime)
            continue;
        prevTime = theTime;

        tempPositions[0] = stod(tokens[1]);
        tempPositions[1] = stod(tokens[2]);
        tempPositions[2] = stod(tokens[3]);

        tempVelocities[0] = stod(tokens[4]);
        tempVelocities[1] = stod(tokens[5]);
        tempVelocities[2] = stod(tokens[6]);

        positions.push_back(tempPositions);
        velocities.push_back(tempVelocities);
        recordTimes.push_back(theTime - l1aEpoch);

        ephStr.clear();
        lineStream.clear();
    }
    inStream.close();

    if (recordTimes.size() == 0) {
        cerr << "-E- No records found in definitive ephemeris file\n";
        exit(EXIT_FAILURE);
    }
}

void interpolateMissingScanTimes(vector<double> &scanStartTimes, vector<short> &sfl) {
    vector<uint32_t> invalidTimeIndices, validTimeIndices;

    for (size_t i = 0; i < scanStartTimes.size(); i++) {
        if (scanStartTimes[i] == -999 || scanStartTimes[i] == -32767)
            invalidTimeIndices.push_back(i);
        else
            validTimeIndices.push_back(i);
    }

    size_t numInvalidTimes = invalidTimeIndices.size();
    if (numInvalidTimes == 0)
        return;  // No missing scan times

    if (validTimeIndices.size() < 2) {
        printf("-E- %s:%d - need at least 2 good scan times\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    // Interpolate valid times to fill missing values
    for (size_t i = 0; i < invalidTimeIndices.size(); i++) {
        // Check for missing time before valid time

        size_t currBadTimeIndex = invalidTimeIndices[i];

        if (currBadTimeIndex < validTimeIndices[0]) {  // if there are no good previous times

            double firstValidTimeIndex = validTimeIndices[0];
            double secondValidTimeIndex = validTimeIndices[1];

            scanStartTimes[currBadTimeIndex] =
                scanStartTimes[firstValidTimeIndex] -
                (scanStartTimes[secondValidTimeIndex] - scanStartTimes[firstValidTimeIndex]) /
                    (secondValidTimeIndex - firstValidTimeIndex) * (firstValidTimeIndex - currBadTimeIndex);

        } else if (currBadTimeIndex > validTimeIndices.back()) {  // if there are no good next times

            size_t i = validTimeIndices.size() - 1;
            double firstValidTimeIndex = validTimeIndices[i - 1];
            double secondValidTimeIndex = validTimeIndices[i];

            scanStartTimes[currBadTimeIndex] =
                scanStartTimes[secondValidTimeIndex] +
                (scanStartTimes[secondValidTimeIndex] - scanStartTimes[firstValidTimeIndex]) /
                    (secondValidTimeIndex - firstValidTimeIndex) * (currBadTimeIndex - secondValidTimeIndex);

        } else {  // Interpolation is easy

            size_t beforeValidTimeIndex = 0;
            for (int j = validTimeIndices.size() - 1; j >= 0; j--) {
                if (validTimeIndices[j] < currBadTimeIndex) {
                    beforeValidTimeIndex = j;
                    break;
                }
            }

            size_t beforeTimeIndex = validTimeIndices[beforeValidTimeIndex];
            size_t afterTimeIndex = validTimeIndices[beforeValidTimeIndex + 1];
            scanStartTimes[currBadTimeIndex] =
                scanStartTimes[beforeTimeIndex] +
                (scanStartTimes[afterTimeIndex] - scanStartTimes[beforeTimeIndex]) /
                    (afterTimeIndex - beforeTimeIndex) * (currBadTimeIndex - beforeTimeIndex);
        }
        sfl[currBadTimeIndex] |= 2;  // set missing time flag for all nf(ailed)
    }
}

void getThetaCorrections(const size_t numPix, const GeoData &geoData, const GeoData::GeoLut &geoLut,
                         const vector<double> &timeOffsets, int mceSpinId, size_t currScan,
                         vector<double> &thetaCorrections) {
    size_t pixelIndex = 0;
    size_t encoderValueIndex;

    while (pixelIndex < numPix) {
        bool encoderInterpErrorFound = false;
        double constexpr SAMPLE_DELTA_T = 0.001;  // Time between encoder samples
        double sampleTime;                        // Encoder sample time for this pixel
        size_t numInterpPix = 0;  // Incremented only when pix are interpolated against encoder data

        for (size_t i = 0; i < geoData.mceTelem.encoderSampleCounts[mceSpinId]; i++) {
            double encoderSampleTime = i * SAMPLE_DELTA_T;
            double pixelTime = timeOffsets[pixelIndex];

            if (encoderSampleTime <= pixelTime && pixelTime < (encoderSampleTime + SAMPLE_DELTA_T)) {
                encoderValueIndex = i;
                sampleTime = encoderSampleTime;
                encoderInterpErrorFound = true;
                break;
            }
        }

        if (!encoderInterpErrorFound) {
            cout << "Encoder interpolation error" << endl;
            exit(-1);
        }

        for (size_t i = 0; i < numPix; i++) {
            if (sampleTime + SAMPLE_DELTA_T <= timeOffsets[i] || timeOffsets[i] < sampleTime)
                continue;

            double timeFraction = (timeOffsets[i] - sampleTime) / SAMPLE_DELTA_T;
            double rtaInterp =
                (1 - timeFraction) * geoData.mceTelem.rtaEncData[mceSpinId][encoderValueIndex] +
                timeFraction * geoData.mceTelem.rtaEncData[mceSpinId][encoderValueIndex + 1];
            double hamInterp =
                (1 - timeFraction) * geoData.mceTelem.hamEncData[mceSpinId][encoderValueIndex] +
                timeFraction * geoData.mceTelem.hamEncData[mceSpinId][encoderValueIndex + 1] +
                geoLut.hamCrossTrackAngles[geoData.hamSides[currScan]];

            thetaCorrections[i] = rtaInterp - hamInterp * 0.236;
            numInterpPix++;
        }
        pixelIndex += numInterpPix;
    }
}

int getEarthViewVectors(GeoData &geoData, const GeoData::GeoLut &geoLut, const int device,
                        const double earthViewTimeOffset, const vector<double> &delt, size_t currScan,
                        vector<array<float, 3>> &vectors) {
    size_t numPixels;
    vector<double> &scanAngles = (device == CCD) ? geoData.ccdScanAngles : geoData.swirScanAngles;
    if (device == CCD)
        numPixels = geoData.numCcdPix;
    else if (device == SWIR)
        numPixels = geoData.numSwirPix;
    else
        throw invalid_argument("-E- Unknown device type passed to getEarthViewVectors");

    size_t pixelOffset = currScan * numPixels;

    const int32_t MAX_ENCODER_COUNT = 0x20000;  // 2^17
    double constexpr RADIANS_TO_ARCSECONDS = OEL_RADEG * 3600;
    double constexpr TAU = 2 * OEL_PI;
    vector<double> timeOffsets(numPixels);  // Use in computing ideal scan angles for science pixels

    float pprAngle = TAU * (geoData.mceTelem.pprOffset - geoLut.rtaNadir[geoData.mceTelem.mceBoardId % 2]) /
                     MAX_ENCODER_COUNT;
    if (pprAngle > OEL_PI)
        pprAngle -= TAU;

    for (size_t i = 0; i < numPixels; i++) {
        timeOffsets[i] = delt[i] + earthViewTimeOffset;
        scanAngles[pixelOffset + i] = pprAngle + 2 * OEL_PI * geoData.mceTelem.comRotRate * timeOffsets[i];
    }

    int mceSpinId = -1;  // The MCE spin ID for this spin number
    for (size_t i = 0; i < geoData.mceTelem.numMceScans; i++) {
        if (geoData.mceTelem.mceSpinIds[i] == geoData.spinIds[currScan]) {
            mceSpinId = (int)i;
            break;
        }
    }

    if (mceSpinId == -1) {
        cout << "No MCE encoder data for spin: " << geoData.spinIds[currScan] << endl;
        return EXIT_FAILURE;
    }

    vector<double> thetaCorrections(numPixels);

    getThetaCorrections(numPixels, geoData, geoLut, timeOffsets, mceSpinId, currScan, thetaCorrections);

    // Calculate planarity deviations and view vectors
    vector<double> alongScanPlanDev(numPixels);
    vector<double> alongTrackPlanDev(numPixels);
    for (size_t i = 0; i < numPixels; i++) {
        scanAngles[pixelOffset + i] -= thetaCorrections[i] / RADIANS_TO_ARCSECONDS;

        if (vectors.size() == 0)  // Only want scan angles
            continue;

        alongScanPlanDev[i] = geoLut.alongTrackPlanarity[0];
        alongTrackPlanDev[i] = geoLut.acrossTrackPlanarity[0];
        for (size_t k = 1; k < 5; k++) {
            alongScanPlanDev[i] += geoLut.alongTrackPlanarity[k] * pow(scanAngles[pixelOffset + i], k);
            alongTrackPlanDev[i] += geoLut.acrossTrackPlanarity[k] * pow(scanAngles[pixelOffset + i], k);
        }

        vectors[i][0] = -sin(alongTrackPlanDev[i] / RADIANS_TO_ARCSECONDS);
        vectors[i][1] = sin(scanAngles[pixelOffset + i] - alongScanPlanDev[i] / RADIANS_TO_ARCSECONDS);
        vectors[i][2] = cos(scanAngles[pixelOffset + i] - alongScanPlanDev[i] / RADIANS_TO_ARCSECONDS);
    }
    return EXIT_SUCCESS;
}

int getEarthView(double comRotRate, int16_t *dataTypes, int16_t *spatialZoneLines,
                 int16_t *spatialAggregation, uint16_t &numHyperSciPix, uint16_t &numSwirPixels,
                 double &earthViewTimeOffset, float *scienceLines, float *swirLines,
                 double *sciencePixeloffsets, double *swirPixelOffsets, bool isDark) {
    const uint8_t ON = 1;
    const uint8_t OFF = 0;

    // Find end of no-data zone
    int16_t indexZone = 0, firstDataLine = 0;
    int returnStatus = -1;

    while (dataTypes[indexZone] == NO_DATA) {
        firstDataLine += spatialZoneLines[indexZone];
        indexZone++;
    }

    if (indexZone == 10)
        return 0;

    // Find number of pixels in Earth views
    numHyperSciPix = 0;
    numSwirPixels = 0;
    int16_t currDataLine = firstDataLine;

    // Is this data type collected on orbit? See enum DataType for definitions
    uint8_t collectOnOrbit[] = {OFF, ON, OFF, ON, ON, ON, ON, ON, ON, ON, OFF, OFF, OFF, OFF, OFF};

    if (isDark)
        collectOnOrbit[DARK_CALIBRATION] = ON;

    for (size_t i = indexZone; i < 9; i++) {
        // Check for not dark or no-data
        if (collectOnOrbit[dataTypes[i]] != ON) {
            continue;
        }

        uint16_t numAggregatedPix = spatialZoneLines[i] / spatialAggregation[i];
        for (size_t j = 0; j < numAggregatedPix; j++) {
            scienceLines[numHyperSciPix + j] =
                currDataLine + j * spatialAggregation[i] + 0.5 * spatialAggregation[i] - 64;
        }
        numHyperSciPix += numAggregatedPix;
        uint16_t ns = spatialZoneLines[i] / 8;
        for (size_t j = 0; j < ns; j++) {
            swirLines[numSwirPixels + j] = currDataLine + j * 8 + 4 - 64;
        }
        numSwirPixels += ns;
        returnStatus = 0;

        currDataLine += spatialZoneLines[i];
    }

    // Calculate times
    for (size_t i = 0; i < (size_t)numHyperSciPix; i++)
        sciencePixeloffsets[i] = comRotRate * scienceLines[i];
    earthViewTimeOffset = 0.5 * (sciencePixeloffsets[0] + sciencePixeloffsets[numHyperSciPix - 1]);

    for (size_t i = 0; i < (size_t)numHyperSciPix; i++)
        sciencePixeloffsets[i] -= earthViewTimeOffset;

    for (size_t i = 0; i < (size_t)numSwirPixels; i++)
        swirPixelOffsets[i] = comRotRate * swirLines[i] - earthViewTimeOffset;

    return returnStatus;
}

MceTlm readMceTelemetry(netCDF::NcFile *l1aFile, GeoData::GeoLut &geoLut, netCDF::NcGroup egrData) {
    MceTlm output;
    output.numMceScans = l1aFile->getDim("number_of_mce_scans").getSize();
    size_t numMceScans = output.numMceScans;
    uint32_t numEncoderChannels = l1aFile->getDim("encoder_samples").getSize();
    uint32_t mceBlockSize = l1aFile->getDim("MCE_block").getSize();
    uint32_t numDdcTelemetryPackets = l1aFile->getDim("DDC_tlm").getSize();
    uint32_t numTelemetryPackets = l1aFile->getDim("tlm_packets").getSize();

    boost::multi_array<uint8_t, 2> mceTelemetry(boost::extents[numMceScans][mceBlockSize]);
    boost::multi_array<int16_t, 2> encoderData(boost::extents[numMceScans][4 * numEncoderChannels]);
    output.mceSpinIds = vector<int32_t>(output.numMceScans);
    boost::multi_array<uint8_t, 2> ddcTelemetry(boost::extents[numTelemetryPackets][numDdcTelemetryPackets]);

    egrData.getVar("MCE_telemetry").getVar(&mceTelemetry[0][0]);
    egrData.getVar("MCE_encoder_data").getVar(&encoderData[0][0]);
    egrData.getVar("MCE_spin_ID").getVar(output.mceSpinIds.data());
    egrData.getVar("DDC_telemetry").getVar(&ddcTelemetry[0][0]);

    // Check for missing MCE telemetry
    bool noMCE = true;
    for (size_t i = 0; i < numMceScans; i++) {
        if (output.mceSpinIds[i] >= 0) {
            noMCE = false;
            break;
        }
    }

    if (noMCE) {
        cout << "No MCE telemetry in L1A file" << endl;
        exit(EXIT_FAILURE);
    }

    // Get MCE board ID
    output.mceBoardId = (int16_t)mceTelemetry[0][322] / 16;

    // Get ref_pulse_divider and compute commanded rotation rate
    uint32_t buf;
    uint32_t refPulseDiv[2];
    swapc_bytes2((char *)&mceTelemetry[0][0], (char *)&buf, 4, 1);
    refPulseDiv[0] = buf % 0x1000000;  // 2^24
    swapc_bytes2((char *)&mceTelemetry[0][4], (char *)&buf, 4, 1);
    refPulseDiv[1] = buf % 0x1000000;  // 2^24

    int32_t ref_pulse_sel = mceTelemetry[0][9] / 128;

    ref_pulse_sel == 0
        ? output.comRotRate =
              geoLut.masterClock / 2 / MAX_ENC_COUNT / (refPulseDiv[ref_pulse_sel] / 256.0 + 1)
        : output.comRotRate = geoLut.mceClock / 2 / MAX_ENC_COUNT / (refPulseDiv[ref_pulse_sel] / 256.0 + 1);

    // Check for static or reverse scan
    int32_t averageStepSpeed = mceTelemetry[0][428] * 256 + mceTelemetry[0][429];
    if (abs(averageStepSpeed) < 1000) {
        // Static mode
        output.comRotRate = 0.0;
        cout << "OCI static mode" << endl;
    } else if (averageStepSpeed < 0) {
        // Reverse spin
        output.comRotRate = -output.comRotRate;
        cout << "OCI reverse scan" << endl;
    }

    // Get PPR offset
    swapc_bytes2((char *)&mceTelemetry[0][8], (char *)&buf, 4, 1);
    output.pprOffset = buf % MAX_ENC_COUNT;

    // Get TDI time and compute time increment per line
    int16_t tdiTime;
    swapc_bytes2((char *)&ddcTelemetry[0][346], (char *)&tdiTime, 2, 1);
    output.lineRate = (tdiTime + 1) / geoLut.masterClock;

    // Get valid encoder count, HAM and RTA encoder data
    output.hamEncData = new double *[numMceScans];
    output.rtaEncData = new double *[numMceScans];
    for (size_t i = 0; i < numMceScans; i++) {
        output.encoderSampleCounts.push_back(mceTelemetry[i][475]);
        output.hamEncData[i] = new double[mceBlockSize];
        output.rtaEncData[i] = new double[mceBlockSize];
        for (size_t j = 0; j < numEncoderChannels; j++) {
            output.hamEncData[i][j] = encoderData[i][4 * j + 0] * geoLut.hamEncoderScale;
            output.rtaEncData[i][j] = encoderData[i][4 * j + 1] * geoLut.rtaEncoderScale;
        }
    }

    return output;
}

GeoData locateOci(NcFile *l1aFile, Level1bFile &outfile, const oel::L1bOptions &options,
                  LocatingContext locatingContext, oel::ProcessingTracker &processingTracker) {
    GeoData geoData;
    NcGroup l1aScanLineAttributes = l1aFile->getGroup("scan_line_attributes");
    NcGroup l1aSpatialSpectralModes = l1aFile->getGroup("spatial_spectral_modes");
    NcGroup l1aEngineeringData = l1aFile->getGroup("engineering_data");
    NcGroup l1aNavigationData = l1aFile->getGroup("navigation_data");
    NcGroup l1aOnboardCalData = l1aFile->getGroup("onboard_calibration_data");
    NcGroup l1aScienceData = l1aFile->getGroup("science_data");

    NcGroupAtt attributeBuffer;
    NcVar variableBuffer;

    // Append call sequence to existing history
    string history = get_history(l1aFile);

    // Get date (change this when year and day are added to time field)
    string timeCoverageStart, timeCoverageEnd;
    attributeBuffer = l1aFile->getAtt("time_coverage_start");
    attributeBuffer.getValues(timeCoverageStart);
    geoData.unixTimeStart = isodate2unix(timeCoverageStart.c_str());
    attributeBuffer = l1aFile->getAtt("time_coverage_end");
    attributeBuffer.getValues(timeCoverageEnd);
    geoData.unixTimeEnd = isodate2unix(timeCoverageEnd.c_str());

    size_t pos = timeCoverageStart.find("T");
    if (pos == string::npos) {
        cout << "-E- L1A file time_coverage_start, bad format\n";
        exit(EXIT_FAILURE);
    }
    string granuleDate = timeCoverageStart.substr(0, pos);

    // Vector of the positions of OCI for 32 hours starting at timeCoverageStart with a resolution of 1 minute
    vector<vector<double>> ephPosVec;
    // Vector of the velocities of OCI for 32 hours starting at timeCoverageStart with a resolution of 1
    // minute
    vector<vector<double>> ephVelVec;
    vector<double> ephOTime;

    if (!options.ephFile.empty()) {
        geoData.unixTimeStart = isodate2unix(timeCoverageStart.c_str());
        geoData.unixTimeEnd = isodate2unix(timeCoverageEnd.c_str());

        int16_t year, month, day;
        double seconds;
        unix2ymds(geoData.unixTimeStart, &year, &month, &day, &seconds);
        double l1aEpoch = ymds2unix(year, month, day, 0.0);

        readDefinitiveEphemeris(options.ephFile, ephOTime, ephPosVec, ephVelVec, l1aEpoch);

        // Ensure ephem file covers L1A file time
        if ((geoData.unixTimeStart < (ephOTime.front() + l1aEpoch)) ||
            (geoData.unixTimeEnd > (ephOTime.back() + l1aEpoch))) {
            cerr << "-E- Definitive ephemeris file does not cover L1A file time" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // Scan time, spin ID and HAM side
    NcDim scansDim = l1aFile->getDim("number_of_scans");
    uint32_t numberOfScans = options.endingLine == -1 ? scansDim.getSize() : options.endingLine;
    vector<size_t> scansStart(1, 0);
    vector<size_t> scansEnd(1, numberOfScans);

    // Number of swir bands
    uint32_t numSwirBands = l1aFile->getDim("SWIR_bands").getSize();

    geoData.scanStartTimes = vector<double>(numberOfScans);
    l1aScanLineAttributes.getVar("scan_start_time")
        .getVar(scansStart, scansEnd, geoData.scanStartTimes.data());
    geoData.spinIds = vector<int32_t>(numberOfScans);
    l1aScanLineAttributes.getVar("spin_ID").getVar(scansStart, scansEnd, geoData.spinIds.data());
    geoData.hamSides = vector<uint8_t>(numberOfScans);
    l1aScanLineAttributes.getVar("HAM_side").getVar(scansStart, scansEnd, geoData.hamSides.data());

    uint32_t numValidScans = 0;
    for (size_t i = 0; i < numberOfScans; i++) {
        if (geoData.spinIds[i] > 0) {
            geoData.scanStartTimes[numValidScans] = geoData.scanStartTimes[i];
            geoData.hamSides[numValidScans] = geoData.hamSides[i];
            numValidScans++;
        }
    }
    geoData.numGoodScans = numValidScans;

    // Check for and fill in missing scan times
    vector<short> scanQualityFlags(numValidScans, 0);
    interpolateMissingScanTimes(geoData.scanStartTimes, scanQualityFlags);

    uint32_t numAttitudeRecords = l1aFile->getDim("att_records").getSize();
    vector<double> attTime(numAttitudeRecords);
    l1aNavigationData.getVar("att_time").getVar(attTime.data());

    uint32_t numQuaternionElements = l1aFile->getDim("quaternion_elements").getSize();
    ERROR_EXIT(numQuaternionElements);

    float **attitudeQuaternions = allocate2d_float(numAttitudeRecords, numQuaternionElements);

    l1aNavigationData.getVar("att_quat").getVar(&attitudeQuaternions[0][0]);

    uint32_t numOrbitRecords = l1aFile->getDim("orb_records").getSize();
    vector<double> orbitTimes(numOrbitRecords);
    if (options.ephFile.empty()) {
        l1aNavigationData.getVar("orb_time").getVar(orbitTimes.data());
    }

    orbArray *orbitPositions = new orbArray[numOrbitRecords];  // 3d position at each orbit record
    l1aNavigationData.getVar("orb_pos").getVar(&orbitPositions[0][0]);

    orbArray *orbitVelocities = new orbArray[numOrbitRecords];  // 3d velocity at each orbit record
    l1aNavigationData.getVar("orb_vel").getVar(&orbitVelocities[0][0]);
    uint32_t numTiltSamples = l1aFile->getDim("tilt_samples").getSize();

    vector<float> tiltAngle(numTiltSamples);
    l1aNavigationData.getVar("tilt").getVar(tiltAngle.data());

    vector<double> tiltTimes(numTiltSamples);  // Seconds of day
    l1aNavigationData.getVar("tilt_time").getVar(tiltTimes.data());

    // **************************** //
    // *** Read geolocation LUT *** //
    // **************************** //

    geoData.geoLut.file = new NcFile(options.geolocationLutFilename, NcFile::read);
    NcGroup timeParameters, coordinateTranslationParams, rtaHamParameters, planarityParams;

    timeParameters = geoData.geoLut.file->getGroup("time_params");
    timeParameters.getVar("master_clock")
        .getVar(&geoData.geoLut.masterClock);                             // Freq of OCI master clock in Hz
    timeParameters.getVar("MCE_clock").getVar(&geoData.geoLut.mceClock);  // Freq of OCI MCE clock in Hz

    coordinateTranslationParams = geoData.geoLut.file->getGroup("coord_trans");
    coordinateTranslationParams.getVar("sc_to_tilt").getVar(&geoData.geoLut.craftToTilt);
    coordinateTranslationParams.getVar("tilt_axis").getVar(&geoData.geoLut.tiltAxis);
    // Tilt angles at fixed positions (aft, forward)
    coordinateTranslationParams.getVar("tilt_angles").getVar(&geoData.geoLut.tiltAngles);
    coordinateTranslationParams.getVar("tilt_home").getVar(&geoData.geoLut.tiltHome);
    // Tilt platform to OCI mechanical transformation
    coordinateTranslationParams.getVar("tilt_to_oci_mech").getVar(&geoData.geoLut.tiltToOciMech);
    // OCI mechanical to optical transformation
    coordinateTranslationParams.getVar("oci_mech_to_oci_opt").getVar(&geoData.geoLut.ociMechToOciOpt);

    rtaHamParameters = geoData.geoLut.file->getGroup("RTA_HAM_params");
    rtaHamParameters.getVar("RTA_axis")
        .getVar(&geoData.geoLut.rtaAxis);  // Rotating Telescope Assembly rotation axis
    rtaHamParameters.getVar("HAM_axis").getVar(&geoData.geoLut.hamAxis);  // Half Angle Mirror rotation axis
    // Along-track mirror-to-axis angles
    rtaHamParameters.getVar("HAM_AT_angles").getVar(geoData.geoLut.hamAlongTrackAngles);
    // Cross-track mirror-to-axis angles
    rtaHamParameters.getVar("HAM_CT_angles").getVar(geoData.geoLut.hamCrossTrackAngles);
    // RTA encoder conversion to arcseconds
    rtaHamParameters.getVar("RTA_enc_scale").getVar(&geoData.geoLut.rtaEncoderScale);
    // HAM encoder conversion to arcseconds
    rtaHamParameters.getVar("HAM_enc_scale").getVar(&geoData.geoLut.hamEncoderScale);
    // Pulse per revolution offset from RTA nadir angle measured in encoder counts
    rtaHamParameters.getVar("RTA_nadir").getVar(geoData.geoLut.rtaNadir);

    planarityParams = geoData.geoLut.file->getGroup("planarity");
    planarityParams.getVar("along_scan_planarity").getVar(&geoData.geoLut.alongTrackPlanarity);
    planarityParams.getVar("along_track_planarity")
        .getVar(&geoData.geoLut.acrossTrackPlanarity);  // PACE prograde

    geoData.geoLut.file->close();

    geoData.mceTelem = readMceTelemetry(l1aFile, geoData.geoLut, l1aEngineeringData);

    int16_t year, dayOfYear;
    double seconds;
    unix2yds(isodate2unix(timeCoverageStart.c_str()), &year, &dayOfYear, &seconds);
    int intYear = year;
    int intDay = dayOfYear;

    double ecrMatrix[3][3];
    quaternion *quaternions = new quaternion[numAttitudeRecords];
    vector<orbArray> positions;   // Array of 1x3 vectors
    vector<orbArray> velocities;  // Array of 1x3 vectors

    // Transform orbit (if needed) and attitude from J2000 to ECR
    if (!options.ephFile.empty()) {
        size_t ephOTimeLBIndex = -1;
        size_t ephOTimeUBIndex = -1;

        for (size_t i = 0; i < ephOTime.size(); i++) {  // Find bracketing ephemeris records
            double value = ephOTime[i];

            if (value <= geoData.scanStartTimes[0]) {
                ephOTimeLBIndex = i;
            }
            if (value >= geoData.scanStartTimes[numberOfScans - 1]) {
                ephOTimeUBIndex = i;
                break;
            }
        }

        // pad the first and last record due to earth view time offset
        // to guarantee time coverage
        if (ephOTimeLBIndex > 0)
            ephOTimeLBIndex--;
        if (ephOTimeUBIndex < (ephOTime.size() - 2))
            ephOTimeUBIndex++;

        double omegaE = 7.29211585494e-5;
        size_t numEphRecords = ephOTimeUBIndex - ephOTimeLBIndex + 1;

        numOrbitRecords = numEphRecords;
        positions.resize(numOrbitRecords);
        velocities.resize(numOrbitRecords);

        for (size_t i = 0; i < positions.size(); i++) {
            for (size_t j = 0; j < positions[i].size(); j++) {
                positions[i][j] = 0;
                velocities[i][j] = 0;
            }
        }

        orbitTimes.resize(numOrbitRecords);

        for (size_t idxRec = 0; idxRec < numEphRecords; idxRec++) {
            size_t ephIndex = idxRec + ephOTimeLBIndex;
            orbitTimes[idxRec] = ephOTime[ephIndex];

            j2000ToEcr(year, dayOfYear, ephOTime[ephIndex], ecrMatrix);

            // Matrix multiplication. positions = ecrMatrix*ephPosVec, velocities = ecrMatrix*ephVelVec
            for (size_t idxLine = 0; idxLine < 3; idxLine++) {
                for (size_t idxValue = 0; idxValue < 3; idxValue++) {
                    positions[idxRec][idxLine] +=
                        ecrMatrix[idxLine][idxValue] * ephPosVec[ephIndex][idxValue];
                    velocities[idxRec][idxLine] +=
                        ecrMatrix[idxLine][idxValue] * ephVelVec[ephIndex][idxValue];
                }
            }

            // Correction for angular momentum along x and y
            velocities[idxRec][0] = velocities[idxRec][0] + positions[idxRec][1] * omegaE;
            velocities[idxRec][1] = velocities[idxRec][1] - positions[idxRec][0] * omegaE;
        }
    } else {
        positions.resize(numOrbitRecords);
        velocities.resize(numOrbitRecords);
        // Orbit transformation
        for (size_t i = 0; i < numOrbitRecords; i++) {
            for (size_t j = 0; j < 3; j++) {
                positions[i][j] = orbitPositions[i][j] / 1000;
                velocities[i][j] = orbitVelocities[i][j] / 1000;
            }
        }  // i loop
    }

    delete[] orbitPositions;
    delete[] orbitVelocities;

    // Attitude
    for (size_t i = 0; i < numAttitudeRecords; i++) {
        double ecrQuaternion[4], resultQuaternion[4];
        float attitudeQuaternion[4];
        j2000ToEcr(year, dayOfYear, attTime[i], ecrMatrix);
        matrixToQuaternion(ecrMatrix, ecrQuaternion);

        memcpy(attitudeQuaternion, &attitudeQuaternions[i][0], 3 * sizeof(float));
        attitudeQuaternion[3] = attitudeQuaternions[i][3];

        multiplyQuaternions(ecrQuaternion, attitudeQuaternion, resultQuaternion);
        for (size_t j = 0; j < 4; j++)
            quaternions[i][j] = resultQuaternion[j];
    }  // i loop

    // ******************************************** //
    // *** Get spatial and spectral aggregation *** //
    // ******************************************** //
    NcDim spatialZones = l1aFile->getDim("spatial_zones");
    uint32_t numSpatialZones = spatialZones.getSize();

    vector<int16_t> spatialZoneDataType(numSpatialZones);
    l1aSpatialSpectralModes.getVar("spatial_zone_data_type").getVar(spatialZoneDataType.data());

    vector<int16_t> spatialZoneLines(numSpatialZones);
    l1aSpatialSpectralModes.getVar("spatial_zone_lines").getVar(spatialZoneLines.data());

    vector<int16_t> spatialAggregation(numSpatialZones);
    l1aSpatialSpectralModes.getVar("spatial_aggregation").getVar(spatialAggregation.data());

    vector<float> scienceLines(32400), swirLines(4050);
    uint16_t numHyperSciPix, numSwirPix;
    double earthViewTimeOffset;
    vector<double> sciencePixelOffset(32400), swirPixelOffset(4050);
    bool isDark = false;
    int earthViewStatus = getEarthView(geoData.mceTelem.lineRate, spatialZoneDataType.data(),
                                       spatialZoneLines.data(), spatialAggregation.data(), numHyperSciPix,
                                       numSwirPix, earthViewTimeOffset, scienceLines.data(), swirLines.data(),
                                       sciencePixelOffset.data(), swirPixelOffset.data(), isDark);
    if (earthViewStatus < 0) {
        cout << "No Earth view in file: " << l1aFile->getName() << endl;
        l1aFile->close();
        return geoData;
    }
    geoData.numCcdPix = numHyperSciPix;
    geoData.numSwirPix = numSwirPix;
    geoData.earthViewTimeOffset = earthViewTimeOffset;
    geoData.isDark = isDark;
    geoData.sciPixOffset = sciencePixelOffset;
    geoData.swirPixOffset = swirPixelOffset;

    vector<double> earthViewTimes(numValidScans);  // Times in which Earth was in view of the RTA
    for (size_t i = 0; i < numValidScans; i++)
        if (geoData.scanStartTimes[i] == BAD_FLT)
            earthViewTimes[i] = geoData.scanStartTimes[i];
        else
            earthViewTimes[i] = geoData.scanStartTimes[i] + earthViewTimeOffset;
    geoData.earthViewTimes = earthViewTimes;

    // Interpolate orbit, attitude and tilt to scan times
    vector<orbArray> interpolatedPos(numValidScans);
    vector<orbArray> interpolatedVel(numValidScans);
    vector<orbArray> attitudeAngles(numValidScans);  // Roll, pitch, yaw
    interpolateOrbitVectors(numOrbitRecords, numValidScans, orbitTimes, positions, velocities, earthViewTimes,
                            interpolatedPos, interpolatedVel);

    quaternion *interpolatedQuats = new quaternion[numValidScans]();
    interpolateAttitudeForScanTimes(numAttitudeRecords, numValidScans, attTime.data(), quaternions,
                                    earthViewTimes.data(), interpolatedQuats);

    vector<float> tiltPositions(numValidScans);
    interpolateTilt(numTiltSamples, numValidScans, tiltTimes.data(), tiltAngle.data(), earthViewTimes.data(),
                    tiltPositions.data());

    vector<orbArray> sunVectors(numValidScans);
    vector<float> earthSunDist(numValidScans, 0.0);
    geoData.qualityFlag = vector<uint8_t>(numValidScans * numHyperSciPix, 0);

    for (size_t i = 0; i < numValidScans; i++) {
        tiltPositions[i] += geoData.geoLut.tiltHome;  // Add tilt home position to angles
        tiltPositions[i] =
            tiltPositions[i] < geoData.geoLut.tiltAngles[0] ? geoData.geoLut.tiltAngles[0] : tiltPositions[i];
        tiltPositions[i] =
            tiltPositions[i] > geoData.geoLut.tiltAngles[1] ? geoData.geoLut.tiltAngles[1] : tiltPositions[i];
        // for any 2 adjecent and = tilts clear tilt flag (ie, bit AND with 110b)
        if ((i > 0) && (tiltPositions[i - 1] == tiltPositions[i])) {
            scanQualityFlags[i - 1] &= 0b110;
            scanQualityFlags[i] &= 0b110;
        } else
            scanQualityFlags[i] |= 1;

        // Setting up geoData.auCorrection
        l_sun_(&intYear, &intDay, &earthViewTimes.data()[i], sunVectors[i].data(), &earthSunDist[i]);
    }

    double earthSunDistCorrection = pow(earthSunDist[numValidScans / 2], 2);
    geoData.auCorrection = earthSunDistCorrection;

    // Declaring these out here saves memory

    double tiltMatrix[3][3];  // Tilt rotation
    float tiltQuaternion[4];  // Tile rotation
    gsl_matrix_view sourceMatrix;
    gsl_matrix_view transformMatrix;
    double *sensorOrientationMatrix;
    double qmat[3][3];  // Temporary variable
    bool eclipseMessageShown = false;
    NcFile *watermaskFilePointer = new NcFile(options.demFile, NcFile::read);
    const size_t WATERMASK_LAT_POINTS = watermaskFilePointer->getDim("lat").getSize();
    const size_t WATERMASK_LON_POINTS = watermaskFilePointer->getDim("lon").getSize();
    double constexpr DELTA_LAT = 180.0;  // Degrees of change in latitude over whole Earth
    double constexpr DELTA_LON = 360.0;  // Degrees of change in longitude over whole Earth
    vector<array<float, 3>> pointingVector(numHyperSciPix);
    double earthMoonDist;
    const double ASTRONOMICAL_UNIT = 149597870.7;
    uint8_t *eclipsedScans = new uint8_t[numValidScans];
    bool fileIsEclipsed = false;  // File has any of its scanlines eclipsed
    geoData.pixelLongtitudes = vector<float>(numValidScans * numHyperSciPix, BAD_FLT);
    geoData.pixelLatitudes = vector<float>(numValidScans * numHyperSciPix, BAD_FLT);
    geoData.solarZeniths = vector<float>(numValidScans * numHyperSciPix, BAD_FLT);
    geoData.solarAzimuths = vector<float>(numValidScans * numHyperSciPix, BAD_FLT);
    geoData.sensorZeniths = vector<float>(numValidScans * numHyperSciPix, BAD_FLT);
    geoData.sensorAzimuths = vector<float>(numValidScans * numHyperSciPix, BAD_FLT);
    geoData.height = vector<short>(numValidScans * numHyperSciPix, BAD_FLT);
    vector<uint8_t> watermask(numValidScans * numHyperSciPix);
    geoData.ccdScanAngles = vector<double>(numValidScans * numHyperSciPix, BAD_FLT);
    geoData.swirScanAngles = vector<double>(numValidScans * numSwirPix, BAD_FLT);
    Gring gringHelper;

    gsl_matrix *craftToTilt = gsl_matrix_alloc(3, 3);      // From spacecraft to tilt platform
    gsl_matrix *craftToOciMount = gsl_matrix_alloc(3, 3);  // From spacecraft to OCI mount
    gsl_matrix *craftToOci = gsl_matrix_alloc(3, 3);       // From spacecraft to OCI instrument
    gsl_matrix *earthTransform = gsl_matrix_alloc(3, 3);
    gsl_vector *positionRelativeToMoon = gsl_vector_alloc(3);
    gsl_vector *velocityRelativeToMoon = gsl_vector_alloc(3);
    gsl_vector *moonPosRotating = gsl_vector_alloc(3);            // Moon's position in ECR
    gsl_matrix *selenographicTransform = gsl_matrix_alloc(3, 3);  // Becomes sensor orientation matrix

    //////////////////////////////
    //  Locate each scan line   //
    //////////////////////////////

    switch (locatingContext) {
        case SELENO:
            cout << "Selenolocating file" << endl;
            break;
        case HELIO:
            cout << "Processing sensor and solar angles" << endl;
            break;
        case GEO:
        default:
            cout << "Geolocating file" << endl;
            break;
    }

    try {
        processingTracker = oel::ProcessingTracker(options.startingLine, numValidScans, 0.10);
    } catch (const invalid_argument &e) {
        cerr << "-E- geolocate_oci.cpp: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    size_t startingLine = processingTracker.getStartingLine();
    size_t endingLine = processingTracker.getEndingLine();
    for (size_t currScan = startingLine; currScan < endingLine; currScan++) {
        processingTracker.update();

        if (geoData.scanStartTimes[currScan] == -999.0) {
            scanQualityFlags[currScan] |= 2;
            continue;
        }

        orbArray selectedPosition;  // Relative to either Earth or Moon
        orbArray selectedVelocity;  // Relative to either Earth or Moon

        gsl_matrix *ecrToSelenographic = gsl_matrix_alloc(3, 3);

        // Prepare for eclipse calculation or selenographic transformations
        orbArray moonVector;
        getMoonVector(year, dayOfYear, geoData.scanStartTimes[currScan], moonVector,
                      earthMoonDist);  // l_moon

        if (locatingContext == GEO) {
            double sunMoonOverlapFactor = sunVectors[currScan][0] * moonVector[0] +
                                          sunVectors[currScan][1] * moonVector[1] +
                                          sunVectors[currScan][2] * moonVector[2];

            if (sunMoonOverlapFactor > 0.9995) {
                fileIsEclipsed = true;
                eclipsedScans[currScan] = 1;
            } else {
                eclipsedScans[currScan] = 0;
            }

            if (fileIsEclipsed && !eclipseMessageShown) {
                cout << "Eclipse detected" << endl;
                eclipseMessageShown = true;
            }
        }
        tiltQuaternion[0] = geoData.geoLut.tiltAxis[0] * sin(tiltPositions[currScan] * OEL_DEGRAD / 2);
        tiltQuaternion[1] = geoData.geoLut.tiltAxis[1] * sin(tiltPositions[currScan] * OEL_DEGRAD / 2);
        tiltQuaternion[2] = geoData.geoLut.tiltAxis[2] * sin(tiltPositions[currScan] * OEL_DEGRAD / 2);
        tiltQuaternion[3] = cos(tiltPositions[currScan] * OEL_DEGRAD / 2);

        quaternionToMatrix(tiltQuaternion, tiltMatrix);

        // Combine tilt and alignments

        // sc2tiltp = tiltm#geo_lut.craftToTilt
        sourceMatrix = gsl_matrix_view_array((double *)geoData.geoLut.craftToTilt, 3, 3);
        transformMatrix = gsl_matrix_view_array((double *)tiltMatrix, 3, 3);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &sourceMatrix.matrix, &transformMatrix.matrix, 0.0,
                       craftToTilt);
        sensorOrientationMatrix = gsl_matrix_ptr(craftToTilt, 0, 0);

        // sc2ocim = geo_lut.tilt_to_oci_mech#sc2tiltp
        transformMatrix = gsl_matrix_view_array((double *)geoData.geoLut.tiltToOciMech, 3, 3);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, craftToTilt, &transformMatrix.matrix, 0.0,
                       craftToOciMount);
        sensorOrientationMatrix = gsl_matrix_ptr(craftToOciMount, 0, 0);

        // sc_to_oci = geo_lut.oci_mech_to_oci_opt#sc2ocim
        transformMatrix = gsl_matrix_view_array((double *)geoData.geoLut.ociMechToOciOpt, 3, 3);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, craftToOciMount, &transformMatrix.matrix, 0.0,
                       craftToOci);
        sensorOrientationMatrix = gsl_matrix_ptr(craftToOci, 0, 0);

        // Convert quaternion to matrix
        quaternionToMatrix(interpolatedQuats[currScan], qmat);

        // smat = sc_to_oci#qmat
        transformMatrix = gsl_matrix_view_array(&qmat[0][0], 3, 3);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, craftToOci, &transformMatrix.matrix, 0.0,
                       earthTransform);

        // Compute attitude angles (informational only)
        getAttitudeAngles(interpolatedPos[currScan], interpolatedVel[currScan], qmat,
                          attitudeAngles[currScan]);

        // Get scan ellipse coefficients
        double scanEllipseCoefs[10];

        if (locatingContext == GEO) {
            sensorOrientationMatrix = gsl_matrix_ptr(earthTransform, 0, 0);
            getEllipsoidScanTrackCoefs(interpolatedPos[currScan], (double (*)[3])sensorOrientationMatrix,
                                       locatingContext, scanEllipseCoefs);

            selectedPosition = interpolatedPos[currScan];
            selectedVelocity = interpolatedVel[currScan];

        } else if (locatingContext == SELENO) {
            // Compute selenographic transformation
            ecrToSeleno(year, dayOfYear, seconds, ecrToSelenographic);  // ecr_to_seleno

            // smat = sc_to_oci#qmat#transpose(selmat[*,*,iscn])
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, earthTransform, ecrToSelenographic, 0.0,
                           selenographicTransform);
            sensorOrientationMatrix = gsl_matrix_ptr(selenographicTransform, 0, 0);

            //   posmr = posi[*,iscn] - xmoon[*,iscn]*rmoon[iscn]
            gsl_vector *scaledMoonVec = getGslVector(moonVector, earthMoonDist);  // xmoon[*,iscn]*rmoon[iscn]
            gsl_vector *interpolatedPosGsl = getGslVector(interpolatedPos[currScan]);
            gsl_vector_sub(interpolatedPosGsl, scaledMoonVec);

            // posm = selmat[*,*,iscn]#posmr
            gsl_blas_dgemv(CblasNoTrans, 1.0, ecrToSelenographic, interpolatedPosGsl, 0.0,
                           positionRelativeToMoon);

            // velm = selmat[*,*,iscn]*veli[*,iscn]
            gsl_vector *interpolatedVelGsl = getGslVector(interpolatedVel[currScan]);
            gsl_blas_dgemv(CblasNoTrans, 1.0, ecrToSelenographic, interpolatedVelGsl, 0.0,
                           velocityRelativeToMoon);

            orbArray velRelMoonAsOrbArray;
            velRelMoonAsOrbArray[0] = static_cast<float>(velocityRelativeToMoon->data[0]);
            velRelMoonAsOrbArray[1] = static_cast<float>(velocityRelativeToMoon->data[1]);
            velRelMoonAsOrbArray[2] = static_cast<float>(velocityRelativeToMoon->data[2]);

            orbArray posRelMoonAsOrbArray;
            posRelMoonAsOrbArray[0] = static_cast<float>(positionRelativeToMoon->data[0]);
            posRelMoonAsOrbArray[1] = static_cast<float>(positionRelativeToMoon->data[1]);
            posRelMoonAsOrbArray[2] = static_cast<float>(positionRelativeToMoon->data[2]);

            getEllipsoidScanTrackCoefs(posRelMoonAsOrbArray, (double (*)[3])sensorOrientationMatrix,
                                       locatingContext, scanEllipseCoefs);
            selectedPosition = posRelMoonAsOrbArray;
            selectedVelocity = velRelMoonAsOrbArray;
        }

        // Generate pointing vector and relative time arrays in instrument frame

        /* It's possible to refactor getEarthViewVectors so we only have one call here, but it might be less
           readable in total */

        try {
            int returnStatus = getEarthViewVectors(geoData, geoData.geoLut, CCD, earthViewTimeOffset,
                                                   sciencePixelOffset, currScan, pointingVector);

            scanQualityFlags[currScan] |= 4 * returnStatus;

            returnStatus = getEarthViewVectors(geoData, geoData.geoLut, SWIR, earthViewTimeOffset,
                                               sciencePixelOffset, currScan, pointingVector);
            scanQualityFlags[currScan] |= 4 * returnStatus;
        } catch (const invalid_argument &e) {
            cout << "-E- getEarthViewVectors failed: " << e.what() << endl;
            exit(EXIT_FAILURE);
        }

        // Geolocate pixels
        size_t pixelIndex = currScan * numHyperSciPix;  // Index of first pixel in this scan
        if (!options.disableGeolocation) {
            orbArray moonVecRotating;  // Vector from Earth to Moon in ECR coords
            if (fileIsEclipsed) {
                moonVecRotating[0] = moonVector[0] * earthMoonDist;
                moonVecRotating[1] = moonVector[1] * earthMoonDist;
                moonVecRotating[2] = moonVector[2] * earthMoonDist;
            } else {
                moonVecRotating[0] = moonVector[0];
                moonVecRotating[1] = moonVector[1];
                moonVecRotating[2] = moonVector[2];
            }

            orbArray selectedSunVector;

            if (locatingContext == SELENO) {
                gsl_vector *sunM = gsl_vector_alloc(3);
                for (size_t i = 0; i < sunM->size; i++) {
                    gsl_vector_set(sunM, i,
                                   sunVectors[currScan][i] * earthSunDist[currScan] * ASTRONOMICAL_UNIT -
                                       moonVector[i] * gsl_vector_get(moonPosRotating, i));
                }

                gsl_vector *sunVecSeleno = gsl_vector_alloc(3);
                gsl_blas_dgemv(CblasNoTrans, 1.0, selenographicTransform, sunM, 0.0, sunVecSeleno);
                gsl_vector_mul(sunM, sunM);
                double totalSunMSquared = gsl_vector_sum(sunM);
                double sqrtTotalSunMSquared = sqrt(totalSunMSquared);
                gsl_vector_scale(sunVecSeleno, 1 / sqrtTotalSunMSquared);
                orbArray temp;
                temp[0] = static_cast<float>(gsl_vector_get(sunVecSeleno, 0));
                temp[1] = static_cast<float>(gsl_vector_get(sunVecSeleno, 0));
                temp[2] = static_cast<float>(gsl_vector_get(sunVecSeleno, 0));
                selectedSunVector = temp;
            } else {
                selectedSunVector = sunVectors[currScan];
            }
            locatePixelsOci(options.demFile.c_str(), selectedPosition, selectedVelocity,
                            (double (*)[3])sensorOrientationMatrix, scanEllipseCoefs, locatingContext,
                            selectedSunVector, fileIsEclipsed, moonVector, earthSunDist[currScan],
                            pointingVector, numHyperSciPix, swirPixelOffset.data(), geoData.qualityFlag,
                            currScan, geoData.pixelLatitudes.data() + pixelIndex,
                            geoData.pixelLongtitudes.data() + pixelIndex, &geoData.solarZeniths[pixelIndex],
                            geoData.solarAzimuths.data() + pixelIndex,
                            geoData.sensorZeniths.data() + pixelIndex,
                            geoData.sensorAzimuths.data() + pixelIndex, geoData.height.data() + pixelIndex);
        }

        for (size_t i = 0; i < numHyperSciPix; i++) {
            const double lat = geoData.pixelLatitudes[pixelIndex + i];
            const double lon = geoData.pixelLongtitudes[pixelIndex + i];

            if (lat == BAD_FLT || lon == BAD_FLT) {
                watermask[currScan * numHyperSciPix + i] = BAD_UBYTE;
                continue;
            }

            // Putting lat and lon into a 360 degree space makes some calculations easier
            size_t latIndex = static_cast<size_t>((lat + 90.0) / DELTA_LAT * (WATERMASK_LAT_POINTS - 1));
            size_t lonIndex = static_cast<size_t>((lon + 180.0) / DELTA_LON * (WATERMASK_LON_POINTS - 1));

            latIndex = std::min(latIndex, WATERMASK_LAT_POINTS - 1);
            lonIndex = std::min(lonIndex, WATERMASK_LON_POINTS - 1);

            // Only grabbing one point because the watermask probably has a different interval
            watermaskFilePointer->getVar("watermask")
                .getVar({latIndex, lonIndex}, {1, 1}, &watermask[pixelIndex + i]);
        }

        gringHelper.processScan(geoData.pixelLatitudes.data() + pixelIndex,
                                geoData.pixelLongtitudes.data() + pixelIndex, geoData.numCcdPix, currScan,
                                geoData.numGoodScans);
    }  // scan loop

    // Gring wants the very last scan but the stuff above doesn't
    int alsoPixelIndex = (numValidScans - 1) * numHyperSciPix;  // Final scan
    gringHelper.processScan(geoData.pixelLatitudes.data() + alsoPixelIndex,
                            geoData.pixelLongtitudes.data() + alsoPixelIndex, geoData.numCcdPix,
                            numValidScans, geoData.numGoodScans);

    switch (locatingContext) {
        case SELENO:
            cout << "Selenolocation finished, cleaning up and writing data" << endl;
            break;
        case HELIO:
            cout << "Angles processing finished, cleaning up and writing data" << endl;
            break;
        case GEO:
        default:
            cout << "Geolocation finished, cleaning up and writing data" << endl;
            break;
    }

    gsl_matrix_free(earthTransform);
    gsl_matrix_free(craftToTilt);
    gsl_matrix_free(craftToOciMount);
    gsl_matrix_free(craftToOci);

    // Get number of bands
    NcDim numTapsDim = l1aFile->getDim("number_of_taps");
    const size_t numTaps = numTapsDim.getSize();

    vector<int16_t> blueBands(numTaps);
    l1aSpatialSpectralModes.getVar("blue_spectral_mode").getVar(blueBands.data());

    vector<int16_t> redBands(numTaps);
    l1aSpatialSpectralModes.getVar("red_spectral_mode").getVar(redBands.data());

    vector<size_t> tapAggFactors(numTaps);
    array<size_t, 16> binCounts;
    size_t numBlueInsBands, numRedInsBands;
    size_t numBlueBands, numRedBands;
    aggregateBands(numTaps, tapAggFactors.data(), binCounts.data(), blueBands.data(), numBlueInsBands,
                   numBlueBands);
    aggregateBands(numTaps, tapAggFactors.data(), binCounts.data(), redBands.data(), numRedInsBands,
                   numRedBands);

    vector<size_t> start(3, 0);
    vector<size_t> count(1, numValidScans);

    outfile.createFile(numValidScans, numBlueBands, numRedBands, numHyperSciPix, numSwirPix, numSwirBands,
                       geoData.geoLut.rtaNadir, options.radianceGenerationEnabled, locatingContext,
                       options.deflate);

    variableBuffer = outfile.scanLineAttributes.getVar("time");
    variableBuffer.putVar(start, count, earthViewTimes.data());
    variableBuffer.putAtt("units", "seconds since " + granuleDate);

    outfile.scanLineAttributes.getVar("HAM_side").putVar(start, count, geoData.hamSides.data());
    outfile.scanLineAttributes.getVar("scan_quality_flags").putVar(start, count, scanQualityFlags.data());
    outfile.navigationData.getVar("tilt_angle").putVar(start, count, tiltPositions.data());
    count.push_back(4);
    outfile.navigationData.getVar("att_quat").putVar(start, count, interpolatedQuats);

    count[1] = 3;

    outfile.navigationData.getVar("att_ang").putVar(start, count, attitudeAngles.data());
    outfile.navigationData.getVar("orb_pos").putVar(start, count, interpolatedPos.data());
    outfile.navigationData.getVar("orb_vel").putVar(start, count, interpolatedVel.data());
    outfile.navigationData.getVar("sun_ref").putVar(start, count, sunVectors.data());

    count[1] = numHyperSciPix;

    outfile.geolocationData.getVar("latitude").putVar(start, count, geoData.pixelLatitudes.data());
    outfile.geolocationData.getVar("longitude").putVar(start, count, geoData.pixelLongtitudes.data());
    outfile.geolocationData.getVar("quality_flag").putVar(start, count, geoData.qualityFlag.data());
    {  // Delete this when we change the L1B product to float. Don't forget to delete algorithm unless used
       // elsewhere
        const auto scaleBy100 = [](float angle) {
            if (angle != BAD_FLT)  // Wouldn't want to change the value if it's fill already
                return angle * 100.0f;
            else
                return angle;
        };
        // &: captures local variables by reference
        auto scaleAndWrite = [&](const string &varName, const vector<float> &data) {  // lambdas are cool
            vector<float> scaledData(data.size());
            transform(data.begin(), data.end(), scaledData.begin(), scaleBy100);
            outfile.geolocationData.getVar(varName).putVar(start, count, scaledData.data());
        };

        scaleAndWrite("sensor_azimuth", geoData.sensorAzimuths);
        scaleAndWrite("sensor_zenith", geoData.sensorZeniths);
        scaleAndWrite("solar_azimuth", geoData.solarAzimuths);
        scaleAndWrite("solar_zenith", geoData.solarZeniths);
    }
    // Uncomment these when we change the L1B product to float
    // outfile.geolocationData.getVar("sensor_azimuth").putVar(start, count, geoData.sensorAzimuths.data());
    // outfile.geolocationData.getVar("sensor_zenith").putVar(start, count, geoData.sensorZeniths.data());
    // outfile.geolocationData.getVar("solar_azimuth").putVar(start, count, geoData.solarAzimuths.data());
    // outfile.geolocationData.getVar("solar_zenith").putVar(start, count, geoData.solarZeniths.data());
    outfile.geolocationData.getVar("height").putVar(start, count, geoData.height.data());
    outfile.geolocationData.getVar("watermask").putVar(start, count, watermask.data());

    // write global attributes, including history and date_created
    set_global_attrs(outfile.l1bFile, history, options.digitalObjectId, options.processingVersion);

    outfile.l1bFile->putAtt("earth_sun_distance_correction", NC_DOUBLE, 1, &geoData.auCorrection);

    array<float, 4> geospatialExtremes{90, -90, 180, -180};
    if (locatingContext == GEO || locatingContext == SELENO) {  // Location is possible for these contexts
        string latitudeCsvList;
        string longitudeCsvList;
        string sequence;
        gringHelper.getGringStrings(longitudeCsvList, latitudeCsvList, sequence);

        string polygon = gringHelper.getWktString(latitudeCsvList, longitudeCsvList, sequence);
        outfile.l1bFile->putAtt("geospatial_bounds_crs", "EPSG:4326");
        outfile.l1bFile->putAtt("geospatial_bounds", polygon);

        geospatialExtremes = gringHelper.getGeospatialExtremes();
    }

    outfile.l1bFile->putAtt("geospatial_lat_min", NC_FLOAT, 1, &geospatialExtremes[0]);
    outfile.l1bFile->putAtt("geospatial_lat_max", NC_FLOAT, 1, &geospatialExtremes[1]);
    outfile.l1bFile->putAtt("geospatial_lon_min", NC_FLOAT, 1, &geospatialExtremes[2]);
    outfile.l1bFile->putAtt("geospatial_lon_max", NC_FLOAT, 1, &geospatialExtremes[3]);

    delete[] (quaternions);

    free2d_float(attitudeQuaternions);
    processingTracker.reset();

    return geoData;
}

int j2000ToEcr(int32_t year, int32_t dayOfYear, double secondOfDay, double ecrMatrix[3][3]) {
    double j2000ToMeanOfDateMatrix[3][3];
    double nutationMatrix[3][3], ut1MinusUtc;

    // Get transformation from J2000 to mean-of-date inertial
    j2000ToMod(year, dayOfYear, secondOfDay, j2000ToMeanOfDateMatrix);
    // Get nutation and UT1-UTC (once per run)
    getNutation(year, dayOfYear, nutationMatrix);
    getUt1ToUtc(year, dayOfYear, ut1MinusUtc);

    // Compute Greenwich hour angle for time of day
    double fractionalDay = dayOfYear + (secondOfDay + ut1MinusUtc) / 86400;
    double greenwichHourAngle, greenwichHourAngleMatrix[3][3];
    gha2000(year, fractionalDay, &greenwichHourAngle);

    double cosGHA = cos(greenwichHourAngle / OEL_RADEG);
    double sinGHA = sin(greenwichHourAngle / OEL_RADEG);

    greenwichHourAngleMatrix[0][0] = cosGHA;
    greenwichHourAngleMatrix[1][1] = cosGHA;
    greenwichHourAngleMatrix[2][2] = 1;
    greenwichHourAngleMatrix[0][1] = sinGHA;
    greenwichHourAngleMatrix[1][0] = -sinGHA;

    greenwichHourAngleMatrix[0][2] = 0;
    greenwichHourAngleMatrix[2][0] = 0;
    greenwichHourAngleMatrix[1][2] = 0;
    greenwichHourAngleMatrix[2][1] = 0;
    // Combine all transformations
    gsl_matrix_view ghaMatrixView = gsl_matrix_view_array(&greenwichHourAngleMatrix[0][0], 3, 3);
    gsl_matrix_view nutationMatrixView = gsl_matrix_view_array(&nutationMatrix[0][0], 3, 3);
    gsl_matrix *tempMatrix1 = gsl_matrix_alloc(3, 3);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &ghaMatrixView.matrix, &nutationMatrixView.matrix, 0.0,
                   tempMatrix1);

    gsl_matrix_view j2000ToMeanOfDateMatrixView = gsl_matrix_view_array(&j2000ToMeanOfDateMatrix[0][0], 3, 3);
    gsl_matrix *tempMatrix2 = gsl_matrix_alloc(3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempMatrix1, &j2000ToMeanOfDateMatrixView.matrix, 0.0,
                   tempMatrix2);
    double *finalMatrixPtr = gsl_matrix_ptr(tempMatrix2, 0, 0);

    memcpy(ecrMatrix, finalMatrixPtr, 9 * sizeof(double));

    gsl_matrix_free(tempMatrix1);
    gsl_matrix_free(tempMatrix2);

    return 0;
}

int j2000ToMod(int32_t year, int32_t idy, double sec, double j2mod[3][3]) {
    int16_t iyear = year;
    int16_t day = idy;

    double t = jday(iyear, 1, day) - J_2000_MODIFIED + sec / 86400;
    t /= 36525;

    double zeta0 = t * (2306.2181 + 0.302 * t + 0.018 * t * t) / OEL_RADEG / 3600;
    double thetap = t * (2004.3109 - 0.4266 * t - 0.04160 * t * t) / OEL_RADEG / 3600;
    double xip = t * (2306.2181 + 1.095 * t + 0.018 * t * t) / OEL_RADEG / 3600;

    double sinZeta0 = sin(zeta0);
    double cosZeta0 = cos(zeta0);
    double sinXip = sin(xip);
    double cosXip = cos(xip);
    double sinThetaP = sin(thetap);
    double cosThetaP = cos(thetap);

    j2mod[0][0] = -sinZeta0 * sinXip + cosZeta0 * cosXip * cosThetaP;
    j2mod[0][1] = -cosZeta0 * sinXip - sinZeta0 * cosXip * cosThetaP;
    j2mod[0][2] = -cosXip * sinThetaP;
    j2mod[1][0] = sinZeta0 * cosXip + cosZeta0 * sinXip * cosThetaP;
    j2mod[1][1] = cosZeta0 * cosXip - sinZeta0 * sinXip * cosThetaP;
    j2mod[1][2] = -sinXip * sinThetaP;
    j2mod[2][0] = cosZeta0 * sinThetaP;
    j2mod[2][1] = -sinZeta0 * sinThetaP;
    j2mod[2][2] = cosThetaP;

    return 0;
}

int getNutation(int32_t year, int32_t idy, double xnut[3][3]) {
    int16_t iyear = year;
    int16_t day = idy;

    double t = jday(iyear, 1, day) - J_2000_MODIFIED;

    double xls, gs, xlm, omega;
    double dpsi, eps, epsm;
    ephparms(t, &xls, &gs, &xlm, &omega);
    nutate(t, xls, gs, xlm, omega, &dpsi, &eps);
    epsm = (double)23.439291 - ((double)3.560e-7) * t;

    double cos_dpsi_over_radeg = cos(dpsi / OEL_RADEG);
    double sin_dpsm_over_radeg = sin(dpsi / OEL_RADEG);
    double cos_epsm_over_radeg = cos(epsm / OEL_RADEG);
    double sin_epsm_over_radeg = sin(epsm / OEL_RADEG);
    double cos_eps_over_radeg = cos(eps / OEL_RADEG);
    double sin_eps_over_radeg = sin(eps / OEL_RADEG);

    xnut[0][0] = cos_dpsi_over_radeg;
    xnut[1][0] = -sin_dpsm_over_radeg * cos_epsm_over_radeg;
    xnut[2][0] = -sin_dpsm_over_radeg * sin_epsm_over_radeg;
    xnut[0][1] = sin_dpsm_over_radeg * cos_eps_over_radeg;
    xnut[1][1] = cos_dpsi_over_radeg * cos_eps_over_radeg * cos_epsm_over_radeg +
                 sin_eps_over_radeg * sin_epsm_over_radeg;
    xnut[2][1] = cos_dpsi_over_radeg * cos_eps_over_radeg * sin_epsm_over_radeg -
                 sin_eps_over_radeg * cos_epsm_over_radeg;
    xnut[0][2] = sin_dpsm_over_radeg * sin_eps_over_radeg;
    xnut[1][2] = cos_dpsi_over_radeg * sin_eps_over_radeg * cos_epsm_over_radeg -
                 cos_eps_over_radeg * sin_epsm_over_radeg;
    xnut[2][2] = cos_dpsi_over_radeg * sin_eps_over_radeg * sin_epsm_over_radeg +
                 cos_eps_over_radeg * cos_epsm_over_radeg;

    return 0;
}

int getUt1ToUtc(int32_t year, int32_t dayOfYear, double &ut1ToUtc) {
    int16_t iyear = year;
    int16_t day = dayOfYear;
    static int32_t julianDay[25000];
    static double ut1[25000];
    string utcDataFile = "$OCVARROOT/modis/utcpole.dat";
    static bool firstCall = true;
    int k = 0;

    if (firstCall) {
        expandEnvVar(&utcDataFile);

        ifstream utcpfile(utcDataFile.c_str());
        if (!utcpfile.is_open()) {
            cout << utcDataFile.c_str() << " not found" << endl;
            exit(1);
        }

        istringstream istr;
        string line;

        getline(utcpfile, line);
        getline(utcpfile, line);
        while (getline(utcpfile, line)) {
            istr.clear();
            istr.str(line.substr(0, 5));
            istr >> julianDay[k];
            istr.clear();
            istr.str(line.substr(44, 9));
            istr >> ut1[k];
            k++;
        }
        julianDay[k] = 0;
        utcpfile.close();
        firstCall = false;

    }  // if (first)

    k = 0;
    int daysSinceEpoch = jday(iyear, 1, day) - 2400000;
    while (julianDay[k] > 0) {
        if (daysSinceEpoch == julianDay[k]) {
            ut1ToUtc = ut1[k];
            break;
        }
        k++;
    }

    return 0;
}

int matrixToQuaternion(double rotationMatrix[3][3], double quaternion[4]) {
    //  Convert direction cosine matrix to equivalent quaternion

    double eulerAxis[3];

    //  Compute Euler angle
    double phi;  // Euler angle
    double cphi = (rotationMatrix[0][0] + rotationMatrix[1][1] + rotationMatrix[2][2] - 1) / 2;
    if (fabs(cphi) < 0.98) {
        phi = acos(cphi);
    } else {
        // Squared sine of Euler angle (called phi)
        double ssphi = (pow(rotationMatrix[0][1] - rotationMatrix[1][0], 2) +
                        pow(rotationMatrix[2][0] - rotationMatrix[0][2], 2) +
                        pow(rotationMatrix[1][2] - rotationMatrix[2][1], 2)) /
                       4;

        phi = cphi < 0 ? OEL_PI - asin(sqrt(ssphi)) : asin(sqrt(ssphi));
    }

    //  Compute Euler axis
    double sinPhi2x = sin(phi) * 2;
    eulerAxis[0] = (rotationMatrix[2][1] - rotationMatrix[1][2]) / sinPhi2x;
    eulerAxis[1] = (rotationMatrix[0][2] - rotationMatrix[2][0]) / sinPhi2x;
    eulerAxis[2] = (rotationMatrix[1][0] - rotationMatrix[0][1]) / sinPhi2x;

    double norm =
        sqrt(eulerAxis[0] * eulerAxis[0] + eulerAxis[1] * eulerAxis[1] + eulerAxis[2] * eulerAxis[2]);

    eulerAxis[0] /= norm;
    eulerAxis[1] /= norm;
    eulerAxis[2] /= norm;

    //  Compute quaternion
    double sinPhiOver2 = sin(phi / 2);
    quaternion[0] = eulerAxis[0] * sinPhiOver2;
    quaternion[1] = eulerAxis[1] * sinPhiOver2;
    quaternion[2] = eulerAxis[2] * sinPhiOver2;
    quaternion[3] = cos(phi / 2);

    return 0;
}

int multiplyQuaternions(double q1[4], float q2[4], double q3[4]) {
    // Compute the product of two quaternions q3 = q1*q2

    q3[0] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
    q3[1] = -q1[0] * q2[2] + q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
    q3[2] = q1[0] * q2[1] - q1[1] * q2[0] + q1[2] * q2[3] + q1[3] * q2[2];
    q3[3] = -q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] + q1[3] * q2[3];

    return 0;
}

int multiplyQuaternions(float q1[4], float q2[4], float q3[4]) {
    // Compute the product of two quaternions q3 = q1*q2

    q3[0] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
    q3[1] = -q1[0] * q2[2] + q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
    q3[2] = q1[0] * q2[1] - q1[1] * q2[0] + q1[2] * q2[3] + q1[3] * q2[2];
    q3[3] = -q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] + q1[3] * q2[3];

    return 0;
}

void interpolateOrbitVectors(size_t numOrbRec, size_t numScans, const vector<double> &orbitTimes,
                             const vector<orbArray> &positions, const vector<orbArray> &velocities,
                             const vector<double> &scanTimes, vector<orbArray> &posI,
                             vector<orbArray> &velI) {
    double constantTerm[3], linearTerm[3], quadraticTerm[3], cubicTerm[3];

    // initialize the arrays
    for (int i = 0; i < 3; i++) {
        constantTerm[i] = 0.0;
        linearTerm[i] = 0.0;
        quadraticTerm[i] = 0.0;
        cubicTerm[i] = 0.0;
    }

    // orbit time
    double startOrbTime = orbitTimes[0];
    double endOrbTime = orbitTimes[numOrbRec - 1];
    const double NAV_FILL = -9999999.0;

    size_t ind = 0;
    for (size_t i = 0; i < numScans; i++) {
        // scantime outside valid orbit time, use fill value
        if (scanTimes[i] < startOrbTime || scanTimes[i] > endOrbTime) {
            for (int vector = 0; vector < 3; vector++) {
                posI[i][vector] = NAV_FILL;
                velI[i][vector] = BAD_FLT;
            }
            // skip calulations below
            continue;
        }

        //  Find input orbit vectors bracketing scan
        for (size_t j = ind; j < numOrbRec; j++) {
            if (scanTimes[i] > orbitTimes[j] && scanTimes[i] <= orbitTimes[j + 1]) {
                ind = j;
                break;
            }
        }

        //  Set up cubic interpolation
        double dt = orbitTimes[ind + 1] - orbitTimes[ind];
        for (size_t j = 0; j < 3; j++) {
            constantTerm[j] = positions[ind][j];
            linearTerm[j] = velocities[ind][j] * dt;
            if (dt >= 3.5) {
                quadraticTerm[j] = 3 * positions[ind + 1][j] - 3 * positions[ind][j] -
                                   2 * velocities[ind][j] * dt - velocities[ind + 1][j] * dt;
                cubicTerm[j] = 2 * positions[ind][j] - 2 * positions[ind + 1][j] + velocities[ind][j] * dt +
                               velocities[ind + 1][j] * dt;
            } else {
                quadraticTerm[j] = (velocities[ind + 1][j] - velocities[ind][j]) / 2;
                cubicTerm[j] = 0.0;
            }
        }

        //  Interpolate orbit position and velocity components to scan line time
        double x = (scanTimes[i] - orbitTimes[ind]) / dt;
        double x2 = x * x;
        double x3 = x2 * x;
        for (size_t j = 0; j < 3; j++) {
            posI[i][j] = constantTerm[j] + linearTerm[j] * x + quadraticTerm[j] * x2 + cubicTerm[j] * x3;
            velI[i][j] = (linearTerm[j] + 2 * quadraticTerm[j] * x + 3 * cubicTerm[j] * x2) / dt;
        }
    }  // i-loop
}

int interpolateAttitudeForScanTimes(size_t numAttRec, size_t numScans, double *quatTimes, quaternion *quats,
                                    double *scanTimes, quaternion *quatsInterpolated) {
    // Attitude
    double startAttTime = quatTimes[0];
    double endAttTime = quatTimes[numAttRec - 1];
    size_t ind = 0;
    for (size_t i = 0; i < numScans; i++) {
        // if scan time is outside the valid att time, use fill value
        if (scanTimes[i] < startAttTime || scanTimes[i] > endAttTime) {
            for (int element = 0; element < 4; element++) {
                quatsInterpolated[i][element] = BAD_FLT;
            }
            continue;
        }

        //  Find input attitude vectors bracketing scan
        for (size_t j = ind; j < numAttRec; j++) {
            if (scanTimes[i] > quatTimes[j] && scanTimes[i] <= quatTimes[j + 1]) {
                ind = j;
                break;
            }
        }

        //  Set up quaternion interpolation
        double deltaT = quatTimes[ind + 1] - quatTimes[ind];
        double quatToInterpolate[4];
        quatToInterpolate[0] = -quats[ind][0];
        quatToInterpolate[1] = -quats[ind][1];
        quatToInterpolate[2] = -quats[ind][2];
        quatToInterpolate[3] = quats[ind][3];

        double relativeRotationAxis[3], resultQuat[4];
        multiplyQuaternions(quatToInterpolate, quats[ind + 1], resultQuat);
        memcpy(relativeRotationAxis, resultQuat, 3 * sizeof(double));
        double rotationMagnitude = sqrt(relativeRotationAxis[0] * relativeRotationAxis[0] +
                                        relativeRotationAxis[1] * relativeRotationAxis[1] +
                                        relativeRotationAxis[2] * relativeRotationAxis[2]);
        for (size_t j = 0; j < 3; j++)
            relativeRotationAxis[j] /= rotationMagnitude;

        // Interpolate quaternion to scan times
        double interpolationFactor = (scanTimes[i] - quatTimes[ind]) / deltaT;
        float relativeRotationQuat[4], interpolatedQuat[4];
        for (size_t j = 0; j < 3; j++)
            relativeRotationQuat[j] = relativeRotationAxis[j] * rotationMagnitude * interpolationFactor;
        relativeRotationQuat[3] = 1.0;
        multiplyQuaternions(quats[ind], relativeRotationQuat, interpolatedQuat);
        memcpy(quatsInterpolated[i], interpolatedQuat, 4 * sizeof(float));
    }

    return 0;
}

int interpolateTilt(size_t numTilts, size_t numScans, double *tiltTimes, float *tiltin, double *scanTimes,
                    float *tiltInterpolated) {
    double startTiltTime = tiltTimes[0];
    double endTiltTime = tiltTimes[numTilts - 1];

    size_t ind = 0;
    for (size_t i = 0; i < numScans; i++) {
        if (scanTimes[i] < startTiltTime || scanTimes[i] > endTiltTime) {
            tiltInterpolated[i] = BAD_FLT;
            continue;
        }

        //  Find input tilt vectors bracketing scan
        for (size_t j = ind; j < numTilts; j++) {
            if (scanTimes[i] > tiltTimes[j] && scanTimes[i] <= tiltTimes[j + 1]) {
                ind = j;
                break;
            }
        }

        double interpolationFactor = (scanTimes[i] - tiltTimes[ind]) / (tiltTimes[ind + 1] - tiltTimes[ind]);
        tiltInterpolated[i] = (1 - interpolationFactor) * tiltin[ind] + interpolationFactor * tiltin[ind + 1];
    }

    return 0;
}

double calculateGreenwichHourAngle(int32_t year, int32_t day, double seconds) {
    double timeOfDay = day + seconds / SECONDS_IN_DAY;
    double greenwichHourAngle;
    gha2000(year, timeOfDay, &greenwichHourAngle);
    return greenwichHourAngle * OEL_RADEG;
}

int getEcrSunVector(size_t numScans, int32_t year, int32_t day, double *sec, orbArray *sunr, double *rs) {
    //  Get unit Sun vector in geocentric inertial coordinates
    getInertialSunVector(numScans, year, day, sec, sunr, rs);

    //  Get Greenwich mean sidereal angle
    for (size_t i = 0; i < numScans; i++) {
        double greenwichHourAngleDegrees = calculateGreenwichHourAngle(year, day, sec[i]);

        //  Transform Sun vector into geocentric rotating frame
        sunr[i][0] =
            sunr[i][0] * cos(greenwichHourAngleDegrees) + sunr[i][1] * sin(greenwichHourAngleDegrees);
        sunr[i][1] =
            sunr[i][1] * cos(greenwichHourAngleDegrees) - sunr[i][0] * sin(greenwichHourAngleDegrees);
    }

    return EXIT_SUCCESS;
}

int getInertialSunVector(size_t numScans, int32_t year, int32_t dayOfYear, double *secondsOfDay,
                         orbArray *sunVector, double *sunEarthDistance) {
    float constexpr aberrationConstant = 0.0056932;  // Constant of aberration

    for (size_t i = 0; i < numScans; i++) {
        //   Compute floating point days since Jan 1.5, 2000
        //    Note that the Julian day starts at noon on the specified date
        double daysSinceJ2000 = jday((int16_t)year, 1, (int16_t)dayOfYear) - J_2000 +
                                (secondsOfDay[i] - (SECONDS_IN_DAY / 2)) / SECONDS_IN_DAY;
        double meanSolarLon, meanAnomalySun, meanLunarLon, ascendingNodeLon;
        double nutationInLongitude, obliquityOfEcliptic;

        //  Compute solar ephemeris parameters
        ephparms(daysSinceJ2000, &meanSolarLon, &meanAnomalySun, &meanLunarLon, &ascendingNodeLon);

        // Compute nutation corrections for this day
        nutate(daysSinceJ2000, meanSolarLon, meanAnomalySun, meanLunarLon, ascendingNodeLon,
               &nutationInLongitude, &obliquityOfEcliptic);

        //  Compute planet mean anomalies. These are the points at which the planets find themselves relative
        //  to their respective perihelions.
        //   Venus Mean Anomaly
        double venusMeanAnomaly = 50.40828 + 1.60213022 * daysSinceJ2000;
        venusMeanAnomaly = fmod(venusMeanAnomaly, (double)360);

        //   Mars Mean Anomaly
        double marsMeanAnomaly = 19.38816 + 0.52402078 * daysSinceJ2000;
        marsMeanAnomaly = fmod(marsMeanAnomaly, (double)360);

        //  Jupiter Mean Anomaly
        double jupiterMeanAnomaly = 20.35116 + 0.08309121 * daysSinceJ2000;
        jupiterMeanAnomaly = fmod(jupiterMeanAnomaly, (double)360);

        //  Compute solar distance (AU)
        sunEarthDistance[i] = 1.00014 - 0.01671 * cos(meanAnomalySun / OEL_RADEG) -
                              0.00014 * cos(2. * meanAnomalySun / OEL_RADEG);

        //  Compute Geometric Solar Longitude
        double geometricSolarLonCorrection =
            (6893. - 4.6543463e-4 * daysSinceJ2000) * sin(meanAnomalySun / OEL_RADEG) +
            72. * sin(2. * meanAnomalySun / OEL_RADEG) -
            7. * cos((meanAnomalySun - jupiterMeanAnomaly) / OEL_RADEG) +
            6. * sin((meanLunarLon - meanSolarLon) / OEL_RADEG) +
            5. * sin((4. * meanAnomalySun - 8. * marsMeanAnomaly + 3. * jupiterMeanAnomaly) / OEL_RADEG) -
            5. * cos((2. * meanAnomalySun - 2. * venusMeanAnomaly) / OEL_RADEG) -
            4. * sin((meanAnomalySun - venusMeanAnomaly) / OEL_RADEG) +
            4. * cos((4. * meanAnomalySun - 8. * marsMeanAnomaly + 3. * jupiterMeanAnomaly) / OEL_RADEG) +
            3. * sin((2. * meanAnomalySun - 2. * venusMeanAnomaly) / OEL_RADEG) -
            3. * sin(jupiterMeanAnomaly / OEL_RADEG) -
            3. * sin((2. * meanAnomalySun - 2. * jupiterMeanAnomaly) / OEL_RADEG);

        double geometricSolarLon = meanSolarLon + geometricSolarLonCorrection / 3600;

        //  Compute Apparent Solar Longitude// includes corrections for nutation
        //  in longitude and velocity aberration
        double apparentSolarLon =
            geometricSolarLon + nutationInLongitude - aberrationConstant / sunEarthDistance[i];

        //   Compute unit Sun vector
        sunVector[i][0] = cos(apparentSolarLon / OEL_RADEG);
        sunVector[i][1] = sin(apparentSolarLon / OEL_RADEG) * cos(obliquityOfEcliptic / OEL_RADEG);
        sunVector[i][2] = sin(apparentSolarLon / OEL_RADEG) * sin(obliquityOfEcliptic / OEL_RADEG);
    }  // i-loop

    return 0;
}

int quaternionToMatrix(float quaternion[4], double rotationMatrix[3][3]) {
    // Convert quaternion to equivalent direction cosine matrix

    rotationMatrix[0][0] = quaternion[0] * quaternion[0] - quaternion[1] * quaternion[1] -
                           quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3];
    rotationMatrix[0][1] = 2 * (quaternion[0] * quaternion[1] + quaternion[2] * quaternion[3]);
    rotationMatrix[0][2] = 2 * (quaternion[0] * quaternion[2] - quaternion[1] * quaternion[3]);

    rotationMatrix[1][0] = 2 * (quaternion[0] * quaternion[1] - quaternion[2] * quaternion[3]);
    rotationMatrix[1][1] = -quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1] -
                           quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3];
    rotationMatrix[1][2] = 2 * (quaternion[1] * quaternion[2] + quaternion[0] * quaternion[3]);

    rotationMatrix[2][0] = 2 * (quaternion[0] * quaternion[2] + quaternion[1] * quaternion[3]);
    rotationMatrix[2][1] = 2 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3]);
    rotationMatrix[2][2] = -quaternion[0] * quaternion[0] - quaternion[1] * quaternion[1] +
                           quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3];

    return EXIT_SUCCESS;
}

int getEllipsoidScanTrackCoefs(orbArray craftPos, double sensorOrientMatrix[3][3],
                               LocatingContext locatingContext, double scanTrackCoefs[10]) {
    //  Compute constants for navigation model using Earth radius values
    double reciprocalFlatDiff = 1 / EARTH_OMF2;
    double radius = EARTH_RADIUS;

    if (locatingContext == SELENO) {
        reciprocalFlatDiff = 1 / MOON_OMF2;
        radius = MOON_RADIUS;
    }

    //  Compute coefficients of intersection ellipse in scan plane
    scanTrackCoefs[0] = 1 + (reciprocalFlatDiff - 1) * sensorOrientMatrix[0][2] * sensorOrientMatrix[0][2];
    scanTrackCoefs[1] = 1 + (reciprocalFlatDiff - 1) * sensorOrientMatrix[1][2] * sensorOrientMatrix[1][2];
    scanTrackCoefs[2] = 1 + (reciprocalFlatDiff - 1) * sensorOrientMatrix[2][2] * sensorOrientMatrix[2][2];
    scanTrackCoefs[3] = (reciprocalFlatDiff - 1) * sensorOrientMatrix[0][2] * sensorOrientMatrix[1][2] * 2;
    scanTrackCoefs[4] = (reciprocalFlatDiff - 1) * sensorOrientMatrix[0][2] * sensorOrientMatrix[2][2] * 2;
    scanTrackCoefs[5] = (reciprocalFlatDiff - 1) * sensorOrientMatrix[1][2] * sensorOrientMatrix[2][2] * 2;
    scanTrackCoefs[6] = (sensorOrientMatrix[0][0] * craftPos[0] + sensorOrientMatrix[0][1] * craftPos[1] +
                         sensorOrientMatrix[0][2] * craftPos[2] * reciprocalFlatDiff) *
                        2;
    scanTrackCoefs[7] = (sensorOrientMatrix[1][0] * craftPos[0] + sensorOrientMatrix[1][1] * craftPos[1] +
                         sensorOrientMatrix[1][2] * craftPos[2] * reciprocalFlatDiff) *
                        2;
    scanTrackCoefs[8] = (sensorOrientMatrix[2][0] * craftPos[0] + sensorOrientMatrix[2][1] * craftPos[1] +
                         sensorOrientMatrix[2][2] * craftPos[2] * reciprocalFlatDiff) *
                        2;
    scanTrackCoefs[9] = craftPos[0] * craftPos[0] + craftPos[1] * craftPos[1] +
                        craftPos[2] * craftPos[2] * reciprocalFlatDiff - radius * radius;

    return 0;
}

void convertVelVecToGroundSpeed(float position[3], float velocity[3], double sensorOrientationMatrix[3][3],
                                double pixelVector[3], double *deltaT, float geoPosition[3],
                                float sensorToPixelVector[3], gsl_vector *C, size_t pixel) {
    double positionMagnitude =
        sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]);
    double geocentricLatitudeCosine =
        sqrt(position[0] * position[0] + position[1] * position[1]) / positionMagnitude;
    double localEarthRadius =
        EARTH_RADIUS * (1. - EARTH_OMF2) /
        sqrt(1. - (2. - EARTH_OMF2) * EARTH_OMF2 * geocentricLatitudeCosine * geocentricLatitudeCosine);

    double groundVelocity[3];
    groundVelocity[0] = velocity[0] * localEarthRadius / positionMagnitude;
    groundVelocity[1] = velocity[1] * localEarthRadius / positionMagnitude;
    groundVelocity[2] = velocity[2] * localEarthRadius / positionMagnitude;
    //  Transform vector from sensor to geocentric frame
    gsl_matrix_view transformMatrix = gsl_matrix_view_array((double *)sensorOrientationMatrix, 3, 3);
    gsl_vector_view pixelVectorView = gsl_vector_view_array(pixelVector, 3);

    gsl_blas_dgemv(CblasTrans, 1.0, &transformMatrix.matrix, &pixelVectorView.vector, 0.0, C);

    double *transformedVector = gsl_vector_ptr(C, 0);
    for (size_t j = 0; j < 3; j++) {
        sensorToPixelVector[j] = transformedVector[j];
        geoPosition[j] = position[j] + sensorToPixelVector[j] + groundVelocity[j] * deltaT[pixel];
    }
}

void calculateLocalCoordinateSystem(const float geoPosition[3], float localVerticalVector[3],
                                    float eastVector[3], float northVector[3], float &horizontalDistance) {
    horizontalDistance = geoPosition[0] * geoPosition[0] + geoPosition[1] * geoPosition[1];
    float verticalDistance =
        sqrt(geoPosition[2] * geoPosition[2] + EARTH_OMF2 * EARTH_OMF2 * horizontalDistance);

    localVerticalVector[0] = EARTH_OMF2 * geoPosition[0] / verticalDistance;
    localVerticalVector[1] = EARTH_OMF2 * geoPosition[1] / verticalDistance;
    localVerticalVector[2] = geoPosition[2] / verticalDistance;
    float horizontalProjectionLength = sqrt(localVerticalVector[0] * localVerticalVector[0] +
                                            localVerticalVector[1] * localVerticalVector[1]);

    eastVector[0] = -localVerticalVector[1] / horizontalProjectionLength;
    eastVector[1] = localVerticalVector[0] / horizontalProjectionLength;
    eastVector[2] = 0.0;

    // no = crossp(up,ea)

    northVector[0] = -localVerticalVector[2] * eastVector[1];
    northVector[1] = localVerticalVector[2] * eastVector[0];
    northVector[2] = localVerticalVector[0] * eastVector[1] - localVerticalVector[1] * eastVector[0];
}

void computeSolarZeniths(float eastVector[3], float northVector[3], float localVerticalVector[3],
                         float sensorToPixelVector[3], float sunUnitVector[3], float pixelToSpacecraftVec[3],
                         float sunToEarthVector[3], short *solarZeniths, int i) {
    // Transform the pixel-to-spacecraft and Sun vectors into local frame
    pixelToSpacecraftVec[0] = -eastVector[0] * sensorToPixelVector[0] -
                              eastVector[1] * sensorToPixelVector[1] - eastVector[2] * sensorToPixelVector[2];
    pixelToSpacecraftVec[1] = -northVector[0] * sensorToPixelVector[0] -
                              northVector[1] * sensorToPixelVector[1] -
                              northVector[2] * sensorToPixelVector[2];
    pixelToSpacecraftVec[2] = -localVerticalVector[0] * sensorToPixelVector[0] -
                              localVerticalVector[1] * sensorToPixelVector[1] -
                              localVerticalVector[2] * sensorToPixelVector[2];

    sunToEarthVector[0] = sunUnitVector[0] * eastVector[0] + sunUnitVector[1] * eastVector[1] +
                          sunUnitVector[2] * eastVector[2];
    sunToEarthVector[1] = sunUnitVector[0] * northVector[0] + sunUnitVector[1] * northVector[1] +
                          sunUnitVector[2] * northVector[2];
    sunToEarthVector[2] = sunUnitVector[0] * localVerticalVector[0] +
                          sunUnitVector[1] * localVerticalVector[1] +
                          sunUnitVector[2] * localVerticalVector[2];

    //  Compute the solar zenith and azimuth
    solarZeniths[i] = (short)(100 * OEL_RADEG *
                              atan2(sqrt(sunToEarthVector[0] * sunToEarthVector[0] +
                                         sunToEarthVector[1] * sunToEarthVector[1]),
                                    sunToEarthVector[2]));
}

int locatePixelsOci(const char *demFile, orbArray pos, orbArray vel, double smat[3][3],
                    double scanPathCoefs[10], LocatingContext locatingContext, orbArray sunUnitVector,
                    bool fileIsEclipsed, orbArray moonVector, double earthSunDist,
                    vector<array<float, 3>> &sensorViews, size_t numPix, double *delT,
                    vector<uint8_t> &qualityFlags, size_t currScan, float *latitudes, float *longitudes,
                    float *solarZeniths, float *solarAzimuths, float *sensorZeniths, float *sensorAzimuths,
                    short *terrainHeights) {
    // Ellipsoid parameters
    float FLAT;  // How squished the ellipsoid is

    switch (locatingContext) {
        case SELENO:
            FLAT = MOON_FLATTENING;
        case HELIO:  // There's probably something useful to do for HELIO
        case GEO:
        default:
            FLAT = 1 / 298.257;
    }
    const float ECCENTRICITY = (1 - FLAT) * (1 - FLAT);  // AKA omf2

    gsl_vector *C = gsl_vector_alloc(3);

    for (size_t i = 0; i < numPix; i++) {
        // Compute sensor-to-surface vectors for all scan angles
        // Compute terms for quadratic equation
        double o = scanPathCoefs[0] * sensorViews[i][0] * sensorViews[i][0] +
                   scanPathCoefs[1] * sensorViews[i][1] * sensorViews[i][1] +
                   scanPathCoefs[2] * sensorViews[i][2] * sensorViews[i][2] +
                   scanPathCoefs[3] * sensorViews[i][0] * sensorViews[i][1] +
                   scanPathCoefs[4] * sensorViews[i][0] * sensorViews[i][2] +
                   scanPathCoefs[5] * sensorViews[i][1] * sensorViews[i][2];
        double p = scanPathCoefs[6] * sensorViews[i][0] + scanPathCoefs[7] * sensorViews[i][1] +
                   scanPathCoefs[8] * sensorViews[i][2];
        double q = scanPathCoefs[9];
        double discriminant = p * p - 4 * q * o;  // Line from sensor to pixel

        // If discriminant < 0, this pixel views above the horizon. If o == 0, there would be a divide by 0
        if (discriminant < 0 || o == 0) {
            latitudes[i] = BAD_FLT;
            longitudes[i] = BAD_FLT;
            solarZeniths[i] = BAD_FLT;
            solarAzimuths[i] = BAD_FLT;
            sensorZeniths[i] = BAD_FLT;
            sensorAzimuths[i] = BAD_FLT;
            terrainHeights[i] = BAD_INT;
            continue;
        }

        //  Solve for magnitude of sensor-to-pixel vector and compute components
        double d = (-p - sqrt(discriminant)) / (2 * o);
        double x1[3];
        for (size_t j = 0; j < 3; j++)
            x1[j] = d * sensorViews[i][j];

        //  Convert velocity vector to ground speed
        double v[3];
        float re = 6378.137;
        double pm = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
        double clatg = sqrt(pos[0] * pos[0] + pos[1] * pos[1]) / pm;
        double rg = re * (1. - FLAT) / sqrt(1. - (2. - FLAT) * FLAT * clatg * clatg);
        if (locatingContext == GEO) {
            v[0] = vel[0] * rg / pm;
            v[1] = vel[1] * rg / pm;
            v[2] = vel[2] * rg / pm;
        } else {
            v[0] = vel[0];
            v[1] = vel[1];
            v[2] = vel[2];
            double rgOverPm = rg / pm;
            scaleElements(v, rgOverPm);
        }

        //  Transform vector from sensor to geocentric frame
        gsl_matrix_view A = gsl_matrix_view_array((double *)smat, 3, 3);
        gsl_vector_view B = gsl_vector_view_array(x1, 3);

        gsl_blas_dgemv(CblasTrans, 1.0, &A.matrix, &B.vector, 0.0, C);

        float rh[3], geovec[3];
        double *ptr_C = gsl_vector_ptr(C, 0);
        for (size_t j = 0; j < 3; j++) {
            rh[j] = ptr_C[j];
            if (locatingContext == GEO)
                geovec[j] = pos[j] + rh[j] + v[j] * delT[i];
            else if (locatingContext == SELENO)
                geovec[j] = pos[j] + rh[j];
        }

        // Compute the local vertical, East and North unit vectors
        float uxy = geovec[0] * geovec[0] + geovec[1] * geovec[1];
        float temp = sqrt(geovec[2] * geovec[2] + ECCENTRICITY * ECCENTRICITY * uxy);

        float up[3];
        up[0] = ECCENTRICITY * geovec[0] / temp;
        up[1] = ECCENTRICITY * geovec[1] / temp;
        up[2] = geovec[2] / temp;
        float upxy = sqrt(up[0] * up[0] + up[1] * up[1]);

        float ea[3];
        ea[0] = -up[1] / upxy;
        ea[1] = up[0] / upxy;
        ea[2] = 0.0;

        // no = crossp(up,ea)
        float no[3];
        no[0] = -up[2] * ea[1];
        no[1] = up[2] * ea[0];
        no[2] = up[0] * ea[1] - up[1] * ea[0];

        //  Compute geodetic latitude and longitude
        float xlat_ = OEL_RADEG * asin(up[2]);
        float xlon_ = OEL_RADEG * atan2(up[1], up[0]);
        // check if it is in the range
        if (std::abs(xlat_) > 90.0e0 || std::abs(xlon_) > 180.0e0) {
            continue;
        }

        // Transform the pixel-to-spacecraft and Sun vectors into local frame
        float rl[3], sl[3];
        rl[0] = -ea[0] * rh[0] - ea[1] * rh[1] - ea[2] * rh[2];
        rl[1] = -no[0] * rh[0] - no[1] * rh[1] - no[2] * rh[2];
        rl[2] = -up[0] * rh[0] - up[1] * rh[1] - up[2] * rh[2];

        sl[0] = sunUnitVector[0] * ea[0] + sunUnitVector[1] * ea[1] + sunUnitVector[2] * ea[2];
        sl[1] = sunUnitVector[0] * no[0] + sunUnitVector[1] * no[1] + sunUnitVector[2] * no[2];
        sl[2] = sunUnitVector[0] * up[0] + sunUnitVector[1] * up[1] + sunUnitVector[2] * up[2];

        //  Compute the solar zenith and azimuth
        solarZeniths[i] = OEL_RADEG * atan2(sqrt(sl[0] * sl[0] + sl[1] * sl[1]), sl[2]);

        // Check for zenith close to zero
        if (solarZeniths[i] > 0.01)  // As a zenith approaches 0, azimuth becomes meaningless
            solarAzimuths[i] = OEL_RADEG * atan2(sl[0], sl[1]);
        else
            solarAzimuths[i] = 0.0;

        // Compute the sensor zenith and azimuth
        sensorZeniths[i] = OEL_RADEG * atan2(sqrt(rl[0] * rl[0] + rl[1] * rl[1]), rl[2]);

        // Check for zenith close to zero
        if (sensorZeniths[i] > 0.01)  // As a zenith approaches 0, azimuth becomes meaningless
            sensorAzimuths[i] = OEL_RADEG * atan2(rl[0], rl[1]);
        else
            sensorAzimuths[i] = 0.0;

        if (fileIsEclipsed) {
            // Reference Wertz, App. l, Table L-4
            double sunAngularRadius = 0.53313 / earthSunDist / OEL_RADEG / 2;
            double moonDist;

            orbArray earthToMoon;
            for (size_t j = 0; j < 2; j++) {
                earthToMoon[j] = moonVector[j] - geovec[j];
            }
            moonDist = sqrt(pow(earthToMoon[0], 2.0) + pow(earthToMoon[1], 2.0) + pow(earthToMoon[2], 2.0));

            orbArray sunMoonDiff;
            for (size_t j = 0; j < 2; j++) {
                sunMoonDiff[j] = earthToMoon[j] / moonDist - sunUnitVector[j];
            }
            // The angle described between the Sun and the Moon with Earth as the apex
            double sunEarthMoon =
                sqrt(pow(sunMoonDiff[0], 2.0) + pow(sunMoonDiff[1], 2.0) + pow(sunMoonDiff[2], 2.0));

            double moonAngularRadius = 1738.2 / moonDist;

            // Set eclipse flag where angle is less than combined Sun and Moon radii
            if (sunEarthMoon < (sunAngularRadius + moonAngularRadius)) {
                qualityFlags[currScan * numPix + i] |= SOLAR_ECLIPSE;
            }
        }

        float tempHeight = 0;
        get_nc_height(demFile, &xlon_, &xlat_, &sensorZeniths[i], &sensorAzimuths[i], &tempHeight);

        // there is also a bug in terrain correction
        bool badLatLons = abs(xlat_) > 90.0 || abs(xlon_) > 180.0;
        if (badLatLons) {
            latitudes[i] = BAD_FLT;
            longitudes[i] = BAD_FLT;
            solarZeniths[i] = BAD_FLT;
            solarAzimuths[i] = BAD_FLT;
            sensorZeniths[i] = BAD_FLT;
            sensorAzimuths[i] = BAD_FLT;
            continue;
        } else {
            latitudes[i] = xlat_;
            longitudes[i] = xlon_;
        }

        terrainHeights[i] = (short)tempHeight;

    }  // pixel loop

    gsl_vector_free(C);

    return EXIT_SUCCESS;
}

int getAttitudeAngles(orbArray pos, orbArray vel, double smat[3][3], orbArray &attitudeAngles) {
    double rem = 6371;               // TODO: This is suspiciously close to Earth radius in meters...
    double f = 1 / (double)298.257;  // TODO: This is FLAT
    double omf2 = (1 - f) * (1 - f);

    // Determine local vertical reference axes
    double p[3], v[3];
    for (size_t j = 0; j < 3; j++) {
        p[j] = (double)pos[j];
        v[j] = (double)vel[j];
    }
    v[0] -= p[1] * omegaE;
    v[1] += p[0] * omegaE;

    //  Compute Z axis as local nadir vector
    double pm = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    double omf2p = (omf2 * rem + pm - rem) / pm;
    double pxy = p[0] * p[0] + p[1] * p[1];
    double temp = sqrt(p[2] * p[2] + omf2p * omf2p * pxy);

    double z[3];
    z[0] = -omf2p * p[0] / temp;
    z[1] = -omf2p * p[1] / temp;
    z[2] = -p[2] / temp;

    // Compute Y axis along negative orbit normal
    double on[3];
    on[0] = v[1] * z[2] - v[2] * z[1];
    on[1] = v[2] * z[0] - v[0] * z[2];
    on[2] = v[0] * z[1] - v[1] * z[0];

    double onm = sqrt(on[0] * on[0] + on[1] * on[1] + on[2] * on[2]);

    double y[3];
    for (size_t j = 0; j < 3; j++)
        y[j] = -on[j] / onm;

    // Compute X axis to complete orthonormal triad (velocity direction)
    double x[3];
    x[0] = y[1] * z[2] - y[2] * z[1];
    x[1] = y[2] * z[0] - y[0] * z[2];
    x[2] = y[0] * z[1] - y[1] * z[0];

    // Store local vertical reference vectors in matrix
    double om[3][3];
    memcpy(&om[0][0], &x, 3 * sizeof(double));
    memcpy(&om[1][0], &y, 3 * sizeof(double));
    memcpy(&om[2][0], &z, 3 * sizeof(double));

    // Compute orbital-to-spacecraft matrix
    double rm[3][3];
    gsl_matrix_view rm_mat = gsl_matrix_view_array(&rm[0][0], 3, 3);

    int s;  // Only used to fill a parameter

    gsl_permutation *perm = gsl_permutation_alloc(3);

    // Compute the LU decomposition of this matrix
    gsl_matrix_view B = gsl_matrix_view_array(&om[0][0], 3, 3);
    gsl_linalg_LU_decomp(&B.matrix, perm, &s);

    // Compute the  inverse of the LU decomposition
    double inv[3][3];
    gsl_matrix_view inv_mat = gsl_matrix_view_array(&inv[0][0], 3, 3);

    gsl_linalg_LU_invert(&B.matrix, perm, &inv_mat.matrix);

    gsl_matrix_view A = gsl_matrix_view_array(&smat[0][0], 3, 3);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &inv_mat.matrix, 0.0, &rm_mat.matrix);

    gsl_permutation_free(perm);

    // Compute attitude angles
    attitudeAngles[0] = OEL_RADEG * atan(-rm[2][1] / rm[2][2]);
    double cosp = sqrt(rm[2][1] * rm[2][1] + rm[2][2] * rm[2][2]);
    if (rm[2][2] < 0)
        cosp *= -1;
    attitudeAngles[1] = OEL_RADEG * atan2(rm[2][0], cosp);
    attitudeAngles[2] = OEL_RADEG * atan(-rm[1][0] / rm[0][0]);

    return 0;
}
