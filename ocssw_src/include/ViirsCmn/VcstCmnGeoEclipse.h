/*****************************************************************************
 *
 * NAME:  VcstCmnGeoEclipse
 *
 * DESCRIPTION: Common Geolocation Eclipse object class.
 *
 * Adapted directly from ProSdrCmnGeoEclipse.h developed by Raytheon Company.
 *
 * REFERENCES:
 * Meeus, Jean, "Astronomical Algorithms, 2nd Edition," Willman-Bell Inc,
 * Richmond VA, 1998, 477 pp.
 *
 *
 *****************************************************************************/

#ifndef VcstCmnGeoEclipse_h
#define VcstCmnGeoEclipse_h

#include <iomanip>
#include <iostream>

/**
 * VcstCmnGeoEclipse serves as a cache for eclipse parameters computed
 * by VcstCmnGeo.  It stores the key parameters for a given eclipse,
 * including the lunar cycle number, the time of the full or new moon,
 * whether an eclipse is occuring at the particular lunar phase and if
 * so, the begin and end times of the (penumbral) eclipse.
 */
class VcstCmnGeoEclipse {
public:

    /**
     * Default constructor.
     */
    VcstCmnGeoEclipse();

    /**
     * Copy Constructor
     */
    VcstCmnGeoEclipse(const VcstCmnGeoEclipse &data);

    /**
     * Assignment operator
     */
    VcstCmnGeoEclipse& operator=(const VcstCmnGeoEclipse &data);

    /**
     * Destructor
     */
    virtual ~VcstCmnGeoEclipse();

    /**
     * Calculate lunar eclipse parameters for the Full Moon specified by
     * Meeus' k-factor.  Upon successful return, the eclipse paramaters
     * for this object will be set to the eclipse parameters for the given
     * moon cycle/phase and may be retrieved via the 'getParams' method.
     *
     * @param kFactor  input "k" factor, which specifies the lunar phase
     *                 cycle since the New Moon of January 6, 2000.  For
     *                 Full Moons, k is an integer increased by 0.5.
     *                 See Meeus pgs. 349 and 350.
     */
    void calcLunarEclipse(const double kFactor);

    /**
     * Calculate solar eclipse parameters for the New Moon specified by
     * Meeus' k-factor.  Upon successful return, the eclipse paramaters
     * for this object will be set to the eclipse parameters for the given
     * moon cycle/phase and may be retrieved via the 'getParams' method.
     *
     * @param kFactor  input "k" factor, which specifies the lunar phase
     *                 cycle since the New Moon of January 6, 2000.  For
     *                 New Moons, k is a integer value.  See Meeus pgs.
     *                 349 and 350.
     */
    void calcSolarEclipse(const double kFactor);

    /**
     * Get the time of the eclipse data, in TAI. If an eclipse is occurring,
     * this time will be the time of maximum eclipse; if no eclipse is
     * occurring, the time will be the time of the mean New or Full Moon.
     */
    /* inline */
    double getTAI() const;

    /**
     * true indicates eclipse is occurring for this new/full moon,
     * false indicates it is not.
     */
    bool getEclipseFlag() const;

    /**
     * TAI time of the beginning of eclipse, if one is occurring.
     */
    double getEclipseBegin() const;

    /**
     * TAI time of the end of eclipse, if one is occurring.
     */
    double getEclipseEnd() const;

    /**
     * Give output operator access to member data.
     */
    friend std::ostream& operator<<(std::ostream &strm,
            const VcstCmnGeoEclipse &eclipse);
protected:

    /**
     * Get the time in Julian centuries since the epoch 2000.  Meeus
     * refers to this parameter as "T".  See eq (49.3), pg 350.
     *
     * @param kFactor  input "k" factor, which specifies the lunar phase
     *                 cycle since the New Moon of January 6, 2000.
     *                 See Meeus, pg 349-350.
     * @return         time in Julian centuries since epoch 2000.
     */
    double calcT(const double kFactor);

    /**
     * Get the time of the mean phase of the Moon (in TDT Julian Days)
     * for the specific lunar cycle and phase represented by the
     * previously computed "k" (lunar cycles) and "T" (Julian centuries)
     * factors.
     * See Meeus, eq (49.1), pg. 349.
     *
     * @return         time of the mean phase of the moon in TDT Julian days
     */
    double calcTjdPhase();

    /**
     * Get the Moon's argument of latitude for the time represented by the
     * previously computed "k" (lunar cycles) and "T" (Julian centuries)
     * factors.  See Meeus parameter "F", eq (49.6), pg. 350.
     *
     * @return         the Moon's Argument of Latitude, in radians and
     *                 normalized with the range 0 to 2 PI.
     */
    double calcMoonArgLat();

    /**
     * Get the Sun's mean anomaly for the time represented by the
     * previously computed "k" (lunar cycles) and "T" (Julian centuries)
     * factors.  See Meeus parameter "M", eq (49.4), pg. 350.
     *
     * @return         the Sun's mean anomaly, in radians and
     *                 normalized with the range 0 to 2 PI.
     */
    double calcSunMeanAnom();

    /**
     * Get the Moon's mean anomaly for the time represented by the
     * previously computed "k" (lunar cycles) and "T" (Julian centuries)
     * factors.  See Meeus parameter "M'", eq (49.5), pg. 350.
     *
     * @return         the Moon's mean anomaly, in radians and
     *                 normalized with the range 0 to 2 PI.
     */
    double calcMoonMeanAnom();

    /**
     * Get the longitude of the ascending node of the lunar orbit for the
     * time represented by the previously computed "k" (lunar cycles) and
     * "T" (Julian centuries) factors.
     * See Meeus parameter "Omega", eq (49.7), pg. 350.
     *
     * @return         the Moon's longitude of ascending node, in radians
     *                 and normalized within the range 0 to 2 PI.
     */
    double calcMoonLongAscNode();

    /**
     * Calculate the eccentricity of the Earth's orbit around the Sun for the
     * time represented by the previously computed "k" (lunar cycles) and
     * "T" (Julian centuries) factors.  See Meeus eq. (47.6), pg. 338.
     *
     * @return   Eccentricity of the Earth's orbit
     */
    double calcEarthEccen();

    /**
     * Calculates the Meeus parameter F1, which is used to
     * calculate the time of maximum eclipse.  See Meeus, pg. 380.
     *
     * @param moonArgLat the Moon's Argument of Latitude, in radians and
     *                   normalized with the range 0 to 2 PI.
     * @param moonLAN    the Moon's longitude of ascending node, in radians
     *                   and normalized within the range 0 to 2 PI.
     *
     * @return "F1" factor, in radians
     */
    double calcF1(const double moonArgLat, const double moonLAN);

    /**
     * Calculates the Meeus parameter A1, which is used to
     * calculate the time of maximum eclipse.  See Meeus, pg. 380.
     *
     * @return "A1" factor, in radians
     */
    double calcA1();

    /**
     * Calculate Meeus' "P" factor; see Meeus, pg. 381.
     *
     * @param earthEccen   Eccentricity of the Earth's orbit, "E"
     * @param sunMeanAnom  the Sun's mean anomaly, "M", in radians and
     *                     normalized with the range 0 to 2 PI.
     * @param moonMeanAnom the Moon's mean anomaly, "M'" in radians and
     *                     normalized with the range 0 to 2 PI.
     * @param F1           Meeus' parameter "F1" from pg. 380
     *
     * @return "P", in units of the equatorial radius of the Earth.
     */
    double calcP(const double earthEccen, const double sunMeanAnom,
            const double moonMeanAnom, const double F1);

    /**
     * Calculate Meeus' "Q" factor; see Meeus, pg. 381.
     *
     * @param earthEccen   Eccentricity of the Earth's orbit, "E"
     * @param sunMeanAnom  the Sun's mean anomaly, "M", in radians and
     *                     normalized with the range 0 to 2 PI.
     * @param moonMeanAnom the Moon's mean anomaly, "M'" in radians and
     *                     normalized with the range 0 to 2 PI.
     *
     * @return "Q", in units of the equatorial radius of the Earth.
     */
    double calcQ(const double earthEccen, const double sunMeanAnom,
            const double moonMeanAnom);

    /**
     * Calculate Meeus' "gamma", which has slightly different meanings
     * depending upon whether you're computing solar or lunar eclipses.
     * For solar eclipses, gamma represents the least distance from the
     * axis of the Moon's shadow to the center of the Earth.  For lunar
     * eclipses, gamma represents the least distance from the center of
     * the Moon to the axis of the Earth's shadow.  In both cases, units
     * are in Earth equatorial radii.  Gamma may be positive or negative.
     * See equation for gamma in Meeus, pg. 381.
     *
     * @param P  Meeus' "P" parameter; see Meeus, pg. 381
     * @param Q  Meeus' "Q" parameter; see Meeus, pg. 381
     * @param F1 Meeus' "F1" parameter, see Meeus, pg. 380.
     *
     * @return "gamma", in units of Earth's equatorial radius
     */
    double calcGamma(const double P, const double Q, const double F1);

    /**
     * Calculate Meeus' "u", which denotes the radius of the Moon's
     * umbral cone in the fundamental plane, which is the plane through
     * the center of the Earth and perpendicular to the axis of the Moon's
     * shadow.  See equation for u in Meeus, pg. 381.
     *
     * @param earthEccen   Eccentricity of the Earth's orbit, "E"
     * @param sunMeanAnom  the Sun's mean anomaly, "M", in radians and
     *                     normalized with the range 0 to 2 PI.
     * @param moonMeanAnom the Moon's mean anomaly, "M'" in radians and
     *                     normalized with the range 0 to 2 PI.
     *
     * @return "u", in units of Earth's equatorial radius
     */
    double calcU(const double earthEccen, const double sunMeanAnom,
            const double moonMeanAnom);

    /**
     * Take the given angle, in radians, and normalize so that
     * 0 <= result <= 2PI.
     *
     * @param angle Angle in radians
     * @return normalized angle in radians
     */
    double normalizeAngle(const double angle) const;

    /**
     * Convert the given time from Terrestrial Dynamic Time in Julian days
     * (Meeus' JDE) to TAI.
     *
     * @param tjd Time in TDT Julian Days (JDE)
     * @return TAI equivalent
     */
    double convertTJDtoTAI(const double tjd) const;

    /**
     * Convert the given time from TAI to Terrestrial Dynamic Time in
     * Julian days (Meeus' JDE).
     *
     * @param inTAI Time in TAI
     * @return tjd equivalent
     */
    double convertTAItoTJD(const double inTAI) const;

private:

    /**
     * TAI time of full or new moon for which this eclipse data is valid.
     * If an eclipse is occurring, the time will actually be that of
     * maximum eclipse.
     */
    double taiPhase_;

    /**
     * true indicates eclipse is occurring for this new/full moon, false
     * indicates it is not.
     */
    bool eclipseFlag_;

    /**
     * TAI time of the beginning of eclipse, if one is occurring.
     */
    double eclipseBegin_;

    /**
     * TAI time of end of eclipse, if one is occurring.
     */
    double eclipseEnd_;

    /**
     * "k" factor from Meeus pgs 349-350.  Indicates moon phase and
     * cycle since New Moon of 2000 January 6.
     */
    double kFactor_;

    /**
     * Time of full or new moon, in Terrestrial Dynamic Time Julian Days
     * (which is equivalent to Meeus parameter "JDE", eq. (49.1), pg. 349).
     * If an eclipse is occurring, the time will actually be that of
     * maximum eclipse.
     */
    double tjdPhase_;

    /**
     * Julian centuries since epoch 2000
     */
    double tjdCent_;

    /**
     * (Julian centuries since epoch 2000) ^ 2
     */
    double tjdCentSqd_;

    /**
     * (Julian centuries since epoch 2000) ^ 3
     */
    double tjdCentCube_;

    /**
     * (Julian centuries since epoch 2000) ^ 4
     */
    double tjdCentTo4_;

    //------------------------------------------------------------------------
    // Define the constants used in the Meeus algorithms for calculating
    // eclipses.  These constants come from Meeus (see reference above),
    // Ch. 49, pages 349-350, and Ch. 54, pages 379-388.
    //------------------------------------------------------------------------

    // Eq (49.3), pg 350
    static const double LUNATIONS_PER_CENT;

    // Constants for the calculation of mean lunar phase, Eq (49.1), pg 349
    static const double TJD_PHASE_C1;
    static const double TJD_PHASE_C2;
    static const double TJD_PHASE_C3;

    // Constants for Moon's argument of latitude, Eq (49.6), pg 350
    static const double MOON_ARG_LAT_C1;
    static const double MOON_ARG_LAT_C2;
    static const double MOON_ARG_LAT_C3;
    static const double MOON_ARG_LAT_C4;
    static const double MOON_ARG_LAT_C5;

    // Constants for Sun's mean anomaly, Eq (49.4), pg 350
    static const double SUN_MEAN_ANOM_C1;
    static const double SUN_MEAN_ANOM_C2;
    static const double SUN_MEAN_ANOM_C3;
    static const double SUN_MEAN_ANOM_C4;

    // Constants for Moon's mean anomaly, Eq (49.5), pg 350
    static const double MOON_MEAN_ANOM_C1;
    static const double MOON_MEAN_ANOM_C2;
    static const double MOON_MEAN_ANOM_C3;
    static const double MOON_MEAN_ANOM_C4;
    static const double MOON_MEAN_ANOM_C5;

    // Constants for Moon's longitude of ascending node, Eq (49.7), pg. 350
    static const double MOON_LAN_C1;
    static const double MOON_LAN_C2;
    static const double MOON_LAN_C3;
    static const double MOON_LAN_C4;

    // Critical value of Moon's argument of latitude (F) for eclipse to exist
    // Meeus pg. 380
    static const double CRIT_SIN_F;

    // Constants for Earth's orbital eccentricity, Eq (47.6), pg. 338
    static const double EARTH_ECCEN_C1;
    static const double EARTH_ECCEN_C2;

    // Constants used for calculation of F1 and A1 parameters, pg. 380.
    static const double F1_ADJUST;
    static const double A1_C1;
    static const double A1_C2;
    static const double A1_C3;

    // Constants used in the calculation of "P" and "Q", pg. 381.
    static const double P_C1;
    static const double P_C2;
    static const double P_C3;
    static const double P_C4;
    static const double P_C5;
    static const double P_C6;
    static const double P_C7;
    static const double Q_C1;
    static const double Q_C2;
    static const double Q_C3;
    static const double Q_C4;
    static const double Q_C5;
    static const double Q_C6;

    // Constants used in the calculation of "gamma" and "u", pg. 381
    static const double GAMMA_C1;
    static const double U_C1;
    static const double U_C2;
    static const double U_C3;
    static const double U_C4;
    static const double U_C5;

    // Constants used to calculate magnitude of eclipses
    static const double LUNAR_MAG_C1;
    static const double LUNAR_MAG_C2;
    static const double SOLAR_MAG_C1;

    // Constants used to calculate time of maximum eclipse
    static const double MAX_ECL_C1_LUNAR;
    static const double MAX_ECL_C2_LUNAR;
    static const double MAX_ECL_C1_SOLAR;
    static const double MAX_ECL_C2_SOLAR;
    static const double MAX_ECL_C3;
    static const double MAX_ECL_C4;
    static const double MAX_ECL_C5;
    static const double MAX_ECL_C6;
    static const double MAX_ECL_C7;
    static const double MAX_ECL_C8;
    static const double MAX_ECL_C9;
    static const double MAX_ECL_C10;
    static const double MAX_ECL_C11;
    static const double MAX_ECL_C12;
    static const double MAX_ECL_C13;
    static const double MAX_ECL_C14;
    static const double MAX_ECL_C15;
    static const double MAX_ECL_C16;

    // Constants used to calculate semidurations of eclipse
    static const double LUNAR_SD_C1;
    static const double LUNAR_SD_C2;
    static const double LUNAR_SD_C3;
    static const double SOLAR_WORST_SD;
};

//---------------------------------------------------------------------------
// inline functions
//---------------------------------------------------------------------------

inline double VcstCmnGeoEclipse::getTAI() const {
    return taiPhase_;
}

/**
 * Output operator
 */
inline std::ostream& operator<<(std::ostream &strm,
        const VcstCmnGeoEclipse &eclipse) {
    strm << "k=" << eclipse.kFactor_ << std::endl << std::fixed
            << std::setprecision(5) << "  tjdPhase = " << eclipse.tjdPhase_
            << ", " << "  taiPhase = " << eclipse.taiPhase_ << std::endl
            << "  eclipseFlag = " << eclipse.eclipseFlag_ << std::endl
            << "  begin = " << eclipse.eclipseBegin_ << ", " << "  end = "
            << eclipse.eclipseEnd_ << std::endl;
    return strm;
}

#endif
