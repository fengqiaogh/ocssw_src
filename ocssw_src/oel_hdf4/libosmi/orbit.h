/*
 *----------------------------------------------------------------------
 *  @(#) orbit.h	1.03	5 May 95	<shc>
 *  Copyright (c) 1993, CSIRO Division of Oceanography.
 *----------------------------------------------------------------------
 *
 *  Various structures, constants and macros used in the orbit
 *  package routines.
 *
 * History:
 *	1.0	18 Oct 93	<shc>
 *	    Initial version.
 *
 *	1.01	12 Dec 94	<shc>
 *	    Added errstr stuff for integration with Python
 *
 *	1.02	27 Jan 95	<shc>
 *	    Added SGP4 stuff.
 *
 *	1.03	5 May 95	<shc>
 *	    Added TLE_VALUES struct and prototypes for tle.c
 *
 *	1.04	24 Jan 96	<shc>
 *	    Modified SGP4 struct so secular rates are more accessible,
 *	    changed rates from rad/min to rad/sec for conformance with
 *	    package standards.
 *
 *----------------------------------------------------------------------
 */

#ifndef __ORBIT__
#define __ORBIT__

#include <stdint.h>
#include <math.h>
#include <time.h>
#ifndef NULL
#include <stdio.h>
#endif

#include "earth.h"
#include "lunsol.h"

#define DPI  3.14159265358979323846
#define D2PI  6.28318530717958647692
#define D3PI  9.42477796076937971538

#define FP_ERRVAL -1.0e-20 /* Error value from f.p. routines */

/* Degrees to radians */
#define DTOR(X) ((X) * DPI/180.0)

/* Radians to degrees */
#define RTOD(X) ((X) * 180.0/DPI)

/* Fold an angle to the range [0..2 Pi] radians */
#define AN2PI(X) fmod(fmod((X), D2PI) + D2PI, D2PI)

/* Fold an angle to the range [-Pi .. Pi] radians */
#define ANPI(X) (fmod(fmod((X), D2PI) + D3PI, D2PI) - DPI)

/* Fold an angle to the range [0..360] degrees */
#define AN360(X) fmod(fmod((X), 360.0) + 360.0, 360.0)

/* Arc-seconds to radians */
#define ASTOR(X) ((X) * DPI/(180.0 * 3600.0))

/* Arc-seconds per century to radians per day */
#define ASPCTORPD(X) ((X) * DPI/(180.0 * 3600.0 * 36525.0))

/* Epoch(s) in my time scheme (JD % 2400000) */
#define J2000 (51545.0)

/*
 * Struct for defining the attitude (orientation) of a satellite
 * or satellite scanner vector.  Angles are in radians.
 */

struct ATTITUDE {
    double pitch;
    double roll;
    double yaw;
};

/*
 * Structure used by the Brouwer routines.  The first 6 structure
 * members, which are set by the user, are Brouwer mean elements.
 * The remaining elements are initialised by blinit() for subsequent use by
 * blosc().  Of these, only "darp", "dran" and "dman", which are the time
 * derivatives of "arp", "ran" and "man", are independently useful.
 */

struct BROUWER {
    double sma, /* Semi-major axis (metres) */
    ecc, /* Eccentricity */
    inc, /* Inclination (radians) */
    arp, /* Argument of perigee (radians) */
    ran, /* Right ascension of the ascending node (radians) */
    man; /* Mean anomaly (radians) */

    double darp, /* Secular rate of "arp" (radians/sec) */
    dran, /* Secular rate of "ran" (radians/sec) */
    dman; /* Secular rate of "man" (radians/sec) */

    /* You really don't want to know about these. */
    double a20, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10,
    b11, b12, b13, b14, b15, edp2, eta, eta2, eta3,
    eta6, gm2, gmp2, sidp, sidph, cidph, theta, theta2;
};

/*
 * Structure for specifying geodetic coordinates (e.g. an observer's
 * location).
 */

struct GEODETIC {
    double lat; /* Geodetic latitude (radians) */
    double lon; /* Longitude (radians) */
    double height; /* Height above the geoid (metres) */
};

/*
 * Structure for transfer of Keplerian orbital elements.
 */

struct KEPLER {
    double sma, /* Semi-major axis (metres) */
    ecc, /* Eccentricity */
    inc, /* Inclination (radians) */
    arp, /* Argument of perigee (radians) */
    ran, /* Right ascension of the ascending node (radians) */
    man; /* Mean anomaly (radians) */
};

/*
 * Structure for specifying nutation angles (radians).
 */
struct NUTATION {
    double mood; /* Mean obliquity of date */
    double nutlon; /* Nutation in longitude */
    double nutobl; /* Nutation in obliquity */
};

/*
 * Structure for specifying precession angles (radians).
 */
struct PRECESSION {
    double zeta;
    double z;
    double theta;
};

/*
 * Structure for specifying right ascension and declination (radians).
 */
struct RADEC {
    double ra;
    double dec;
};

/*
 * Structure used by the SGP4 routines.  The first 7 members are the
 * orbital elements at epoch + the drag coefficient.  The remaining
 * members are initialised by sgp4_init() for subsequent use by
 * sgp4_pred().
 */

struct SGP4 {
    double xno, /* Mean motion (radians/sec) */
    eo, /* Eccentricity */
    xincl, /* Inclination (radians) */
    omegao, /* Argument of perigee (radians) */
    xnodeo, /* Longitude of ascending node (radians) */
    xmo, /* Mean anomaly (radians) */
    bstar; /* Drag coefficient */

    double omgdot, /* Secular rate of omega (radians/sec) */
    xnodot, /* Secular rate of xnode (radians/sec) */
    xmdot; /* Secular rate of xm    (radians/sec) */

    /* You probably don't want to know about these. */
    int isimp;
    double aodp, aycof, c1, c1sq, c4, c5, cosio, d2, d3, d4, delmo,
    eta, omgcof, sinio, sinmo, t2cof, t3cof, t4cof, t5cof,
    x1mth2, x3thm1, x7thm1, xhdot1, xlcof, xmcof, xnodcf, xnodp;
};

/*
 * Two-line-element values.  The function tle_parse() can be used to
 * read the values from an element set into this structure.
 */

struct TLE_VALUES {
    int satellite_number; /* NORAD satellite number */
    char designator[9]; /* International designator */
    int year; /* Last two digits of year */
    double dayofyr; /* Day of year (1 origin) */
    double mmdt; /* First derivative of mean motion */
    double mmdt2; /* Second derivative of mean motion */
    double bstar; /* B* drag term */
    int ephem_type;
    int element_number;
    double incl; /* Inclination (deg) */
    double raan; /* Right ascension of asc node (deg) */
    double ecc; /* Eccentricity */
    double argp; /* Argument of perigee (deg) */
    double manom; /* Mean anomaly (deg) */
    double mmotion; /* Mean motion (revs/day) */
    int revolution; /* Epoch orbit number */
};

/*
 * Structure for specifying topocentric coordinates.
 */

struct TOPOCENTRIC {
    double az; /* Azimuth (radians) */
    double el; /* Elevation (radians) */
    double range; /* Range (metres) */
};



/*
 * Function prototypes.
 */

/* brouwer.c */
int blinit(struct BROUWER *bp, char **errstr);
int blosc(struct BROUWER *bp, double delt, struct KEPLER *kp,
        char **errstr);

/* ctogd.c */
void ctogd(double r[3], double gha, struct GEODETIC *geod,
        double dgeod[3][3]);

/* ctok.c */
int ctok(double r[6], struct KEPLER *keplr, char **errmsg);

/* ctotc.c */
void ctotc(double r[3], struct GEODETIC *obs, double gha,
        struct TOPOCENTRIC *top, double daer[3][3]);

/* eanom.c */
double eanom(double manom, double ecc, char **errstr);

/* efitoeci.c */
int efitoeci(double r[6], double gha, double rp[6]);

/* efrtoeci.c */
int efrtoeci(double r[6], double gha, double rp[6]);

/* ecitoefi.c */
int ecitoefi(double r[6], double gha, double rp[6]);

/* ecitoefr.c */
int ecitoefr(double r[6], double gha, double rp[6]);

/* gmha.c */
double gmha(double tjd);
double gmhadot(double tjd);

/* julian.c */
void caldat(int32_t julian, int *year, int *month, int *day);
double itojul(int32_t date, int32_t time);
int32_t julday(int year, int month, int day);
void jultoi(double tjd, int32_t *date, int32_t *time);
void jultotm(double tjd, struct tm *time);
double tmtojul(struct tm *time);

/* ktoc.c */
int ktoc(struct KEPLER *kp, double r[6], char **errstr);

/* ktosgp4.c */
int ktosgp4(struct KEPLER *kepler, struct SGP4 *sgp4, char **errmsg);

/* locate.c */
int locate(double r[6], struct ATTITUDE *sat,
        struct ATTITUDE *scanner, double loc[3]);

/* moonpos.c */
void moonpos(double tjd, double r[3]);

/* nutate.c */
void nut_angles(double tjd, struct NUTATION *angles, struct NUTATION *rates);
void nut_matrix(double tjd, double a[3][3]);

/* obliq.c */
void obliq(double tjd, double *mood, double *dmood);

/* precess.c */
void pre_angles(double tjd0, double tjd1, struct PRECESSION *angles);
void pre_reduce(double tjd0, struct RADEC *rad0,
        double tjd1, struct RADEC *rad1);

/* rdtotc.c */
void rdtotc(struct RADEC *rad, struct GEODETIC *obs,
        double gmha, struct TOPOCENTRIC *top);

/* refract.c */
double refract(double alt, double temp, double pressure);
double unrefract(double alt, double temp, double pressure);

/* sgp4.c */
int sgp4_init(struct SGP4 *sgp4, char **errstr);
int sgp4_pred(struct SGP4 *sgp4, double tsince, double rect[6],
        char **errstr);
double sgp4_orbit(struct SGP4 *sgp4, struct TLE_VALUES *tv, double tsince);
int sgp4_tle(double *epoch, struct SGP4 *sgp4, char *line1,
        char *line2, char **errstr);

/* sunpos.c */
int sunpos(double tjd, double r[3], char **errstr);

/* tconv.c */
double tconv(double tjd, char conv[], char **errstr);
int tcset(double x, char what[], char **errstr);

/* tle.c */
int tle_parse(char *line1, char *line2, struct TLE_VALUES *tv,
        char **errstr);

#endif /* __ORBIT__ */
