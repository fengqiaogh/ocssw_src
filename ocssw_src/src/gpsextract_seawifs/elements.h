#ifndef _ELEMENTS_H
#define _ELEMENTS_H

typedef unsigned char BYTE;
typedef short int INT16;
typedef int32_t INT32;
typedef float FLOAT32;
typedef double FLOAT64;

typedef struct hdr_struct {
    INT32 nrecs;
    BYTE fill[124];
} hdrstr;

typedef struct orb_struct {
    FLOAT64 sma; /* Semi-Major Axis (km)           */
    FLOAT64 eccen; /* Eccentricity                   */
    FLOAT64 incl; /* Inclination (deg)              */
    FLOAT64 ra; /* Lon. of Asc. Node (deg)        */
    FLOAT64 ap; /* Arg of Periapsis (deg)         */
    FLOAT64 ma; /* Mean anomaly (deg)             */
} orbstr;

typedef struct elements_struct {
    INT32 year; /* Year (e.g. 1996)               */
    INT32 day; /* Day of year                    */
    FLOAT64 sec; /* Seconds of day                 */
    orbstr orb1; /* Initial Orbit                  */
    orbstr orb2; /* Final (fitted) Orbit           */
    FLOAT64 cdrg; /* Drag Coefficient               */
    INT32 type; /* 0=init, 1=final (fitted)       */
    INT32 fill; /* spare                          */
} elmstr;

#endif
