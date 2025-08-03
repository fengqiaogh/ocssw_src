// Compute geodetic longitude, latitude and altitude 
//  from Cartesian orbit vector in ECR coordinates

// Uses geodetic approximation from Patt and Gregg, IJRS, 1993.

//	Arguments
//     
//     	Name    Type 	I/O 	Description
//     	----	---- 	--- 	-----------
//       nRecords   int   		I     number of records
//       orb      	double    	I     3 x nRecords array of orbit vectors
//       lon    	double   	O     nRecords array of longitude
//       lat     	double    	O     nRecords array of latitude
//       alt     	double    	O     nRecords array of altitude

// Liang Hong, Jan 14, 2016, Ported from IDL
// V0.1, Feb 1, 2016

#include <math.h>
#include <genutils.h>

int orb2lla(int nRecords, double orb[], double lon[], double lat[], double alt[]) {
    // Define constants
    double re = 6378.137; // Earth equatorial radius (km)
    double rem = 6371; // Earth mean radius (km)
    double f = 1.0 / 298.257; // Earth flattening factor

    int i;
    double omf2, omf2p, rad, pxy, temp, rl;
    omf2 = (1.0 - f)*(1.0 - f);

    for (i = 0; i < nRecords; i++) {
        // Compute longitude
        lon[i] = (180.0 / OEL_PI) * atan2(orb[1 + i * 6], orb[0 + i * 6]);

        // Compute geodetic latitude 
        rad = sqrt(orb[0 + i * 6] * orb[0 + i * 6] + orb[1 + i * 6] * orb[1 + i * 6] + orb[2 + i * 6] * orb[2 + i * 6]);
        omf2p = (omf2 * rem + rad - rem) / rad;
        pxy = orb[0 + i * 6] * orb[0 + i * 6] + orb[1 + i * 6] * orb[1 + i * 6];
        temp = sqrt(orb[2 + i * 6] * orb[2 + i * 6] + omf2p * omf2p * pxy);
        lat[i] = (180.0 / OEL_PI) * asin(orb[2 + i * 6] / temp);

        // Compute altitude
        rl = re * (1.0 - f) / sqrt(1.0 - (2.0 - f) * f * pxy / (rad * rad));
        alt[i] = rad - rl;
    }
    return 0;
}
