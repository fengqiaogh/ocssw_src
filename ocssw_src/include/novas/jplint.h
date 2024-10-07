

#ifndef _JPLINT_
   #define _JPLINT_

/* structure JPLEphemType: contains JPL Planetary Ephemeris data for
                           major bodies in space

      tjdStart:     Start time of ephemeris data in tjd format

      tjdEnd:       End time of ephemeris data in tjd format

      span:         Span of each record

      au:           Astronomical unit

      emrat:        Earth Moon ratio

      ipt:          Start index, # of coeffs, # of sets for each body type

      coeffs:       Ephemeris data
*/
#define MAX_COEFFS 460

typedef struct JPLEphemType
{
  double tjdStart;
  double tjdEnd;
  double span;
  double au;
  double emrat;
  int    ipt[13][3];
  char   pad[4];
  double coeffs[MAX_COEFFS][826];
}JPL_EPHEM_TYPE;


/* structure JPLIntUtilType: contains fields formerly placed in 'COMMON'
                            section of JPL Fortran code

      pvsun:    Array containing barycentric position and velocity of the sun

      list:     Array specifying what interpolation is wanted for each body

      saveBary: Variable saving value of bary field

      bary:     Defines output center.  Default is 'False'
                    True(1):  Center is solar-system barycenter
                    False(0): Center is sun

      km:       Defines physical units of the output.  Default is 'False'
                    True(1):  KM and KM/sec
                    False(0): AU and AU/day

 
*/

typedef struct JPLIntUtilType
{
  double pvsun[6];
  int list[12];
  int saveBary;
  int bary;
  int km;  
}JPLINT_UTIL_TYPE;



/* Function prototypes
*/

  void jplint(double tjdTDB, int target, int center,
              double *outPosVel, int *error);

  int state(double tjdTDB, JPLIntUtilType *util, double posvel[13][6],
            double *pnut);

  void interp(double *ephemPtr, int startLoc, double *inTime, int numCoefs,
              int numCom, int numSets, int velFlag, double *posvel);

  int setJPLPtr(JPLEphemType *inPtr);

#endif

