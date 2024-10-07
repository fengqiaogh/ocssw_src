#ifndef _IAS_STRUCTURES_H_
#define _IAS_STRUCTURES_H_


/* A data structure for decomposed date/time information */
typedef struct IAS_DATETIME
{
    int    year;
    int    month;
    int    day_of_month;
    int    day_of_year;
    int    hour;
    int    minute;
    double second;            /* Allows for fractional seconds */
}   IAS_DATETIME;


typedef struct IAS_DBL_XY
    {
    double x;           /* X value                           */
    double y;           /* Y value                           */
    } IAS_DBL_XY;

typedef struct IAS_VECTOR
    {
    double x;           /* Vector X component */
    double y;           /* Vector Y component */
    double z;           /* Vector Z component */
    } IAS_VECTOR;

typedef struct IAS_FLOAT_VECTOR
{
    float x;
    float y;
    float z;
}IAS_FLOAT_VECTOR;

typedef struct IAS_QUATERNION
{
    IAS_VECTOR vector;
    double scalar;
}IAS_QUATERNION;

typedef struct IAS_COMPLEX
   {
   double re;
   double im;
   } IAS_COMPLEX;

typedef struct IAS_LNG_XY
{
  int x;                  /* X value                           */
  int y;                  /* Y value                           */
} IAS_LNG_XY;

typedef struct IAS_LNG_HW
{
  int height;             /* Height value                      */
  int width;              /* Width value                       */
} IAS_LNG_HW;

typedef struct IAS_LNG_LS
{
  int line;               /* Line value                        */
  int samp;               /* Sample value                      */
} IAS_LNG_LS;

typedef struct IAS_DBL_LS
{
  double line;             /* Line value                        */
  double samp;             /* Sample value                      */
} IAS_DBL_LS;

typedef struct IAS_DBL_LAT_LONG
{
  double lat;              /* Latitude value                    */
  double lng;              /* Longitude value                   */
} IAS_DBL_LAT_LONG;

#define COEFS_SIZE 4
typedef struct IAS_COEFFICIENTS
{
  double a[COEFS_SIZE];    /* Array of a coefficients           */
  double b[COEFS_SIZE];    /* Array of b coefficients           */
} IAS_COEFFICIENTS;

typedef struct IAS_CORNERS
{
  IAS_DBL_XY upleft;    /* X/Y value of the upper left corner  */
  IAS_DBL_XY upright;   /* X/Y value of the upper right corner */
  IAS_DBL_XY loleft;    /* X/Y value of the lower left corner  */
  IAS_DBL_XY loright;   /* X/Y value of the lower right corner */
} IAS_CORNERS;

typedef struct IAS_IMAGE
{
  int data_type;
  int nl;
  int ns;
  int band_number;
  short *data; 
  double pixel_size_x;
  double pixel_size_y;
  IAS_CORNERS corners;
} IAS_IMAGE;

#endif
