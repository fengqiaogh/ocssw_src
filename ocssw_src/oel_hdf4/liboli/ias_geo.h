#ifndef IAS_GEO_H
#define IAS_GEO_H

#include "ias_const.h"
#include "ias_structures.h"

/* Define a zone number for all projections except State Plane and UTM. */
#define NULLZONE 62

#define WGS84_SPHEROID 12

/* Type defines for projection related structures */
typedef struct ias_geo_proj_transformation IAS_GEO_PROJ_TRANSFORMATION;
/* The ias_projection structure matches the gctp_projection structure
   definition.  The gctp_projection structure is not included here to prevent
   needing to modify the build to find gctp.h everywhere ias_geo.h is used. */
typedef struct ias_projection
{
    int proj_code;      /* Projection code */
    int zone;           /* Projection zone number - only has meaning for
                           projections like UTM and stateplane */
    int units;          /* Units of coordinates */
    int spheroid;       /* Spheroid code for the projection */
    double parameters[IAS_PROJ_PARAM_SIZE];
                        /* Array of projection parameters */
} IAS_PROJECTION;

int ias_geo_convert_deg2dms
(
    double deg,         /* I: Angle in seconds, minutes, or degrees */
    double *dms,        /* O: Angle converted to DMS */
    const char *check   /* I: Angle usage type (LAT, LON, or DEGREES) */
);

int ias_geo_convert_dms2deg
(
    double angle_dms,     /* I: Angle in DMS (DDDMMMSSS) format */
    double *angle_degrees,/* O: Angle in decimal degrees */
    const char *type      /* I: Angle usage type (LAT, LON, or DEGREES) */
);

void ias_geo_convert_geod2cart
(
    double latitude,    /* I: Lat of geodetic coordinates in radians */
    double longitude,   /* I: Long of geodetic coordinates in radians*/
    double height,      /* I: Height (elevation) of geodetic coord in meters */
    double semimajor,   /* I: Reference ellipsoid semi-major axis in meters */
    double flattening,  /* I: Flattening of the ellipsoid 
                           (semimajor-semiminor)/semimajor */
    IAS_VECTOR *cart    /* O: Cartesian vector for the coord */
);

void ias_geo_find_deg
(
    double angle,  /* I: Angle in total degrees */
    int  *degree   /* O: Degree portion of the angle */
);


void ias_geo_find_min 
(
    double angle,       /* I: Angle in total degrees */
    int  *minute        /* O: Minute portion of the angle */
);

void ias_geo_find_sec 
(
    double angle,       /* I: Angle in total degrees */
    double *second      /* O: Second portion of the angle */
);

int ias_geo_get_units 
(
    const char *unit_name,   /* I: Units name */
    int *unit_num            /* O: Units number */
);

int ias_geo_does_cross_180
(
    int unit,             /* I: The angular unit of the angles passed in */
    double const corner_longitudes[4]
                          /* I: The longitude for each corner of the scene */
);

int ias_geo_add_once_around
(
    int unit,        /* I: The angular unit to use to interpret lon */
    double *lon      /* I/O: The angle to be adjusted */
);


void ias_geo_lagrange_interpolate
( 
    const double *seconds_from_ref, /* I: Array of n_pts reference times */
    const IAS_VECTOR *position, /* I: Array of n_pts position vectors */
    const IAS_VECTOR *velocity, /* I: Array of n_pts velocity vectors */
    int n_pts,           /* I: Number of points to use in interpolation */
    double delta_time,   /* I: Delta time from the reference time */
    IAS_VECTOR *interpolated_position,
                         /* O: New satellite position at delta_time */
    IAS_VECTOR *interpolated_velocity
                         /* O: New satellite velocity at delta_time */
);

int ias_geo_report_proj_err 
(
    int err              /* I: Error returned from the GCTP call */
);

int ias_geo_transform_projection 
(
    int inproj,             /* I: Input projection code */
    int inunit,             /* I: Input projection units code */
    int inzone,             /* I: Input projection zone code */
    const double *inparm,   /* I: Array of 15 projection parameters--input */
    int inspheroid,         /* I: Input spheroid code */
    int outproj,            /* I: Output projection code */
    int outunit,            /* I: Output projection units code */
    int outzone,            /* I: Output projection zone code */
    const double *outparm,  /* I: Array of 15 projection parameters--output */
    int outspheroid,        /* I: Output spheroid code */
    double inx,             /* I: Input X projection coordinate */
    double iny,             /* I: Input Y projection coordinate */
    double *outx,           /* O: Output X projection coordinate */
    double *outy            /* O: Output Y projection coordinate */
);

IAS_GEO_PROJ_TRANSFORMATION *ias_geo_create_proj_transformation
(
    const IAS_PROJECTION *source_projection, /* I: source projection */
    const IAS_PROJECTION *target_projection  /* I: target projection */
);

void ias_geo_destroy_proj_transformation
(
    IAS_GEO_PROJ_TRANSFORMATION *trans
);

void ias_geo_only_allow_threadsafe_transforms();

int ias_geo_transform_coordinate
(
    const IAS_GEO_PROJ_TRANSFORMATION *trans, /* I: transformation to use */
    double inx,             /* I: Input X projection coordinate */
    double iny,             /* I: Input Y projection coordinate */
    double *outx,           /* O: Output X projection coordinate */
    double *outy            /* O: Output Y projection coordinate */
);

void ias_geo_set_projection
(
    int proj_code,          /* I: input projection code */
    int zone,               /* I: input zone */
    int units,              /* I: input units */
    int spheroid,           /* I: input spheroid */
    const double *parms,    /* I: input projection parameters */
    IAS_PROJECTION *proj    /* I: target projection structure */
);

#endif
