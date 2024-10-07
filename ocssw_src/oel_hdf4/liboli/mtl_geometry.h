/*
 * Name:	mtl_geometry.h
 *
 * Purpose:	Include file defining structures and function prototypes
 *		for MTL file geometric processing, including scene gridding,
 *              for angle generation.
 */
#ifndef	_MTL_GEOMETRY_
#define	_MTL_GEOMETRY_

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "smeta.h"
#include "smeta_exploit.h"
#include "smeta_geometry.h"

#define	MIN(a,b)	(a<b?a:b)
#define	MAX(a,b)	(a>b?a:b)

typedef struct WRS2
{
  double	semimajor;		/* Earth semi-major axis in meters */
  double	semiminor;		/* Earth semi-minor axis in meters */
  double	inertialvel;		/* Earth inertial angular velocity in radians/second */
  double	orbitincl;		/* Orbital inclination in degrees */
  double	orbitrad;		/* Orbital radius in meters */
  double	path1long;		/* Longitude of path 1 descending node in degrees */
  int		numpath;		/* Number of WRS paths (233) */
  int		numrow;			/* Number of WRS rows (248) */
  int		cycle;			/* Number of days in one WRS2 cycle (16) */
  int		row0lat;		/* Equator crossing row (60) */
} WRS2;

typedef struct VECTOR
{
  double	x;			/* X coordinate */
  double	y;			/* Y coordinate */
  double	z;			/* Z coordinate */
} VECTOR;

int frame_scene(			/* Fits the nominal WRS to the MTL scene frame */
  WRS2 wrsparms,			/* I: Structure of WRS-2 parameters */
  SMETA smeta,				/* I: Input scene metadata structure */
  double *dellon,			/* O: WRS longitude offset */
  double *delrow );			/* O: Scene length in fractional WRS rows */

/* different variant for older LS algo */
int frame_scene_ls(			/* Fits the nominal WRS to the MTL scene frame */
  WRS2 wrsparms,			/* I: Structure of WRS-2 parameters */
  SMETA smeta,				/* I: Input scene metadata structure */
  double *dellon,			/* O: WRS longitude offset */
  double *delrow );			/* O: Scene length in fractional WRS rows */


int project_point(			/* Projects a LOS vector to the Earth */
  WRS2 wrsparms,			/* I: Structure of WRS-2 parameters */
  SMETA smeta,				/* I: Input scene metadata structure */
  double dellon,			/* I: WRS longitude offset */
  int band,				/* I: Band number (user band) */
  int sca,				/* I: SCA number (1 to nsca) */
  double frow,				/* I: Fractional WRS row coordinate */
  double ndet,				/* I: Normalized detector coordinate */
  double *l1t_line,			/* O: Projected line coordinate */
  double *l1t_samp );			/* O: Projected sample coordinate */

/* Slightly different variant for older LS */
int project_point_ls(			/* Projects a LOS vector to the Earth */
  WRS2 wrsparms,			/* I: Structure of WRS-2 parameters */
  SMETA smeta,				/* I: Input scene metadata structure */
  double dellon,			/* I: WRS longitude offset */
  int band,				/* I: Band number (user band) */
  double frow,				/* I: Fractional WRS row coordinate */
  double nsamp,				/* I: Sample coordinate */
  double nline,				/* I: Line coordinate */
  double *l1t_line,			/* O: Projected line coordinate */
  double *l1t_samp );			/* O: Projected sample coordinate */

void pathrow_to_latlon(			/* Calculates nominal WRS scene center */
  WRS2  parms,				/* I: WRS2 system parameters */
  int   path,				/* I: WRS2 path */
  double wrsrow,			/* I: WRS2 fractional row */
  double *lat,				/* O: Nominal latitude */
  double *lon );			/* O: Nominal longitude */

double dotprod(				/* Returns the dot product of two vectors */
  VECTOR a,				/* I: First input vector */
  VECTOR b );				/* I: Second input vector */

void crossprod(				/* Computes vector cross product */
  VECTOR a,				/* I: First input vector */
  VECTOR b,				/* I: Second input vector */
  VECTOR *c );				/* O: Output vector */

double vecnorm(				/* Returns vector magnitude */
  VECTOR a );				/* I: Input vector */

void unitvec(				/* Normalizes input vector */
  VECTOR a,				/* I: Input vector */
  VECTOR *b);				/* O: Normalized output vector */

void rotatevec(
  double mat[3][3],			/* I: 3-by-3 rotation matrix */
  VECTOR a,				/* I: Original vector */
  VECTOR *b );				/* O: Rotated vector */

int calc_yaw_steering_gp(		/* Returns 0 for success and -1 for failure */
  WRS2 parms,				/* I: WRS2 system parameters */
  VECTOR pos,				/* I: Position vector */
  VECTOR vel,				/* I: Velocity vector */
  VECTOR i_los,				/* I: Instrument line-of-sight unit vector */
  double roll,				/* I: Off-nadir roll angle (in degrees) */
  double *gp_lat,			/* O: Ground point geodetic latitude (in degrees) */
  double *gp_lon);			/* O: Ground point geodetic longitude (in degrees) */

/* variant of calc_yaw_steering_gp for other LS */
int calc_yaw_steering_gp_ls(		/* Returns 0 for success and -1 for failure */
  WRS2 parms,				/* I: WRS2 system parameters */
  VECTOR pos,				/* I: Position vector */
  VECTOR vel,				/* I: Velocity vector */
  VECTOR i_los,				/* I: Instrument line-of-sight unit vector */
  double roll,				/* I: Off-nadir roll angle (in degrees) */
  double *gp_lat,			/* O: Ground point geodetic latitude (in degrees) */
  double *gp_lon);			/* O: Ground point geodetic longitude (in degrees) */

int calc_los(
  SMETA_BAND band_smeta,		/* I: Metadata band structure */
  int sca,				/* I: SCA number */
  double ndet,				/* I: Normalized detector coordinate */
  VECTOR *los );			/* O: Line of sight unit vector */

/* Slightly different variant for older LS */
int calc_los_ls(
  SMETA_BAND_LS band_smeta,		/* I: Metadata band structure */
  double frow,                          /* I: Fractional row */
  double npix,				/* I: Sample number */
  double ndet,				/* I: Line number */
  VECTOR *los );			/* O: Line of sight unit vector */

void pathrow_to_posvel(
  WRS2  parms,				/* I: WRS2 system parameters */
  int   path,				/* I: WRS2 path */
  double dellon,                        /* I: Path longitude adjustment (in degrees) */
  double wrsrow,			/* I: WRS2 fractional row */
  VECTOR *pos,				/* O: Nominal position vector */
  VECTOR *vel );			/* O: Nominal velocity vector */

int get_ecfsun( 
  WRS2 wrsparms,			/* I: Structure of WRS-2 parameters */
  SMETA smeta,				/* I: Input scene metadata structure */
  double dellon,			/* I: WRS longitude offset */
  VECTOR *ecfsun );			/* O: Sun vector in ECEF coordinates */

/* Slightly different variant for older LS */
int get_ecfsun_ls(
  WRS2 wrsparms,			/* I: Structure of WRS-2 parameters */
  SMETA smeta,				/* I: Input scene metadata structure */
  double dellon,			/* I: WRS longitude offset */
  VECTOR *ecfsun );			/* O: Sun vector in ECEF coordinates */
#endif
