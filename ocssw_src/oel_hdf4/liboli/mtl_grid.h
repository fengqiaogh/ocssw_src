/*
 * Name:	mtl_grid.h
 *
 * Purpose:	Include file defining structures and function prototypes
 *		for scene gridding for angle generation.
 */
#ifndef	_MTL_GRID_
#define	_MTL_GRID_

#include "mtl_geometry.h"
#include "smeta_exploit.h"
#include "smeta_geometry.h"

#define	NUM_GRID_COL	9
#define	NUM_GRID_ROW	15
#define	NUM_GRID_ROW_LS	13
#define	GRID_COL_INC	0.25
#define	GRID_ROW_INC	0.1
#define	GRID_ROW_BASE	-0.7
#define	GRID_ROW_BASE_LS -0.6
#define	GRID_COL_BASE	-1.0

typedef struct MTL_GRID_PT
{
  double	npix;			/* Normalized pixel coordinate - added for older LS */
  double	ndet;			/* Normalized detector coordinate */
  double	frow;			/* Fractional row offset */
  double	l1t_line;		/* L1T line coordinate */
  double	l1t_samp;		/* L1T sample coordinate */
  VECTOR	V;			/* LSR view vector */
  VECTOR	S;			/* LSR sun vector */
} MTL_GRID_PT;

typedef struct MTL_GRID_SCA
{
  int		sca;					/* SCA number */
  double	min_line;				/* Minimum line coordinate */
  double	max_line;				/* Maximum line coordinate */
  double	min_samp;				/* Minimum sample coordinate */
  double	max_samp;				/* Maximum sample coordinate */
  MTL_GRID_PT	pgrid[NUM_GRID_ROW][NUM_GRID_COL];	/* SCA grid points */
  double	ls_to_ndet[4];				/* Line/samp to ndet coefficients */
  double	ls_to_frow[4];				/* Line/samp to frow coefficients */
} MTL_GRID_SCA;

typedef struct MTL_GRID_SCA_LS
{
  int		sca;					/* SCA number */
  double	min_line;				/* Minimum line coordinate */
  double	max_line;				/* Maximum line coordinate */
  double	min_samp;				/* Minimum sample coordinate */
  double	max_samp;				/* Maximum sample coordinate */
  MTL_GRID_PT	pgrid[NUM_GRID_ROW_LS][NUM_GRID_COL];	/* SCA grid points */
  double	ls_to_ndet[4];				/* Line/samp to ndet coefficients */
  double	ls_to_frow[4];				/* Line/samp to frow coefficients */
} MTL_GRID_SCA_LS;

typedef struct MTL_GRID_BAND
{
  int		band;			/* Band number */
  int		nsca;			/* Number of SCAs */
   /* Sudipta new addition for OLI to extract SCA and Det numbers */
  int		num_detectors;				/* Number of detectors per SCA */
  /* End Sudipta new addition for OLI */
  MTL_GRID_SCA	*sgrid;			/* Array of SCA grids */
} MTL_GRID_BAND;

typedef struct MTL_GRID_BAND_LS
{
  int		band;			/* Band number */
  int		nsca;			/* Number of SCAs */
  MTL_GRID_SCA_LS	*sgrid;			/* Array of SCA grids */
} MTL_GRID_BAND_LS;

int calc_band_grid(
  WRS2 parms,				/* I: WRS2 system parameters */
  SMETA smeta,				/* I: Scene metadata structure */
  double dellon,               	        /* I: Path longitude adjustment (in degrees) */
  VECTOR ecfsun,			/* I: ECF sun vector */
  int band_index,			/* I: Index for current band */
  MTL_GRID_BAND *bgrid );		/* O: Band grid structure */

/* variant for older LS */
int calc_band_grid_ls(
  WRS2 parms,				/* I: WRS2 system parameters */
  SMETA smeta,				/* I: Scene metadata structure */
  double dellon,               	        /* I: Path longitude adjustment (in degrees) */
  VECTOR ecfsun,			/* I: ECF sun vector */
  int band_index,			/* I: Index for current band */
  MTL_GRID_BAND_LS *bgrid );		/* O: Band grid structure */

int calc_grid_point(
  WRS2 parms,				/* I: WRS2 system parameters */
  SMETA smeta,				/* I: Scene metadata structure */
  double dellon,                        /* I: Path longitude adjustment (in degrees) */
  int band_index,			/* I: Index for current band */
  int sca,				/* I: SCA number (1 to nsca) */
  VECTOR ecfsun,			/* I: ECF sun vector */
  MTL_GRID_PT *gpt );			/* I/O: Grid point with only ndet and frow on input */


int calc_grid_fit(
  MTL_GRID_SCA *sgrid );		/* I/O: SCA grid array to compute polynomial fit for */

/* variant of calc_grid_fit for older LS */
int calc_grid_fit_ls(
  MTL_GRID_SCA_LS *sgrid );		/* I/O: SCA grid array to compute polynomial fit for */

int simple_inverse(
  int n,                                /* I: Dimension of matrix to invert */
  double a[4][4] );                        /* I/O: Matrix to invert */

int smeta_angles(
  double line,                  /* I: L1T line coordinate */
  double samp,                  /* I: L1T sample coordinate */
  MTL_GRID_BAND bgrid,          /* I: MTL grid structure for current band */
  /* Sudipta new addition for OLI to extract SCA and Det numbers */
  short *sca_num,		/* O: SCA containing current point */
  short *det_num,		/* O: Detector generating current point */
  /* End Sudipta new addition for OLI */
  double *sat_zn,                /* O: Viewing zenith angle scaled to 0.01 degrees */
  double *sat_az,                /* O: Viewing azimuth angle scaled to 0.01 degrees */
  double *sun_zn,                /* O: Solar zenith angle scaled to 0.01 degrees */
  double *sun_az );              /* O: Solar azimuth angle scaled to 0.01 degrees */


int smeta_angles_ls(
  double line,                  /* I: L1T line coordinate */
  double samp,                  /* I: L1T sample coordinate */
  MTL_GRID_BAND_LS bgrid,          /* I: MTL grid structure for current band */
  short *sat_zn,                /* O: Viewing zenith angle scaled to 0.01 degrees */
  short *sat_az,                /* O: Viewing azimuth angle scaled to 0.01 degrees */
  short *sun_zn,                /* O: Solar zenith angle scaled to 0.01 degrees */
  short *sun_az );              /* O: Solar azimuth angle scaled to 0.01 degrees */
#endif
