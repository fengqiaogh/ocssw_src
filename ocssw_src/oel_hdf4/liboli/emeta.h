/*
 *	Name:	emeta.h
 *
 *	Purpose:
 *	Include file defining the enhanced metadata structure that is used to capture the
 *	scene viewing geometry information needed to calculate per-pixel view angles. This
 *      information includes linkages between L1T line/sample pixel coordinates and the
 *	corresponding time of observation, so it can also be used to generate per-pixel
 *	sun angles.
 */
#ifndef	_EMETA_H_
#define	_EMETA_H_

/* Defines the maximum number of normal bands present on any of the supported
   satellites.  Note that this should only be used if absolutely necessary.
   It is preferred that arrays be sized dynamically. */
#define IAS_MAX_NBANDS 11

/* Defines the maximum number of SCAs present on any of the supported 
   sensors.   Note that this should only be used if absolutely necessary.
   It is preferred that arrays be sized dynamically. */
#define IAS_MAX_NSCAS 14

#include "ias_const.h"
#include "ias_structures.h"			/* IAS_VECTOR and IAS_CORNERS                */

#define		NUM_RPC_COEF	5		/* Number of image RPC coefficents in model  */
#define		ANG_RPC_COEF	10		/* Number of angle RPC coefficents in model  */

typedef struct EMETA_IMG_RPC
{
    int		sca_num;			/* Number of current SCA                     */
    double	mean_l1r[2];			/* L1R mean line, sample offset              */
    double	mean_l1t[2];			/* L1T mean line, sample offset              */
    double	mean_hgt;			/* Mean height offset                        */
    double	line_num[NUM_RPC_COEF];		/* L1R line RPC numerator coefficients       */
    double	line_den[NUM_RPC_COEF-1];	/* L1R line RPC denominator coefficients     */
    double	samp_num[NUM_RPC_COEF];		/* L1R sample RPC numerator coefficients     */
    double	samp_den[NUM_RPC_COEF-1];	/* L1R sample RPC denominator coefficients   */
} EMETA_IMG_RPC;

typedef struct EMETA_ANG_RPC
{
    double	mean_l1t[2];			/* L1T mean line, sample offset              */
    double	mean_hgt;			/* Mean height offset                        */
    double	mean_l1r[2];			/* L1R mean line, sample offset              */
    double	mean_vec[3];			/* Mean vector component offsets             */
    double	x_num[ANG_RPC_COEF];		/* X axis RPC numerator coefficients         */
    double	x_den[ANG_RPC_COEF-1];		/* X axis RPC denominator coefficients       */
    double	y_num[ANG_RPC_COEF];		/* Y axis RPC numerator coefficients         */
    double	y_den[ANG_RPC_COEF-1];		/* Y axis RPC denominator coefficients       */
    double	z_num[ANG_RPC_COEF];		/* Z axis RPC numerator coefficients         */
    double	z_den[ANG_RPC_COEF-1];		/* Z axis RPC denominator coefficients       */
} EMETA_ANG_RPC;

typedef struct EMETA_SCENE_PROJ
{
    /* ==== WGS84 ellipsoid parameters ====                                                  */
    double	wgs84_major_axis;		/* WGS 84 ellipsoid semi-major axis (meters) */
    double	wgs84_minor_axis;		/* WGS 84 ellipsoid semi-minor axis (meters) */
    /* ==== Projection parameters and info. ====                                             */
    char	units[IAS_UNITS_SIZE];		/* Projection units string                   */
    int		code;				/* Projection code for the output space image.
						   Values for this field are defined in the
						   "gctp.h" include file.                    */
    char	datum[IAS_DATUM_SIZE];		/* Projection datum string                   */
    //int		spheroid;			/* Projection spheroid code                  */
    int		zone;				/* Projection zone code for UTM or
						   State Plane projections.                  */
    double	projprms[IAS_PROJ_PARAM_SIZE];
		/* Array of 15 projection coefficients as required by the projection 
		   transformation package.  Refer to the projection package documentation 
		   for a description of each field for a given projection.                   */
    /* ==== Grid corners. ====                                                               */
    struct	IAS_CORNERS corners;		/* Projection coordinates of the resulting 
						   output image's four corners.              */
} EMETA_SCENE_PROJ;

typedef struct EMETA_BAND
{
    int			band;			/* User band number                          */
    int			nsca;           	/* Number of SCAs in the band                */
    int			l1t_lines;		/* Number of lines in the L1T image.         */
    int			l1t_samps;		/* Number of samples in the L1T image.       */
    int			l1r_lines;		/* Number of lines in the L1R image.         */
    int			l1r_samps;		/* Number of samples in one L1R SCA.         */
    double		pixsize;		/* Projection distance per pixel in meters.  */
    double		img_start_time;		/* Image start time, offset in seconds from
						   ephemeris epoch                           */
    double		line_time;		/* L1R line time increment in seconds        */
    EMETA_ANG_RPC	sat;			/* Satellite view vector RPC structure       */
    EMETA_ANG_RPC	sun;			/* Solar illumination vector RPC structure   */
    EMETA_IMG_RPC	rpc[IAS_MAX_NSCAS];	/* Per SCA RPC coefficient structures        */
} EMETA_BAND;

typedef struct EMETA_EPHEM
{
    double	eph_sample_time;		/* Sample time from ephemeris epoch          */
    IAS_VECTOR	ecef_position;			/* ECEF positions                            */
} EMETA_EPHEM;

typedef struct EMETA
{
    //int			wrs_path;		/* WRS path number (target)                  */
    //int			wrs_row;		/* WRS row number (target)                   */
    EMETA_SCENE_PROJ	projection;		/* Ground reference for this L1T scene       */
    int			num_band;		/* Number of bands in the metadata structure */
    EMETA_BAND		band_emeta[IAS_MAX_NBANDS];	/* Enhanced metadata for each band   */
    double		eph_epoch_time[3];	/* UTC time ([0]=year,
						             [1]=day of year,
						             [2]=second of day )
						   of ephemeris epoch                        */
    int			ephem_count;		/* Number of ephemeris samples               */
    EMETA_EPHEM		*ephemeris;		/* Ephemeris data for scene                  */
    EMETA_EPHEM		*sunvector;		/* ECEF solar vectors for scene              */
} EMETA;

int emeta_allocate_ephemeris(
    EMETA	*emeta,				/* Enhanced metadata structure               */
    int		nephem );			/* Number of ephemeirs points required       */

void emeta_free_ephemeris(
    EMETA	*emeta );			/* Enhanced metadata structure               */

#endif
