/*
 *	Name:	smeta.h
 *
 *	Purpose:
 *	Include file defining the product metadata structure that is used to generate the
 *	scene viewing geometry information needed to calculate per-pixel view angles.
 */
#ifndef	_SMETA_H_
#define	_SMETA_H_

/* Defines the maximum number of normal bands present on any of the supported
   satellites.  Note that this should only be used if absolutely necessary.
   It is preferred that arrays be sized dynamically. */
#define IAS_MAX_NBANDS 11
#define IAS_MAX_NBANDS_LS 8
/* Defines the maximum number of SCAs present on any of the supported 
   sensors.   Note that this should only be used if absolutely necessary.
   It is preferred that arrays be sized dynamically. */
#define IAS_MAX_NSCAS 14

#define IAS_MAX_NSCANS 375

#include "ias_const.h"
#include "ias_structures.h"			/* IAS_VECTOR and IAS_CORNERS                */

typedef enum LANDSAT_SPACECRAFT
{
    unknown_spacecraft_type = -1,
    IAS_L1=1,
    IAS_L2,
    IAS_L3,
    IAS_L4,
    IAS_L5,
    IAS_L7,
    IAS_L8
} LANDSAT_SPACECRAFT;

typedef enum LANDSAT_SENSOR
{
    unknown_sensor_type = -1,
    IAS_MSS=1,
    IAS_TM,
    IAS_ETM,
    IAS_OLI_TIRS
} LANDSAT_SENSOR;

typedef enum SENSOR_MODE
{
    unknown_mode = -1,
    MSS=1,
    TM_SAM,
    TM_BUMPER,
    OLI
} SENSOR_MODE;


typedef struct SMETA_SCENE_PROJ
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
    int		spheroid;			/* Projection spheroid code                  */
    int		zone;				/* Projection zone code for UTM or
						   State Plane projections.                  */
    double	projprms[IAS_PROJ_PARAM_SIZE];
		/* Array of 15 projection coefficients as required by the projection 
		   transformation package.  Refer to the projection package documentation 
		   for a description of each field for a given projection.                   */
    /* ==== Grid corners. ====                                                               */
    struct	IAS_CORNERS corners;		/* Projection coordinates of the resulting 
						   output image's four corners.              */
} SMETA_SCENE_PROJ;

typedef struct SMETA_BAND
{
    int			band;			/* User band number                          */
    int			l1t_lines;		/* Number of lines in the L1T image.         */
    int			l1t_samps;		/* Number of samples in the L1T image.       */
    int			nsca;			/* Number of SCAs in this band               */
    /* Sudipta new addition for OLI to extract SCA and Det numbers */
    int			num_detectors;		/* Number of detectors per SCA for this band */
     /* End Sudipta new addition for OLI */
    double		pixsize;		/* Projection distance per pixel in meters.  */
    double		align[3][3];		/* Instrument to spacecraft alignment        */
    double		legendre[IAS_MAX_NSCAS][2][4];   /* Per SCA Legendre coefficients    */
} SMETA_BAND;

/* for older LS */
typedef struct SMETA_BAND_LS
{
    int			band;			/* User band number                          */
    int			l1t_lines;		/* Number of lines in the L1T image.         */
    int			l1t_samps;		/* Number of samples in the L1T image.       */
    int                 nsca;                   /* Number of SCAs                            */
    int			nscans;			/* Number of scans in this band              */
    int			nlines;			/* Number of lines within a scans in this band */
    int			nsamps;			/* Number of samples within a scans in this band */
    double		pixsize;		/* Projection distance per pixel in meters.  */
    double		align[3][3];		/* Instrument to spacecraft alignment        */
    double		fascan;                 /* Forward along scan field of view          */
    double		rascan;                 /* Reverse along scan field of view          */
    double		xscan;                  /* Across scan field of view                 */
    double		afov;                   /* Along detector field of view              */
    double		xfov;                   /* Across detector field of view             */
    double		aoffset;                /* Along scan band offset                    */
    double		xoffset;                /* Across scan band offset                   */
    double              center_sample;          /* MSS Center sample location                */
    double              sample_slope;           /* MSS detectors / radian scale              */
} SMETA_BAND_LS;


typedef struct SMETA
{
    char		scene_id[22];		/* Scene ID as Lx8ppprrryyyydddLGNnn         */
    LANDSAT_SPACECRAFT  spacecraft_id;          /* Landsat spacecraft ID          */
    LANDSAT_SENSOR      sensor_id;              /* Landsat sensor ID                         */
    SENSOR_MODE		sensor_mode;		/* Sensor operating mode (for TM/ETM+)       */
    int			wrs_path;		/* WRS path number (target)                  */
    int			wrs_row;		/* WRS row number (target)                   */
    double		roll_angle;		/* Off-nadir roll angle in degrees           */
    int			acq_date[3];		/* Acquisition year, month, day              */
    double		sun_azimuth;		/* Scene center sun azimuth in degrees       */
    double		sun_elevation;		/* Scene center sun elevation in degrees     */
    SMETA_SCENE_PROJ	projection;		/* Ground reference for this L1T scene       */
    int			num_band;		/* Number of bands in the metadata structure */
    SMETA_BAND		band_smeta[IAS_MAX_NBANDS];	/* Metadata for each band            */
    SMETA_BAND_LS	band_smeta_ls[IAS_MAX_NBANDS_LS];	/* Metadata for each band for older LS */
} SMETA;

#endif
