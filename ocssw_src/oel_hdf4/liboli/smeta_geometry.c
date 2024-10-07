/*
 *	Name:	smeta_geometry.c
 *
 *	Purpose:
 *	Source file containing the functions used to relate Level 1T pixel line/sample
 *	coordinates to the corresponding viewing and solar illumination angles. These functions
 *	operate on the SMETA data type defined in "smeta.h".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "gctp.h"				/* Projection constants                      */
#include "ias_logging.h"
#include "ias_math.h"
#include "ias_geo.h"
#include "ias_const.h"
#include "smeta_geometry.h"			/* Definition of SMETA structure and
						   function prototypes.                      */
#define	NUM_INTERP	4			/* Number of points to use in Lagrange interpolation */
#define	PIX_BUFFER	3			/* Resampling kernel reduction in effective SCA size */

static double	wgs84_a;			/* WGS84 semimajor axis in meters            */
static double	wgs84_flat;			/* WGS84 flattening                          */
static int	wgs84_init = 0;			/* Initialization flag                       */
static IAS_GEO_PROJ_TRANSFORMATION *trans=NULL;	/* Map projection transformation structure   */

int smeta_geodetic_to_ecef(
    double	lat,				/* Input latitude in radians                 */
    double	lon,				/* Input longitude in radians                */
    double	hgt,				/* Input height in meters                    */
    IAS_VECTOR	*ecef )				/* Output ECEF vector                        */
{
    /* Make sure ellipsoid constants were initialized */
    if ( wgs84_init < 1 )
    {
        IAS_LOG_ERROR("WGS84 ellipsoid constants not initialized.");
        return(ERROR);
    }

    /* Use library routine to perform the computations */
    ias_geo_convert_geod2cart( lat, lon, hgt, wgs84_a, wgs84_flat, ecef );

    return(SUCCESS);
}

int smeta_geodetic_to_ecf2lsr(
    double      lat,				/* Input latitude in radians                 */
    double      lon,				/* Input longitude in radians                */
    double      ecf2lsr[3][3] )			/* ECEF to LSR rotation matrix               */
{
    double	clat, slat;			/* cosine and sine of latitude               */
    double	clon, slon;			/* cosine and sine of longitude              */

    clat = cos( lat );
    slat = sin( lat );
    clon = cos( lon );
    slon = sin( lon );
    ecf2lsr[0][0] = -slon;
    ecf2lsr[0][1] =  clon;
    ecf2lsr[0][2] =  0.0;
    ecf2lsr[1][0] = -slat*clon;
    ecf2lsr[1][1] = -slat*slon;
    ecf2lsr[1][2] =  clat;
    ecf2lsr[2][0] =  clat*clon;
    ecf2lsr[2][1] =  clat*slon;
    ecf2lsr[2][2] =  slat;

    return(SUCCESS);
}

int smeta_proj_to_l1t(
    SMETA	*smeta,				/* Input enhanced metadata info              */
    int		band,				/* Input band number                         */
    double	proj_x,				/* Input projection X coordinate             */
    double	proj_y, 			/* Input projection Y coordinate             */
    double	*l1t_line,			/* Output L1T line number                    */
    double	*l1t_samp )			/* Output L1T sample number                  */
{
    int		band_index;			/* Index corresponding to input band number  */

    /* Look for the input band number in the metadata band list */
    *l1t_line = 0.0;
    *l1t_samp = 0.0;
    band_index = smeta_band_number_to_index( smeta, band );
    /* See if we found the band */
    if ( band_index < 0 )
    {
        IAS_LOG_ERROR("Band %d not found in metadata.", band);
        return(ERROR);
    }
    
    /* Convert projection X/Y to L1T line/sample */
    /* Note, this only works for north-up images */
    if ( smeta->spacecraft_id == IAS_L8 )
    {
        *l1t_samp = (proj_x - smeta->projection.corners.upleft.x)/smeta->band_smeta[band_index].pixsize;
        *l1t_line = (smeta->projection.corners.upleft.y - proj_y)/smeta->band_smeta[band_index].pixsize;
    }
    else
    {
        *l1t_samp = (proj_x - smeta->projection.corners.upleft.x)/smeta->band_smeta_ls[band_index].pixsize;
        *l1t_line = (smeta->projection.corners.upleft.y - proj_y)/smeta->band_smeta_ls[band_index].pixsize;
    }

    return(SUCCESS);
}

int smeta_init_projection(
    SMETA_SCENE_PROJ    proj )			/* Input enhanced metadata projection info   */
{
    int		proj_units = 2;			/* GCTP numeric code for input units         */
    int		geo_units = 0;			/* Output angular units of radians           */
    int		geo_code = 0;			/* Code for geographic output                */
    int		geo_zone = NULLZONE;		/* Null zone number for geographic output    */
    double	geo_parms[IAS_PROJ_PARAM_SIZE];	/* Zero parameter array for geographic output*/
    int		i;				/* Loop index                                */
    IAS_PROJECTION in_proj;	/* Source projection parameters */
    IAS_PROJECTION out_proj;	/* Destination projection parameters */

    /* Initialized the geographic projection parameters */
    for ( i=0; i<IAS_PROJ_PARAM_SIZE; i++ ) geo_parms[i] = 0.0;

    /* Set the input units code */
    if (!strcmp (proj.units, "RADIANS"))
        proj_units = 0;
    else if (!strcmp (proj.units, "FEET"))
        proj_units = 1;
    else if (!strcmp (proj.units, "METERS"))
        proj_units = 2;
    else if (!strcmp (proj.units, "SECONDS"))
        proj_units = 3;
    else if (!strcmp (proj.units, "DEGREES"))
        proj_units = 4;
    else if (!strcmp (proj.units, "DMS"))
        proj_units = 5;

    /* Load the projection structures */
    ias_geo_set_projection( geo_code, geo_zone, geo_units, proj.spheroid, geo_parms, &in_proj );
    ias_geo_set_projection( proj.code,  proj.zone,  proj_units,  proj.spheroid,  proj.projprms,  &out_proj );

    /* Create transformation */
    if ( (trans = ias_geo_create_proj_transformation( &in_proj, &out_proj )) == NULL )
    {
        IAS_LOG_ERROR("Initializing map projection transformation.");
        return(ERROR);
    }

    /* Store ellipsoid constants */
    wgs84_a = proj.wgs84_major_axis;
    wgs84_flat = (proj.wgs84_major_axis - proj.wgs84_minor_axis) / proj.wgs84_major_axis;
    wgs84_init = 1;

    return(SUCCESS);
}

int smeta_geodetic_to_proj(
    SMETA_SCENE_PROJ	proj,			/* Input enhanced metadata projection info   */
    double	lat,				/* Input latitude (radians)                  */
    double	lon, 				/* Input longitude (radians)                 */
    double	*proj_x,			/* Output projection X coordinate            */
    double	*proj_y )			/* Output projection Y coordinate            */
{
    /* Make sure the projection transformation has been initialized */
    if ( trans == NULL )
    {
        IAS_LOG_ERROR("Map projection transformation not initialized.");
        return(ERROR);
    }

    /* Convert the projection coordinates to latitude/longitude */
    if ( smeta_transform_projection( lon, lat, proj_x, proj_y ) != SUCCESS )
    {
        IAS_LOG_ERROR("Converting projection X/Y to lat/long.");
        return(ERROR);
    }

    return(SUCCESS);
}

int smeta_transform_projection(
    double		inx,		/* I: Input X projection coordinate            */
    double		iny,		/* I: Input Y projection coordinate            */
    double		*outx,		/* O: Output X projection coordinate           */
    double		*outy )		/* O: Output Y projection coordinate           */
{
    int			status;		/* Projection package return code              */

    /* Apply transformation */
    status = ias_geo_transform_coordinate( trans, inx, iny, outx, outy );

    return( status );
}

void smeta_release_projection()
{
    /* Release the transformation */
    ias_geo_destroy_proj_transformation( trans );

    trans = NULL;
    return;
}

int smeta_band_number_to_index(
    SMETA	*smeta,				/* Input metadata info                       */
    int		band )				/* Input band number                         */
{
    int		i;				/* Loop index                                */
    int		band_index;			/* Index corresponding to input band number  */

    /* Look for the input band number in the metadata band list */
    band_index = -1;
    for ( i=0; i<smeta->num_band; i++ )
    {
        if ( smeta->spacecraft_id == IAS_L8 )
        {
            if ( band == smeta->band_smeta[i].band )
            {
                band_index = i;
                break;
            }
        }
        else
        {
            if ( band == smeta->band_smeta_ls[i].band )
            {
                band_index = i;
                /* printf("in smeta_band_number_to_index: band = %d, band_index = %d\n", band, band_index);*/
                break;
            }
        }
    }

    return(band_index);
}

