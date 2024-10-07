/*
 *	Name:	smeta_geometry.c
 *
 *	Purpose:
 *	Source file containing the functions used to relate Level 1T pixel line/sample
 *	coordinates to the corresponding viewing and solar illumination angles. These functions
 *	operate on the SMETA data type defined in "smeta.h".
 */
#ifndef _SMETA_GEOMETRY_H_
#define _SMETA_GEOMETRY_H_

#include "ias_math.h"				/* IAS_VECTOR declaration                    */
#include "smeta.h"				/* Definition of SMETA structure             */

int smeta_geodetic_to_ecef(
    double	lat,				/* Input latitude in radians                 */
    double	lon,				/* Input longitude in radians                */
    double	hgt,				/* Input height in meters                    */
    IAS_VECTOR	*ecef );			/* Output ECEF vector                        */

int smeta_geodetic_to_ecf2lsr(
    double	lat,				/* Input latitude in radians                 */
    double	lon,				/* Input longitude in radians                */
    double	ecf2lsr[3][3] );	 	/* ECEF to LSR rotation matrix               */

int smeta_proj_to_l1t(
    SMETA	*smeta,				/* Input enhanced metadata info              */
    int		band,				/* Input band number                         */
    double	proj_x,				/* Input projection X coordinate             */
    double	proj_y,  			/* Input projection Y coordinate             */
    double	*l1t_line,			/* Output L1T line number                    */
    double	*l1t_samp );			/* Output L1T sample number                  */

int smeta_geodetic_to_proj(
    SMETA_SCENE_PROJ	proj,			/* Input enhanced metadata projection info   */
    double	lat,				/* Input latitude (radians)                  */
    double	lon,  				/* Input longitude (radians)                 */
    double	*proj_x,			/* Output projection X coordinate            */
    double	*proj_y );			/* Output projection Y coordinate            */

int smeta_init_projection(
    SMETA_SCENE_PROJ    proj );			/* Input enhanced metadata projection info   */

int smeta_transform_projection(
    double		inx,			/* I: Input X projection coordinate            */
    double		iny,			/* I: Input Y projection coordinate            */
    double		*outx,			/* O: Output X projection coordinate           */
    double		*outy );		/* O: Output Y projection coordinate           */

void smeta_release_projection();

int smeta_band_number_to_index(
    SMETA		*smeta,			/* Input enhanced metadata info                */
    int			band );			/* Input band number                           */

int smeta_calc_vectors(
    SMETA		*smeta,			/* Enhanced metadata                           */
    int			band,			/* Current band number                         */
    double		l1t_line,		/* Current L1T line number                     */
    double		l1t_samp,		/* Current L1T sample number                   */
    IAS_VECTOR		*satvec,		/* Satellite view vector                       */
    IAS_VECTOR		*sunvec );		/* Solar illumination vector                   */

#endif
