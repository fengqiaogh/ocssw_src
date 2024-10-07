/*
 *	Name:	emeta_geometry.c
 *
 *	Purpose:
 *	Source file containing the functions used to relate Level 1T pixel line/sample
 *	coordinates to the corresponding viewing and solar illumination angles. These functions
 *	operate on the EMETA data type defined in "emeta.h".
 */
#ifndef _EMETA_GEOMETRY_H_
#define _EMETA_GEOMETRY_H_

#include "emeta.h"				/* Definition of EMETA structure             */

int emeta_l1t_to_l1r(
    EMETA_BAND	*emeta_band,			/* Enhanced metadata band structure          */
    double	l1t_line,			/* Input L1T line coordinate                 */
    double	l1t_samp,			/* Input L1T sample coordinate               */
    double	height,				/* Input L1T height coordinate               */
    int		*sca_list,			/* Array of SCA numbers containing point     */
    double	*l1r_line,			/* Array of L1R line coordinates by SCA      */
    double	*l1r_samp );			/* Array of L1R sample coordinates by SCA    */

int emeta_band_number_to_index(
    EMETA		*emeta,			/* Input enhanced metadata info                */
    int			band );			/* Input band number                           */
	
#endif
