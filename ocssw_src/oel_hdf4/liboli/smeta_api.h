#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "ias_math.h"
#include "ias_logging.h"          /* IAS logging library                        */
#include "mtl_geometry.h"	  /* MTL file geometric operations		*/
#include "mtl_grid.h"	          /* MTL grid generation functions		*/
#include "smeta_exploit.h"        /* EMETA file read routine                    */
#include "smeta_geometry.h"       /* Geometric operations using EMETA structure */

#ifndef _SMETA_API_
#define _SMETA_API_

/* Initialize Metadata Access */
int smeta_init( 
    char	*smeta_filename );	/* Product metadata file name       */

/* Connect to scene ID text string */
char *smeta_file_name();

/* Retrieve scene framing information from the SMETA structure */
int smeta_get_frame(
    int		band,                   /* Input band number                */
    META_FRAME	*frame );               /* Output image frame info          */


/* Calculate view and sun angles for a L1T pixel for L5 and L7 satellites*/
int smeta_calc_ls(
    int		in_line,		/* L1T line coordinate              */
    int		in_samp,		/* L1T sample coordinate            */
    int		in_band,		/* Spectral band number             */
    short	*sat_zn,		/* Satellite zenith angle           */
    short	*sat_az,		/* Satellite azimuth angle          */
    short	*sun_zn,		/* Solar zenith angle               */
    short	*sun_az );		/* Solar azimuth angle              */

/* Calculate view and sun angles for a L1T pixel for OLI */
int smeta_calc(
    int		in_line,		/* L1T line coordinate              */
    int		in_samp,		/* L1T sample coordinate            */
    int		in_band,		/* Spectral band number             */
    /* Sudipta new addition for OLI to extract SCA and Det numbers */
    short	*sca_num, 		/* SCA number                       */
    short	*det_num,		/* Detector number                  */
    /* End Sudipta new addition for OLI */
    double	*sat_zn,		/* Satellite zenith angle           */
    double	*sat_az,		/* Satellite azimuth angle          */
    double	*sun_zn,		/* Solar zenith angle               */
    double	*sun_az );		/* Solar azimuth angle              */

/* Disconnect from the MTL interface */
void smeta_close();

#endif
