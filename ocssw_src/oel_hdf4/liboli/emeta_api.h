#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "ias_logging.h"          /* IAS logging library                        */
#include "emeta_exploit.h"        /* EMETA file read routine                    */
#include "emeta_geometry.h"       /* Geometric operations using EMETA structure */

#ifndef _EMETA_API_
#define _EMETA_API_

/* Initialize Enhanced Metadata Access */
int emeta_init( 
    char	*emeta_filename );	/* Enhanced metadata file name      */

/* Retrieve scene framing information from the EMETA structure */
int emeta_get_frame(
    int		band,                   /* Input band number                */
    EMETA_FRAME	*frame );               /* Output image frame info          */


/* Calculate view and sun angles for a L1T pixel */
int emeta_calc(
    int		in_line,		/* L1T line coordinate              */
    int		in_samp,		/* L1T sample coordinate            */
    int		in_band,		/* Spectral band number             */
    short	*sat_zn,		/* Satellite zenith angle           */
    short	*sat_az,		/* Satellite azimuth angle          */
    short	*sun_zn,		/* Solar zenith angle               */
    short	*sun_az );		/* Solar azimuth angle              */

/* Disconnect from the EMD interface */
void emeta_close();

#endif
