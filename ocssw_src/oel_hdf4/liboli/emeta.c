/*
 *	Name:	emeta.c
 *
 *	Purpose:
 *	Source file containing the functions used to allocate and release memory in the
 *	enhanced metadata structure.
 */
#include <stdlib.h>
#include "ias_logging.h"
#include "emeta.h"

int emeta_allocate_ephemeris(
    EMETA	*emeta,				/* Enhanced metadata structure               */
    int		nephem )			/* Number of ephemeris points required       */
{
    EMETA_EPHEM	*local_ephem = NULL;		/* Local pointer to allocated buffer         */

    /* Check the number of points requested */
    if ( nephem < 1 )
    {
        IAS_LOG_ERROR("Invalid number of ephemeris samples requested.");
        return(ERROR);
    }

    /* Allocate the ephemeris array */
    local_ephem = (EMETA_EPHEM *)malloc( nephem * sizeof( EMETA_EPHEM ) );
    if ( local_ephem == NULL )
    {
        IAS_LOG_ERROR("Unable to allocate enhanced metadata ephemeris structure.");
        return(ERROR);
    }

    /* Set the ephemeris pointer in the EMETA structure */
    emeta->ephemeris = local_ephem;

    /* Allocate the sunvector array */
    local_ephem = (EMETA_EPHEM *)malloc( nephem * sizeof( EMETA_EPHEM ) );
    if ( local_ephem == NULL )
    {
        IAS_LOG_ERROR("Unable to allocate enhanced metadata solar ECEF vector structure.");
        free( emeta->ephemeris );
        return(ERROR);
    }

    /* Set the sunvector pointer in the EMETA structure */
    emeta->sunvector = local_ephem;

    return(SUCCESS);
}

void emeta_free_ephemeris(
    EMETA	*emeta )			/* Enhanced metadata structure               */
{
    free( emeta->ephemeris );
    free( emeta->sunvector );
    return;
}

