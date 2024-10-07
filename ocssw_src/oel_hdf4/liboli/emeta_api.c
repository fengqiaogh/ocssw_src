#include "emeta_api.h"        		/* EMETA structure and routines   */

static EMETA		emeta;		/* Enhanced metadata structure    */

int emeta_init( 			/* Returns num bands in EMETA     */
    char	*emeta_filename )	/* Enhanced metadata file name    */
{
    /* Read the enhanced metadata file */
    if ( emeta_read( &emeta, emeta_filename ) != SUCCESS )
    {
        IAS_LOG_ERROR("Reading the enhanced metadata file %s.", emeta_filename);
        return(0);
    }

    return(emeta.num_band);
}

int emeta_get_frame(
    int         eband,			/* Input band index        */
    EMETA_FRAME *frame )		/* Output image frame info */
{

    if ( eband > emeta.num_band-1 )
    {
        IAS_LOG_ERROR("Invalid band requested");
        return(ERROR);
    }

    /* Establish the output image file frame */
    frame->band = emeta.band_emeta[eband].band;
    frame->nlines = emeta.band_emeta[eband].l1t_lines;
    frame->nsamps = emeta.band_emeta[eband].l1t_samps;
    frame->code = emeta.projection.code;
    frame->zone = emeta.projection.zone;
    frame->pixsize = emeta.band_emeta[eband].pixsize;
    frame->ul_x = emeta.projection.corners.upleft.x;
    frame->ul_y = emeta.projection.corners.upleft.y;

    return(SUCCESS);
}

int emeta_calc(
    int		in_line,		/* L1T line coordinate              */
    int		in_samp,		/* L1T sample coordinate            */
    int		in_band,		/* Spectral band number             */
    short	*sat_zn,		/* Satellite zenith angle           */
    short	*sat_az,		/* Satellite azimuth angle          */
    short	*sun_zn,		/* Solar zenith angle               */
    short	*sun_az )		/* Solar azimuth angle              */
{
    double	l1t_line;		/* Float line coordinate            */
    double	l1t_samp;		/* Float sample coordinate          */
    double	satang[2];		/* Satellite viewing angles         */
    double	sunang[2];		/* Solar viewing angles             */
    double	r2d;			/* Radian to degree conversion      */
    int		eband;			/* Band index                       */
    int		band;			/* Band number                      */
    int		status;			/* Status return flag               */

    /* Initialize conversion factor from radians to hundredths of degrees */
    r2d = 4500.0 / atan(1.0);
    status = ERROR;

    /* Process the angles for each band */
    for ( eband=0; eband<emeta.num_band; eband++ )
    {
        band = emeta.band_emeta[eband].band;
        if ( band != in_band ) continue;
        status = SUCCESS;

        /* Calculate the angles for the current line/sample */
        l1t_line = (double)in_line;
        l1t_samp = (double)in_samp;

	/* Calculate the viewing geometry using the RPC method */
        status = emeta_angles_rpc( l1t_line, l1t_samp, &emeta, band, satang, sunang );
        if ( status != SUCCESS ) IAS_LOG_ERROR("Evaluating angles in band %d.", band);

        /* Quantize the satellite angles */
        *sat_zn = (short)floor( r2d*satang[0] + 0.5 );
        *sat_az = (short)floor( r2d*satang[1] + 0.5 );
        /* Quantize the solar angles */
        *sun_zn = (short)floor( r2d*sunang[0] + 0.5 );
        *sun_az = (short)floor( r2d*sunang[1] + 0.5 );
    }

    /* Return status */
    return( status );
}

void emeta_close()
{
    /* Free the ephemeris structure */
    emeta_free_ephemeris( &emeta );

    return;
}
