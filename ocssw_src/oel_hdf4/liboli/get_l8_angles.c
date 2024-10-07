#include "emeta_api.h"       /* Enhanced metadata functionality */
#define MAXBAND 7
/**
 * Get the Landsat OLI angles from an Enhanced Meta Data (EMD) file
 * @param in input EMD file
 * @param in input number of pixels per line
 * @param in input number of lines
 * @param in input scan line to work on
 * @param out output satellite zenith angle
 * @param out output satellite azimuth angle
 * @param out output solar zenith angle
 * @param out output solar azimuth angle
 */
int get_oli_angles( char *emeta_filename, int32_t npix, int32_t nscan, int32_t iscan,
        float *solz, float *sola, float *senz, float *sena,
        int32_t isGeomPerBand)
{
    static int  firstCall = 1;
    static int  nband;                     /* Number of bands                */
    static EMETA_FRAME       frame[MAXBAND];         /* Output image frame info        */
    int         band;                       /* User band number               */
    int         band_index;                 /* Band index                     */
    int         is;                         /* Sample (pixel) index           */
    int         status;                     /* Status return flag             */
    short       *sat_zn;                    /* Satellite zenith angle         */
    short       *sat_az;                    /* Satellite azimuth angle        */
    short       *sun_zn;                    /* Solar zenith angle             */
    short       *sun_az;                    /* Solar azimuth angle            */
    int         isb;                        /* index counter for band based arrays */

    if (firstCall == 1) {

        /* Initialize the enhanced metadata package */
        if ( (nband = emeta_init( emeta_filename )) < 1 )
        {
            IAS_LOG_ERROR("Initializing the enhanced metadata file %s.", emeta_filename);
            fprintf(stderr,
                    "-E- %s line %d: Initializing the metadata file %s.\n",
                    __FILE__, __LINE__, emeta_filename);
            exit (EXIT_FAILURE);
        }
    }
    /* Allocate the intermediary buffers for calculating average over MAXBAND*/
    if ( (sat_zn = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        emeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate satellite elevation angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    if ( (sat_az = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        free(sat_zn);
        emeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate satellite azimuth angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    if ( (sun_zn = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        free(sat_zn);
        free(sat_az);
        emeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate solar elevation angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    if ( (sun_az = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        free(sat_zn);
        free(sat_az);
        free(sun_zn);
        emeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate solar azimuth angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    /* Initialize nominal values */
    for ( is=0; is<npix; is++ )
    {
        senz[is] = 0;
        sena[is] = 0;
        solz[is] = 0;
        sola[is] = 0;
    }

    /* Process the angles for each band */
    for ( band_index=0; band_index<MAXBAND; band_index++ )
    {
        /* Get framing information for this band */
        if ( firstCall == 1 && emeta_get_frame( band_index, &frame[band_index] ) != SUCCESS )
        {
            fprintf(stderr,
                    "-E- %s line %d: Unable to get emeta frame.\n",
                    __FILE__, __LINE__);
            exit (EXIT_FAILURE);
        }
        band = frame[band_index].band;

        /* Loop through the L1T  samples */
        for ( is=0; is<npix; is++ )
        {
            isb = is*MAXBAND+band_index;

            status = emeta_calc( (int)iscan, is, band, &sat_zn[is], &sat_az[is], &sun_zn[is], &sun_az[is] );
            if ( status != SUCCESS )
            {
                fprintf(stderr,
                        "-E- %s line %d: Error Evaluating angles in band %d.\n",
                        __FILE__, __LINE__,band);
                exit (EXIT_FAILURE);
            }
            // Populate the nominal values as band averages
            /* Average the bands and convert to decimal degrees */
            if (isGeomPerBand)
            {
                 senz[isb] = (float)sat_zn[is]/100.;
                 sena[isb] = (float)sat_az[is]/100.;
                 solz[isb] = (float)sun_zn[is]/100.;
                 sola[isb] = (float)sun_az[is]/100.;
            } else {
                senz[is] += (float) sat_zn[is] / (100 * MAXBAND);
                sena[is] += (float) sat_az[is] / (100 * MAXBAND);
                solz[is] += (float) sun_zn[is] / (100 * MAXBAND);
                sola[is] += (float) sun_az[is] / (100 * MAXBAND);
            }
        }
     }

    /* Free the angle arrays */
    free( sat_zn );
    free( sat_az );
    free( sun_zn );
    free( sun_az );

    /* Release the enhanced metadata */
    //emeta_close(); /* Probably should close this properly somewhere after processing */
    firstCall = 0;
    return (SUCCESS);
}
