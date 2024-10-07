#include "smeta_api.h"       /* Product metadata functionality */
#define MAXBAND 7
/**
 * Get the Landsat Nominal OLI angles from a MTL file
 * @param in input MTL file
 * @param in input number of pixels per line
 * @param in input number of lines
 * @param in input scan line to work on
 * @param out output satellite zenith angle
 * @param out output satellite azimuth angle
 * @param out output solar zenith angle
 * @param out output solar azimuth angle
 */
int get_l57tm_nom_angles( char *meta_filename, int32_t npix, int32_t nscan, int32_t iscan,
        float *solz, float *sola, float *senz, float *sena)
{
    static int  firstCall = 1;
    static int  nband;                     /* Number of bands                */
    static META_FRAME   frame[MAXBAND];    /* Output image frame info        */
    int         band;	  			        /* User band number               */
    int			band_index;			        /* Band index                     */
    int         is;                         /* Sample (pixel) index           */
    int			status;				        /* Status return flag             */
    short		*sat_zn;			        /* Satellite zenith angle         */
    short		*sat_az;			        /* Satellite azimuth angle        */
    short		*sun_zn;			        /* Solar zenith angle             */
    short		*sun_az;			        /* Solar azimuth angle            */

    if (firstCall == 1) {
        /* Initialize the metadata package */
        if ( (nband = smeta_init( meta_filename )) < 1 )
        {
            IAS_LOG_ERROR("Initializing the metadata file %s.", meta_filename);
            fprintf(stderr,
                    "-E- %s line %d: Initializing the metadata file %s.\n",
                    __FILE__, __LINE__, meta_filename);
            exit (EXIT_FAILURE);
        }
    }
    
    
    /* printf("nband = %d\n", nband); we get nband = 7 */
    
    /* Allocate the intermediary buffers for calculating average over MAXBAND*/
    if ( (sat_zn = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        smeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate satellite elevation angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    if ( (sat_az = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        free(sat_zn);
        smeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate satellite azimuth angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }
    if ( (sun_zn = (short *)malloc( npix*sizeof(short) )) == NULL )
    {
        free(sat_zn);
        free(sat_az);
        smeta_close();
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
        smeta_close();
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate solar azimuth angle array.\n",
                __FILE__, __LINE__);
        exit (EXIT_FAILURE);
    }

    /* Initialize  */
    for ( is=0; is<npix; is++ )
    {
        senz[is] = 0;
        sena[is] = 0;
        solz[is] = 0;
        sola[is] = 0;
    }

    /* Process the angles for each band */
    //printf("nbnd in get_l57tm_nom_angles = %d\n", nband);
    for ( band_index=0; band_index<nband; band_index++ )
    {
        /* Get framing information for this band */
        /* for Landsat 5 we need to skip band 6
           for labndsat 7 we need to skip band 6 and 8 */
        if ((band_index == 5) || (band_index == 7))
            continue;

        if ( firstCall == 1 && smeta_get_frame( band_index, &frame[band_index] ) != SUCCESS )
        {
            fprintf(stderr,
                    "-E- %s line %d: Unable to get smeta frame.\n",
                    __FILE__, __LINE__);
            exit (EXIT_FAILURE);
        }
        band = frame[band_index].band;
        /* } */
        
        
        //printf("band_index = %d, band = %d\n", band_index, band);
        
        //if ( (iscan) % 500 == 0 ) printf("GET_L57TM_NOM_ANGLES: Band/Line %d/%d\n", band, iscan );
        /* Loop through the L1T lines and samples */
            for ( is=0; is<npix; is++ )
            {
                status = smeta_calc_ls( (int) iscan, is, band, &sat_zn[is], &sat_az[is], &sun_zn[is], &sun_az[is] );
                if ( status != SUCCESS )
                {
                    fprintf(stderr,
                            "-E- %s line %d: Error Evaluating angles in band %d.\n",
                            __FILE__, __LINE__,band);
                    exit (EXIT_FAILURE);
                }
                /* Average the bands and convert to decimal degrees */
                //printf("band=%d sat_zn[%d]=%d, sun_zn = %d\n",band, is,sat_zn[is], sun_zn[is]);
                senz[is] += (float)sat_zn[is]/(100*(MAXBAND-1));
                sena[is] += (float)sat_az[is]/(100*(MAXBAND-1));
                solz[is] += (float)sun_zn[is]/(100*(MAXBAND-1));
                sola[is] += (float)sun_az[is]/(100*(MAXBAND-1));
            }
    }
        
    /* for ( is=0; is<npix; is+=10)
    {
        printf("senz[%d]=%f\n",is,senz[is]);
        printf("sena[%d]=%f\n",is,sena[is]);
        printf("solz[%d]=%f\n",is,solz[is]);
        printf("sola[%d]=%f\n",is,sola[is]);
    } */
        
    /* Free the angle arrays */
    free( sat_zn );
    free( sat_az );
    free( sun_zn );
    free( sun_az );

    /* Release the metadata */
    //smeta_close();  /* Probably should close this properly somewhere after processing */
    firstCall = 0;
    return (SUCCESS);
}
