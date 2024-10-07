/*
 *	Name:	emeta_geometry.c
 *
 *	Purpose:
 *	Source file containing the functions used to relate Level 1T pixel line/sample
 *	coordinates to the corresponding viewing and solar illumination angles. These functions
 *	operate on the EMETA data type defined in "emeta.h".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "ias_logging.h"
#include "ias_const.h"
#include "emeta_geometry.h"			/* Definition of EMETA structure and
						   function prototypes.                      */
#define	PIX_BUFFER	3			/* Resampling kernel reduction in effective SCA size */

int emeta_l1t_to_l1r(
    EMETA_BAND	*emeta_band,			/* Enhanced metadata band structure          */
    double	l1t_line,			/* Input L1T line coordinate                 */
    double	l1t_samp,			/* Input L1T sample coordinate               */
    double	height,				/* Input L1T height coordinate               */
    int		*sca_list,			/* Array of SCA numbers containing point     */
    double	*l1r_line,			/* Array of L1R line coordinates by SCA      */
    double	*l1r_samp )			/* Array of L1R sample coordinates by SCA    */
{
    int		nsca;				/* Number of SCAs covering the point (0-2)   */
    int		isca;				/* SCA index                                 */
    double	l1t_l;				/* Offset L1T line coordinate                */
    double	l1t_s;				/* Offset L1T sample coordinate              */
    double 	l1t_h;				/* Offset L1T height coordinate              */
    double	l1r_l;				/* RPC value of current L1R line             */
    double	l1r_s;				/* RPC value of current L1R sample           */

    /* Loop through the SCAs */
    nsca = 0;
    for ( isca=0; isca<emeta_band->nsca; isca++ )
    {
        /* Evaluate the RPCs for the current SCA */
        l1t_l = l1t_line - emeta_band->rpc[isca].mean_l1t[0];
        l1t_s = l1t_samp - emeta_band->rpc[isca].mean_l1t[1];
        l1t_h = height - emeta_band->rpc[isca].mean_hgt;
        l1r_l = (emeta_band->rpc[isca].line_num[0] + emeta_band->rpc[isca].line_num[1]*l1t_l
              + emeta_band->rpc[isca].line_num[2]*l1t_s + emeta_band->rpc[isca].line_num[3]*l1t_h
              + emeta_band->rpc[isca].line_num[4]*l1t_l*l1t_s)
              / (1.0 + emeta_band->rpc[isca].line_den[0]*l1t_l + emeta_band->rpc[isca].line_den[1]*l1t_s
              +  emeta_band->rpc[isca].line_den[2]*l1t_h + emeta_band->rpc[isca].line_den[3]*l1t_l*l1t_s)
              + emeta_band->rpc[isca].mean_l1r[0];
        l1r_s = (emeta_band->rpc[isca].samp_num[0] + emeta_band->rpc[isca].samp_num[1]*l1t_l
              + emeta_band->rpc[isca].samp_num[2]*l1t_s + emeta_band->rpc[isca].samp_num[3]*l1t_h
              + emeta_band->rpc[isca].samp_num[4]*l1t_l*l1t_s)
              / (1.0 + emeta_band->rpc[isca].samp_den[0]*l1t_l + emeta_band->rpc[isca].samp_den[1]*l1t_s
              +  emeta_band->rpc[isca].samp_den[2]*l1t_h + emeta_band->rpc[isca].samp_den[3]*l1t_l*l1t_s)
              + emeta_band->rpc[isca].mean_l1r[1];
        /* See if this is a valid point for the current SCA */
        if ( l1r_l >= PIX_BUFFER && l1r_l <= (emeta_band->l1r_lines-1-PIX_BUFFER) &&
             l1r_s >= PIX_BUFFER && l1r_s <= (emeta_band->l1r_samps-1-PIX_BUFFER) )
        {
              if ( nsca > 1 )
              {
                  IAS_LOG_ERROR("Too many valid SCAs found.");
                  return(0);
              }
              l1r_line[nsca] = l1r_l;
              l1r_samp[nsca] = l1r_s;
              sca_list[nsca] = emeta_band->rpc[isca].sca_num;
              nsca++;
        }
    }

    return(nsca);
}

int emeta_band_number_to_index(
    EMETA	*emeta,				/* Input enhanced metadata info              */
    int		band )				/* Input band number                         */
{
    int		i;				/* Loop index                                */
    int		band_index;			/* Index corresponding to input band number  */

    /* Look for the input band number in the metadata band list */
    band_index = -1;
    for ( i=0; i<emeta->num_band; i++ )
    {
        if ( band == emeta->band_emeta[i].band )
        {
            band_index = i;
            break;
        }
    }

    return(band_index);
}

