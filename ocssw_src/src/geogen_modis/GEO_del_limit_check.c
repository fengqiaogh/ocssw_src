/* file: GEO_del_limit_check.c */

#include "PGS_MODIS_35251.h"
#include "smfio.h"
#include "GEO_validation.h"

int GEO_del_limit_check(
	double data_samples[],
	int const number_of_samples,
	double const delta_limit,
	int sample_flags[]
	)
/*
!C****************************************************************************

!Description:   
	Routine in Input group of the Level-1A geolocation software to validate
	a set of data samples by comparing differences between successive
	samples against a limit.  Previously flagged values are not checked.
	Flags set to BAD_DATA for all samples not within limits.             

!Input Parameters:
	data_samples		array of input data samples
	number_of_samples	number of samples to check
	delta_limit		difference limit

!Output Parameters:

!Input/Output Parameters:
	sample_flags		array of validation flags with previously
				invalidated samples flagged 

Return Values:
	SUCCESS - if it can find at least one valid pair of values.
	FAIL	- otherwise

Externally Defined:
	BAD_DATA		"GEO_validation.h"
	FAIL			"GEO_basic.h"
	MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
	MODIS_E_ONE_UNFLAG	"PGS_MODIS_35251.h"
	SUCCESS			"GEO_basic.h"

Called by:
	GEO_prepare_mirr_data
	GEO_validate_attit_data

Routines Called:
	GEO_find_next_flag
	modsmf

!Revision History:
$Log: GEO_del_limit_check.c,v $
Revision 4.2  2003/08/28 16:01:45  kuyper
Corrected prolog.

Revision 4.1  2003/02/21 22:53:21  kuyper
Corrected to use void* pointers with %p format code.

Revision 3.2  2001/03/16 20:55:05  kuyper
Changed to not return FAIL for number_of_samples = 0

Revision 3.1  2001/03/07  19:33:25  kuyper
Major rewrite, to mark both sides of a sudden jump as bad.

James.R.Kuyper.1@gsfc.nasa.gov

Requirements:
	PR03-F-2.4-1
	PR03-F-2.4-2
	PR03-F-2.5-1
	PR03-F-2.5-2
	PR03-F-2.5-3
	PR03-I-1
	PR03-I-2
	PR03-S-1

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes

!END*************************************************************************
 */
{

    /* Local variable declarations */

    double del = 0.0; /* absolute value of difference */
    double limit = 0.0; /* limit value */
    int num_valid = 0; /* number of valid data samples found */
    int f1 = -1; /* array index */
    int f2 = -1; /* array index */
    char msgbuf[64];

    if(data_samples==NULL || sample_flags==NULL)
    {
	sprintf(msgbuf, "data_samples:%p sample_flags:%p",
	    (void*)data_samples, (void*)sample_flags);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, __FILE__ ", GEO_del_limit_check");

	return FAIL;
    }
    /* Begin program logic */
    while( (f2=GEO_find_next_flag(sample_flags,number_of_samples,f2+1))>=0)
    {
	if(f1 < 0) /* First time	*/
	    num_valid++;
	else
	{
	    del = fabs(data_samples[f2] - data_samples[f1]);
	    limit = delta_limit * (double)(f2 - f1);
	    if (del > limit) {
		if(sample_flags[f1] != BAD_DATA)
		{
		    num_valid--;	/* We counted it as valid earlier.*/
		    sample_flags[f1] = BAD_DATA;
		}

		sample_flags[f2] = BAD_DATA;
	    }
	    else
	      num_valid++;
	}

	f1 = f2;
    }

    if (num_valid ==0 && number_of_samples > 0)
	return FAIL;
      
    if (num_valid == 1)
	modsmf(MODIS_E_ONE_UNFLAG, "", __FILE__ ", GEO_del_limit_check");

    return SUCCESS;
}

