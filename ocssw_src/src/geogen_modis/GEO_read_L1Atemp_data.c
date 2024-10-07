#include "PGS_MODIS_35251.h"
#include "GEO_input.h"
#include "GEO_geo.h"
#include "GEO_parameters.h"
#include "mapiL1A.h"
#include "mapi.h"
#include "smfio.h"

PGSt_SMF_status GEO_read_L1Atemp_data(
	const GEO_param_struct	* const geo_params,
	MODFILE			* const l1a_file,
	int			const num_scans,
	float32			average_temp[NUM_TEMPS])

/******************************************************************************
!C

!Description:   
	Subtroutine in the input group of the Level-1A geolocation software to
	read in thermometer readings, convert them to temperature values, and
	average them over the period of the granule.

!Input Parameters:
	geo_params		Geolocation parameters struct
	l1a_file		M-API file handle for the L1A file
	num_scans		The number of scans in the granule

!Output Parameters:
	average_temp		The average over the granule of each temperature

Return Values:
	MODIS_E_BAD_INPUT_ARG	If any input argument was invalid.
	MODIS_E_GEO		If any subroutine failed.
	PGS_S_SUCCESS		Otherwise

Externally Defined:
	M01CR_CP_A_ON_M		"mapiL1A.h"
	M01CR_CP_B_ON_M		"mapiL1A.h"
	M01CR_RC_CFPA_T1SET	"mapiL1A.h"
	M01CR_RC_CFPA_T3SET	"mapiL1A.h"
	M01MAJCYC1OF7		"mapiL1A.h"
	M01MAJCYC5AOF7		"mapiL1A.h"
	M01MAJCYCALL1		"mapiL1A.h"
	M01MAJCYC0OF63		"mapiL1A.h"
	M01MAJCYC3OF63		"mapiL1A.h"
	M01MAJCYC4OF63		"mapiL1A.h"
	M01MAJCYC21OF63		"mapiL1A.h"
	M01MAJCYC23OF63		"mapiL1A.h"
	M01TA_RC_SMIR_CFPA	"mapiL1A.h"
	M01TP_AO_SMIR_OBJ	"mapiL1A.h"
	M01TP_MF_CALBKHD_SR	"mapiL1A.h"
	M01TP_MF_Z_BKHD_BB	"mapiL1A.h"
	M01TP_SA_RCT1_MIR	"mapiL1A.h"
	M01TP_SR_SNOUT		"mapiL1A.h"
	MAPIOK			"mapi.h"
	MAX_SCAN_NUMBER		"GEO_geo.h"
	MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
	MODIS_E_GEO		"PGS_MODIS_35251.h"
	NUM_TEMPS		"GEO_parameters.h"
	PGS_S_SUCCESS		"PGS_SMF.h"
	TEMP_FVALUE		"GEO_geo.h"
	TEMP_ORDER		"GEO_parameters.h"

Called by:
	GEO_prepare_l1a_data

Routines Called:
	getMODIStable		"mapi.h"
	modsmf			"smfio.h"

!Revision History:
        $Log: GEO_read_L1Atemp_data.c,v $
        Revision 4.2  2003/03/17 18:53:23  vlin
        Updated after code walkthrough

        Revision 4.1  2003/02/11 15:15:12  vlin
        vlin@saicmodis.com

Requirements:      None

!Team-unique Header:

	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

!END
******************************************************************************/

#define NUM_CONFIGS 4                   /* number of configuration fields */

{
  uint16    readings[NUM_TEMPS+NUM_CONFIGS][MAX_SCAN_NUMBER];
                                        /* the thermometer readings */
  int       eqn;                        /* used to identify which equation in 
                                           each section should be used */
  const int base[NUM_TEMPS] = {0, 3, 5, 5, 3, 5};
                                        /* used to identify the index for the 
                                           first equation in the equation list 
                                           that applies for each temperature */
  int       field, scan, temp, term;
  double    temperature, total[NUM_TEMPS] = {0.0};
  int       count[NUM_TEMPS] = {0};
  long int  buffsize, tstart = 0L;
  char      msg[128];
  char      filefunc[] = __FILE__ ", GEO_read_L1Atemp_data";

  struct field{
       char *vdata_name;
       char *field_name;
  } field_list[] = {
     {M01MAJCYC1OF7,  M01CR_CP_A_ON_M},
     {M01MAJCYC1OF7,  M01CR_CP_B_ON_M},
     {M01MAJCYC5AOF7, M01CR_RC_CFPA_T1SET},
     {M01MAJCYC5AOF7, M01CR_RC_CFPA_T3SET},
     {M01MAJCYCALL1,  M01TA_RC_SMIR_CFPA},
     {M01MAJCYC0OF63, M01TP_AO_SMIR_OBJ},
     {M01MAJCYC3OF63, M01TP_MF_CALBKHD_SR},
     {M01MAJCYC4OF63, M01TP_MF_Z_BKHD_BB},
     {M01MAJCYC21OF63,M01TP_SA_RCT1_MIR},
     {M01MAJCYC23OF63,M01TP_SR_SNOUT}
  };

  if (l1a_file == NULL || geo_params == NULL || average_temp == NULL ||
      num_scans < 0 || num_scans > MAX_SCAN_NUMBER) {
      sprintf(msg,"l1a_file: %p, geo_params: %p, num_scans: %d, average_temp: "
              "%p", (void *)l1a_file, (void *)geo_params, num_scans, 
             (void *)average_temp);
      modsmf(MODIS_E_BAD_INPUT_ARG, msg, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  for (temp = 0; temp < NUM_TEMPS; temp++) 
      average_temp[temp] = TEMP_FVALUE;

  if (num_scans == 0)
      return PGS_S_SUCCESS;

  for (field = 0; field < NUM_CONFIGS+NUM_TEMPS; field++){
     buffsize = (long)sizeof(readings[field]);
     if (getMODIStable(l1a_file, field_list[field].vdata_name, NULL,
         field_list[field].field_name, tstart, (long)num_scans, &buffsize,
         (unsigned char *)readings[field]) != MAPIOK) {
         sprintf(msg,"getMODIStable(%s, %s)",field_list[field].vdata_name,
                 field_list[field].field_name);
         modsmf(MODIS_E_GEO, msg, filefunc);
         return MODIS_E_GEO;
     }
  }

  for (scan = 0; scan < num_scans ; scan++) {
    /* TA_RC_SMIR_CFPA: an active sensor, with different calibrations
                 depending upon which heaters have been turned on. */ 
     if (readings[2][scan] == 1 && readings[3][scan] == 0)
	eqn = 0;              /* 83 K settings: TA07p83 or TA07f83 */	
     else if (readings[2][scan] == 0 && readings[3][scan] == 0)
        eqn = 1;              /* 85 K settings: TA07p85 or TA07f85 */
     else if (readings[2][scan] == 0 && readings[3][scan] == 1)
        eqn = 2;              /* 88 K settings: TA07p88 or TA07f88 */
     else
        eqn = -1;             /* Invalid combination. */

    /* NUM_CONFIGS is the offset of the first non-configuration reading. */
     if (eqn == -1 || readings[NUM_CONFIGS][scan] == 0xFFFF)
	 temperature = TEMP_FVALUE;
     else {
         temperature = geo_params->temp_coeff[base[0]+eqn][TEMP_ORDER-1];
         for (term = TEMP_ORDER-2; term >= 0; term--)
              temperature = temperature*readings[NUM_CONFIGS][scan] +
                            geo_params->temp_coeff[base[0]+eqn][term];
         total[0] += temperature;
	 count[0]++;
     }

    /* The same equations are used for both electronics sides. */
     if (strcmp(geo_params->spacecraft_ID, "Terra") == 0)
         eqn = 0;              /* TP01p or TP02p */
     else if (readings[0][scan] == 0 && readings[1][scan] == 1)
         eqn = 0;              /* TP01fa or TP02fa */
     else if (readings[0][scan] == 1 && readings[1][scan] == 0)
         eqn = 1;              /* TP01fb or TP02fb */
     else
         eqn = -1;

     if (eqn == -1)
         continue;

     for (temp = 1; temp < NUM_TEMPS; temp++) {
         if (readings[temp][scan] == 0xFFFF)
             continue;
         temperature = geo_params->temp_coeff[base[temp]+eqn][TEMP_ORDER-1];
         for (term = TEMP_ORDER-2; term >= 0; term--)
              temperature = temperature*readings[NUM_CONFIGS+temp][scan] +
                            geo_params->temp_coeff[base[temp]+eqn][term];
         total[temp] += temperature;
         count[temp]++;
     }
  }

  for (temp = 0; temp < NUM_TEMPS; temp++)
      if (count[temp] > 0)
          average_temp[temp] = total[temp]/count[temp];

  return PGS_S_SUCCESS;
}
