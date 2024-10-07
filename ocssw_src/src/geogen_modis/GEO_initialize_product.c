#include "PGS_EPH.h"
#include "GEO_geo.h"
#include "GEO_output.h"
#include "GEO_product.h"
#include "PGS_MODIS_35251.h"

PGSt_SMF_status GEO_initialize_product(
	int				const number_of_scans,
	fill_values_struct const	* fill_values,
	MODFILE				* const geo_file,
	int32				const swfid,
	GEO_param_struct const		* GEO_param
)
/*
!C*****************************************************************************
!Description:   
		Routine in output group of the Level-1A geolocation
                software to initialize the HDF objects in the output 
		product for subsequent scan-by-scan output.  It calls 
		MAPI routines to initialize the scan metadata, to
                initialize the scan data in the geolocation file, and
                assigns a name to the SDS dimension.  

!Input Parameters:
       number_of_scans   the number of scans.
       fill_values       fill values for scan level data
       geo_file          M-API structure for the geolocation file.
       swfid             Swath (HDF-EOS) filehandle.
       GEO_param         struct containing runtime parameters.

!Output Parameters:      None

Return Values:
       MODIS_E_BAD_INPUT_ARG     If any argument is invalid
       MODIS_E_GEO               If any subroutine fails
       PGS_S_SUCCESS             Otherwise

Externally Defined:
		ATT_QUALITY		"GEO_product.h"
                ATTIT_ANG               "GEO_product.h"
                AVERAGE_TEMPERATURES    "GEO_product.h"
                DATATYPELENMAX          "mapi.h"
                DFNT_CHAR               "hntdefs.h"
                DFNT_FLOAT32            "hntdefs.h"
                EV_FRAMES               "GEO_product.h"
                EVTIME                  "GEO_product.h"
                GEO_QUALITY             "GEO_product.h"
                GFLAGS                  "GEO_product.h"
                GFLAGS_FVALUE           "GEO_geo.h"
                _HDF_VDATA              "vg.h"  
                HEIGHT                  "GEO_product.h"
                HGHT_FVALUE             "GEO_geo.h"
                IMPULSE_ENC             "GEO_product.h"
                IMPULSE_TIME            "GEO_product.h"
                I16                     "mapi.h"
                I32                     "mapi.h"
                I8                      "mapi.h"
                L_SMASK_FVALUE          "GEO_geo.h"
                L1_QUALITY              "GEO_product.h"
                L1A_ENGINEERING_GRP     "GEO_product.h"
                LAND_SEAMASK            "GEO_product.h"
                LAT_FVALUE              "GEO_geo.h"
                LATITUDE                "GEO_product.h"
                LONG_FVALUE             "GEO_geo.h"
                LONGITUDE               "GEO_product.h"
                M03PITCH_ELEM           "mapiL1Bgeo.h"
                M03ROLL_ELEM            "mapiL1Bgeo.h"
                M03YAW_ELEM             "mapiL1Bgeo.h"
                MAPIOK                  "mapi.h"
                MAX_IMPULSE_NUMBER      "GEO_parameters.h"
                MAX_SCAN_NUMBER         "GEO_geo.h"
                MDATA_RANGE             "mapi.h" 
                MFILL_VALUE             "mapi.h"
                MODIS_E_GEO             "PGS_MODIS_35251.h"
                MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
                MOON_VECTOR             "GEO_product.h"
                MSIDE                   "GEO_product.h"
                MSLOPE                  "mapi.h"
                MUNITS                  "mapi.h"
                NSCANS                  "GEO_product.h"
                NUM_IMPULSE             "GEO_product.h"
                NUM_L1A_QUALITY_FLAGS   "GEO_geo.h"
                NUM_TEMPS               "GEO_parameters.h"
                NUMENC                  "GEO_product.h"
                NUMQUAL                 "GEO_product.h"
                ORB_POS                 "GEO_product.h"
                ORB_VEL                 "GEO_product.h"
                PGS_S_SUCCESS           "PGS_SMF.h"
                RANGE                   "GEO_product.h"
                RANGE_FVALUE            "GEO_geo.h"
                R32                     "mapi.h"
                R64                     "mapi.h"
                S_NUM                   "GEO_product.h"
                S_TYPE                  "GEO_product.h"
                SCAN_GRP                "GEO_product.h"
                SCAN_META_GRP           "GEO_product.h" 
                SCAN_TYPE_LEN           "GEO_parameters.h"
                SCTIME                  "GEO_product.h"
                SD_FRAMES               "GEO_product.h"
                SDTIME                  "GEO_product.h"
                SEN_AZIMUTH             "GEO_product.h"
                SEN_ZENITH              "GEO_product.h"
                SENSORAZIM_FVALUE       "GEO_geo.h"
                SENSORZEN_FVALUE        "GEO_geo.h"
                SOLARAZIM_FVALUE        "GEO_geo.h"
                SOLARZEN_FVALUE         "GEO_geo.h"
                SUCCEED                 "hdf.h"
                SUN_AZIMUTH             "GEO_product.h"
                SUN_REF                 "GEO_product.h"
                SUN_ZENITH              "GEO_product.h"
                SV_FRAMES               "GEO_product.h"
                T_INST2ECR              "GEO_product.h"
                TA_RC_SMIR_CFPA         "GEO_product.h"
                TEMP_FVALUE             "GEO_geo.h"
                TP_AO_SMIR_OBJ          "GEO_product.h"
                TP_MF_CALBKHD_SR        "GEO_product.h"
                TP_MF_Z_BKHD_BB         "GEO_product.h"
                TP_SA_RCT1_MIR          "GEO_product.h"
                TP_SR_SNOUT             "GEO_product.h"
                TXT                     "mapi.h"
                UI16                    "mapi.h"
                UI8                     "mapi.h"
                VECDIM                  "GEO_product.h"
		WATER_PRESENT		"GEO_product.h"

Called by:
                GEO_locate_one_granule()
		
Routines called:
                createMODISarray        "mapi.h"
                GEO_create_swath        "GEO_output.h"
                modsmf                  "SDST.h"
                putMODISarinfo          "mapi.h"
                putMODISarray           "mapi.h"
                putMODISdimname         "mapi.h"
                VSattach                "hproto.h"
                VSdetach                "hproto.h"
                VSfdefine               "hproto.h"
                VSsetattr               "hproto.h"
                VSsetfields             "hproto.h"
                VSsetname               "hproto.h"

!Revision History:
     $Log: GEO_initialize_product.c,v $
     Revision 6.5  2011/02/14 21:32:12  kuyper
     Corrected const-qualification of parameters that point at input data.
     Added waterpresent SDS.
     Added valid_range attribute for landsea mask.

     Revision 6.4  2010/06/25 19:16:37  kuyper
     Change to make attitude and ephemeris quality flags not be Swath data fields.
     Added roll/pitch/yaw index SDS attributes to the Thermal Correction SDS.

     Revision 6.3  2010/05/19 20:23:08  kuyper
     Removed work-around for HDF bug that has since been fixed.
     Added Thermal Correction, Attitude and Ephemeris Quality SDSs.

     James Kuyper Jr.		James.R.Kuyper@NASA.gov

     Revision 6.2  2009/05/30 22:44:07  xgeng
     All strings in error messages are quoted.

     Revision 6.1  2009/05/30 19:28:44  xgeng
     Added hi-resolution output SDSs, if N_samp is greater than 1.
     Changed to expect a status code from GEO_create_swath().

     Xu Geng (xu.geng@saic.ocm)
 
     Revision 5.1  2004/08/18 19:29:51  vlin
     1. create the Average Temperatures vdata when number_of_scans is 0.
     2. input argument "num_detectors" is replaced by "fill_values".
     3. Add _FillValue attributes for members of fill_values_struct.

     Revision 4.3  2003/08/20 14:50:36  kuyper
     Reinstated redundant num_detectors parameter; interface is frozen for now.

     Revision 4.2  2003/07/28 20:33:44  vlin
     updated after code walkthrough.

     Revision 4.1  2003/07/03 17:26:55  vlin
     Added temperature vdata, Dropped num_detectors argument,
     Changed some error messages to MODIS_E_GEO

!Team-unique Header:

     This software is developed by the MODIS Science Data Support
     Team for the National Aeronautics and Space Administration,
     Goddard Space Flight Center, under contract NAS5-32373.

     HDF portions developed at the National Center for Supercomputing
     Applications at the University of Illinois at Urbana-Champaign.

!END**************************************************************************/

{
  enum dimtype{
    nscans,
    mframes,
    vecdim,
    numqual,
    numenc,
    stypelen,
    dimcount,
    none
  };

 static struct {
    char *name;
 long int  value;
 } dimensions[dimcount] = { 
   {NSCANS, 0},                      /* nscans  */
   {MFRAMES, MAX_FRAMES},	     /* mframes */
   {VECDIM, 3},                      /* vecdim  */
   {NUMQUAL, NUM_L1A_QUALITY_FLAGS}, /* numqual */
   {NUMENC, MAX_IMPULSE_NUMBER},     /* numenc  */
   {S_TYPE, SCAN_TYPE_LEN}              /* stypelen */
 }; 

 static struct {
    char *name;
    int   value;
 } element[] = { 
   {M03ROLL_ELEM,  0},
   {M03PITCH_ELEM, 1},
   {M03YAW_ELEM,   2},
 };

  /* The following arrays contain the min, and max pointed to
     by the initializers for the SDS_object structure.	*/

  static uint8 ui8value[][3]={
    {(uint8)0, (uint8)MAX_IMPULSE_NUMBER},	/* num_impulse  */
    {0, 0, (uint8)GFLAGS_FVALUE},		/* gflags fill value  */
    {(uint8)SHALLOW_OCEAN, (uint8)DEEP_OCEAN, (uint8)L_SMASK_FVALUE},
    /* land sea mask range, fill value  */
    {0, 8, (uint8)GFLAGS_FVALUE}	/* water present range, fill value */
  };

  static uint16 ui16value[][3]={
    {(uint16)0, (uint16)1},			/* Mirror side	*/
    {(uint16)27000, (uint16)65535, (uint16)RANGE_FVALUE}, /* Range	*/
    {(uint16)0, (uint16)1400},			/* EV frames	*/
  };

  static uint32 ui32value[][3]={
    {0, PGSd_INTERPOLATED_POINT | EPH_QUALITY_UNAVAIL | EPH_REPAIRED |
	EPH_LONG_PRECEED | EPH_SHORT_FOLLOW | EPH_REDHI | EPH_YELLOWLO |
	EPH_DATA | EPH_OVERALL, PGSd_NO_DATA},	/* Attitude Quality */
    {0, PGSd_INTERPOLATED_POINT | EPH_REPAIRED | EPH_LONG_PRECEED |
	EPH_SHORT_FOLLOW | EPH_REDHI | EPH_YELLOWLO | EPH_DATA | EPH_OVERALL,
	PGSd_NO_DATA},				/* Ephemeris Quality */
  };

  static int16 i16value[][3]={
    {(int16)-400, (int16)10000, (int16)HGHT_FVALUE},         /* height	 */
    {(int16)0, (int16)(18000), (int16)SENSORZEN_FVALUE},     /* SensorZenith */
    {(int16)-18000, (int16)18000, (int16)SENSORAZIM_FVALUE}, /* SensorAzimuth */
    {(int16)0, (int16)18000, (int16)SOLARZEN_FVALUE},	     /* SolarZenith */
    {(int16)-18000, (int16)18000, (int16)SOLARAZIM_FVALUE},  /* SolarAzimuth */
  };

  static float32 r32value[][3]={
    {(float32)-180.0, (float32)180.0, (float32)LONG_FVALUE},   /* longitude */
    {(float32)(-90.0), (float32)(90.0), (float32)LAT_FVALUE},  /* latitude  */
    {(float32)(-1.0), (float32)(1.0)},                         /* sun_ref   */
    /* Thermal attitude correction: */
    {(float32)-180.0, (float32)180.0, (float32)THERMCORR_FVALUE},
  };

  static float64 r64value[][2]={
    {(float64)-7200000.0, (float64)7200000.0},		/* orb_pos      */
    {(float64)-7600.0, (float64)7600.0},		/* orb_vel	*/
    {(float64)0.0, (float64)16384.0},			/* impulse_enc	*/
    {(float64)0.0, (float64)1.5},			/* impulse_time	*/
    {(float64)-1.0, (float64)1.0}			/* T_inst2ECR  */
  };

  static int8 i8value[3]={(int8)-127, (int8)127, HIRES_FVALUE};

  static struct{
    char const		*name;		/* SDS object name.    */
    char const		*group;		/* HDF vgroup.         */
    long		rank;		/* Array rank.         */
    enum dimtype	dims[3];	/* Array dimensions    */
    char		data_type[DATATYPELENMAX];	/* M-API data type. */
    unsigned char	not_swath;	/* flag defining whether an SDS needs to
                                           be created and dimensions named.
                                           If true, the SDS must be created and
                                           the dimension named.         */
    char const		*units;		/* units of measurement. 	*/
    void const		*range;		/* min/max valid data pair.	*/
    void const		*fillvalue;	/* Value for unset array cells.	*/
    void const		*scale_val;	/* SDS attribute scale value.   */
  } SDS_object[]={
    {S_NUM,       SCAN_META_GRP, 1, {nscans},  I16, 1},	/* 0 */
    {EV_FRAMES,	  SCAN_META_GRP, 1, {nscans}, UI16, 1, NULL, ui16value[2]},/*1*/
    {SD_FRAMES,	  SCAN_META_GRP, 1, {nscans}, UI16, 1},	/* 2 */
    {SV_FRAMES,	  SCAN_META_GRP, 1, {nscans}, UI16, 1}, /* 3 */
    {EVTIME,      SCAN_META_GRP, 1, {nscans},  R64, 1, "seconds"},	/* 4 */
    {SDTIME,      SCAN_META_GRP, 1, {nscans},  R64, 1, "seconds"},	/* 5 */
    {SVTIME,      SCAN_META_GRP, 1, {nscans},  R64, 1, "seconds"},	/* 6 */
    {SCTIME,      SCAN_META_GRP, 1, {nscans},  R64, 1, "seconds"},	/* 7 */
    {MSIDE,       SCAN_META_GRP, 1, {nscans},  UI16,1, NULL, ui16value[0]},/*8*/
    {SUN_ZENITH,  SCAN_META_GRP, 1, {nscans},  R32, 1, "radians"},	/* 9 */
    {SUN_AZIMUTH, SCAN_META_GRP, 1, {nscans},  R32, 1, "radians"},	/* 10 */
    {MOON_VECTOR, SCAN_META_GRP, 2, {nscans, 
                                     vecdim}, R32, 1},			/* 11 */
    {L1_QUALITY,  SCAN_META_GRP, 2, {nscans, 
                                     numqual}, I32, 1},			/* 12 */
    {GEO_QUALITY, SCAN_META_GRP, 2, {nscans, 
                                     numqual}, I8, 1},			/* 13 */
    {ORB_POS,     SCAN_META_GRP, 2, {nscans, 				/* 14 */
                                     vecdim},  R64, 1, "meters", r64value[0]},
    {ORB_VEL,     SCAN_META_GRP, 2, {nscans, 				/* 15 */
                                     vecdim},  R64, 1, "meters per second", 
                                                        r64value[1]},
    {T_INST2ECR,  SCAN_META_GRP, 3, {nscans, 
                                     vecdim, 				/* 16 */
                                     vecdim},  R64, 1, NULL, r64value[4]},
    {ATTIT_ANG,	  SCAN_META_GRP, 2, {nscans, 
                                     vecdim},  R64, 1, "radians"},	/* 17 */
    {SUN_REF,	  SCAN_META_GRP, 2, {nscans, 				/* 18 */
                                     vecdim}, R32, 1, NULL, r32value[2]},
    {NUM_IMPULSE, SCAN_META_GRP, 1, {nscans},  UI8, 1, NULL, ui8value[0]},/*19*/
    {IMPULSE_ENC, SCAN_META_GRP, 2, {nscans, numenc},  R64, 1, "encoder pulses",
	r64value[2]},							/* 20 */
    {IMPULSE_TIME,SCAN_META_GRP, 2, {nscans, numenc},  R64, 1, "seconds",
	r64value[3]},							/* 21 */
    {LONGITUDE,	  SCAN_GRP,	 0, {none},    R32, 0, "degrees", r32value[0],
	&r32value[0][2]},						/* 22 */
    {LATITUDE,	  SCAN_GRP,      0, {none},    R32, 0, "degrees", r32value[1],
	&r32value[1][2]},						/* 23 */
    {HEIGHT, 	  SCAN_GRP,	 0, {none},    I16, 0, "meters",  i16value[0],
	&i16value[0][2]},						/* 24 */
    {SEN_ZENITH,  SCAN_GRP,	 0, {none},    I16, 0, "degrees", i16value[1],
	&i16value[1][2]},						/* 25 */
    {SEN_AZIMUTH, SCAN_GRP,	 0, {none},    I16, 0, "degrees", i16value[2],
	&i16value[2][2]},						/* 26 */
    {RANGE,       SCAN_GRP,	 0, {none},   UI16, 0, "meters",  ui16value[1],
	&ui16value[1][2]},						/* 27 */
    {SOL_ZENITH,  SCAN_GRP,      0, {none},    I16, 0, "degrees", i16value[3],
	&i16value[3][2]},						/* 28 */
    {SOL_AZIMUTH, SCAN_GRP,	 0, {none},    I16, 0, "degrees", i16value[4],
	&i16value[4][2]},						/* 29 */
    {LAND_SEAMASK, SCAN_GRP,	 0, {none},    UI8, 0, NULL, ui8value[2],
	&ui8value[2][2]},						/* 30 */
    {WATER_PRESENT, SCAN_GRP,	 0, {none},    UI8, 0, NULL, ui8value[3],
	&ui8value[3][2]},						/* 31 */
    {GFLAGS,	   SCAN_GRP,     0, {none},    UI8, 0, NULL, NULL,
	&ui8value[1][2]},						/* 32 */
    {S_TYPE,      SCAN_META_GRP, 2, {nscans,stypelen},  TXT, 1},	/* 33 */
    {THERMCORR, SCAN_META_GRP, 2, {nscans, vecdim}, R32, 1, "degrees",
	r32value[3], &r32value[3][2]},					/* 34 */
    {ATT_QUALITY, SCAN_META_GRP, 2, {nscans, mframes}, UI32, 1, NULL,	
	ui32value[0], &ui32value[0][2]},				/* 35 */
    {EPH_QUALITY, SCAN_META_GRP, 2, {nscans, mframes}, UI32, 1, NULL,
	ui32value[1], &ui32value[0][2]},				/* 36 */
    /* The next objects must be the last ones in the array. */
    {SCAN_OFFSET, SCAN_GRP, 0, {none}, I8, 0, "km IFOV", i8value, &i8value[2]}, 
    {TRACK_OFFSET, SCAN_GRP, 0, {none}, I8, 0, "km IFOV", i8value, &i8value[2]},
    {HEIGHT_OFFSET, SCAN_GRP, 0, {none}, I8, 0, "km", i8value, &i8value[2]}
  };

  enum eqn_value{
       TA07, 
       TP01, 
       TP02, 
       EQN_TYPES
  };
  const float32 temp_limits[EQN_TYPES][2] = { {50.0, 118.0}, 
                               {-22.0, 87.0}, {-23.0, 87.0}};
  float32 fill_temp[1];
  const int eqn_type[6] = {TA07, TP02, TP01, TP01, TP02, TP02};
  const char *temp_units[] = {"Kelvin ", "Celsius", "Celsius", 
                              "Celsius", "Celsius", "Celsius"};
  char  *temperature_field[6] =  {TA_RC_SMIR_CFPA, TP_AO_SMIR_OBJ, 
             TP_MF_CALBKHD_SR, TP_MF_Z_BKHD_BB, TP_SA_RCT1_MIR, TP_SR_SNOUT};
  char  TEMPERATURE_FIELDS[] =   TA_RC_SMIR_CFPA","TP_AO_SMIR_OBJ","
             TP_MF_CALBKHD_SR","TP_MF_Z_BKHD_BB","TP_SA_RCT1_MIR","TP_SR_SNOUT; 
  static struct {
      char*	name;
      uint32	value;
  } quality_attr[] = {
      {QFL_OVERALL,	EPH_OVERALL},
      {QFL_DATASUM,	EPH_DATA},
      {QFL_RED_LO,	EPH_REDLO},
      {QFL_YEL_LO,	EPH_YELLOWLO},
      {QFL_YEL_HI,	EPH_YELLOWHI},
      {QFL_RED_HI,	EPH_REDHI},
      {QFL_LONG_FOLL,	EPH_LONG_FOLLOW},
      {QFL_SHRT_FOLL,	EPH_SHORT_FOLLOW},
      {QFL_SHRT_PREC,	EPH_SHORT_PRECEED},
      {QFL_LONG_PREC,	EPH_LONG_PRECEED},
      {QFL_REPAIRED,	EPH_REPAIRED},
      {QFL_QFL_PROB,	EPH_QUALITY_UNAVAIL},
      {QFL_INTERP,	PGSd_INTERPOLATED_POINT}
  };
  long int          dim_sizes[3];
  PGSt_SMF_status   ret_val = PGS_S_SUCCESS;
  int i, obj, dim, temp;
  int num_objs = (int)(sizeof SDS_object /sizeof SDS_object[0]); 
  int32 vdata_id, status_32;
  intn  status_n;
  char  msgbuf[PGS_SMF_MAX_MSGBUF_SIZE] = "";
  char  filefunc[] = __FILE__ ", GEO_initialize_product";


  if (number_of_scans < 0 || number_of_scans > MAX_SCAN_NUMBER ||
      fill_values == NULL || geo_file == NULL || swfid == FAIL || 
      GEO_param == NULL) {
     sprintf(msgbuf,"geo_file = %p, number_of_scans = %d, swfid = %ld\n"
             "fill_values = %p, GEO_param = %p", (void *)geo_file, 
             number_of_scans, (long)swfid, (void *)fill_values, 
             (void *)GEO_param);
     modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
     return MODIS_E_BAD_INPUT_ARG;
  }

/* C5 feature: 
 *    Create the Average Temperatures vdata even if number_of_scans is 0.
 */

  vdata_id = VSattach((int32)geo_file->hdf_id, _HDF_VDATA, "w");
  if (vdata_id == FAIL) {
      modsmf(MODIS_E_GEO, "VSattach()", filefunc); 
      return MODIS_E_GEO;
  }

  status_32 = VSsetname(vdata_id, AVERAGE_TEMPERATURES);
  if (status_32 != SUCCEED) {
      modsmf(MODIS_E_GEO, "VSsetname()", filefunc);
      ret_val = MODIS_E_GEO; 
  }

  for (temp = 0; temp < NUM_TEMPS; temp++) {
       status_n  = VSfdefine(vdata_id, temperature_field[temp], 
                    DFNT_FLOAT32, 1);
       if (status_n != SUCCEED) {
           modsmf(MODIS_E_GEO, "VSfdefine()", filefunc);
           ret_val = MODIS_E_GEO;
       } 
  }

  status_n  = VSsetfields(vdata_id, TEMPERATURE_FIELDS);
  if (status_n != SUCCEED) {
      modsmf(MODIS_E_GEO, "VSsetfields()", filefunc);
      ret_val = MODIS_E_GEO;
  }

  fill_temp[0] = TEMP_FVALUE;
  if (VSsetattr(vdata_id, _HDF_VDATA, MFILL_VALUE, DFNT_FLOAT32, 1, 
                fill_temp) != SUCCEED) {
      sprintf(msgbuf, "VSsetattr(\"%s\", \"%s\")", AVERAGE_TEMPERATURES, MFILL_VALUE);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      ret_val = MODIS_E_GEO;
  } 

  for (temp = 0; temp < NUM_TEMPS; temp++) {
      if (VSsetattr(vdata_id, temp, MUNITS, DFNT_CHAR, 7, 
          temp_units[temp]) != SUCCEED) {
          sprintf(msgbuf, "VSsetattr(\"%s\",%d,\"%s\",\"%s\")", AVERAGE_TEMPERATURES,
                  temp, MUNITS, temp_units[temp]);
          modsmf(MODIS_E_GEO, msgbuf, filefunc);
          ret_val = MODIS_E_GEO;
      }
      if (VSsetattr(vdata_id, temp, MDATA_RANGE, DFNT_FLOAT32, 2, 
                    temp_limits[eqn_type[temp]]) != SUCCEED) {
          sprintf(msgbuf, "VSsetattr(\"%s\", %d, \"%s\")", AVERAGE_TEMPERATURES,
                  temp, MDATA_RANGE);
          modsmf(MODIS_E_GEO, msgbuf, filefunc);
          ret_val = MODIS_E_GEO;
      }
  }

  status_32 = VSdetach(vdata_id);
  if (status_32 != SUCCEED) {
      modsmf(MODIS_E_GEO, "VSdetach()", filefunc);
      ret_val = MODIS_E_GEO;
  }

  if (number_of_scans == 0)
     return PGS_S_SUCCESS;

  dimensions[nscans].value = (long)number_of_scans;
  SDS_object[25].scale_val = &GEO_param->angle_scale;  /* SensorZenith  */
  SDS_object[26].scale_val = &GEO_param->angle_scale;  /* SensorAzimuth */
  SDS_object[27].scale_val = &GEO_param->range_scale;  /* Range         */
  SDS_object[28].scale_val = &GEO_param->angle_scale;  /* SolarZenith   */
  SDS_object[29].scale_val = &GEO_param->angle_scale;  /* SolarAzimuth  */
  SDS_object[37].scale_val = &GEO_param->hires_scale;  /* scan offset  */
  SDS_object[38].scale_val = &GEO_param->hires_scale;  /* track offset  */
  SDS_object[39].scale_val = &GEO_param->hires_scale;  /* height offset  */
  SDS_object[0].fillvalue = &fill_values->scan_number;
  SDS_object[4].fillvalue = &fill_values->EV_start_time;
  SDS_object[5].fillvalue = &fill_values->SD_start_time;
  SDS_object[6].fillvalue = &fill_values->SV_start_time;
  SDS_object[7].fillvalue = &fill_values->EV_center_time;
  SDS_object[8].fillvalue = &fill_values->mirr_side;
  SDS_object[9].fillvalue = &fill_values->SD_sun_zenith;
  SDS_object[10].fillvalue = &fill_values->SD_sun_azimuth;
  SDS_object[11].fillvalue = &fill_values->moon_vector;
  SDS_object[12].fillvalue = &fill_values->L1_scan_quality;
  SDS_object[13].fillvalue = &fill_values->geo_scan_quality;
  SDS_object[14].fillvalue = &fill_values->orb_pos;
  SDS_object[15].fillvalue = &fill_values->orb_vel;
  SDS_object[16].fillvalue = &fill_values->T_inst2ECR;
  SDS_object[17].fillvalue = &fill_values->attitude_angels;
  SDS_object[18].fillvalue = &fill_values->sun_ref;
  SDS_object[20].fillvalue = &fill_values->impulse_enc;
  SDS_object[21].fillvalue = &fill_values->impulse_time;

  if (GEO_create_swath(number_of_scans, GEO_param->num_detectors, swfid) != 
      PGS_S_SUCCESS)
  {
     ret_val = MODIS_E_GEO;
     modsmf(MODIS_E_GEO, "GEO_create_swath()", filefunc);
  }

  /* For this to work, the offsets must be the last objects in the array. */
  if(GEO_param->geometry_params.N_samp[GEO_param->geometry_params.band_number]
    == 1)
    num_objs -= 3; 
   
  /* Create the SDS arrays, and write their attributes.	*/
  for (obj=0; obj<num_objs; obj++) {

    if ( SDS_object[obj].not_swath )
    {
       memset (dim_sizes, 0, sizeof(dim_sizes));
       for (i=0; i<(int)SDS_object[obj].rank; i++)
           dim_sizes[i] = dimensions[SDS_object[obj].dims[i]].value;

       if (createMODISarray(geo_file, SDS_object[obj].name,
          SDS_object[obj].group, SDS_object[obj].data_type,
          SDS_object[obj].rank, dim_sizes) != MAPIOK) {
          sprintf(msgbuf,"createMODISarray(\"%s\",\n \"%s\", \"%s\")", geo_file->filename,
                  SDS_object[obj].name, SDS_object[obj].data_type);
          modsmf(MODIS_E_GEO, msgbuf, filefunc);
          ret_val = MODIS_E_GEO;
       } 
       else {
          for (dim=0; dim<(int)SDS_object[obj].rank; dim++){
             if (putMODISdimname(geo_file, SDS_object[obj].name,
                 SDS_object[obj].group, (long)dim,
                 dimensions[SDS_object[obj].dims[dim]].name) != MAPIOK){
                 sprintf(msgbuf, "putMODISdimname(\"%s\",\n \"%s\", \"%s\", \"%s\")",
                 geo_file->filename, SDS_object[obj].name, SDS_object[obj].group, 
                 dimensions[SDS_object[obj].dims[dim]].name);
                 modsmf(MODIS_E_GEO, msgbuf, filefunc);
                 ret_val = MODIS_E_GEO;
             } 
          }
       }
     }
    
     if ((SDS_object[obj].units) && putMODISarinfo(geo_file, 
        SDS_object[obj].name, SDS_object[obj].group, MUNITS, TXT, 
        (long)strlen(SDS_object[obj].units), SDS_object[obj].units) != MAPIOK){
        sprintf(msgbuf,
	    "putMODISarinfo(\"%s\",\n \"%s\", \"%s\", \"%s\", \"%s\")",
	    geo_file->filename, SDS_object[obj].name, SDS_object[obj].group,
	    MUNITS, SDS_object[obj].units);
        modsmf(MODIS_E_GEO, msgbuf, filefunc);
        ret_val = MODIS_E_GEO;
     }

     if ((SDS_object[obj].range) && putMODISarinfo(geo_file, 
        SDS_object[obj].name, SDS_object[obj].group, MDATA_RANGE, 
        SDS_object[obj].data_type, 2L, SDS_object[obj].range) != MAPIOK) {
        sprintf(msgbuf, "putMODISarinfo(\"%s\",\n \"%s\", \"%s\", \"%s\", \"%s\")",
                geo_file->filename, SDS_object[obj].name,
                SDS_object[obj].group, MDATA_RANGE, SDS_object[obj].data_type);
        modsmf(MODIS_E_GEO, msgbuf, filefunc);
        ret_val = MODIS_E_GEO;
     }

     if (SDS_object[obj].fillvalue) {
	dim_sizes[0] = dim_sizes[1] = dim_sizes[2] = 1;
        if (putMODISarinfo(geo_file, SDS_object[obj].name, 
            SDS_object[obj].group, MFILL_VALUE, SDS_object[obj].data_type, 1L, 
            SDS_object[obj].fillvalue) != MAPIOK) {
            sprintf(msgbuf, "putMODISarinfo(\"%s\",\n \"%s\", \"%s\", \"%s\", \"%s\")",
              geo_file->filename, SDS_object[obj].name, SDS_object[obj].group,
              MFILL_VALUE, SDS_object[obj].data_type);
            modsmf(MODIS_E_GEO, msgbuf, filefunc);
            ret_val = MODIS_E_GEO;
        }
     }

     if ((SDS_object[obj].scale_val) && putMODISarinfo(geo_file, 
        SDS_object[obj].name, SDS_object[obj].group, MSLOPE, R64, 1L, 
        SDS_object[obj].scale_val) != MAPIOK) {
        sprintf(msgbuf,"putMODISarinfo(\"%s\"\n\"%s\", \"%s\", \"%s\", \"%s\")",
           geo_file->filename, SDS_object[obj].name, SDS_object[obj].group,
           MSLOPE, SDS_object[obj].data_type);
        modsmf(MODIS_E_GEO, msgbuf, filefunc);
        ret_val = MODIS_E_GEO;
     }
  }

  for (obj=0; obj<(int)(sizeof element / sizeof element[0]); obj++) {
    if (putMODISarinfo(geo_file, ATTIT_ANG, SCAN_META_GRP,
        element[obj].name, UI32, 1L, &element[obj].value) != MAPIOK)
    {
        sprintf(msgbuf, "putMODISarinfo(\"%s\",\n\"%s\", \"%s\")",
	    geo_file->filename, ATTIT_ANG, element[obj].name);
        modsmf(MODIS_E_GEO, msgbuf, filefunc);
        ret_val = MODIS_E_GEO;
    }
    if (putMODISarinfo(geo_file, THERMCORR, SCAN_META_GRP,
        element[obj].name, UI32, 1L, &element[obj].value) != MAPIOK)
    {
        sprintf(msgbuf, "putMODISarinfo(\"%s\",\n\"%s\", \"%s\")",
	    geo_file->filename, THERMCORR, element[obj].name);
        modsmf(MODIS_E_GEO, msgbuf, filefunc);
        ret_val = MODIS_E_GEO;
    }
  }

  for (obj=0; obj<(int)(sizeof(quality_attr)/sizeof(quality_attr[0])); obj++)
  {
      if (putMODISarinfo(geo_file, ATT_QUALITY, SCAN_META_GRP,
	  quality_attr[obj].name, UI32, 1L, &quality_attr[obj].value) != MAPIOK)
      {
	  sprintf(msgbuf, "putMODISarinfo(\"%s\",\n\"%s\", \"%s\")",
	      geo_file->filename, ATT_QUALITY, quality_attr[obj].name);
	  modsmf(MODIS_E_GEO, msgbuf, filefunc);
	  ret_val = MODIS_E_GEO;
      }

      if (putMODISarinfo(geo_file, EPH_QUALITY, SCAN_META_GRP,
	  quality_attr[obj].name, UI32, 1L, &quality_attr[obj].value) != MAPIOK)
      {
	  sprintf(msgbuf, "putMODISarinfo(\"%s\",\n\"%s\", \"%s\")",
	      geo_file->filename, EPH_QUALITY, quality_attr[obj].name);
	  modsmf(MODIS_E_GEO, msgbuf, filefunc);
	  ret_val = MODIS_E_GEO;
      }
  }

  return ret_val;
}
