#include    "Preprocess.h"   /* contains externals, includes L1B_Tables.h */
#include    "PreprocessP.h"  /* contains local definitions and prototypes */
#include    "HDF_Lib.h"
#include    "PGS_PC.h"
#include    "PGS_MET.h"
#include    "PGS_Error_Codes.h"
#include    "FNames.h"
#include    <math.h>
#include    "granule_extract_metadata.h"

/***************************************************************************
Developer's note: (Jim Rogers)

In the the various functions, there is a return statement after each
L1BErrorMsg function call.  These are actually not necessary if the last
argument to L1BErrorMsg is "True". These returns make it easier to test
error out conditions (with SMF_ERROR used in UNIT_TEST_MODE_ONLY -- it does
not actually exit).
****************************************************************************/

int32 SUBSAMPLES_AT_RES[NUM_L1A_RESOLUTIONS] = {4, 2, 1, 1};

int16 TARGET_1km_FRAMES[NUM_TARGETS] = { 50, 10, 50, 50, 1354 };

/* Enclosed the rows in the initializer of the 2-D array with braces for
   ANSI-C compliance (Razor Issue 173). */
char *L1A_TARGET_SDS_NAME[NUM_TARGETS][NUM_L1A_RESOLUTIONS] = {
 {"SD_250m",   "SD_500m",   "SD_1km_day",   "SD_1km_night"},
 {"SRCA_250m", "SRCA_500m", "SRCA_1km_day", "SRCA_1km_night"},
 {"BB_250m",   "BB_500m",   "BB_1km_day",   "BB_1km_night"},
 {"SV_250m",   "SV_500m",   "SV_1km_day",   "SV_1km_night"},
 {"EV_250m",   "EV_500m",   "EV_1km_day",   "EV_1km_night"}
};

/*--------------------------------------------------------------
The following
L1A_BAND_DIM_NAME, L1A_TRACK_DIM_NAME, & L1A_FRAME_DIM_NAME
are used only in Preprocess.c to write OBC SDS's, and
may not be the same as in L1A input file. --zh
--------------------------------------------------------------*/
char *L1A_BAND_DIM_NAME[NUM_L1A_RESOLUTIONS] = {
  "Band_250m", "Band_500m", "Band_1km_day", "Band_1km_night"
};

char *L1A_TRACK_DIM_NAME[NUM_L1A_RESOLUTIONS] = {
  "40*nscans", "20*nscans", "10*nscans", "10*nscans"
};

char *L1A_FRAME_DIM_NAME[NUM_TARGETS-1][NUM_L1A_RESOLUTIONS] = {
 {"4*SD_frames",   "2*SD_frames",   "SD_frames",   "SD_frames"},
 {"4*SRCA_frames", "2*SRCA_frames", "SRCA_frames", "SRCA_frames"},
 {"4*BB_frames",   "2*BB_frames",   "BB_frames",   "BB_frames"},
 {"4*SV_frames",   "2*SV_frames",   "SV_frames",   "SV_frames"}
};

  /* January 4, 2002 (Razor Issue #154)
   * Terra- and Aqua- specific temperature coefficients have been moved
   * to PreprocessP.h.
   * Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)
   */

extern int16 RFLAG;
int32 extract_pixel_offset;
int32 extract_pixel_count;
int32 extract_line_offset;
int32 extract_line_count;

PGSt_SMF_status Calculate_Temp_QA(int32              num_scans,
                                  Temperatures_t     *temps,
                                  QA_tables_t        *tables,
			       	                    QA_Data_t          *QA)
/*
!C**************************************************************************
!Description:   Calculate the ratio of temperature variation to the value
                stored in look up table. The temperatures being calculated
                are BB thermistor temperatures, averaged BB thermistor
                temperature, lwir, mwir, nir, vis focal plane temperatures,
                mirror side 1 temperature, mirror side 2 temperature, averaged
                mirror side temperature, lwir adc temperature, mwir adc
                temperature, instrument temperature, cavity temperature.

!Input Parameters:
     int32                 num_scans    number of data for each kind of
                                        temperature values
     Temperatures_t        *temps       temperature values
     QA_tables_t           *tables,     contains the predetermined variation
                                        of temperatures

!Output Parameters:
     QA_Data_t        *QA               contains the scaled integer, which is
                                        ratio being mapped to 0 - 255.

!Revision History:
 $Log: Preprocess.c,v $
 Revision 1.35  2015-03-17 10:48:12-04  xgeng
 Fix the anomaly before and after the sector rotation due to the timing mismatch between the command performed and the status updated.

 Revision 1.34  2011-07-28 15:44:33-04  xgeng
 Default b1 is derived using real-time LWIR FPA temperature.

 Revision 1.32  2011-04-07 14:39:08-04  xgeng
 1. RSB &TEB uncertainty algorithm update; 2. The quadratic RSB RVS changed to 4th order.

 Revision 1.31  2010-11-15 11:40:04-05  xgeng
 Added a new function called Get_Electronics_index_special() for PCLW check; and no calibration is performed if both sides of PCLW electronics are on.

 Revision 1.30  2009/11/27 19:18:15  xgeng
 Fixed a bug in Revision 1.29

 Revision 1.29  2009/07/24 19:54:29  xgeng
 fixed the sector rotation abnormal

 Revision 1.27  2008/01/24 15:50:36  ltan
 Changed to relax the RVS correction limit range from [0.8, 1.2] to [0.4, 2.4].

 Revision 1.24  2005/01/18 21:57:39  ltan
 Added new file attributes prefixed with HDFEOS_FractionalOffset


 Revision 02.12, September 30, 1999
 Added assignment of num_thermistor_outliers.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.11 May 12, 1999
 Eliminated duplication of code by adding function.  Function also checks
 for a bad value of temperature (<= 0 is bad).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Initial development
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

!Team-unique Header:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32         i              = 0;
  int32         S              = 0;
  int32         n              = 0;
  float32       avg            = 0;
  float32       var            = 0;
  float32       ratio          = 0;
  float32       tbbavg[MAX_NUM_SCANS];

#define CALC_ONE_TEMP_QA(t,tablevar,output)                                 \
  returnStatus = Get_Temp_Avg_And_Variance(num_scans, t, &avg, &var);       \
  if (returnStatus != MODIS_S_OK)                                           \
    SMF_ERROR(MODIS_F_NOK,                                                  \
      "Get_Temp_Avg_And_Variance failed, Calculate_Temp_QA, Preprocess.c"); \
  ratio = var / tablevar;                                                   \
  CONVERT_TO_UINT8(ratio, output);

      /* Initialize the array holding the granule averages */

  for (i = 0; i < MAX_NUM_GRAN_AVERAGES; i++)
    QA->QA_common.granule_averages[i] = 0;

      /* First compute average BB temperature on each scan with
       * bad value rejection but not outlier rejection (the value in
       * temps->bb_avg has outlier rejection applied).
       * Also, assign number of outliers on each scan to QA_emiss.
       */

  for (S = 0; S < num_scans; S++)
  {
    avg = 0;
    n = 0;
    for (i = 0; i < NUM_BB_THERMISTORS; i++)
    {
      if (temps->bb[i][S] > 0)
      {
        avg += temps->bb[i][S];
        n++;
      }
    }
    if (n > 0)
      tbbavg[S] = avg / n;
    else
      tbbavg[S] = -1;
    QA->QA_emiss.num_thermistor_outliers[S] =
        temps->num_thermistor_outliers[S];
  }

      /* Now calculate noise, for all */

  for (i = 0; i < NUM_BB_THERMISTORS; i++)
  {
    CALC_ONE_TEMP_QA(temps->bb[i],
                     tables->emiss_QA_tables.var_T_bb[i],
                     QA->QA_emiss.noise_T_bb[i]);
    QA->QA_common.granule_averages[i] = avg;
  }

  CALC_ONE_TEMP_QA(tbbavg,
                   tables->emiss_QA_tables.var_T_bb_avg,
                   QA->QA_emiss.noise_T_bb_avg);

  CALC_ONE_TEMP_QA(temps->fp[3],
                   tables->emiss_QA_tables.var_T_lwir,
                   QA->QA_emiss.noise_T_lwir);
  QA->QA_common.granule_averages[GRAN_AVG_TA_RC_LWIR_CFPAE] = avg;

  CALC_ONE_TEMP_QA(temps->fp[2],
                   tables->emiss_QA_tables.var_T_mwir,
                   QA->QA_emiss.noise_T_mwir);
  QA->QA_common.granule_averages[GRAN_AVG_TA_RC_SMIR_CFPAE] = avg;

  CALC_ONE_TEMP_QA(temps->fp[1],
                   tables->refl_QA_tables.var_NIR_FPA,
                   QA->QA_refl.noise_T_nir);
  QA->QA_common.granule_averages[GRAN_AVG_TA_AO_NIR_FPAE] = avg;

  CALC_ONE_TEMP_QA(temps->fp[0],
                   tables->refl_QA_tables.var_visual_FPA,
                   QA->QA_refl.noise_T_vis);
  QA->QA_common.granule_averages[GRAN_AVG_TA_AO_VIS_FPAE] = avg;

  CALC_ONE_TEMP_QA(temps->mir1,
                   tables->emiss_QA_tables.var_T_mir1,
                   QA->QA_emiss.noise_T_mir1);
  QA->QA_common.granule_averages[GRAN_AVG_TP_SA_RCT1_MIRE] = avg;

  CALC_ONE_TEMP_QA(temps->mir2,
                   tables->emiss_QA_tables.var_T_mir2,
                   QA->QA_emiss.noise_T_mir2);
  QA->QA_common.granule_averages[GRAN_AVG_TP_SA_RCT2_MIRE] = avg;

  CALC_ONE_TEMP_QA(temps->mir_avg,
                   tables->emiss_QA_tables.var_T_mir_avg,
                   QA->QA_emiss.noise_T_mir_avg);

  CALC_ONE_TEMP_QA(temps->ins,
                   tables->emiss_QA_tables.var_T_ins,
                   QA->QA_emiss.noise_T_ins);

  CALC_ONE_TEMP_QA(temps->cav,
                   tables->emiss_QA_tables.var_T_cav,
                   QA->QA_emiss.noise_T_cav);

  if (Granule_Average_Temperature(num_scans, temps->scn_mtr,
      &QA->QA_common.granule_averages[GRAN_AVG_TP_SA_A_MTR]) != MODIS_S_OK)
  {
    SMF_ERROR(MODIS_F_NOK,
        "avg cavity temperature failed in Calculate_Temp_QA");
  }

  for (n = 0, i = GRAN_AVG_TP_MF_CALBKHD_SR; n < NUM_T_CAV_THERMISTORS
            && i <= GRAN_AVG_TP_MF_CVR_OP_SR; n++, i++)
  {
    if (Granule_Average_Temperature(num_scans, temps->cav_temp[n],
                  &QA->QA_common.granule_averages[i]) != MODIS_S_OK)
    {
      SMF_ERROR(MODIS_F_NOK,
          "avg cavity temperature failed in Calculate_Temp_QA");
    }
  }
  for (n = 0, i = GRAN_AVG_TP_AO_SMIR_OBJ; n < NUM_T_INS_THERMISTORS
            && i <= GRAN_AVG_TP_AO_LWIR_LENS; n++, i++)
  {
    if (Granule_Average_Temperature(num_scans, temps->ins_temp[n],
                  &QA->QA_common.granule_averages[i]) != MODIS_S_OK)
    {
      SMF_ERROR(MODIS_F_NOK,
          "avg inst. temperature failed in Calculate_Temp_QA");
    }
  }
  for (n = 0, i = GRAN_AVG_TA_RC_CS; n < NUM_T_RC_VALUES
            && i <= GRAN_AVG_TA_RC_OS_OG; n++, i++)
  {
    if (Granule_Average_Temperature(num_scans, temps->rc_temp[n],
                  &QA->QA_common.granule_averages[i]) != MODIS_S_OK)
    {
      SMF_ERROR(MODIS_F_NOK,
          "avg RC temperature failed in Calculate_Temp_QA");
    }
  }
  avg = 0;
  n = 0;
  for (S = 0; S < num_scans; S++) {
    if (temps->vr_lwir[S] != VOLTAGE_BAD_VALUE) {
      n++;
      avg += temps->vr_lwir[S];
    }
  }
  if (n > 0)
    QA->QA_common.granule_averages[GRAN_AVG_VR_RC_LW_FPA_HTR] = avg / n;
  else
    QA->QA_common.granule_averages[GRAN_AVG_VR_RC_LW_FPA_HTR] =
        VOLTAGE_BAD_VALUE;

  return returnStatus;
}

PGSt_SMF_status Calculate_Planck (float32 *RSR,
                                  float32 *wl,
                                  int16   size,
                                  float32 T,
                                  float32 *planck)
/*
!C************************************************************************
!Description:  Calculate Planck function integration over RSR

!Input Parameters:
    float32  *RSR   * RSR values for a particular channel *
    int16    size   * number of valid points for that RSR *
    float32  T      * temperature *
    float32  *wl    * wavelength with same array size*

!Output Parameters:
    float32 planck  * the planck function integration *

!Revision History:
 Revision 02.20 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Initial development
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16   i      = 0;
  float32 sum    = 0;
  float32 xL     = 0;
  float32 lambda = 0;
  float32 dw     = 0;

  *planck = 0;            /* initialize for error state on return */

  if (size <= 0 || T < TOLERANCE)    /* check for bad values */
    return (returnStatus);

  if (size == 1)                      /* check for single value */
  {
    *planck = L(wl[0], T);
    return (returnStatus);
  }

  lambda = wl[0];
  dw = 0.5 * ( wl[1] - wl[0] );
  xL = RSR[0] * L(lambda, T) * dw;
  sum = RSR[0] * dw;

  for (i = 1; i< size - 1; i++)
  {
    lambda = wl[i];
    dw = 0.5 * ( wl[i+1] - wl[i-1] );
    xL += RSR[i] * L(lambda, T) * dw;
    sum += RSR[i] * dw;
  }

  lambda = wl[size-1];
  dw = 0.5 * ( wl[size-1] - wl[size-2] );
  xL += RSR[size-1] * L(lambda, T) * dw;
  sum += RSR[size-1] * dw;

  if (sum > 0)
    *planck = xL / sum;

  return (returnStatus);
}

PGSt_SMF_status Get_Emiss_Coeff_Per_Scan
                          (int8               moon_in_SV_KOB,
                           int16              *BB,
                           int16              *SV,
                           int16              B,
                           int16              D,
                           int16              D_emiss,
                           int16              MS,
                           Emiss_Cal_Coeff_t  *RVS_Coeff,
                           emiss_tables_t     *tables,
                           float32            T_bb,
                           float32            T_mir,
                           float32            T_cav,
                           float32            T_ins,
                           float32            T_fp_lwir,
                           float32            *Xdn_bb_31,
                           float32            *Xb1,
                           float32            *xdLbb,
                           float32            *xdnbb,
                           float32            *xLbb,
                           float32            *xLcav,
                           float32            *xdnsv,
                           float32            *xdnsv_var,
                           uint32             *sv_omask,
                           float32            *SNR,
                           int32              satellite_ID)
/*
!C****************************************************************************
!Description:
 For a given emissive band and detector in the band, this function computes
 Emissive Calibration Coefficients for one scan. The quantities calculated
 include an averaged space-view (SV) digital number (DN_sv), the band-averaged
 radiance difference between the black-body and the space view (dL_bb),
 the linear response coefficient (b1), and the signal to noise ratio (SNR).

!Input Parameters:

 int8              moon_in_SV_KOB  Flag denoting if the moon is inside
                                   the SV KOB
 int16 *           BB           Black-body DNs (one set of frames)
 int16 *           SV           Space-view DNs (one set of frames)
 int16             B            Band index within the resolution.
 int16             D            Detector index within one emissive band.
 int16             D_emiss      Detector index within the set of emissive
                                detectors.
 int16             MS           Mirror side index
 Emiss_Cal_Coeff_t *RVS_Coeff   Emissive RVS Correction terms
 emiss_tables_t *  tables       emissive lookup tables, contains ADC
                                coefficients and the number frames for the BB
                                and SV averages, and the ADC and PCX correction
                                switches
 float32           T_bb         black body average temperature.
 float32           T_mir        average mirror temperature
 float32           T_cav        average cavity temperature.
 float32           T_ins        instrument temperature.
 float32           *Xdn_bb_31   the address of the blackbody dn value for the
                                detector D in band 31 and for current scan. It
                                is an input for bands 32-36 and an output for
                                band 31.

 float32 *         SNR          Address of the SNR. If entry value is > 0,
                                then calculate SNR. Otherwise, do not calculate.

 int32            satellite_ID 0 = TERRA, 1 = AQUA, -1 = INVALID_SATELLITE_ID


!Output Parameters:

 float32 *        Xb1          Address of one element of b1 of the middle
                               granule.
 float32 *        xdLbb        Address of one element of delta radiance of
                               black body for middle granule.
 float32 *        xdnbb        Address of one element of average blackbody
                               dn for middle granule.
 float32 *        xLbb         Address of one element of radiance of black
                               body for middle granule.
 float32 *        xLcav        Address of one element of radiance of cavity
                               for middle granule.
 float32 *        xdnsv        Address of one element of average space-view
                               DN for middle granule.
 float32 *        xdnsv_var    Address of one element of average space-view
                               variance for middle granule.
 uint32 *         sv_omask     array of two elements to return outlier mask.
 float32 *        SNR          Address of the SNR. may or may not be calculated
                               depending on input value.

!Revision History:
 Revision 02.31 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.30, January 29, 2002   Razor Issue #175
 Added satellite_ID to inputs.  Added a check on average BB temperature for Aqua
 bands (currently 33, 35, and 36):  If average BB temperature is above a
 band-specific cutoff value, use a value of b1 specified in the
 BB_T_sat_default_b1_aqua LUT if the BB_T_sat_switch_aqua switch is on.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)


 Revision 02.21  March 10, 2002   Razor Issue #174
 Added passed parameter Emiss_Cal_Coeff_t
   and changed use of emissive lookup tables to these arrays for
   RVS correction
 Alice Isaacman  SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.20, Jan 29, 2001
 Improved efficiency of calculation of a0 and a2.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.19 Sep 11, 2000
 Made changes for Razor issue 133 (TEB algorithm for calculating <DN-SV> when
 moon is in the SVP).  Added input argument "moon_in_SV_KOB".  Added logic
 for calculating <DN-SV> if the moon is inside the SV KOB.
 Corrected bug in calculating SNR.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.18 August 1999
 Removed xdnbb, one element of average black body DN for middle granule, from
 argument list.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.17 August 1999
 See SDF, "L1B Code Change to meet ADC correction algorithm change"
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.16 May, 1999
 See SDF, "Stage 3: Changes to Emissive Algorithm and Code" notes.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.15  April 12, 1999
 Removed cubic term a3
 Steven Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 02.11  Feb 9, 1999
 Added use of variables ADC_correction_switch and PCX_correction_switch
 set to macro values to remove compiler warnings.  Macros in upper case.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 01.00 Feb. 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

 WARNING: SNR is used as both an input and an output.  As an input, it is a
 flag to determine whether or not to calculate the SNR.  As an output,
 it holds the value of the SNR, if calculated, or -1 if not.
 The only purpose of having the flag is to save a little CPU time when
 processing scans from the leading and trailing granules.

!END**************************************************************************
*/

{
  PGSt_SMF_status returnStatus      = MODIS_S_OK;
  int16           F                 = 0;
  int16           B_xt              = 0;
  int16           B_emiss           = 0;
  float32         DN_bb             = 0;
  float32         DN_sv             = 0;
  float32         dL_bb             = 0;
  float32         dn_bb             = 0;
  float32         dn_bb_31;
  float32         bb_corr           = 0;
  float32         sv_corr           = 0;
  float32         a0                = 0;
  float32         a2                = 0;
  float32         L_sm              = 0;
  float32         L_bb              = 0;
  float32         L_cav             = 0;
  float32         sdev              = 0;
  int8            PCX_correction_switch;
  int8            BB_T_sat_switch_aqua;
  int16           rejects[MAX_1KM_OBC_FRAME_DIM];
  float32         delta_T_bb;
  int16           i;
  float32         BB_temp_threshold = 1000.0;
  boolean         calculate_SNR;

  /*   Determine the emissive band index. */
  if (B <= BAND26)
    B_emiss = B;
  else
    B_emiss = B - 1;

  /*
   * Determine if we will calculate SNR later in this function
   * (do this here so that we can initialize SNR with a return value)
   */

  if (*SNR > 0)
    calculate_SNR = True;
  else
    calculate_SNR = False;

  /*
   * Initialize return values in case of error exits
   * (Xb1 is already initialized)
   */

  *xdLbb     = -1;
  *xdnbb     = -1;
  *xLbb      = -1;
  *xLcav     = -1;
  *xdnsv     = -1;
  *xdnsv_var = -1;
  *SNR       = -1;

  PCX_correction_switch = tables->PCX_correction_switch;
  BB_T_sat_switch_aqua = tables->BB_T_sat_switch_aqua;

  /*
   * initialize "rejects" array to zero for all elements.
   */

  for (F = 0; F < MAX_1KM_OBC_FRAME_DIM; F++)
    rejects[F] = 0;

  /*
   * Compute space-view (SV) average DN, standard deviation and rejects.
   * If the moon is in the SV port, use the "lowest N" algorithm, with all
   * 50 SV frames.  If the moon is not in the SV KOB, use the normal
   * algorithm.
   */

  if (moon_in_SV_KOB == MOON_INSIDE_SV_KOB)
  {
    returnStatus =
       Get_DN_Avg_SDev_Rejects_LowN(
                            tables->SV_DN_moon_include_frames,
                            0,
                            SV_1km_FRAMES,
                            1,
                            SV,
                            SATURATED_DN - 1,
                            0,
                            &DN_sv,
                            &sdev,
                            &rejects[0]);

    if (returnStatus != MODIS_S_OK)
      SMF_ERROR(returnStatus,
      "Get_DN_Avg_SDev_Rejects_LowN() in "
          "Get_Emiss_Coeff_Per_Scan(), Preprocess.c");

  }
  else
  {
    returnStatus =
       Get_DN_Avg_SDev_Rejects(tables->SV_DN_first_frame_to_use,
                               tables->SV_DN_number_of_frames_to_use,
                               1,
                               SV,
                               SATURATED_DN - 1,
                               0,
                               &DN_sv,
                               &sdev,
                               &rejects[tables->SV_DN_first_frame_to_use]);

    if (returnStatus != MODIS_S_OK)
      SMF_ERROR(returnStatus,
        "Get_DN_Avg_SDev_Rejects() in "
          "Get_Emiss_Coeff_Per_Scan(), Preprocess.c");
  }

  returnStatus = Pack_Rejects_In_Outlier_Mask(SV_1km_FRAMES,
                                              rejects,
                                              sv_omask);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
      "Pack_Rejects_In_Outlier_Mask() in "
          "Get_Emiss_Coeff_Per_Scan(), Preprocess.c");

  /*
   * Check for a failure to calculate average SV DN.  A failure is not
   * fatal.  It just means that "b1" for this scan will not be calculated.
   * Return the value of MODIS_S_OK to allow the code to keep going.
   */

  if (DN_sv < 0)
    return (MODIS_S_OK);

  *xdnsv     = DN_sv;
  *xdnsv_var = sdev * sdev;

  /*
   * Compute black-body (BB) average DN.
   */

  returnStatus =
     Get_DN_Avg_SDev_Rejects(tables->BB_DN_first_frame_to_use,
                             tables->BB_DN_number_of_frames_to_use,
                             1,
                             BB,
                             SATURATED_DN - 1,
                             0,
                             &DN_bb,
                             &sdev,
                             &rejects[tables->BB_DN_first_frame_to_use]);

  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
      "Get_DN_Avg_SDev_Rejects() in "
      "Get_Emiss_Coeff_Per_Scan(), Preprocess.c");


  /*
   * BB WARMUP SATURATION FIX FOR FM1 (MODIS/AQUA)
   *
   * On FM1, the BB dns saturate prematurely on BB warmup and consequently
   * values of b1 cannot be calculated.  This problem currently (1/30/2002)
   * affects FM1 bands 33, 35, and 36.  The fix is to check the BB temperature
   * for these bands, and if it is above a band-specific threshold, use
   * a LUT value of b1 rather than attempt to calculate b1.
   */

  /*
   * Check to see if satellite is Aqua and the BB_T_sat_switch_aqua
   * correction switch is on.
   */

  if (BB_T_sat_switch_aqua == ON && satellite_ID == AQUA)
  {
    /* See if this is one of the bands which needs attention. */
    for (i = 0; i < NUM_AQUA_BB_SAT_BANDS; i++)
    {
      if (B_emiss == bb_sat_bands_indices_aqua[i])
      {
        /* The Blackbody temperature threshold: */

        BB_temp_threshold = tables->BB_T_sat_aqua[i];

        /* Is the average blackbody temperature over the threshold? */

        if (T_bb > BB_temp_threshold)
        {
           /* If so, use the LUT value of b1 and return. */
           /* 07/21/2011 apply temperature correction to the default b1. by X. Geng */

          *Xb1 = tables->BB_T_sat_default_b1_baseline_aqua[i][D][MS]*
                 (1. + tables->BB_T_sat_default_b1_c1_aqua[i][D][MS]*
                  (T_fp_lwir - tables->BB_T_sat_default_b1_Tlwir_baseline_aqua));

          return (MODIS_S_OK);

        }  /* END if (T_bb >> BB_temp_threshold) */

      }    /* END if (B_emiss == bb_sat_bands_indices_aqua[i]) */

    }      /* END for (i = 0; i < NUM_AQUA_BB_SAT_BANDS; i++) */

  }        /* END if (BB_T_sat_switch_aqua == ON && satellite_ID == AQUA  */

  /*
   * END BB WARMUP SATURATION FIX FOR FM1 (MODIS/AQUA)
   */

  /*
   * Check for a failure to calculate average BB DN.  A failure is not
   * fatal.  It just means that "b1" for this scan will not be calculated.
   * The value of MODIS_S_OK is returned to allow the code to keep going.
   * NOTE:  if DN_sv is OK but DN_bb is not (all were missing), this will
   * show up in the product as the "Could not compute zero point DN"
   * unusable value.
   */

  if (DN_bb < 0)
    return (MODIS_S_OK);

  /*
   * Compute dn_bb
   */

  dn_bb = DN_bb - DN_sv;

  /*
   * PCX correction
   */

  if ( PCX_correction_switch == ON && B >= BAND31)
  {

    /* For band 31, Xdn_bb_31 is an output. So store dn bb in this address.
     * For bands 32-36, Xdn_bb_31 is an input and the value stored
     * in this address is used for PC cross-talk correction.
     */

    if ( B == BAND31 ) *Xdn_bb_31 = dn_bb;
    else if( B >= BAND32 && B <= BAND36 )
    {
      dn_bb_31 = *Xdn_bb_31;
      B_xt = B - BAND32;
      dn_bb += (tables->PC_XT[B_xt][D][3] -
                tables->PC_XT[B_xt][D][0]*0.01*
                dn_bb_31*tables->PC_XT[B_xt][D][2]);
    }
  }

  if (dn_bb < TOLERANCE)     /* We cannot calibrate this scan */
    return returnStatus;

  if (MS == -1)
    return returnStatus;

  *xdnbb = dn_bb;

  /*
   * Calculate a0, a2
   * These are quadratics in instrument temperature (T_ins).
   */

  a0 = tables->A0[0][MS][D_emiss] + T_ins * (tables->A0[1][MS][D_emiss]
       + T_ins * tables->A0[2][MS][D_emiss]);

  a2 = tables->A2[0][MS][D_emiss] + T_ins * (tables->A2[1][MS][D_emiss]
       + T_ins * tables->A2[2][MS][D_emiss]);

  /*
   * Calculate dL_bb
   */

  delta_T_bb = tables->delta_T_bb_beta[D_emiss] * (T_cav - T_bb)
               + tables->delta_T_bb_delta[D_emiss];

  returnStatus = Calculate_Planck(tables->RSR[D_emiss],
                                  tables->wavelength[D_emiss],
                                  tables->NUM_RSR_vs_Lambda[D_emiss],
                                  T_mir, &L_sm);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
              "Calculate_Planck() in "
        "Get_Emiss_Coeff_Per_Scan()");

  returnStatus = Calculate_Planck(tables->RSR[D_emiss],
                                  tables->wavelength[D_emiss],
                                  tables->NUM_RSR_vs_Lambda[D_emiss],
                                  T_bb+delta_T_bb, &L_bb);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
              "Calculate_Planck() in "
        "Get_Emiss_Coeff_Per_SCan()");

  *xLbb = L_bb;

  returnStatus = Calculate_Planck(tables->RSR[D_emiss],
                                  tables->wavelength[D_emiss],
                                  tables->NUM_RSR_vs_Lambda[D_emiss],
                                  T_cav,
                                  &L_cav);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
              "Calculate_Planck() in "
        "Get_Emiss_Coeff_Per_SCan()");

  *xLcav = L_cav;

  bb_corr = RVS_Coeff->RVS_1km_Emiss_BB[D_emiss][MS];
  sv_corr = RVS_Coeff->RVS_1km_Emiss_SV[D_emiss][MS];

  dL_bb = bb_corr * tables->epsilon_bb[D_emiss] * L_bb +
            (sv_corr - bb_corr) * L_sm +
          (1.0 - tables->epsilon_bb[D_emiss]) *
          tables->epsilon_cav[D_emiss] *
            bb_corr * L_cav ;

  *Xb1 = (dL_bb - a0 - a2 * dn_bb * dn_bb) / dn_bb;
  *xdLbb = dL_bb;

  if (calculate_SNR == True)
  {
    /*
     *  The anomaly data should be removed in SNR calculation, 3/24/2011, xgeng
     */
    /*
    for ( x = 0; x < BB_1km_FRAMES; x++ )
    {
      sum1 += (double) (BB[x] - DN_sv);
      sum2 += (double) (BB[x] - DN_sv) * (double) (BB[x] - DN_sv);
    }
    mean = (float32) (sum1 / (double) x);
    sdev = (float32) sqrt(fabs((double) (sum2 - sum1*sum1/(double)x)) /
              (double)(x - 1));
    if ( sdev < 1.0E-19 )
      *SNR = -1.0;
    else
      *SNR = mean / sdev;
    */
    if ( sdev < 1.0E-19 )
      *SNR = -1.0;
    else
      *SNR = dn_bb / sdev;

  }

  return returnStatus;
}



PGSt_SMF_status Preprocess_L1A_Data (lookup_tables_t    *tables,
                                     L1A_granule_t      *L1A_Gran,
                                     L1B_granule_t      *L1B_Gran,
                                     QA_Data_t          *QA,
                                     Preprocess_Data_t  *PP)
/*
!C****************************************************************************
!Description:
   Read OBC/Eng data from L1A files and calculate a0, averaged b1, a2,
   delta_bb, SV counts, dnbb( BB counts subtract SV counts) for every scan in
   the current granule for emissive bands.  Calculate frame-averaged SV
   counts and standard deviation with outlier rejection for reflective
   bands. Write OBC target SDSs and preprocessed SDSs, and copy
   scan_level_metadata sdss, pixel_quality sdss, engineering_memory sdss,
   and engineering vdata from the middle L1A granule into the L1B_OBC file.

!Input Parameters:
   lookup_tables_t    *tables      Look-up table inputs.
   L1A_granule_t      *L1A_Gran    All data for the middle granule should
                                   be filled and file accesses (sd_id, v_id)
                                   should be open.
   L1B_granule_t      *L1B_Gran    Contains space for RVS correction terms
                                   which will be filled
!Output Parameters:
   QA_Data_t          *QA          All QAs
   Preprocess_Data_t  *PP          OBC/Eng averages/statistics

   Output file: L1B OBC/Eng file.

!Revision History:

   Revision 02.2  March 11, 2002   Razor Issue #174
   Passed L1B_Gran to Preprocess.c, calculating RVS correction terms before
     passing them to Emissive coefficient procedures.
   Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.12, Jan 26, 2001, Razor issue 152
   Moved Read Overlap lines inside function Process_OBCEng_Emiss.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.11 Sep 5, 2000
   Made several changes for Razor issue 130 (errors computing "b1" when
   near a sector rotation or if Ecal ends in the last 40 scans of the
   leading granule or starts in the 1st 40 scans of the trailing granule).
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   ... (many revisions not logged) ...

   Revision 02.10 Mar. 1998
   Added the QA_Data_t in the interface
   updated with the new Preprocess_Data_t
   Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

   Revision 01.00 Jan. 1997
   Initial development
   Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END**************************************************************************
*/
{
  PGSt_SMF_status       returnStatus      = MODIS_S_OK;
  L1A_granule_OBCEng_t  *L1A_OBCEng;
  L1A_OBCEng = (L1A_granule_OBCEng_t*) malloc(sizeof(L1A_granule_OBCEng_t));
  Moon_arrays_t         *moon_arrays;
  moon_arrays = (Moon_arrays_t*) malloc(sizeof(Moon_arrays_t));

  /* Read Extract Metadata */
  read_extract_metadata(L1A_Gran->sd_id,
                        &extract_pixel_offset,
                        &extract_pixel_count,
                        &extract_line_offset,
                        &extract_line_count);

  /*
   * Calculate the RVS Correction terms for use by Emissive
   * Coefficient routines and later in L1B.
   */
  returnStatus = Calculate_RVS_Correction (tables, L1B_Gran);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
               "Calculate_RVS_Correction() in Preprocess_L1A_Data()");


  /*
   * Read OBC & Engineering data from the middle granule.
   * These are data needed only within the Preprocess module.
   */

  returnStatus = Read_L1A_OBCEng (&tables->emiss, L1A_Gran, L1A_OBCEng);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
               "Read_L1A_OBCEng() in Preprocess_L1A_Data()");

  /*
   * If the NAD is closed, the DN values will be small and the dn_stars
   * could be small negative values due to noise. So, set the minimum
   * dn_star values to be -40 for all reflective bands if NAD is closed
   * for any scan within the granule.
   */

  returnStatus = Adjust_dn_star_Min(tables->refl.dn_star_Min,
                                    L1A_Gran->v_id,
                                    L1A_OBCEng->num_scans);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
               "Adjust_dn_star_Min() in Preprocess_L1A_Data()");

  /*
   * Check if moon is in SV KOB.  This reads data from geolocation file.
   */

  returnStatus = Check_For_Moon_in_SV_KOB(L1A_Gran->num_scans,
                                          &tables->QA.common_QA_tables,
                                          QA, moon_arrays);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
               "Check_For_Moon_in_SV_KOB() in Preprocess_L1A_Data()");

  /*
   * Process OBCEng data for Reflective_Cal.
   */

  returnStatus = Process_OBCEng_Refl (L1A_Gran->sd_id,
                                      L1A_Gran->num_scans,
                                      L1A_OBCEng->temps.ins,
                                      L1A_OBCEng->MirrorSide,
                                      &tables->refl,
                                      moon_arrays,
                                      L1A_OBCEng,
                                      &PP->DN_OBC_Avg,
                                      &QA->QA_refl,
                                      &PP->PP_Refl);

  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
               "Process_OBCEng_Refl() in "
               "Preprocess_L1A_Data(), Preprocess.c");

  /*
   * Process OBCEng data for Emissive_Cal.
   */


  returnStatus = Process_OBCEng_Emiss
                         (L1A_Gran,
                         L1A_OBCEng,
                         &L1B_Gran->Emiss_Cal_Coeff,
                         &tables->emiss,
                         &tables->QA,
                         moon_arrays,
                         QA,
                         &PP->PP_Emiss,
                         &PP->DN_OBC_Avg);

  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
               "Process_OBCEng_Emiss"
               "Preprocess_L1A_Data(), Preprocess.c");

  /*
   * Write most of the OBC file here.  Some metadata are added later in other
   * parts of the L1B code.
   */

  returnStatus = Write_L1B_OBCEng (L1A_Gran, tables,
                                   moon_arrays,
                                   &PP->DN_OBC_Avg);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
               "Write_L1B_OBCEng"
               "Preprocess_L1A_Data(), Preprocess.c");

  /*
   * Done with the Preprocessor.
   */

  free(L1A_OBCEng);
  free(moon_arrays);
  
  return returnStatus;
}


PGSt_SMF_status Process_OBCEng_Emiss
                       (L1A_granule_t          *L1A_Gran,
                        L1A_granule_OBCEng_t   *L1A_OBCEng,
                        Emiss_Cal_Coeff_t      *RVS_Coeff,
                        emiss_tables_t         *tables,
                        QA_tables_t            *QA_tables,
                        Moon_arrays_t          *moon_arrays,
                        QA_Data_t              *QA,
                        Preprocess_Emiss_t     *PP,
                        DN_OBC_Avg_t           *DN_OBC_Avg)
/*
!C****************************************************************************
!Description:
 Preprocess OBC and engineering data for Emissive Calibration. Read
 temperature data field, calculated temperatures averaged over a certain
 amount of scans, and the ratio of the averaged temperature variance to
 prelaunch values (Temperature QA), and the emissive calibration coefficients
 (a0, b1, a2, SV dn).

!Input Parameters:
 Overlap_OBCEng_t *     Leading_OBCEng  Contains number of scans, mirror side,
                                        BB and SV DNs and engineering data for
                                        the leading L1A granule.
 Overlap_OBCEng_t *     Trailing_OBCEng Contains number of scans, mirror side,
                                        BB and SV DNs and engineering data for
                                        the trailing L1A granule.
 L1A_granule_OBCEng_t * L1A_OBCEng      Contains number of scans, mirror side,
                                        BB and SV DNs and engineering data for
                                        the middle L1A granule.
 Emiss_Cal_Coeff_t     *RVS_Coeff       Emissive RVS Correction terms
                                        Array of RVS Correction terms for SV
 emiss_tables_t *       tables          Contains all the emissive calibration
                                        LUTs.
 QA_tables_t *          QA_tables       Contains all the emissive quality
                                        assurance LUTs.
 Moon_arrays_t *        moon_arrays     Contains the array of flags for the
                                        moon being in the SV keep-out box (KOB).

!Output Parameters:
 QA_Emiss_t *           QA              Structure containing quality assurance
                                        data to be written to the L1B products.
 Preprocess_Emiss_t *   PP              Structure containing emissive
                                        calibration coefficients and running
                                        averages of temperatures.
 DN_OBC_Avg_t *         DN_OBC_Avg      Structure containing data to fill in
                                        this routine.  The "1km_night" members
                                        are filled (except for band 26, which
                                        was filled in Process_OBCEng_Refl).

!Revision History:
 Revision 02.22, April 29, 2002  Razor Issue #180
 Moved initialization of Xb1 and XMS to inside detector loop.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.21, March 20, 2002  Razor Issue #178
 Removed ADC Correction 
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.20, March 10, 2002  Razor Issue #174
 Added passed parameter Emiss_Cal_Coeff_t
   in Get_Emiss_Coeff_Per_Scan
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.19, February 7, 2002  Razor Issue #180
 Moved initialization of Xb1 and XMS to inside band loop.
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.18, January 29, 2002  Razor Issue 175
 Create satellite ID and pass to
   Get_{Leading,Middle,Trailing}_Gran_Emiss_Coeff
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.17, Jan 26, 2001, Razor issue 152
 Moved Read Overlap lines inside here from parent function.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.16 Sep 5, 2000
 Made changes for Razor issue 130.  Removed Ecal_On from argument list
 and from calling argument lists of subordinate functions.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 ... (many revisions not logged) ...

 Revision 02.15 Oct. 22 1999
 Added checking for Ecal on
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.14 August 1999
 Removed variables and calculations for interpolation correction.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.13 August 1999
 See SDF, "L1B Code Change to meet ADC correction algorithm change"
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.12 May, 1999
 See SDF, "Stage 3: Changes to Emissive Algorithm and Code" notes.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.11  Feb 9, 1999
 Added use of variable ADC_correction_switch set to macro values to
 remove compiler warnings.  Macros in upper case.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 Apr. 1998
 Changed the dynamic to static allocation for the XGranule arrays
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Remove the D_inv, because it has already been taken care of in L1A
 Implemented V2.1 DN vs. L algorithm data preprocessing.  Implemented the
 different number of scans XGranule averaging for DN and Temperatures
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 01.00 Feb. 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
 For all cross-granule arrays, a value of -1 is assigned initially as a
 value that represents a missing or invalid value.  There are two ways that
 elements may remain at this invalid value: (1) if the leading or trailing
 granule is missing, or (2) if scans are missing from the data.
 Additionally, a value of -1 is assigned to the mirror side arrays, also to
 denote missing data.

 Variables with "X" at the beginning of the name span the overlap regions of
 the leading and trailing granules as well as the middle granule.

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16           R            = INDEX_1000M_EMISS;
  int16           B_emiss      = 0;
  int16           D_emiss      = 0;
  int16           B            = 0;
  int16           D            = 0;
  int16           *XMS;
  XMS = (int16*)malloc(MAX_TOTAL_XGRAN_SCANS * sizeof(int16));
  float32         *Xb1;
  Xb1 = (float32*)malloc(MAX_TOTAL_XGRAN_SCANS * sizeof(float32));
  typedef float32 Xdn_bb_31_t[DETECTORS_PER_1KM_BAND];
  Xdn_bb_31_t     *Xdn_bb_31;
  Xdn_bb_31 = (Xdn_bb_31_t*) malloc(MAX_TOTAL_XGRAN_SCANS*sizeof(Xdn_bb_31_t));
  int32           s, fp, n;
  int32           satellite_ID = INVALID_SATELLITE_ID;
                         /* 0 = TERRA (PFM);   1 = AQUA (FM1) */

  boolean         dropped_scans = False;

  Overlap_OBCEng_t      *Leading_OBCEng;
  Leading_OBCEng = (Overlap_OBCEng_t*) malloc(sizeof(Overlap_OBCEng_t));
  Overlap_OBCEng_t      *Trailing_OBCEng;
  Trailing_OBCEng = (Overlap_OBCEng_t*) malloc(sizeof(Overlap_OBCEng_t));

  /* Determine the satellite platform. */
  satellite_ID = L1A_Gran->satellite_id;

  /*
   * Check the macro MAX_TOTAL_XGRAN_SCANS
   */

  n = L1A_OBCEng->num_scans + 2 * tables->num_overlap_scans_b1;
  if (n > MAX_TOTAL_XGRAN_SCANS)
    SMF_ERROR(MODIS_F_NOK,
      "MAX_TOTAL_XGRAN_SCANS is too small in "
      "Process_OBCEng_Emiss(), Preprocess.c");

  /*
   * Read from overlapping granules SV & BB(Emissive only), etc.
   */

  dropped_scans = False;    /* Initialize for extra insurance */

  returnStatus = Read_Overlap_OBCEng (L1A_Gran,
                                     tables,
                                     LEADING_L1A_GRANULE,
                                     Leading_OBCEng,
                                     &dropped_scans);

  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
      "Read_Overlap_OBCEng() (Leading granule) in "
      "Process_OBCEng_Emiss(), Preprocess.c");

  /*
   * If any scans were dropped between the leading and middle granule, set
   *     a QA flag:
   */

  QA->QA_common.leading_granule_scan_gap = dropped_scans;


  dropped_scans = False;    /* Initialize for extra insurance */


  returnStatus = Read_Overlap_OBCEng (L1A_Gran,
                                     tables,
                                     TRAILING_L1A_GRANULE,
                                     Trailing_OBCEng,
                                     &dropped_scans);

  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,
      "Read_Overlap_OBCEng() (Trailing granule) in "
      "Process_OBCEng_Emiss(), Preprocess.c");

  /*
   * If any scans were dropped between the leading and middle granule, set
   *     a QA flag:
   */

  QA->QA_common.trailing_granule_scan_gap = dropped_scans;

  /*
   * Set the QA flags for missing leading and trailing granules.
   */

  if (Leading_OBCEng->num_scans == 0)
    QA->QA_common.missing_leading_granule = True;
  else
    QA->QA_common.missing_leading_granule = False;
  if (Trailing_OBCEng->num_scans == 0)
    QA->QA_common.missing_trailing_granule = True;
  else
    QA->QA_common.missing_trailing_granule = False;

  /*
   * Assign quantities to the Preprocess_Emiss_t structure from the per-scan
   * values in the middle granule.  Also, calculate the granule average
   * temperatures which are written to ECS core metadata (in Metadata.c).
   */


  for (s = 0; s < L1A_OBCEng->num_scans; s++) {
    PP->fp_set_point_state[s] = L1A_OBCEng->temps.fp_set_point_state[s];
    PP->T_bb[s] = L1A_OBCEng->temps.bb_avg[s];
    PP->T_mir[s] = L1A_OBCEng->temps.mir_avg[s];
    PP->T_ins[s] = L1A_OBCEng->temps.ins[s];
    for (fp = 0; fp < NUM_FOCAL_PLANES; fp++)
      PP->T_fp[fp][s] = L1A_OBCEng->temps.fp[fp][s];
  }
  returnStatus = Granule_Average_Temperature(L1A_OBCEng->num_scans,
                                             L1A_OBCEng->temps.mir_avg,
                                             &PP->T_mir_gran);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Granule_Average_Temperature() in "
        "Get_Middle_Gran_Temperatures()");

  returnStatus = Granule_Average_Temperature(L1A_OBCEng->num_scans,
                                             L1A_OBCEng->temps.bb_avg,
                                             &PP->T_bb_gran);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Granule_Average_Temperature() in "
        "Get_Middle_Gran_Temperatures()");

  returnStatus = Granule_Average_Temperature(L1A_OBCEng->num_scans,
                                             L1A_OBCEng->temps.fp[0],
                                             &PP->T_fp1);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Granule_Average_Temperature() in "
        "Get_Middle_Gran_Temperatures()");

  returnStatus = Granule_Average_Temperature(L1A_OBCEng->num_scans,
                                             L1A_OBCEng->temps.fp[1],
                                             &PP->T_fp2);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Granule_Average_Temperature() in "
        "Get_Middle_Gran_Temperatures()");

  returnStatus = Granule_Average_Temperature(L1A_OBCEng->num_scans,
                                             L1A_OBCEng->temps.fp[2],
                                             &PP->T_fp3);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Granule_Average_Temperature() in "
        "Get_Middle_Gran_Temperatures()");

  returnStatus = Granule_Average_Temperature(L1A_OBCEng->num_scans,
                                             L1A_OBCEng->temps.fp[3],
                                             &PP->T_fp4);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Granule_Average_Temperature() in "
        "Get_Middle_Gran_Temperatures()");

/*-----------------------------------------
   Compute temperature QA over a Granule
  -----------------------------------------*/
  returnStatus = Calculate_Temp_QA(L1A_OBCEng->num_scans,
                                   &L1A_OBCEng->temps, QA_tables, QA);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
              "Calculate_Temp_QA() in Process_OBCEng_Emiss()");

/*-------------------------------
    Calculate  PP->Planck_mir
  -------------------------------*/

  returnStatus = Calculate_PP_Planck_Mir(L1A_OBCEng->num_scans, tables, PP);
  if (returnStatus != MODIS_S_OK)
     SMF_ERROR(returnStatus,
        "Calculate_PP_Planck_Mir() in Process_OBCEng_Emiss()");


/*------------------------------------------------
  Compute PP->b1 through XGranule running average
          PP->a0
          PP->a2
          PP->DELTAbb
          PP->dbb
 ------------------------------------------------*/

  B_emiss = 0;
  D_emiss = 0;
  for (B = 0; B < L1A_BANDS_AT_RES[R]; B++)
  {
    if (B == BAND26)  continue;     /*skip band_26 */

    for (D = 0; D < DETECTORS_PER_1KM_BAND; D++)
    {

      /*
       * Initialize the X arrays to the invalid value (used as a flag to
       * avoid  averaging data from scans that are missing or otherwise could
       * not calculate b1).
       *   XMS = Mirror side flags across three granules
       *   Xb1 = b1 values across three granules (before averaging)
       */

      for (s = 0; s < MAX_TOTAL_XGRAN_SCANS; s++)
      {
        Xb1[s] = -1;
        XMS[s] = -1;
      }

     /* LeadingOverlap */
      returnStatus = Get_Leading_Gran_Emiss_Coeff
                                 (Leading_OBCEng,
                                  B,
                                  D,
                                  D_emiss,
                                  RVS_Coeff,
                                  tables,
                                  moon_arrays,
                                  Xdn_bb_31,
                                  Xb1,
                                  XMS,
                                  satellite_ID);

      if (returnStatus != MODIS_S_OK)
         SMF_ERROR(returnStatus,
            "Get_Leading_Gran_Emiss_Coeff() in "
            "Process_OBCEng_Emiss(), Preprocess.c");

     /* Middle granule
      */

      returnStatus = Get_Middle_Gran_Emiss_Coeff
                                (L1A_OBCEng,
                                 B,
                                 D,
                                 D_emiss,
                                 RVS_Coeff,
                                 tables,
                                 &QA_tables->emiss_QA_tables,
                                 moon_arrays,
                                 Xdn_bb_31,
                                 Xb1,
                                 XMS,
                                 &QA->QA_emiss,
                                 PP,
                                 DN_OBC_Avg,
                                 satellite_ID);
      if (returnStatus != MODIS_S_OK)
         SMF_ERROR(returnStatus,
            "Get_Middle_Gran_Emiss_Coeff() in "
            "Process_OBCEng_Emiss(), Preprocess.c");

      /* TrailingOverlap */

      returnStatus = Get_Trailing_Gran_Emiss_Coeff
                                  (L1A_OBCEng->num_scans,
                                   Trailing_OBCEng,
                                   B,
                                   D,
                                   D_emiss,
                                   RVS_Coeff,
                                   tables,
                                   moon_arrays,
                                   Xdn_bb_31,
                                   Xb1,
                                   XMS,
                                   satellite_ID);

      if (returnStatus != MODIS_S_OK)
         SMF_ERROR(returnStatus,
            "Get_Trailing_Gran_Emiss_Coeff() in "
            "Process_OBCEng_Emiss(), Preprocess.c");

      /* Calculate Xgranule running average of b1: PP->b1
      */
      returnStatus = Get_All_Emiss_Coeff(L1A_OBCEng,
                                         B,
                                         D_emiss,
                                         Xb1,
                                         XMS,
                                         tables,
                                         &QA_tables->emiss_QA_tables,
                                         &QA->QA_emiss,
                                         PP);
      if (returnStatus != MODIS_S_OK)
         SMF_ERROR(returnStatus,
            "Get_All_Emiss_Coeff() in "
            "Process_OBCEng_Emiss(), Preprocess.c");
      D_emiss++;
       }/*End of loop over Detectors per band*/

    B_emiss++;
  }/*End of loop over bands*/

  free(XMS);
  free(Xb1);
  free(Xdn_bb_31);
  free(Leading_OBCEng);
  free(Trailing_OBCEng);

  return returnStatus;
}

PGSt_SMF_status Calculate_PP_Planck_Mir(int32               num_scans,
                                        emiss_tables_t      *tables,
                                        Preprocess_Emiss_t  *PP)

/*
!C************************************************************************
!Description:  Calculate radiance of mirror averaged over wavelengths. The
               value will be used in Emissive_Cal().

!Input Parameters:
    int32                  num_scans
    emiss_tables_t *       tables
!Output Parameters:
    Preprocess_Emiss_t *   PP

!Revision History:
 Revision 01.00 Oct. 1998
 Initial development
 (This code is part of the original function Process_OBCEng_Emiss written
  by Zhidong Hao and modified by Shi-Yue Qiu)
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 This function is part of the original function Process_OBCEng_Emiss
 written by: Zhidong Hao and modified by Shi-Yue Qiu

 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status  returnStatus = MODIS_S_OK;
  int16 B          = 0;
  int16 D          = 0;
  int16 S          = 0;
  int16 D_emiss    = 0;

  for (B = 0; B < NUM_EMISSIVE_BANDS; B++)
  {
    for (D = 0; D < DETECTORS_PER_1KM_BAND; D++)
    {
      for (S = 0; S < num_scans; S++)
      {

        returnStatus = Calculate_Planck(tables->RSR[D_emiss],
                                        tables->wavelength[D_emiss],
                                        tables->NUM_RSR_vs_Lambda[D_emiss],
                                        PP->T_mir[S],
                                        &PP->Planck_mir[D_emiss][S]);
        if (returnStatus != MODIS_S_OK)
          SMF_ERROR(returnStatus,
                    "Calculate_Planck() in Process_OBCEng_Emiss()");
      }
      D_emiss++;
    }
  }
  return (returnStatus);
}

PGSt_SMF_status Get_Leading_Gran_Emiss_Coeff
                (Overlap_OBCEng_t  *Leading_OBCEng,
                 int16             B,
                 int16             D,
                 int16             D_emiss,
                 Emiss_Cal_Coeff_t *RVS_Coeff,
                 emiss_tables_t    *tables,
                 Moon_arrays_t     *moon_arrays,
                 float32           Xdn_bb_31[][DETECTORS_PER_1KM_BAND],
                 float32           *Xb1,
                 int16             *XMS,
                 int32              satellite_ID)

/*
!C****************************************************************************
!Description:
 This function calculates the linear calibration coefficient, b1, for the last
 number of overlap scans of the leading L1A granule.  The values are placed in
 the first set of overlap scans of the return array Xb1.  The remaining scan
 elements of Xb1 are filled within other functions. The mirror side flags for
 the scans of Xb1 are also saved and returned in XMS.  If this is band 31, then
 BB dn values are saved in Xdn_bb_31.  Otherwise, those values are input for
 bands 32-36.  If this is band 33, 35, or 36, the BB temperature is checked
 against a threshold value for the band; if the temperature is over the
 threshold, a default b1 value is used for the b1 value.

!Input Parameters: (Prefix "X" indicates extension cross 3 granules)
 Overlap_OBCEng_t * Leading_OBCEng Contains number of scans in leading granule,
                                   SV and BB DNs, mirror side and engineering
                                   temperatures.
 int16              B              Band index within the resolution.
 int16              D              Detector index within one emissive band.
 int16              B_emiss        Band index within the set of emissive bands.
 int16              D_emiss        Detector index within the set of emissive
                                   detectors.
 Emiss_Cal_Coeff_t *RVS_Coeff      Emissive RVS Correction terms
 emiss_tables_t *   tables         emissive lookup tables, contains ADC
                                   coefficients and the number of overlap scans
                                   for b1.
 Moon_arrays_t *    moon_arrays    Contains the array of flags for the moon
                                   being in the SV keep-out box (KOB) for the middle
                                   granule.
 float32 *          Xdn_bb_31[][]  Array of blackbody dn for each detector in
                                   band 31 and for each scan. It is an input for
                                   bands 32-36 and an output for band 31.

!Output Parameters: (Prefix "X" indicates extension across 3 granules)
 float32 *          Xdn_bb_31[][]  (see above)
 float32 *          Xb1            Per-scan values of b1 in the first number
                                   of overlap scans elements of extended array.
 int16 *            XMS            Mirror side flags across granule, matching the
                                   scans in Xb1.
 int32              satellite_ID   0 = TERRA, 1 = AQUA, -1 = INVALID_SATELLITE_ID

!Revision History:
 Revision 02.20 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.03  March 10, 2002   Razor Issue #174
 Added passed parameter Emiss_Cal_Coeff_t
   in Get_Emiss_Coeff_Per_Scan
 Alice Isaacman  SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.10, January 29, 2002, Razor Issue 175
 Input satellite_ID and passed to Get_Emiss_Coeff_Per_Scan.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.02, Dec 19, 2000, Razor issue 148
 Delete XT arrays, use values in temps structure (per-scan values)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.01 Sep 11, 2000
 Made changes for Razor issue 133 (TEB algorithm for calculating <DN-SV> when
 moon is in the SVP).  Now, pass flag denoting if moon is in the SV KOB into
 Get_Emiss_Coeff_Per_Scan.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.00 Sep 5, 2000
 Made several changes for Razor issue 130.  Removed Ecal_On from argument
 list.  Added check for Ecal_On and Sector_Rotated for the correct scan prior
 to calling Get_Emiss_Coeff_Per_Scan.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 ... (many revisions not logged) ...

 Revision 01.04 Oct. 22, 1999
 Added checking for Ecal on
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 August 1999
 Made corresponding change due to the remove of xdnbb in
   Get_Emiss_Coeff_Per_Scan()
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.02 August 1999
 See SDF, "L1B Code Change to meet ADC correction algorithm change"
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 May, 1999Revision 01.02 August 1999
 See SDF, "L1B Code Change to meet ADC correction algorithm change"
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 See SDF, "Stage 3: Changes to Emissive Algorithm and Code" notes.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 Oct. 1998
 Initial development
 (This code is part of the original function Process_OBCEng_Emiss written
  by Zhidong Hao and modified by Shi-Yue Qiu)
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
 If the number of scans in the leading granule is set to zero, then it is
 assumed that the leading granule is undefined or missing.  This function
 simply returns normally to allowed continued processing.

 On a scan-by-scan basis, we do not compute emissive coefficients if:
 (1) Ecal is on, (2) there is a sector
 rotation on this scan or (3) the mirror side flag is -1 (which usually
 means that all data are missing for the scan).

 We determine if the moon is in the SV KOB and pass a flag into
 Get_Emiss_Coeff_Per_Scan.  For determining if the moon is in the SVP, we use
 the 1st scan of the middle granule for all scans of the leading granule.

!END**************************************************************************
*/
{
  PGSt_SMF_status       returnStatus = MODIS_S_OK;
  int16 S               = 0;
  int16 start_scan      = 0;
  int16 MS              = 0;
  int16 XS              = 0;
  int16 track_index     = 0;
  float32 SNR           = 0;
  float32 xdLbb         = 0;
  float32 xdnbb         = 0;
  float32 xLbb          = 0;
  float32 xLcav         = 0;
  float32 xdnsv         = 0;
  float32 xdnsv_var;
  uint32  sv_omask[2];
  int16   B_38;


  /*
   * If the number of scans in the leading granule is zero, then treat
   * the leading granule as missing.  Simply return and continue
   * processing normally.
   */


  if (Leading_OBCEng->num_scans == 0)
    return returnStatus;

  /*
   * start_scan  = absolute scan index in the leading granule for the
   *                  1st overlap scan.
   * XS          = index in the set of overlap scans, starting at 0
   *                  for the leading granule.
   * S           = absolute scan index in the leading granule that
   *                  corresponds to XS
   * B_38        = index of this band in the set of 38 MODIS bands.
   * track_index = index in the BB and SV arrays that corresponds to
   *                  XS and S (the first track index in those arrays
   *                  corresponds to the start_scan index)
   */

  start_scan = Leading_OBCEng->num_scans - tables->num_overlap_scans_b1;
  B_38 = MODIS_BAND20_INDEX + B;

  for (XS = 0; XS < tables->num_overlap_scans_b1; XS++)
  {
    S  = start_scan + XS;
    MS = Leading_OBCEng->MirrorSide[S];
    XMS[XS] = MS;
    track_index = XS * DETECTORS_PER_1KM_BAND + D;

    SNR = 0.0;
    if ((MS == 0 || MS == 1)
        && Leading_OBCEng->Ecal_On[S][B_38] == False
        && Leading_OBCEng->Sector_Rotated[S] == False 
        && Leading_OBCEng->temps.pclw_adc_index[S] != ELECTRONICS_BOTH)
    {
      returnStatus =
        Get_Emiss_Coeff_Per_Scan
                            (moon_arrays->moon_in_SV_KOB[0][B_38],
                             Leading_OBCEng->BB_1km_night[track_index][B],
                             Leading_OBCEng->SV_1km_night[track_index][B],
                             B,
                             D,
                             D_emiss,
                             MS,
                             RVS_Coeff,
                             tables,
                             Leading_OBCEng->temps.bb_avg[S],
                             Leading_OBCEng->temps.mir_avg[S],
                             Leading_OBCEng->temps.cav[S],
                             Leading_OBCEng->temps.ins[S],
                             Leading_OBCEng->temps.fp[FP_LWIR][S],
                             &Xdn_bb_31 [XS][D],
                             &Xb1       [XS],
                             &xdLbb,
                             &xdnbb,
                             &xLbb,
                             &xLcav,
                             &xdnsv,
                             &xdnsv_var,
                             sv_omask,
                             &SNR,
                             satellite_ID);

      if (returnStatus != MODIS_S_OK)
        SMF_ERROR(returnStatus, "Get_Emiss_Coeff_Per_Scan() in "
                                "Get_Leading_Gran_Emiss_Coeff()");
    }
  }

  return (returnStatus);
}

PGSt_SMF_status Get_Middle_Gran_Emiss_Coeff
                (L1A_granule_OBCEng_t   *L1A_OBCEng,
                 int16                  B,
                 int16                  D,
                 int16                  D_emiss,
                 Emiss_Cal_Coeff_t      *RVS_Coeff,
                 emiss_tables_t         *tables,
                 emiss_QA_tables_t      *QA_tables,
                 Moon_arrays_t          *moon_arrays,
                 float32                Xdn_bb_31[][DETECTORS_PER_1KM_BAND],
                 float32                *Xb1,
                 int16                  *XMS,
                 QA_Emiss_t             *QA,
                 Preprocess_Emiss_t     *PP,
                 DN_OBC_Avg_t           *DN_OBC_Avg,
                 int32                  satellite_ID)

/*
!C****************************************************************************
!Description:
 This function calculates the linear calibration coefficient, b1, for the
 middle L1A granule.  The values are placed in the in Xb1 following the
 overlap scans values filled by Get_Leading_Gran_Emiss_Coeff.  This function
 also calculates averaged SV DN for emissive processing. The mirror side flags for
 the scans of Xb1 are also saved and returned in XMS.  If this is band 31, then
 BB dn values are saved in Xdn_bb_31.  Otherwise, those values are input for
 bands 32-36.


!Input Parameters:
 L1A_granule_OBCEng_t *L1A_OBCEng     Contains number of scans in middle
                                      granule, SV and BB DNs, mirror side and
                                      engineering temperatures.
 int16                 B              Band index within the resolution.
 int16                 D              Detector index within one emissive band.
 int16                 B_emiss        Band index within the set of emissive bands.
 int16                 D_emiss        Detector index within the set of
                                      emissive detectors.
 Emiss_Cal_Coeff_t    *RVS_Coeff      Emissive RVS Correction terms
 emiss_tables_t       *tables         emissive lookup tables, contains ADC
                                      coefficients and the number of overlap
                                      scans for b1.
 emiss_QA_tables      *QA_tables      Contains prelaunch NEdL.
 Moon_arrays_t        *moon_arrays    Contains the array of flags for the moon
                                      being in the SV keep-out box (KOB) for the
                                      middle granule.
 float32 *             Xdn_bb_31[][]  Array of blackbody dn for each detector in
                                      band 31 and for each scan. It is an input
                                      for bands 32-36 and an output for band31.
 int32                 satellite_ID   0 = TERRA, 1 = AQUA, -1 = INVALID_SATELLITE_ID

!Output Parameters:
 float32 *             Xdn_bb_31[][]  (see above)
 float32 *             Xb1            Per-scans values of b1 in the first
                                      number of overlap scans elements of
                                      extended array.
 int16 *               XMS            Mirror side flags across granule
                                      matching the scans in Xb1.
 QA_Emiss_t *          QA             Hold the NEdL calculated here.
 Preprocess_Emiss_t *  PP             members filled include DN_sv.
 DN_OBC_Avg_t *        DN_OBC_Avg     Members filled include the 1km_night
                                      arrays.

!Revision History:
 Revision 02.10 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.03  March 11, 2002   Razor Issue #174
 Added passed parameter Emiss_Cal_Coeff_t
   and added to call to Get_Emiss_Coeffs_Per_Scan
 Alice Isaacman  SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.02, January 29, 2002  Razor Issue #175
 Added satellite_ID as input and passed to Get_Emiss_Coeffs_Per_Scan
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.01 Sep 11, 2000
 Made changes for Razor issue 133 (TEB algorithm for calculating <DN-SV> when
 moon is in the SVP).  Now, pass flag denoting if moon is in the SV KOB into
 Get_Emiss_Coeff_Per_Scan.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.00 Sep 5, 2000
 Made several changes for Razor issue 130.  Removed Ecal_On from argument
 list.  Added check for Ecal_On and Sector_Rotated for the correct scan prior
 to calling Get_Emiss_Coeff_Per_Scan.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 ... (many revisions not logged) ...

 Revision 01.04 Oct. 22, 1999
 Added checking for Ecal on
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 August 1999
 Removed calculation of delta radiance and DN for the black body.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.02 August 1999
 See SDF, "L1B Code Change to meet ADC correction algorithm change"
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 May, 1999
 See SDF, "Stage 3: Changes to Emissive Algorithm and Code" notes.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 Oct. 1998
 Initial development
 (This code is part of the original function Process_OBCEng_Emiss written
  by Zhidong Hao and modified by Shi-Yue Qiu)
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
 Prefix "X" indicates extension cross 3 granules.
 On a scan-by-scan basis, we do not compute emissive coefficients if:
 (1) Ecal is on, (2) there is a sector
 rotation on this scan or (3) the mirror side flag is -1 (which usually
 means that all data are missing for the scan).

 We determine if the moon is in the SV KOB and pass a flag into
 Get_Emiss_Coeff_Per_Scan.

!END**************************************************************************
*/
{
  PGSt_SMF_status   returnStatus = MODIS_S_OK;
  int16 rejects[MAX_1KM_OBC_FRAME_DIM];
  int16 S           = 0;
  int16 MS          = 0;
  int16 XS          = 0;
  int16 track_index = 0;
  float32 SNR       = 0;
  float32 sum_NEdL  = 0;
  float32 xdLbb     = 0;
  float32 xdnbb     = 0;
  float32 xLbb      = 0;
  float32 xLcav     = 0;
  float32 xdnsv     = 0;
  float32 xdnsv_var = 0;
  float32 ratio     = 0;
  int32   n;
  int16   B_38;
  int16   F;

  n = 0;
  B_38 = MODIS_BAND20_INDEX + B;
  for (S = 0; S < L1A_OBCEng->num_scans; S++)
  {
    MS = L1A_OBCEng->MirrorSide[S];
    XS = tables->num_overlap_scans_b1 + S;
    XMS[XS] = MS;
    track_index = S * DETECTORS_PER_1KM_BAND + D;
    SNR = 1.0;
    if (L1A_OBCEng->Ecal_On[S][B_38] == False
        && L1A_OBCEng->Sector_Rotated[S] == False
        && L1A_OBCEng->temps.pclw_adc_index[S] != ELECTRONICS_BOTH
        && (MS == 0 || MS == 1))
    {
        returnStatus =
            Get_Emiss_Coeff_Per_Scan
                        (moon_arrays->moon_in_SV_KOB[S][B_38],
                        L1A_OBCEng->BB_1km_night[track_index][B],
                        L1A_OBCEng->SV_1km_night[track_index][B],
                        B,
                        D,
                        D_emiss,
                        MS,
                        RVS_Coeff,
                        tables,
                        L1A_OBCEng->temps.bb_avg[S],
                        L1A_OBCEng->temps.mir_avg[S],
                        L1A_OBCEng->temps.cav[S],
                        L1A_OBCEng->temps.ins[S],
                        L1A_OBCEng->temps.fp[FP_LWIR][S],
                        &Xdn_bb_31 [XS][D],
                        &Xb1       [XS],
                        &xdLbb,
                        &xdnbb,
                        &xLbb,
                        &xLcav,
                        &xdnsv,
                        &xdnsv_var,
                        DN_OBC_Avg->DN_obc_1km_night_outlier_mask
                                        [track_index][B][0],
                        &SNR,
                        satellite_ID);

      if (returnStatus != MODIS_S_OK)
        SMF_ERROR(returnStatus,
            "Get_Emiss_Coeff_Per_Scan() in "
            "Get_Middle_Gran_Emiss_Coeff()");
    }
    else
    {
      for (F = 0; F < MAX_1KM_OBC_FRAME_DIM; F++)
        rejects[F] = 1;
      returnStatus = Pack_Rejects_In_Outlier_Mask
              (SV_1km_FRAMES,
              rejects,
              DN_OBC_Avg->DN_obc_1km_night_outlier_mask
                  [track_index][B][0]);

      if (returnStatus != MODIS_S_OK)
        SMF_ERROR(returnStatus,
            "Pack_Rejects_In_Outlier_Mask() in "
            "Get_Middle_Gran_Emiss_Coeff()");

      xdLbb = -1;
      xdnbb = -1;
      xLbb  = -1;
      xLcav = -1;
      xdnsv = -1;
      xdnsv_var = -1;
      SNR = -1;
    }

    PP->dn_bb[D_emiss][S] = xdnbb;
    PP->L_bb[D_emiss][S] = xLbb;
    PP->L_cav[D_emiss][S] = xLcav;
    PP->DN_sv[D_emiss][S] = xdnsv;
    DN_OBC_Avg->DN_obc_1km_night_avg[track_index][B][0] = xdnsv;
    DN_OBC_Avg->DN_obc_1km_night_var[track_index][B][0] = xdnsv_var;
    if ( SNR > 0 )
    {
      PP->dn_bb_sdev[D_emiss][S] = PP->dn_bb[D_emiss][S]/SNR;
      PP->NEdL[D_emiss][S] = xdLbb/SNR;
      sum_NEdL += xdLbb/SNR;
      n++;
    }
    else
    {
      PP->dn_bb_sdev[D_emiss][S] = 0.;
      PP->NEdL[D_emiss][S] = 0;
    }
  }/*End of loop over middle scans*/

  if (n == 0)
    ratio = 0.0;
  else
  {
    sum_NEdL /= (float32)n;
    ratio = sum_NEdL / QA_tables->NEdL[D_emiss];
  }
  CONVERT_TO_UINT8(ratio, QA->NEdL[D_emiss]);

  return (returnStatus);
}


PGSt_SMF_status Get_Trailing_Gran_Emiss_Coeff
                (int32              num_scans_middle,
                 Overlap_OBCEng_t   *Trailing_OBCEng,
                 int16              B,
                 int16              D,
                 int16              D_emiss,
                 Emiss_Cal_Coeff_t  *RVS_Coeff,
                 emiss_tables_t     *tables,
                 Moon_arrays_t      *moon_arrays,
                 float32            Xdn_bb_31[][DETECTORS_PER_1KM_BAND],
                 float32            *Xb1,
                 int16              *XMS,
                 int32              satellite_ID)


/*
!C****************************************************************************
!Description:
 This function calculates the linear calibration coefficient, b1, for the first
 number of overlap scans of the trailing L1A granule.  The values are placed
 in the last set of overlap scans of the return array Xb1.  The remaining
 scan elements of Xb1 were previously filled within other functions. The
 mirror side flags for the scans of Xb1 are also saved and returned in XMS.
 If this is band 31, then BB dn values are saved in Xdn_bb_31.  Otherwise,
 those values are input for bands 32-36.


!Input Parameters:
 int32              num_scans_middle Number of scans in the middle granule.
 Overlap_OBCEng_t * Trailing_OBCEng  Contains number of scans in trailing
                                     granule, SV and BB DNs, mirror side and
                                     engineering temperatures.
 int16              B                Band index within the resolution.
 int16              D                Detector index within one emissive band.
 int16              B_emiss          Band index within the set of emissive bands.
 int16              D_emiss          Detector index within the set of emissive
                                     detectors.
 Emiss_Cal_Coeff_t *RVS_Coeff        Emissive RVS Correction terms
 emiss_tables_t *   tables           emissive lookup tables, contains ADC
                                     coefficients and the number of overlap
                                     scans for b1.
 Moon_arrays_t *    moon_arrays      Contains the array of flags for the moon
                                     being in the SV keep-out box (KOB) for the
                                     middle granule.
 float32 *          Xdn_bb_31[][]    Array of blackbody dn for each detector in
                                     band 31 and for each scan. It is an input
                                     for bands 32-36 and an output for band 31.
 int32              satellite_ID     0 = TERRA, 1 = AQUA, -1 = INVALID_SATELLITE_ID
!Output Parameters:
 float32 *          Xdn_bb_31[][]    (see above)
 float32 *          Xb1              Per-scans values of b1 in the last number
                                     of overlap scans elements of extended
                                     array.
 int16 *            XMS              Mirror side flags across granule, matching the
                                     scans in Xb1.
 int32              satellite_ID     0 = TERRA, 1 = AQUA, -1 = INVALID_SATELLITE_ID

!Revision History:
 Revision 02.10 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.03  March 11, 2002   Razor Issue #174
 Added passed parameter RVS_Coeff
   and added to call to Get_Emiss_Coeffs_Per_Scan
 Alice Isaacman  SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.02 January 29, 2002   Razor Issue #175
 Added input parameter satellite_ID and passed to Get_Emiss_Coeff_Per_Scan
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.01 Sep 11, 2000
 Made changes for Razor issue 133 (TEB algorithm for calculating <DN-SV> when
 moon is in the SVP).  Now, pass flag denoting if moon is in the SV KOB into
 Get_Emiss_Coeff_Per_Scan.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.00 Sep 5, 2000
 Made several changes for Razor issue 130.  Removed Ecal_On from argument
 list.  Added check for Ecal_On and Sector_Rotated for the correct scan prior
 to calling Get_Emiss_Coeff_Per_Scan.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 ... (many revisions not logged) ...

 Revision 01.04 Oct. 22, 1999
 Added checking for Ecal on
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 August 1999
 Removed the calculation of delta radiance and DN for the black body.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.02 August 1999
 See SDF, "L1B Code Change to meet ADC correction algorithm change"
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 May, 1999
 See SDF, "Stage 3: Changes to Emissive Algorithm and Code" notes.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 Oct. 1998
 Initial development
 (This code is part of the original function Process_OBCEng_Emiss written
  by Zhidong Hao and modified by Shi-Yue Qiu)
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
 If the number of scans in the trailing granule is set to zero, then it is
 assumed that the trailing granule is undefined or missing.  This function
 simply returns normally to allowed continued processing.  The values in the
 cross-granule arrays remain unchanged from their entry state (the parent
 initializes these to an appropriate value).

 On a scan-by-scan basis, we do not compute emissive coefficients if:
 (1) Ecal is on, (2) there is a sector
 rotation on this scan or (3) the mirror side flag is -1 (which usually
 means that all data are missing for the scan).

 We determine if the moon is in the SV KOB and pass a flag into
 Get_Emiss_Coeff_Per_Scan.  For determining if the moon is in the SVP, we use
 the last scan of the middle granule for all scans of the leading granule.

!END**************************************************************************
*/
{
  PGSt_SMF_status       returnStatus = MODIS_S_OK;
  int16 S               = 0;
  int16 MS              = 0;
  int16 XS              = 0;
  int16 track_index     = 0;
  float32 SNR           = 0;
  float32 xdLbb         = 0;
  float32 xdnbb         = 0;
  float32 xLbb          = 0;
  float32 xLcav         = 0;
  float32 xdnsv         = 0;
  float32 xdnsv_var;
  uint32  sv_omask[2];
  int16 B_38;

  if (Trailing_OBCEng->num_scans == 0)
    return returnStatus;

  B_38 = MODIS_BAND20_INDEX + B;

  for (S = 0; S < tables->num_overlap_scans_b1; S++)
  {
     MS = Trailing_OBCEng->MirrorSide[S];
     XS = tables->num_overlap_scans_b1 + num_scans_middle + S;
     XMS[XS] = MS;
     track_index = S * DETECTORS_PER_1KM_BAND + D;

     SNR = 0.0;
     if (moon_arrays->moon_in_SV_KOB[num_scans_middle-1][B_38] !=
         MOON_INSIDE_SV_KOB
         && (MS == 0 || MS == 1)
         && Trailing_OBCEng->Ecal_On[S][B_38] == False
         && Trailing_OBCEng->Sector_Rotated[S] == False 
         && Trailing_OBCEng->temps.pclw_adc_index[S] != ELECTRONICS_BOTH)
     {
       returnStatus = Get_Emiss_Coeff_Per_Scan
                              (moon_arrays->moon_in_SV_KOB
                                               [num_scans_middle - 1][B_38],
                              Trailing_OBCEng->BB_1km_night[track_index][B],
                              Trailing_OBCEng->SV_1km_night[track_index][B],
                              B,
                              D,
                              D_emiss,
                              MS,
                              RVS_Coeff,
                              tables,
                              Trailing_OBCEng->temps.bb_avg[S],
                              Trailing_OBCEng->temps.mir_avg[S],
                              Trailing_OBCEng->temps.cav[S],
                              Trailing_OBCEng->temps.ins[S],
                              Trailing_OBCEng->temps.fp[FP_LWIR][S],
                              &Xdn_bb_31[XS][D],
                              &Xb1[XS],
                              &xdLbb,
                              &xdnbb,
                              &xLbb,
                              &xLcav,
                              &xdnsv,
                              &xdnsv_var,
                              sv_omask,
                              &SNR,
                              satellite_ID);

       if (returnStatus != MODIS_S_OK)
         SMF_ERROR(returnStatus,
             "Get_Emiss_Coeff_Per_Scan() in Process_OBCEng_Emiss()");
     }
   }/*End of loop over trailing scans*/

   return (returnStatus);
}

PGSt_SMF_status Get_All_Emiss_Coeff
                (L1A_granule_OBCEng_t    *L1A_OBCEng,
                 int16                   B,
                 int16                   D_emiss,
                 float32                 *Xb1,
                 int16                   *Xmir,
                 emiss_tables_t          *tables,
                 emiss_QA_tables_t       *QA_tables,
                 QA_Emiss_t              *QA,
                 Preprocess_Emiss_t      *PP)

/*
!C************************************************************************
!Description:  Calculate averaged b1, a0, a2 for each scan for the current
               granule(corresponding to the middle L1A granule) and the b1
               change since prelaunch.

!Input Parameters:
    L1A_granule_OBCEng_t   *L1A_OBCEng        Containing the number of scans
                                              and the mirror side for each
                                              scan for the middle granule.
    int16                  B                  Band index within the resolution.
    int16                  D_emiss            Detector index within the set of
                                              emissive detectors.
    float32 *              Xb1                Cross-granule per-scan values of
                                              b1 calculated in other functions.
    int16 *                Xmir               Cross-granule per-scan values of mirror
                                              side.
    emiss_tables_t *       tables             Containing b1 for band 21 and the
                                              A0 and A2 values for calculating a0
                                              and a2, respectively.
    emiss_QA_tables_t *    QA_tables          Containing the prelaunch values of b1,
                                              i.e. a1.
!Output Parameters:
    QA_Emiss_t *           QA                 Containing the value change of b1 since
                                              prelaunch.
    Preprocess_Emiss_t *   PP                 Store the coefficients b1, a0 and a2 of
                                              the middle granule.

!Revision History:
 Revision 02.19, October 15, 2004  Razor Issue #201
 Changed the for-loop of getting off-line calculated average b1 of Band 21
 to include the new dimension for Mirror Side specification.
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 02.18, Jan 29, 2001
 Improved efficiency of calculation of a0 and a2.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.17 Sep. 12, 2000
 Corrected the calculation of average b1 in case of missing scans
 (where b1 will be set to -1).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.16 August 23, 1999
 Removed calculation of the delta_bb.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.15 April 12, 1999
 Removed cubic term a3
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.00 Oct. 1998
 Initial development
 (This code is part of the original function Process_OBCEng_Emiss written
  by Zhidong Hao and modified by Shi-Yue Qiu)
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 This function is part of the original function Process_OBCEng_Emiss
 written by: Zhidong Hao and modified by Shi-Yue Qiu

 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

 "b1" must be a positive number to be valid.  In here (and elsewhere),
 a negative value of b1 is used as a flag to denote that it could not
 be calculated for some reason.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16    MS     = 0;
  int16    S      = 0;
  float32  avg_b1 = 0;
  float32  ratio  = 0;
  float32  T_ins  = 0;
  int16    ns_avg = 0;

  /* Calculate Xgranule running average of b1: PP->b1
   */

  if (B == BAND21)  /* b1(B21) is fixed */
  {
    /* use D_emiss - 10 for index to avoid passing in new argument */
    for (S = 0; S < L1A_OBCEng->num_scans; S++){
      MS = L1A_OBCEng->MirrorSide[S];
      PP->b1[D_emiss][S] = tables->Band_21_b1[D_emiss - 10][MS];
    }
  }
  else
  {
    returnStatus = Cross_Granule_Sliding_Average
                          (tables->num_overlap_scans_b1,
                           L1A_OBCEng->num_scans,
                           Xb1,
                           Xmir,
                           PP->b1[D_emiss]);
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR(returnStatus,
    "Cross_Granule_Sliding_Average(b1) in Process_OBCEng_Emiss()");
  }

  /* Calculate QA->change_b1 : b1/pre_launch_a1
   */

  avg_b1 = 0;
  ns_avg = 0;
  for ( S = 0; S < L1A_OBCEng->num_scans; S++ ) {
    if (PP->b1[D_emiss][S] > 0) {
      avg_b1 += PP->b1[D_emiss][S];
      ns_avg++;
    }
  }
  if (ns_avg > 0)
    avg_b1 = avg_b1/ns_avg;

  ratio = avg_b1 / QA_tables->a1[D_emiss];
  CONVERT_TO_UINT8(ratio, QA->change_b1[D_emiss]);

  /* Calculate PP->a0, PP->a2
   */
  for (S = 0; S < L1A_OBCEng->num_scans; S++)
  {
    MS = L1A_OBCEng->MirrorSide[S];

    if(PP->T_ins[S] > TOLERANCE && (MS == 0 || MS == 1))
    {
      T_ins = PP->T_ins[S];
      PP->a0[D_emiss][S] = tables->A0[0][MS][D_emiss]
                            + T_ins * (tables->A0[1][MS][D_emiss]
                            + T_ins * tables->A0[2][MS][D_emiss]);

      PP->a2[D_emiss][S] = tables->A2[0][MS][D_emiss]
                            + T_ins * (tables->A2[1][MS][D_emiss]
                            + T_ins * tables->A2[2][MS][D_emiss]);

      PP->sigma_a0[D_emiss][S] = tables->sigma_a0[0][MS][D_emiss]
                            + T_ins * (tables->sigma_a0[1][MS][D_emiss]
                            + T_ins * tables->sigma_a0[2][MS][D_emiss]);

      PP->sigma_a2[D_emiss][S] = tables->sigma_a2[0][MS][D_emiss]
                            + T_ins * (tables->sigma_a2[1][MS][D_emiss]
                            + T_ins * tables->sigma_a2[2][MS][D_emiss]);
    }
    else {
      PP->a0[D_emiss][S] = 0.;
      PP->a2[D_emiss][S] = 0.;
      PP->sigma_a0[D_emiss][S] = 0.;
      PP->sigma_a2[D_emiss][S] = 0.;
    }
  }

  return (returnStatus);
}

PGSt_SMF_status Get_DN_Avg_SDev_Rejects
                       (int32   start_index,
                        int32   N,
                        int32   index_increment,
                        int16   *DN_array,
                        int16   DN_upper_valid_limit,
                        int16   DN_lower_valid_limit,
                        float32 *mean_DN,
                        float32 *sdev_DN,
                        int16   *rejects)
/*
!C**********************************************************************
!Description:
   Given a set of Reflective Solar Band (RSB) or Thermal Emissive Band
 (TEB) digital numbers (DNs), compute the mean and standard deviation
 of a subset of the set.  The subset is determined by supplying a start
 index, number of values and an index increment.  The DN values are checked
 for validity and an outlier-rejection algorithm is included that may also
 reduce the final set of values for which mean and standard deviation are
 computed.  The indices of the values rejected are also returned.


!Input Parameters:
    int32             start_index      starting index for values to include
    int32             N                number of values to include in the set.
    int32             index_increment  amount to increment when forming set
                                         of values
    int16 *           DN_array         array of digital numbers.  This should
                                       be large enough to extract all data
                                       from the subset.
    int16             DN_upper_valid_limit   Upper limit of the valid range of
                                             DNs (typically 4094).
    int16             DN_lower_valid_limit   Lower limit of the valid range of
                                             DNs (typically 0).
!Output Parameters:
    float32 *         mean_DN          Mean value of final set of data after
                                       all  invalid and outlier values are
                                       removed.  A  negative return value means
                                       that all values  were rejected for some
                                       reason.
    float32 *         sdev_DN          standard deviation of the final set of
                                       values.   A negative return value means
                                       that all values  were rejected for some
                                       reason.
    int16 *           rejects          An array of dimension [N], where 1
                                       denotes a  value was removed and 0
                                       denotes that a value  was used in the
                                       final set.

!Revision History:

 Revision 01.10 March 25, 2002   Razor Issue #178
 Strip out ADC Correction (delta_DN)
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.02 November 7, 2001
 Initialized DN_work to 0 (Razor issue #168)
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.01 August 1999
 Changed interface due to new ADC correction algorithm. See design document for
   detail.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 April 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32   i;
  int16   DN;
  double  sum1;
  double  sum2;
  float32 mean;
  float32 sigma;
  float32 three_sigma;
  int32   ninclude;
  float32 delta;
  float32 DN_work[MAX_250M_OBC_FRAME_DIM];

  if (N <= 0 || N > MAX_250M_OBC_FRAME_DIM)
    return (MODIS_F_NOK);

  for (i = 0; i < MAX_250M_OBC_FRAME_DIM; i++)
    DN_work[i] = 0.;

  *mean_DN = -1;
  *sdev_DN = -1;
  sum1 = 0;
  sum2 = 0;
  ninclude = 0;
  for (i = 0; i < N; i++)
  {
    DN = DN_array[start_index + i * index_increment];
    if (DN >= DN_lower_valid_limit && DN <= DN_upper_valid_limit)
    {
      DN_work[i] = DN;

      sum1 += (double) DN_work[i];
      sum2 += (double) DN_work[i] * (double) DN_work[i];
      ninclude++;
      rejects[i] = 0;
    }
    else
    {
      rejects[i] = 1;
    }
  }

  if (ninclude == 0)
    return returnStatus;
  else if (ninclude == 1) {
    *mean_DN = sum1;
    *sdev_DN = 0;
    return returnStatus;
  }

  mean  = sum1/ninclude;
  sigma = sqrt(fabs((double) (sum2 - sum1 * sum1/(double)ninclude)) /
            (double) (ninclude - 1));
  three_sigma = 3 * sigma;

  for (i = 0; i < N; i++)
  {
    delta = fabs((double) (DN_work[i] - mean));
    if (rejects[i] == 0 && delta > three_sigma)
    {
      sum1 -= (double) DN_work[i];
      sum2 -= (double) DN_work[i] * (double) DN_work[i];
      rejects[i] = 1;
      ninclude--;
    }
  }

  if (ninclude == 0 || ninclude == 1)
    return returnStatus;

  mean  = sum1/ninclude;
  sigma = sqrt(fabs((double) (sum2 - sum1 * sum1/(double)ninclude)) /
            (double) (ninclude - 1));
  *mean_DN = mean;
  *sdev_DN = sigma;

  return returnStatus;
}

PGSt_SMF_status Pack_Rejects_In_Outlier_Mask(int32  N,
                                             int16  *rejects,
                                             uint32 *packed_rejects)
/*
!C**********************************************************************
!Description:
   Store the values of flag array "rejects" into a smaller array
 "packed_rejects". Each value is stored in a bit beginning at the least
 significant bit of the first element in "packed_rejects".

!Input Parameters:
   int32      N           number of elements in the array rejects
   int16      rejects     array of length N,  holding 1's and 0's
!Output Parameters:
   uint32 *   packed_rejects    array to pack the rejects in.  It is
                                assumed that there are enough 32 bit
                                words in this variable to pack all N
                                values in the rejects array into.

!Revision History:
 Revision 01.00 April 1999
 Initial development
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 mask;
  int32 nvars;
  int32 nrem;
  int32 ivar;
  int32 i;

  nvars = N / 32;
  nrem  = N % 32;

  for (ivar = 0; ivar < nvars; ivar++)
  {
    mask = 1;
    for (i = ivar * 32; i < (ivar + 1) * 32; i++)
    {
      if (rejects[i] == 1)
        packed_rejects[ivar] = packed_rejects[ivar] | mask;
      mask *= 2;
    }
  }

  mask = 1;
  for (i = nvars * 32; i < nvars * 32 + nrem; i++)
  {
    if (rejects[i] == 1)
      packed_rejects[ivar] = packed_rejects[ivar] | mask;
    mask *= 2;
  }

  return returnStatus;
}

PGSt_SMF_status Process_OBCEng_Refl
                        (int32                sd_id,
                         int32                num_scans,
                         float32              *T_ins,
                         int16                *MirrorSide,
                         refl_tables_t        *refl_tables,
                         Moon_arrays_t        *moon_arrays,
                         L1A_granule_OBCEng_t *L1A_OBCEng,
                         DN_OBC_Avg_t         *DN_OBC_Avg,
                         QA_Refl_t            *QA_refl,
                         Preprocess_Refl_t    *PP_Refl)
/*
!C**********************************************************************
!Description:
   Compute the DN OBC averages for all reflective solar bands.
 Store in DN_OBC_Avg for writing to the OBC file and store in PP_Refl for
 use in reflective calibration.

!Input Parameters:
   int32            sd_id         SD file access ID for the L1A middle granule.
   int32            num_scans     Number of scans in the middle granule.
   float32          *T_ins        instrument temperature array
   int16            *MirrorSide   mirror side flags
   refl_tables_t    *refl_tables  Structure which contains all reflective LUTs
   Moon_arrays_t    *moon_arrays   Structure containing the moon_in_SV_KOB array.
   L1A_granule_OBCEng_t *L1A_OBCEng  Structure containing nir, vis, mwir
                                     adc board index.
!Output Parameters:
   DN_OBC_Avg_t   DN_OBC_Avg   Structure containing arrays that will
                                     ultimately be written to the OBC file.
   QA_Refl_t            QA_refl      Depending on processing values, some
                                     elements of these may be tripped.
   Preprocess_Refl_t *PP_Refl        Structure containing DN saturation value.

!Revision History:

 Revision 02.01  March 25, 2002  Razor Issue #178
 Stripped out ADC correction.
 Alice Isaacman  SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.00 Sep 5, 2000
 Made changes for Razor issue 130.  Removed Ecal_On from argument list
 and modified the passing of Ecal_On in argument lists of subordinate
 functions.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 ... (many revisions not logged) ...

 Revision 01.03 Oct. 22, 1999
 Added checking for Ecal on
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.02 Sept. 2, 1999
 Implemented changes to meet the requirements of the new SWIR algorithm. See SDF.
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.01 August 1999
 Implemented changes to meet the requirements of the new ADC algorithm. See SDF.
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 April 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 S;
  float32 avg;

  for (S = 0; S < num_scans; S++)
  {
    QA_refl->all_SV_DN_bad[S] = 0;
    QA_refl->all_BB_DN_bad[S] = 0;
  }

  if ((RFLAG & 1) == 0) {
  returnStatus = Fill_250m_DN_OBC_Avg(sd_id,
                                      num_scans,
                                      MirrorSide,
                                      refl_tables,
                                      moon_arrays,
                                      L1A_OBCEng->Ecal_On,
                                      DN_OBC_Avg,
                                      QA_refl);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Fill_250m_DN_OBC_Avg() in Preprocess_OBCEng_Refl()");
  }

  if ((RFLAG & 2) == 0) {
  returnStatus = Fill_500m_DN_OBC_Avg(sd_id,
                                      num_scans,
                                      MirrorSide,
                                      refl_tables,
                                      moon_arrays,
                                      L1A_OBCEng->Ecal_On,
                                      DN_OBC_Avg,
                                      QA_refl);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Fill_500m_DN_OBC_Avg() in Preprocess_OBCEng_Refl()");
  }

  if ((RFLAG & 4) == 0) {
  returnStatus = Fill_1km_day_DN_OBC_Avg(sd_id,
                                         num_scans,
                                         MirrorSide,
                                         refl_tables,
                                         moon_arrays,
                                         L1A_OBCEng->Ecal_On,
                                         DN_OBC_Avg,
                                         QA_refl);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Fill_1km_day_DN_OBC_Avg() in Preprocess_OBCEng_Refl()");
  }

  returnStatus = Fill_Band_26_DN_OBC_Avg(sd_id,
                                         num_scans,
                                         MirrorSide,
                                         refl_tables,
                                         moon_arrays,
                                         L1A_OBCEng->Ecal_On,
                                         DN_OBC_Avg,
                                         QA_refl);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,
        "Fill_Band_26_DN_OBC_Avg() in Preprocess_OBCEng_Refl()");

  returnStatus = Get_Temp_Avg_And_Variance(num_scans,
                                           T_ins,
                                           &avg,
                                           &PP_Refl->var_T_ins);

  return (returnStatus);
}


PGSt_SMF_status Cross_Granule_Sliding_Average(int32   num_overlapScans,
                                              int32   num_scans_middle,
                                              float32 *Xarray,
                                              int16   *Xmir,
                                              float32 *avg)
/*
!C****************************************************************************
!Description:
 For a cross-granule array in the emissive preprocessing code, this function
 generates a sliding, windowed average of the input quantity.  The output
 of the function is an averaged value at each scan of the middle granule.  The
 window interval for averaging is twice the number of overlap scans.  Only
 scans of the same mirror side contribute to the average at any scan of the
 middle granule.  If the number of overlap scans is set to zero, then no
 averaging is accomplished. The values are simply assigned on a per-scan basis.

!Input Parameters:

 int32       num_overlapScans  Number of overlap scans the quantity being
                               averaged.
 int32       num_scans_middle  Number of scans in the middle granule.
 float32 *   Xarray            Cross-granule array holding data to be averaged.
                               This should be dimensioned and organized as:
                               [num_overlapScans + num_scans_middle
                                + num_overlapScans].
 int16 *     Xmir              Array of mirror side flags having the same
                               dimension and organization as Xarray.

!Output Parameters:

 float32 *   avg               Array of averaged values, dimensioned
                               [num_scans_middle].

!Revision History:
 Revision 01.01 May, 1999
 Formerly "Average_Over"Scans"
 See SDF, "Stage 3: Changes to Emissive Algorithm and Code" notes.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 Jan. 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
 The quantity being averaged is assumed to be a non-negative number.  Negative
 numbers indicate that the value was missing or bad for some reason.
 Similarly, the mirror side array should contain values of 0 or 1.  Any other
 value of the mirror side implies that the data for that scan is bad.

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16 XS;
  int16 S;
  int16 n;
  int16 mirror_side;
  float32 sum;

  if (num_scans_middle == 0)
    SMF_ERROR(MODIS_F_NOK,
        "invalid middle granule scan numbers in "
        "Cross_Granule_Sliding_Average()");

  if (num_overlapScans == 0)
  {
    for (S = 0; S < num_scans_middle; S++)
      avg[S] = Xarray[S];
    return returnStatus;
  }

  for (S = 0; S < num_scans_middle; S++)
  {
    sum = 0;
    n   = 0;
    mirror_side = Xmir[S + num_overlapScans];
    if (mirror_side == 0 || mirror_side == 1)
    {
      for (XS = S; XS < S + 2 * num_overlapScans + 1; XS++)
      {
        if (Xmir[XS] == mirror_side && Xarray[XS] >= 0)
        {
          sum = sum + Xarray[XS];
          n++;
        }
      }
    }

    if (n > 0)
      avg[S] = sum/n;
    else
      avg[S] = -1;
  }

  return returnStatus;
}

PGSt_SMF_status Granule_Average_Temperature(int32   num_values,
                                            float32 *array,
                                            float32 *avg)
/*
!C****************************************************************************
!Description:
 For an array of positive temperature values, this function generates and
 returns the average.  A value of 0 is assigned if all values of the array
 are bad (not positive).

!Input Parameters:
 int32      num_values   Number values of the quantity being averaged.
 float32 *  array        Array of temperatures dimension [num_values].

!Output Parameters:
 float32 *  avg          Average value is assigned to this address.

!Revision History:
 Revision 01.02 May, 1999, Initial Development
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
 The quantity being averaged is assumed to be a non-negative number.  Negative
 numbers indicate that the value was missing or bad for some reason.

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 n;
  float32 sum;
  int32 i;

  if (num_values == 0 || array == NULL || avg == NULL)
    SMF_ERROR(MODIS_F_NOK,
        "invalid argument in Granule_Average_Temperature()");

  sum = 0;
  n   = 0;
  for (i = 0; i < num_values; i++)
  {
    if (array[i] > 0)
    {
      sum = sum + array[i];
      n++;
    }
  }

  if (n > 0)
    *avg = sum/n;
  else
    *avg = 0.0;

  return returnStatus;
}

PGSt_SMF_status Get_Temp_Avg_And_Variance(int32   N,
                                          float32 *T,
                                          float32 *avg,
                                          float32 *var)
/*
!C****************************************************************************
!Description:
   Given a set of engineering temperatures, calculate the mean and variance
   over all input values.  Bad temperatures (flagged as being less than or
   equal to zero) are excluded.  If all values become excluded, the mean
   and variance are set to -1, implying bad values.  If all but one value
   is excluded, then the mean is set to the remaining value and the variance
   is set to -1.

!Input Parameters:
   int32             N       number of value in array T
   float32 *         T       array of temperatures

!Output Parameters:
    float32 *        avg     mean of the values
    float32 *        var     variance of the values

!Revision History:
 Revision 01.00 April 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32   i, n;
  double sum1;
  double sum2;

  if (N <= 0 || !T || !var || !avg)
    return (MODIS_F_NOK);

  *avg = -1;
  *var = -1;
  sum1 = 0;
  sum2 = 0;
  n = 0;
  for (i = 0; i < N; i++)
  {
    if (T[i] > 0)
    {
      sum1 += (double) T[i];
      sum2 += (double) T[i] * (double) T[i];
      n++;
    }
  }

  if (n == 0)
    return returnStatus;

  else if (n == 1)
  {
    *avg = (float32) sum1;
    return returnStatus;
  }

  *avg  = (float32) (sum1 / ((double) n));
  *var = (float32) ( (sum2 - sum1 * sum1 / ((double) n)) /
                      ((double) (n - 1)) );
  return returnStatus;
}


PGSt_SMF_status Fill_Invalid_Temp_DNs(int32  num_scans,
                                      uint16 miss_DN,
                                      uint16 sat_DN,
                                      uint16 *DN)
/*
!C****************************************************************************
!Description:
 For a set of engineering temperatures, check for missing or saturated values
 and if found, then assign the nearest valid value.  The only way that any
 elements in the array DN do not leave this routine assigned to a valid value
 is if none of the values in the array were valid.

!Input Parameters:
 int32     num_scans   Number of DN values in the array DN
 uint16    miss_DN     Value representing the missing DN value
 uint16    sat_DN      Value representing the saturated DN value
 uint16    *DN         Array of DN values

!Output Parameters:
 uint16    *DN         DN values with invalid values filled in.


!Revision History:
 Revision 01.00 April 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Science Data Support
 Team for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 last_good_index, next_good_index, i, j, found;

  if (num_scans <= 0 || !DN)
  {
    SMF_ERROR(MODIS_F_INVALID_ARGUMENT,
              "Fill_Invalid_Temp_DNs, Preprocess.c");
  }

  last_good_index = num_scans;
  next_good_index = -1;
  for (i = 0; i < num_scans; )  /* do not automatically increment i */
  {
    if (DN[i] != miss_DN && DN[i] != sat_DN)
    {
      last_good_index = i;
      next_good_index = i;
      i++;
    }
    else if (next_good_index > i && last_good_index < i)
    {

           /* assign closest one */

      if ((next_good_index - i) < (i - last_good_index))
        DN[i] = DN[next_good_index];
      else
        DN[i] = DN[last_good_index];
      i++;
    }
    else if (next_good_index > i)
    {

                   /* assign next_good_index value */

      DN[i] = DN[next_good_index];
      i++;
    }
    else           /* find the next valid one */
    {
      for (found = 0, j = i+1; j < num_scans; j++)
      {
        if (DN[j] != miss_DN && DN[j] != sat_DN)
        {
          next_good_index = j;
          found = 1;
          break;
        }
      }
      if (found == 0)   /* assign the rest, then return */
      {
        if (last_good_index < i)
        {
          for (j = i; j < num_scans; j++)
            DN[j] = DN[last_good_index];
        }
        return returnStatus;
      }
    }
  }
  return returnStatus;
}

PGSt_SMF_status Get_Electronics_index(int32  num_scans,
                                      uint16 *Electronics_pri,
                                      uint16 *Electronics_red,
                                      uint16 *Electronics_index)

/*
!C**************************************************************************
!Description:  This function determines whether electronics side "A" or "B"
               of a system was used on each scan over a set of scans. See
               design notes.
!Input Parameters:
     int32      num_scans       Number of scans defining the dimensions of
                                the input and output arrays(could be all the
                                scans for a granule or a portion of the scans
                                of the granule)
     uint16 *   Electronics_pri the array hold the telemetry values for the
                                primary electronics
     uint16 *   Electronics_red the array hold the telemetry values for the
                                redundant electronics

!Output Parameters:
     uint16 *   Electronics_index  array with each entry indicating an
                                   electronics is primary (side "A") or
                                   redundant (side "B") in each scan

!Revision History:
 Revision 01.00 Aug 1999
 Initial development
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16 S               = 0;
  int16 S2              = 0;
  int16 valid_set_found = 0;

  for (S = 0; S < num_scans; S++)
  {
    if ( Electronics_pri[S] == ON && Electronics_red[S] == OFF )
      Electronics_index[S] = ELECTRONICS_PRIMARY;

    else if ( Electronics_pri[S] == OFF && Electronics_red[S] == ON )
      Electronics_index[S] = ELECTRONICS_REDUNDANT;

    else if (Electronics_pri [S] == OFF && Electronics_red[S] == OFF )
    {
      /* missing eng data --> use another scan's Primary/Redundant flag values.

         Search for a valid Electronics Primary/Redundant flag set to use.
         It is highly unlikely the electronics side will have changed in
         any given granule, so this is not a bad work around.
      */

      valid_set_found = 0;

      for (S2 = 0; S2 < num_scans; S2++)
      {
        if ( Electronics_pri[S2] == ON && Electronics_red[S2] == OFF )
        {
          Electronics_index[S] = ELECTRONICS_PRIMARY;
          valid_set_found = 1;
          break; /* usable flag values found, break out of loop */
        }
        else if ( Electronics_pri[S2] == OFF &&
                    Electronics_red[S2] == ON )
        {
          Electronics_index[S] = ELECTRONICS_REDUNDANT;
          valid_set_found = 1;
          break; /* usable flag values found, break out of loop */
        }

      } /* next S2 */

      if ( !valid_set_found )
        /* no usable ELECTRONICS Primary/Redundant flag data in the big
           chunk of whole granule looked at! */
         return (MODIS_F_NOK);
    }
    else  /* unexpected (both Primary and Redundant are 1) --
             should never occur  */
      return (MODIS_F_NOK);


  }/*Loop over scan*/

  return (returnStatus);

}

PGSt_SMF_status Get_Electronics_index_special(int32  num_scans,
                                      uint16 *Electronics_pri,
                                      uint16 *Electronics_red,
                                      uint16 *Electronics_index)

/*
!C**************************************************************************
!Description:  This function is similar to the function Get_Electronics_index,
               except that it also handle the case of both sides A and B on 
               for PCLW, which is permissible according to the design notes.
               
!Input Parameters:
     int32      num_scans       Number of scans defining the dimensions of
                                the input and output arrays(could be all the
                                scans for a granule or a portion of the scans
                                of the granule)
     uint16 *   Electronics_pri the array hold the telemetry values for the
                                primary electronics
     uint16 *   Electronics_red the array hold the telemetry values for the
                                redundant electronics

!Output Parameters:
     uint16 *   Electronics_index  array with each entry indicating an
                                   electronics is primary (side "A") or
                                   redundant (side "B") in each scan

!Revision History:

!Team-unique Header:
    This software is developed by the MCST L1B Team for the 
    National Aeronautics and Space Administration,
    Goddard Space Flight Center.

!References and Credits:
    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16 S               = 0;
  int16 S2              = 0;
  int16 valid_set_found = 0;

  for (S = 0; S < num_scans; S++)
  {
    if ( Electronics_pri[S] == ON && Electronics_red[S] == OFF )
      Electronics_index[S] = ELECTRONICS_PRIMARY;

    else if ( Electronics_pri[S] == OFF && Electronics_red[S] == ON )
      Electronics_index[S] = ELECTRONICS_REDUNDANT;

    else if ( Electronics_pri[S] == ON && Electronics_red[S] == ON ) 
      Electronics_index[S] = ELECTRONICS_BOTH;

    else
    {
      /* missing eng data --> use another scan's Primary/Redundant flag values.

         Search for a valid Electronics Primary/Redundant flag set to use.
         It is highly unlikely the electronics side will have changed in
         any given granule, so this is not a bad work around.
      */

      valid_set_found = 0;

      for (S2 = 0; S2 < num_scans; S2++)
      {
        if ( Electronics_pri[S2] == ON && Electronics_red[S2] == OFF )
        {
          Electronics_index[S] = ELECTRONICS_PRIMARY;
          valid_set_found = 1;
          break; /* usable flag values found, break out of loop */
        }
        else if ( Electronics_pri[S2] == OFF &&
                    Electronics_red[S2] == ON )
        {
          Electronics_index[S] = ELECTRONICS_REDUNDANT;
          valid_set_found = 1;
          break; /* usable flag values found, break out of loop */
        }

      } /* next S2 */

      if ( !valid_set_found )
        /* no usable ELECTRONICS Primary/Redundant flag data in the big
           chunk of whole granule looked at! */
         return (MODIS_F_NOK);

    }

  }/*Loop over scan*/

  return (returnStatus);

}

PGSt_SMF_status Read_L1A_OBCEng (emiss_tables_t        *emiss_tables,
                                 L1A_granule_t         *L1A_Gran,
                                 L1A_granule_OBCEng_t  *L1A_OBCEng)
/*
!C****************************************************************************
!Description: This function reads BB, SV and temperature DNs and converts the
              temperature DNs into engineering units.

!Input Parameters:
 L1A_granule_t *        L1A_Gran     Contains SD file interface ID, Vdata
                                     file interface ID, number of scans and
                                     mirror side for every scan.

!Output Parameters:
 L1A_granule_OBCEng_t * L1A_OBCEng   Contains BB and SV DNs, number of scans
                                     in the granule, mirror side flags and
                                     engineering temperatures.

!Revision History:
 Revision 01.05 January 4, 2002   Razor Issue #154
 Added satellite_ID to parameters fed to Read_Convert_Temperatures.
  Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.04 November 22, 1999
 Added calls to L1BErrorMsg.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 August, 1999
 Added checking if the data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.02 May, 1999
 See SDF, "Stage 3: Changes to Emissive Algorithm and Code" notes.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 Feb, 1999
 Added read of Number of Day mode scans, EV start time and Scan Type here
 to consolidate read of all members of L1A_Gran in one place and avoid
 reading the same data more than once.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 Jan. 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32           S;
  int32           track_dim;
  int32           count;
  int32           satellite_ID = INVALID_SATELLITE_ID;
                      /* 0 = TERRA (PFM);   1 = AQUA (FM1) */
  char            *location    = "Read_L1A_OBCEng";

  /* Determine the satellite platform. */
  satellite_ID = L1A_Gran->satellite_id;
    /*
     * Assign Engineering coefficients depending on whether instrument
     * is Terra (PFM) or Aqua (FM1).
     * Check the platform name.  Currently, the allowed name for Terra is "AM-1"
     * and possibly "Terra", and the allowed name for Aqua could be "PM-1" or
     * "Aqua".
     */

  if (satellite_ID != TERRA && satellite_ID != AQUA)
      { returnStatus = MODIS_F_OUT_OF_RANGE;
        L1BErrorMsg(location, returnStatus,
                    "*** Satellite ID is incorrect ***\n"
                    "Satellite platform must be Terra or Aqua.",
                    NULL, QA_TABLES_FILE,
                    NULL, True);
        return returnStatus;
      }


  /*
   * Assign number of scans and mirror side to the OBCEng structure for use
   * within the preprocess module.
   */

  L1A_OBCEng->num_scans = L1A_Gran->num_scans;
  for (S = 0; S < L1A_OBCEng->num_scans; S++)
    L1A_OBCEng->MirrorSide[S] = L1A_Gran->MirrorSide[S];

  /*
   * Compute the track dimension of sector DN SDSs based on actual
   * number of scans.
   */

  track_dim =
      L1A_OBCEng->num_scans * DETECT_PER_BAND_AT_RES[INDEX_1000M_NIGHT];

  /*
   * Read black-body sector DNs
   */

  returnStatus = read_sds_rank3(L1A_Gran->sd_id,
                                "BB_1km_night",
                                track_dim,
                                NUM_1000M_NIGHT_BANDS,
                                BB_1km_FRAMES,
                                (void *)L1A_OBCEng->BB_1km_night);
  if ( returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus,
                "Could not read BB_1km_night.",
                "read_sds_rank3", FIRST_L1A_GRANULE,
                Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  /*
   * Check valid range on black-body sector DNs
   */

  count = track_dim * NUM_1000M_NIGHT_BANDS * BB_1km_FRAMES;
  returnStatus = Check_Valid_Range("BB_1km_night",
                                   DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   count,
                                   L1A_OBCEng->BB_1km_night);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*
   * Read space-view sector DNs
   */

  returnStatus = read_sds_rank3(L1A_Gran->sd_id,
                                "SV_1km_night",
                                track_dim,
                                NUM_1000M_NIGHT_BANDS,
                                SV_1km_FRAMES,
                                (void *)L1A_OBCEng->SV_1km_night);
  if ( returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus,
                "Could not read SV_1km_night.",
                "read_sds_rank3",
                FIRST_L1A_GRANULE,
                Invalid_MOD01_Msg,
                True);
    return returnStatus;
  }

  /*
   * Check valid range on space-view sector DNs
   */

  count = track_dim * NUM_1000M_NIGHT_BANDS * SV_1km_FRAMES;
  returnStatus = Check_Valid_Range("SV_1km_night",
                                   DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   count,
                                   L1A_OBCEng->SV_1km_night);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*
   * Determine the scan-by-scan flags Ecal_On for this granule.
   */

  returnStatus = Check_For_Ecal_On(FIRST_L1A_GRANULE,
                                   L1A_Gran->num_scans,
                                   L1A_Gran->v_id,
                                   L1A_OBCEng->Ecal_On);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Check_For_Ecal_On", FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*
   * Determine the scan-by-scan flags Sector_Rotated for this granule.
   */

  returnStatus = Check_If_Sector_Rotated(FIRST_L1A_GRANULE,
                                         L1A_Gran->num_scans,
                                         L1A_Gran->v_id,
                                         L1A_OBCEng->Sector_Rotated); 
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Check_If_Sector_Rotated", FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }


  /*
   * Read temperature DNs and convert them into engineering units
   */

  returnStatus = Read_Convert_Temperatures (emiss_tables,
                                            L1A_Gran->v_id,
                                            0,
                                            L1A_Gran->num_scans,
                                            satellite_ID,
                                            &L1A_OBCEng->temps);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Read_Convert_Temperatures", FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  return returnStatus;
}

PGSt_SMF_status Read_Overlap_OBCEng (L1A_granule_t     *L1A_Gran,
                                     emiss_tables_t    *emiss_tables,
                                     PGSt_PC_Logical   lun,
                                     Overlap_OBCEng_t  *Overlap_OBCEng,
                                     boolean           *scan_gap)

/*
!C****************************************************************************
!Description:
  This function reads data needed for emissive calibration preprocessing
  from a subset of the total number of scans from the leading or trailing
  L1A granule.  These data include the black-body digital numbers (DNs),
  space-view DNs, mirror side flags and engineering temperatures.  For the
  leading granule, we read a set of scans at the end of the granule. For the
  trailing granule, we read a set of scans at the start of the granule.

!Input Parameters:
  L1A_granule_t     *L1A_Gran       Contains the satellite id
  emiss_tables_t    *emiss_tables   Contains the number of overlap scans for
                                    processing of b1 and T BB and for other
                                    temperatures (2 different values).
  int32             satellite_ID    0 = TERRA, 1 = AQUA, -1 = INVALID_SATELLITE_ID
  PGSt_PC_Logical   lun             Logical unit number, must be for leading
                                    or trailing L1A granules only.

!Output Parameters:
  Overlap_OBCEng_t  *Overlap_OBCEng Contains BB and SV DNs, number of scans
                                    in the granule, mirror side flags and
                                    engineering temperatures.
  boolean           *scan_gap       Indicates gap of a few scans
                                         between granules.

!Revision History:

  Revision 01.13 March 2003, Razor Issue #173
  Initialized variables "time_to_first_finish", "first_valid_second_time", 
  and "time_from_second_start" to 0.0 for ANSI-C compliance.
  Liqin Tan, SAIC GSO (ltan@saicmodis.com)

  Revision 01.12 April 16, 2002   Razor Issue #166
  Refined time check on leading/trailing granules to allow processing but
    print a warning message if a small gap of a few scans between granules
    exists.
  Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

  Revision 01.11 January 4, 2002 (Razor Issue #154)
  Passed satellite_ID to Read_Convert_Temperatures.
  Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

  Revision 01.1 November 5, 2001 (Razor Issue #166)
  Added a check on time of leading and trailing granules.
  Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

  Revision 01.00 November 3, 1999
  Initial Development -- combination of two functions: Read_LeadingOverlap
  and Read_TrailingOverlap
  Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
  This software is developed by the MODIS Characterization Support
  Team (MCST)for the National Aeronautics and Space Administration,
  Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
  The previous functions (Read_LeadingOverlap Read_TrailingOverlap) were
  substantially developed by Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and
  Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

  HDF portions developed at the National Center for Supercomputing
  Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
  The macro ALLOW_MISSING_GRANULES is used as a switch to allow the code to
  process if a granule is determined to be missing or if all data are missing
  from the granule.  If the production rules change to allow this, then
  define this macro to turn on the capability.

  The use of "return" statements after the L1BErrorMsg (with a "True" for the
  error-out argument) is to allow easier testing of all paths in unit testing.
  Generally, if the error-out flag is True, we do not need the extra return.
  This assumes that a non-OK return status is fatal, which is the case for
  L1B.

  Regardless of the existence of the macro ALLOW_MISSING_GRANULES, there are
  some instances where the granule will be treated as missing and processing
  will continue:
  a) The number of scans in the MOD01 granule is zero (empty granule)
  b) The number of scans in the MOD01 granule exceeds MAX_NUM_SCANS
  c) All telemetry data for key telemetries (in Read_Convert_Temperatures)
     is missing or invalid.
  d) The leading/trailing MOD01 granule does not immediately precede/follow
     the middle granule in time (respectively).

  The value of Overlap_OBCEng->num_scans = 0 upon return from this function
  means that the rest of the program will treat this granule as MISSING.

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  intn            hdf_return   = FAIL;
  PGSt_integer    Version      = 1;
  char            file_name[PGSd_PC_FILE_PATH_MAX];
  int32           file_id      = FAIL;
  int32           v_id         = FAIL;
  int32           num_scans    = 0;
  int16           T            = 0; /*Target index*/
  int16           R            = 0; /*Resolution index*/
  int16           S            = 0; /*Scan index*/
  int32           track_dim    = 0;
  int32           start0       = 0; /* values used in read_part_sds routines*/
  int32           start1       = 0;
  int32           start2       = 0;
  int32           num_overlap_scans_b1    = 0;
  int16           track;  /* track index */
  int16           band;   /* band index */
  int16           frame;  /* frame index */
  int16           j;      /* minor loop index */
  int32           satellite_ID = INVALID_SATELLITE_ID;
                    /* 0 = TERRA (PFM);   1 = AQUA (FM1) */

    /* start times (to determine data collection time */
  float64   EVStartTime_TAIsecond[MAX_NUM_SCANS];

  /*
   * Last valid time of leading/middle (if LEADING/TRAILING resp.)
   * granule
   */
  float64   last_valid_first_time;

  /*
   * time interval from last_valide_first_time to end of
   * leading/middle (if LEADING/TRAILING resp.) granule
   */
  float64   time_to_first_finish = 0.0;

  /*
   * first valid time of mid/trailing (if LEADING/TRAILING resp.)
   * granule
   */
  float64   first_valid_second_time = 0.0;

  /*
   * time interval from start of mid/trailing
   * (if LEADING/TRAILING resp.) granule
   */
  float64   time_from_second_start = 0.0;

  /*
   * Time difference between granules based on times of valid scans and
   * estimate of time difference between first and second granules based
   * on scan numbersof valid times, respectively
   */
  float64 time_difference, estimated_time_difference;

  char *location = "Read_Overlap_OBCEng";
  int32 count;

    /*    Leading/Trailing granule indicator, determined by LUN */
  boolean leading  = FALSE;
  boolean trailing = FALSE;

  char *leading_missing_granule_msg = "\
The leading MOD01 granule will be treated as missing.\n\
The resulting L1B output will be valid if no other errors occur.";

  char *trailing_missing_granule_msg = "\
The trailing MOD01 granule will be treated as missing.\n\
The resulting L1B output will be valid if no other errors occur.";

#ifdef ALLOW_MISSING_GRANULES

    char *leading_granule_time_msg = "\
The leading MOD01 granule has an undefined data collection time or a\n\
data collection time which does not synchronize with the middle\n\
granule data collection time and will be treated as missing.  \n\
The resulting L1B output will be valid if no other errors occur.";

  char *trailing_granule_time_msg = "\
The trailing MOD01 granule has an undefined data collection time or a\n\
data collection time which does not synchronize with the middle\n\
granule data collection time and will be treated as missing.  \n\
The resulting L1B output will be valid if no other errors occur.";

#endif  /* ALLOW_MISSING_GRANULES */

  char *leading_dropped_scan_msg = "\
At least one and no more than";
  
    char *leading_dropped_scan_msg_end = "\n\
scans were dropped between the leading and middle granules.  \n\
The resulting L1B output will be valid if no other errors occur.";
  
    char *trailing_dropped_scan_msg = "\
At least one and no more than";
  
    char *trailing_dropped_scan_msg_end = "\n\
scans were dropped between the middle and trailing granules.  \n\
The resulting L1B output will be valid if no other errors occur.";

  /*
   * Check input arguments.  The LUN should correspond to leading or trailing
   * L1A granule only.
   */

  if (!(lun == LEADING_L1A_GRANULE || lun == TRAILING_L1A_GRANULE)) {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus,
                "LUN is invalid.  Must be for previous or "
                "following MOD01 granule only",
                NULL, lun, NULL, True);
  }

  if (lun == LEADING_L1A_GRANULE) leading = TRUE;
  if (lun == TRAILING_L1A_GRANULE) trailing = TRUE;

    /*
   * Initialize members of the OBCEng structure to fill values.
   * The most important of these is num_scans.  num_scans = 0 implies that
   * the granule is missing.
   */

  Overlap_OBCEng->num_scans = 0;
  for (track = 0; track < MAX_OVERLAP_TRACK_DIM; track++)
    for (band = 0; band < NUM_1000M_NIGHT_BANDS; band++)
      for (frame = 0; frame < BB_1km_FRAMES; frame++)
        Overlap_OBCEng->BB_1km_night[track][band][frame] =
            MISSING_L1A_FLAG;
  for (track = 0; track < MAX_OVERLAP_TRACK_DIM; track++)
    for (band = 0; band < NUM_1000M_NIGHT_BANDS; band++)
      for (frame = 0; frame < SV_1km_FRAMES; frame++)
        Overlap_OBCEng->SV_1km_night[track][band][frame] =
            MISSING_L1A_FLAG;
  for (S = 0; S < MAX_NUM_SCANS; S++)
  {
    Overlap_OBCEng->MirrorSide[S]                         = -1;
    for (j = 0; j < NUM_FOCAL_PLANES; j++)
      Overlap_OBCEng->temps.fp[j][S]                      = 0;
    Overlap_OBCEng->temps.mir1[S]                         = 0;
    Overlap_OBCEng->temps.mir2[S]                         = 0;
    Overlap_OBCEng->temps.mir_avg[S]                      = 0;
    Overlap_OBCEng->temps.cav[S]                          = 0;
    Overlap_OBCEng->temps.ins[S]                          = 0;
    for (j = 0; j < NUM_BB_THERMISTORS; j++)
      Overlap_OBCEng->temps.bb[j][S]                      = 0;
    Overlap_OBCEng->temps.bb_avg[S]                       = 0;
    Overlap_OBCEng->temps.fp_set_point_state[S]           = 0;
    Overlap_OBCEng->temps.lwir_adc_index[S]               = 0;
    Overlap_OBCEng->temps.mwir_adc_index[S]               = 0;
    Overlap_OBCEng->temps.nir_adc_index[S]                = 0;
    Overlap_OBCEng->temps.vis_adc_index[S]                = 0;
    Overlap_OBCEng->temps.pclw_adc_index[S]               = 0;
    Overlap_OBCEng->temps.num_thermistor_outliers[S]      = 0;
  }

  /*
   * Prepare for opening file: convert logical ID to physical name
   */

  returnStatus = PGS_PC_GetReference (lun, &Version, file_name);
#ifdef ALLOW_MISSING_GRANULES
  if (returnStatus != PGS_S_SUCCESS)
  {                                           /* allow continued processing */
    returnStatus = MODIS_S_OK;
    if (leading) L1BErrorMsg(location, returnStatus,
          "WARNING: Could not retrieve leading granule file name from PCF.",
          "PGS_PC_GetReference", lun, leading_missing_granule_msg, False);
    if (trailing) L1BErrorMsg(location, returnStatus,
          "WARNING: Could not retrieve trailing granule file name from PCF.",
          "PGS_PC_GetReference", lun, trailing_missing_granule_msg, False);
    return returnStatus;
  }
#else
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    if (leading) L1BErrorMsg(location, returnStatus,
                "Could not retrieve leading granule file name from PCF.",
                "PGS_PC_GetReference", lun, NULL, True);
    if (trailing) L1BErrorMsg(location, returnStatus,
                "Could not retrieve trailing granule file name from PCF.",
                "PGS_PC_GetReference", lun, NULL, True);
    return returnStatus;
  }
#endif /* ALLOW_MISSING_GRANULES */

  /*
   * Open file for HDF science data (SD) access.
   */

  file_id = SDstart(file_name, DFACC_RDONLY);
#ifdef ALLOW_MISSING_GRANULES
  if (file_id == FAIL)
  {                                           /* allow continued processing */
    returnStatus = MODIS_S_OK;
    if (leading) L1BErrorMsg(location, returnStatus,
        "WARNING: Could not open leading granule file for SD read access.",
        "SDstart", lun, leading_missing_granule_msg, False);
    if (trailing) L1BErrorMsg(location, returnStatus,
        "WARNING: Could not open trailing granule file for SD read access.",
        "SDstart", lun, trailing_missing_granule_msg, False);
    return returnStatus;
  }
#else
  if (file_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    if (leading) L1BErrorMsg(location, returnStatus,
        "Could not open leading granule file for SD read access.",
        "SDstart", lun,
        "The file may be missing, corrupted or not an HDF-4 file.", True);
    if (trailing) L1BErrorMsg(location, returnStatus,
        "Could not open trailing granule file for SD read access.",
        "SDstart", lun,
        "The file may be missing, corrupted or not an HDF-4 file.", True);
    return returnStatus;
  }
#endif /* ALLOW_MISSING_GRANULES */

  /*
   * Determine the satellite instrument is from the leading or trailing granule
   * and check if it matches that determined from the middle granule.
   */

  returnStatus = Get_Satellite_ID(lun, &satellite_ID);
  if (returnStatus != MODIS_S_OK)
  {
    if (leading) L1BErrorMsg(location, returnStatus,
                "Failed to get satellite id from leading granule.",
                "Get_Satellite_ID", lun,
                "The leading granule is invalid.", True);
    if (trailing) L1BErrorMsg(location, returnStatus,
                "Failed to get satellite id from trailing granule.",
                "Get_Satellite_ID", lun,
                "The trailing granule is invalid.", True);
    return returnStatus;
  }

  if (satellite_ID != L1A_Gran->satellite_id)
  {
    returnStatus = MODIS_F_NOK;
    if (leading) L1BErrorMsg(location, returnStatus,
                "The satellite instrument name in the leading granule\n"
                "does not match that in the middle granule.",
                NULL, lun, NULL, True);
    if (trailing) L1BErrorMsg(location, returnStatus,
                "The satellite instrument name in the trailing granule\n"
                "does not match that in the middle granule.",
                NULL, lun, NULL, True);
    return returnStatus;
  }

  /*
   * Open file for HDF Vdata access and initialize Vdata interface.
   * Call Hopen() and Vstart() in the same place.  Here, do not allow
   * processing with a missing granule since we were able to open the
   * file using SDstart.  If we cannot open it with Hopen and Vstart, then
   * the file must be bad and we should terminate.
   */

  v_id = Hopen(file_name,DFACC_RDONLY,0);
  if (v_id == FAIL)
  {
    SDend (file_id);
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    if (leading) L1BErrorMsg(location, returnStatus,
                "Could not open leading granule file for Vdata read access.",
                "Hopen", lun,
                "The file may be corrupted or not an HDF-4 file.", True);
    if (trailing) L1BErrorMsg(location, returnStatus,
                "Could not open trailing granule file for Vdata read access.",
                "Hopen", lun,
                "The file may be corrupted or not an HDF-4 file.", True);
    return returnStatus;
  }
  hdf_return = Vstart (v_id);
  if (hdf_return == FAIL)
  {
    SDend (file_id); Hclose (v_id);
    returnStatus = MODIS_F_HDF_ERROR;
    if (leading) L1BErrorMsg(location, returnStatus,
        "Could not initialize Vdata interface for leading granule file.",
        "Vstart", lun,
        "The file may be corrupted or not an HDF-4 file.", True);
    if (trailing) L1BErrorMsg(location, returnStatus,
        "Could not initialize Vdata interface for trailing granule file.",
        "Vstart", lun,
        "The file may be corrupted or not an HDF-4 file.", True);
    return returnStatus;
  }

  /*
   * Read & check num_scans.  If the number of scans is zero, assume that the
   * L1A granule is valid, but empty.  Treat it as a missing granule.
   * Other errors (except time errors) will cause error out.
   */

  returnStatus = read_attribute (file_id, "Number of Scans",
                                 DFNT_INT32, (void *)&num_scans);
  if (returnStatus != MODIS_S_OK)
  {
    if (leading)
      L1BErrorMsg(location, returnStatus,
        "Could not read Number of Scans from leading granule file.",
        "read_attribute", lun, Invalid_MOD01_Msg, True);
    if (trailing)
      L1BErrorMsg(location, returnStatus,
        "Could not read Number of Scans from trailing granule file.",
        "read_attribute", lun, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  if (num_scans == 0)
  {                                      /* allow continued processing */
    returnStatus = MODIS_S_OK;
    if (leading)
      L1BErrorMsg(location, returnStatus,
               "WARNING: Number of scans in leading granule is zero.",
                NULL, lun, leading_missing_granule_msg, False);
    if (trailing)
      L1BErrorMsg(location, returnStatus,
               "WARNING: Number of scans in trailing granule is zero.",
                NULL, lun, trailing_missing_granule_msg, False);
    goto Read_Overlap_OBCEng_exit;      /* close accesses and files. */
  }

  if (num_scans < 0)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    if (leading) L1BErrorMsg(location, returnStatus,
                "Number of Scans in the leading granule is negative.",
                NULL, lun, Invalid_MOD01_Msg, True);
    if (trailing) L1BErrorMsg(location, returnStatus,
                "Number of Scans in the trailing granule is negative.",
                NULL, lun, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  if (num_scans > MAX_NUM_SCANS)
  {                                      /* allow continued processing */
    returnStatus = MODIS_S_OK;
    if (leading)
      L1BErrorMsg(location, returnStatus,
        "WARNING: Number of scans in leading granule is > MAX_NUM_SCANS.",
        NULL, lun, leading_missing_granule_msg, False);
    if (trailing)
      L1BErrorMsg(location, returnStatus,
        "WARNING: Number of scans in trailing granule is > MAX_NUM_SCANS.",
        NULL, lun, trailing_missing_granule_msg, False);
    goto Read_Overlap_OBCEng_exit;       /* close accesses and files. */
  }

  /*
   * Assign num_scans now because it is a valid number.
   */

  Overlap_OBCEng->num_scans = num_scans;

  /*
   * Check to see if number of overlap scans is zero.  If so, then there is
   * nothing to read from the file (no data will be used).  Go to function
   * exit point, close up files and continue processing.
   */

  num_overlap_scans_b1 = emiss_tables->num_overlap_scans_b1;
  if (num_overlap_scans_b1 == 0)
  {
    returnStatus = MODIS_S_OK;
    goto Read_Overlap_OBCEng_exit;
  }

  /*
   * Check number of scans against the emissive LUTs
   */

  if (num_scans < num_overlap_scans_b1)
  {
    returnStatus = MODIS_F_NOK;
    if (leading)
      L1BErrorMsg(location, returnStatus,
        "Number of Scans in the leading granule "
          "less than number of overlap scans b1.",
        NULL, lun,
        "There are not enough scans in the leading MOD01 granule to allow\n"
        "processing with the current L1B emissive algorithm.", True);
    if (trailing)
      L1BErrorMsg(location, returnStatus,
        "Number of Scans in the trailing granule "
          "less than number of overlap scans b1.",
        NULL, lun,
        "There are not enough scans in the trailing MOD01 granule to allow\n"
        "processing with the current L1B emissive algorithm.", True);
    return returnStatus;
  }

  /*
   * Check the LUT values against the dimension macros.  These checks are
   * independent of which granule it is.
   */

  if (num_overlap_scans_b1 > MAX_NUM_OVERLAP_SCANS)
  {
    returnStatus = MODIS_F_NOK;
    if (leading)
      L1BErrorMsg(location, returnStatus,
        "num_overlap_scans_b1 in the leading granule is greater than \n"
        "macro MAX_NUM_OVERLAP_SCANS.",
        NULL, 0,
        "The L1B code macro must be increased or the LUT value reduced.",
        True);
    if (trailing)
      L1BErrorMsg(location, returnStatus,
        "num_overlap_scans_b1 in the trailing granule is greater than \n"
        "macro MAX_NUM_OVERLAP_SCANS.",
        NULL, 0,
        "The L1B code macro must be increased or the LUT value reduced.",
        True);
    return returnStatus;
  }

  /* Read the scan times to check against the leading/trailing positions. */

  returnStatus = read_sds_rank1 (file_id, "EV start time",
                                 num_scans, (void *)EVStartTime_TAIsecond);
  if (returnStatus != MODIS_S_OK)
  {
    if (leading) L1BErrorMsg(location, returnStatus,
                "Could not read EV start time for the leading granule.",
                "read_sds_rank1", lun, Invalid_MOD01_Msg, True);
    if (trailing) L1BErrorMsg(location, returnStatus,
                "Could not read EV start time for the trailing granule.",
                "read_sds_rank1", lun, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  last_valid_first_time   = 0;
  first_valid_second_time = 0;

  if (leading)
  /*
   * Determine the times of the first valid scan of the middle granule and
   * the last valid scan of the leading granule.
   */
  {

    for (S = num_scans - 1; S >= 0; S--) /* Count backwards for last time */
    {
      if (EVStartTime_TAIsecond[S] > 0)
      {
        last_valid_first_time = EVStartTime_TAIsecond[S];
        time_to_first_finish =
            SCAN_TIME_INTERVAL*((float64) (num_scans - .5 - S));
        break;
      }
    }

    for (S = 0; S < L1A_Gran->num_scans; S++)
    {
      if (L1A_Gran->EVStartTime_TAIsecond[S] > 0)
      {
        first_valid_second_time = L1A_Gran->EVStartTime_TAIsecond[S];
        time_from_second_start =
            SCAN_TIME_INTERVAL*((float64) (.5 + S));
        break;
      }
    }
  }

  if (trailing)
  /*
   * Determine the times of the first valid scan of the trailing granule and
   * the last valid scan of the middle granule.
   */
  {
    for (S = L1A_Gran->num_scans - 1; S >= 0; S--)
    {
      if (L1A_Gran->EVStartTime_TAIsecond[S] > 0)
      {
        last_valid_first_time = L1A_Gran->EVStartTime_TAIsecond[S];
        time_to_first_finish =
            SCAN_TIME_INTERVAL*((float64) (L1A_Gran->num_scans - .5 - S));
        break;
      }
    }

    for (S = 0; S < num_scans; S++)
    {
      if (EVStartTime_TAIsecond[S] > 0)
      {
        first_valid_second_time = EVStartTime_TAIsecond[S];
        time_from_second_start =
            SCAN_TIME_INTERVAL*((float64) (.5 + S));
          break;
      }
    }
  }

  /*
   * In the following explanation, "first" and "second" granule
   * refer to the leading and middle granule respectively if the
   * LEADING granule is being read, and to the middle and trailing
   * granule respectively if the TRAILING granule is being read.
   *
   * The first valid scan of the second granule must have a time
   * after that of the last valid time of the first granule.  The
   * time difference between the last valid time of the first
   * granule and the first valid time of the second granule should
   * be approximately the same as the time difference estimated from
   * the scan numbers of the last valid scan of the first granule
   * and the first valid scan of the second granule.  If the time
   * difference is not acceptable (less than 0 or greater than a set
   * number of scans) then set the number of scans to 0
   * so that the data will be treated as missing, or stop with an
   * error if D_ALLOW_MISSING_GRANULES is false.  If the time difference
   * indicates at least one scan was dropped, continue with a warning message
   * and set a new bit qa flag.
   */

  time_difference = first_valid_second_time - last_valid_first_time;
  estimated_time_difference =
      time_to_first_finish + time_from_second_start;

  /*
   * If the time difference is more than 1 scan but less than the upper
   * allowed number of scans, print a warning message but allow the
   * granule to be used for calibration purposes.  Set a qa flag bit.
   */

  if ((time_difference >= 0.) &&
    (fabs(estimated_time_difference - time_difference) >=
       INNER_GRANULE_TIMEDIFF_TOLERANCE) &&
    (fabs(estimated_time_difference - time_difference) <=
       OUTER_GRANULE_TIMEDIFF_TOLERANCE))
    {

      *scan_gap = True;

      if (leading)
      {
        char dropped_scan_msg[256];
        sprintf (dropped_scan_msg, "%s %d %s", leading_dropped_scan_msg,
            DROPPED_SCAN_UPPER_LIMIT, leading_dropped_scan_msg_end);

        L1BErrorMsg(location, returnStatus,
          "WARNING: the leading granule does not immediately precede \n"
          "the middle granule.",
           NULL, lun, dropped_scan_msg, False);
      }
      if (trailing)
      {
        char dropped_scan_msg[256];
        sprintf (dropped_scan_msg, "%s %d %s", trailing_dropped_scan_msg,
            DROPPED_SCAN_UPPER_LIMIT, trailing_dropped_scan_msg_end);
        L1BErrorMsg(location, returnStatus,
          "WARNING: the trailing granule does not immediately follow \n"
          "the middle granule.",
          NULL, lun, dropped_scan_msg, False);
      }
    }

  /*
   *   Otherwise check to see whether there is a significant time
   *   discrepancy and do not use the leading/trailing granule if so.
   */

  else
  {
    if ((time_difference < 0.) ||
        (fabs(estimated_time_difference - time_difference) >=
           OUTER_GRANULE_TIMEDIFF_TOLERANCE))
#ifdef ALLOW_MISSING_GRANULES
    {                                    /* allow continued processing */

      Overlap_OBCEng->num_scans = 0;  /* flag this granule as missing */
      returnStatus = MODIS_S_OK;
      if (leading)
        L1BErrorMsg(location, returnStatus,
          "WARNING: the leading granule does not immediately precede \n"
          "the middle granule.",
           NULL, lun, leading_granule_time_msg, False);
      if (trailing)
        L1BErrorMsg(location, returnStatus,
          "WARNING: the trailing granule does not immediately follow \n"
          "the middle granule.",
          NULL, lun, trailing_granule_time_msg, False);

      goto Read_Overlap_OBCEng_exit;  /* close accesses and files. */

    }
#else
    {                                    /* do not allow continued processing */

      Overlap_OBCEng->num_scans = 0;     /* flag this granule as missing */
      returnStatus = MODIS_F_NOK;
      if (leading) L1BErrorMsg(location, returnStatus,
                 "The leading granule does not immediately precede \n"
                 "the middle granule.",
                  NULL, lun, NULL, True);
      if (trailing) L1BErrorMsg(location, returnStatus,
                  "The trailing granule does not immediately follow \n"
                  "the middle granule.",
                  NULL, lun, NULL, True);

      goto Read_Overlap_OBCEng_exit;     /* close accesses and files. */
    }
#endif /* ALLOW_MISSING_GRANULES */
  }

  /*
   * Read the sector DN data for BB and SV for the 1km Emissive bands.  We
   * only need to read these if the number of overlap scans for b1 is > 0.
   * The value of overlap scans for b1 could be zero and yet the value of
   * overlap scans for tempes is non-zero (we have already checked to see
   * if both are zero).
   */

  if (num_overlap_scans_b1 > 0)
  {
    R = INDEX_1000M_EMISS;
    if (leading)
      start0 = (num_scans - num_overlap_scans_b1) *
                   DETECT_PER_BAND_AT_RES[R];
    else
      start0 = 0;

    track_dim = num_overlap_scans_b1 * DETECT_PER_BAND_AT_RES[R];

    /*
     * Read BB data and check that values are not out of bounds.
     */

    T = BB_INDEX;
    returnStatus = read_part_sds_rank3 (file_id,
                                        L1A_TARGET_SDS_NAME[T][R],
                                        start0,
                                        start1,
                                        start2,
                                        track_dim,
                                        L1A_BANDS_AT_RES[R],
                                        BB_1km_FRAMES,
                                        Overlap_OBCEng->BB_1km_night);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
          "read_part_sds_rank3", lun, NULL, True);
      return returnStatus;
    }
    count = track_dim * NUM_1000M_NIGHT_BANDS * BB_1km_FRAMES;
    returnStatus = Check_Valid_Range(L1A_TARGET_SDS_NAME[T][R],
                                     DFNT_INT16,
                                     L1A_DN_SDS_LB,
                                     L1A_DN_SDS_UB,
                                     L1A_DN_SDS_FV,
                                     count,
                                     Overlap_OBCEng->BB_1km_night);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
          "Check_Valid_Range", lun, NULL, True);
      return returnStatus;
    }

    /*
     * Read SV data and check that values are not out of bounds.
     */

    T = SV_INDEX;
    returnStatus = read_part_sds_rank3 (file_id,
                                        L1A_TARGET_SDS_NAME[T][R],
                                        start0,
                                        start1,
                                        start2,
                                        track_dim,
                                        L1A_BANDS_AT_RES[R],
                                        SV_1km_FRAMES,
                                        Overlap_OBCEng->SV_1km_night);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
          "read_part_sds_rank3", lun, NULL, False);
      return returnStatus;
    }
    count = track_dim * NUM_1000M_NIGHT_BANDS * SV_1km_FRAMES;
    returnStatus = Check_Valid_Range(L1A_TARGET_SDS_NAME[T][R],
                                     DFNT_INT16,
                                     L1A_DN_SDS_LB,
                                     L1A_DN_SDS_UB,
                                     L1A_DN_SDS_FV,
                                     count,
                                     Overlap_OBCEng->SV_1km_night);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
          "Check_Valid_Range", lun, NULL, True);
      return returnStatus;
    }

  }                                        /* end of if (num_overlap_scans_b1 > 0) */

  /*
   * Read and check MirrorSide array.
   */

  returnStatus = read_sds_rank1(file_id,
                                "Mirror side",
                                Overlap_OBCEng->num_scans,
                                Overlap_OBCEng->MirrorSide);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
        "read_sds_rank1", lun, NULL, True);
    return returnStatus;
  }
  returnStatus = Check_Valid_Range("Mirror side",
                                   DFNT_INT16,
                                   "-1",
                                   "1",
                                   NULL,
                                   Overlap_OBCEng->num_scans,
                                   Overlap_OBCEng->MirrorSide);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                lun, NULL, True);
    return returnStatus;
  }

  /*
   * Determine the scan-by-scan flags Ecal_On for this granule.
   */

  returnStatus = Check_For_Ecal_On(lun,
                                   Overlap_OBCEng->num_scans,
                                   v_id,
                                   Overlap_OBCEng->Ecal_On);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Check_For_Ecal_On", lun, NULL, True);
    return returnStatus;
  }

  /*
   * Determine the scan-by-scan flags Sector_Rotated for this granule.
   */

  returnStatus = Check_If_Sector_Rotated(lun,
                                         Overlap_OBCEng->num_scans,
                                         v_id,
                                         Overlap_OBCEng->Sector_Rotated);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Check_If_Sector_Rotated", lun, NULL, True);
    return returnStatus;
  }


  /*
   * Read engineering telemetry data (temperatures, ADC indexes, etc.)
   * and convert those that are temperatures into physical units.
   */

  if (lun == LEADING_L1A_GRANULE)
    start0 = num_scans - num_overlap_scans_b1;
  else
    start0 = 0;
  returnStatus = Read_Convert_Temperatures (emiss_tables,
                                            v_id,
                                            start0,
                                            num_overlap_scans_b1,
                                            satellite_ID,
                                            &Overlap_OBCEng->temps);


  /*
   * Error code logic:  If the return status is MODIS_F_OUT_OF_RANGE,
   * terminate execution because there is something seriously wrong with the
   * granule.  Otherwise, assume that telemetry was missing and simply proceed
   * as if the granule were missing.
   */


  if (returnStatus == MODIS_F_OUT_OF_RANGE)
  {
    L1BErrorMsg(location, returnStatus, NULL,
        "Read_Convert_Temperatures", lun, NULL, True);
    return returnStatus;
  }
  else if (returnStatus != MODIS_S_OK)
  {                                           /* allow continued processing */
    Overlap_OBCEng->num_scans = 0;            /* flag this granule as missing */
    returnStatus = MODIS_S_OK;
    if (leading)
      L1BErrorMsg(location, returnStatus,
        "WARNING: the leading granule cannot be used "
        "(possible missing telemetry).",
        "Read_Convert_Temperatures", lun,
          leading_missing_granule_msg, False);
    if (trailing)
      L1BErrorMsg(location, returnStatus,
        "WARNING: the trailing granule cannot be used "
        "(possible missing telemetry).",
        "Read_Convert_Temperatures", lun,
          trailing_missing_granule_msg, False);
    goto Read_Overlap_OBCEng_exit;            /* close accesses and files. */
  }

  /*
   * Successfully completed exit code:
   */

  returnStatus = MODIS_S_OK;

  /*
   * Function exit point.
   */

Read_Overlap_OBCEng_exit:

  /*
   * Close file and all accesses.
   */

  if (v_id != FAIL)
  {
    hdf_return = Vend(v_id);
    if (hdf_return == FAIL)
    {
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, "Vend", lun,
         "Memory or the disk file must have become corrupted.", True);
    }
    hdf_return = Hclose(v_id);
    if (hdf_return == FAIL)
    {
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, "Hclose", lun,
         "Memory or the disk file must have become corrupted.", True);
    }
  }
  if (file_id != FAIL)
  {
    if (SDend(file_id) == FAIL)
    {
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, "SDend", lun,
         "Memory or the disk file must have become corrupted.", True);
    }
  }
  return returnStatus;
}

PGSt_SMF_status Read_Convert_Temperatures (emiss_tables_t *tables,
                                           int32          v_id,
                                           int32          start_scan,
                                           int32          num_scans,
                                           int32          satellite_ID,
                                           Temperatures_t *temps)
/*
!C**************************************************************************
!Description:   Read engineering temperature DNs and convert them into
                temperatures in physical units.

!Input Parameters:
     int32          v_id         file id for vdata interfaces
     int32          start_scan   starting scan to read
     int32          num_scans    number of scans to read
     int32          satellite_ID TERRA = 0, AQUA = 1, INVALID_SATELLITE_ID = -1

!Output Parameters:
     Temperatures_t *temps

!Revision History:
 Revision 02.18 March 2003, Razor Issue #173
 Enclosed the rows in the initializers of array-of-struct "temp" with braces, 
 and changed the type of variables "FP_SPS" and "thermistor" from "int8" to 
 "int16" for ANSI-C compliance.
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 02.17   January 4, 2002  (Razor Issue #154)
 Added input variable satellite_ID and passed it to Compute_BB_Temperature.
 Added structure Eng_Coeffs determined by satellite ID and used it for
 conversion coefficients.  Moved computation of CPA/CPB index to this module.

 Revision 02.16  October 24, 2000
 Added code for granule average QA values as per issue 140.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 ... (many changes not logged) ...

 Revision 02.15 November 22, 1999
 Added calls to L1BErrorMsg.  Changed logic of immediately error-out
 to instead simply write error message and return to parent (so that
 parent can identify LUN of file).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.14 August 28, 1999
 Added checking if the data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.13 August, 1999
 Implemented changes to meet the requirements of the new ADC algorithm. See SDF.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.12 May, 1999
 Implemented changes described in SDF.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.11 Oct 26 1998
 put the BB temperature part into function "Compute_BB_Temperature"
 and put the adc index part into function "Get_ADC_Index"
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 02.10 Apr 20 1998
 Change "TA_RC_SWIR_CFPAE" back to "TA_RC_SMIR_CFPAE",
 and TP_AO_SWIR_OBJ to TP_AO_SMIR_OBJ.
 That's the way it is in the L1A file.
 David Catozzi (cato@ltpmail)

 Revision 02.10 Mar 1998
 Change "TA_RC_SMIR_CFPAE" to "TA_RC_SWIR_CFPAE",
 Add new temperatures for the V2.1 TEB algorithm
 Implement the 1/ln algorithm for T_BB
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 01.00 Aug 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

  If an error occurs in this function, write the error message and then
  return to the parent (so that the parent can identify the LUN of
  the file).

!END********************************************************************
*/
{
  char *location = "Read_Convert_Temperatures";
  PGSt_SMF_status   returnStatus   = MODIS_S_OK;
  enum {
    FP_T1SET,
    FP_T3SET,
    LWHTR_ON,
    SMHTR_ON,
    FP1,
    FP2,
    FP3,
    FP4,
    MIR1,
    MIR2,
    SCN_MTR,
    CAV1,
    CAV2,
    CAV3,
    CAV4,
    ADC_LWIR,
    ADC_MWIR,
    INS1,
    INS2,
    INS3,
    INS4,
    CPA,
    CPB,
    BB1,
    BB2,
    BB3,
    BB4,
    BB5,
    BB6,
    BB7,
    BB8,
    BB9,
    BB10,
    BB11,
    BB12,
    LWIR_ADC_PRI,
    LWIR_ADC_RED,
    MWIR_ADC_PRI,
    MWIR_ADC_RED,
    NIR_ADC_PRI,
    NIR_ADC_RED,
    VIS_ADC_PRI,
    VIS_ADC_RED,
    PCLW_ADC_PRI,
    PCLW_ADC_RED,
    RC1,
    RC2,
    RC3,
    RC4,
    RC5,
    VR_RC_LW_FPA,
    NUM_TEMPS
  } i;

  struct {
    char *vname;     /*vdata_name*/
    char *fname;     /*field_name*/
    uint16 buffer[MAX_NUM_SCANS];  /*data buffer; all eng_vdata are type of uint16*/
    uint16 max_valid_value;
  } temp[NUM_TEMPS] = {
    {"Telemetry Major Cycle 5A of 7" , "CR_RC_CFPA_T1SET", {0},  1   }, /* FP_T1SET     */
    {"Telemetry Major Cycle 5A of 7" , "CR_RC_CFPA_T3SET", {0},  1   }, /* FP_T3SET     */
    {"Telemetry Major Cycle 5A of 7" , "CR_RC_LWHTR_ON"  , {0},  1   }, /* LWHTR_ON     */
    {"Telemetry Major Cycle 5A of 7" , "CR_RC_SMHTR_ON"  , {0},  1   }, /* SMHTR_ON     */
    {"Engineering Temperature data"  , "TA_AO_VIS_FPAE"  , {0},  4095}, /* FP1          */
    {"Engineering Temperature data"  , "TA_AO_NIR_FPAE"  , {0},  4095}, /* FP2          */
    {"Engineering Temperature data"  , "TA_RC_SMIR_CFPAE", {0},  4095}, /* FP3          */
    {"Engineering Temperature data"  , "TA_RC_LWIR_CFPAE", {0},  4095}, /* FP4          */
    {"Engineering Temperature data"  , "TP_SA_RCT1_MIRE" , {0},  4095}, /* MIR1         */
    {"Engineering Temperature data"  , "TP_SA_RCT2_MIRE" , {0},  4095}, /* MIR2         */
    {"Telemetry Major Cycle 20 of 63", "TP_SA_A_MTR"     , {0},  511 }, /* SCN_MTR      */
    {"Telemetry Major Cycle 3 of 63" , "TP_MF_CALBKHD_SR", {0},  511 }, /* CAV1         */
    {"Telemetry Major Cycle 23 of 63", "TP_SR_SNOUT"     , {0},  511 }, /* CAV2         */
    {"Telemetry Major Cycle 4 of 63" , "TP_MF_Z_BKHD_BB" , {0},  511 }, /* CAV3         */
    {"Telemetry Major Cycle 3 of 63" , "TP_MF_CVR_OP_SR" , {0},  511 }, /* CAV4         */
    {"Telemetry Major Cycle 14 of 63", "TA_PVLW_PWB4_10" , {0},  511 }, /* ADC_LWIR     */
    {"Telemetry Major Cycle 14 of 63", "TA_PVSM_PWB6_12" , {0},  511 }, /* ADC_MWIR     */
    {"Telemetry Major Cycle 0 of 63" , "TP_AO_SMIR_OBJ"  , {0},  4095}, /* INS1         */
    {"Telemetry Major Cycle 0 of 63" , "TP_AO_LWIR_OBJ"  , {0},  4095}, /* INS2         */
    {"Telemetry Major Cycle 0 of 63" , "TP_AO_SMIR_LENS" , {0},  4095}, /* INS3         */
    {"Telemetry Major Cycle 0 of 63" , "TP_AO_LWIR_LENS" , {0},  4095}, /* INS4         */
    {"Telemetry Major Cycle 1 of 7"  , "CR_CP_A_ON_M"    , {0},  1   }, /* CPA A on (1553 telemetry) */
    {"Telemetry Major Cycle 1 of 7"  , "CR_CP_B_ON_M"    , {0},  1   }, /* CPB B on (1553 telemetry) */
    {"Engineering BB data"           , "TP_BB_TEMP01"    , {0},  4095}, /* BB1          */
    {"Engineering BB data"           , "TP_BB_TEMP02"    , {0},  4095}, /* BB2          */
    {"Engineering BB data"           , "TP_BB_TEMP03"    , {0},  4095}, /* BB3          */
    {"Engineering BB data"           , "TP_BB_TEMP04"    , {0},  4095}, /* BB4          */
    {"Engineering BB data"           , "TP_BB_TEMP05"    , {0},  4095}, /* BB5          */
    {"Engineering BB data"           , "TP_BB_TEMP06"    , {0},  4095}, /* BB6          */
    {"Engineering BB data"           , "TP_BB_TEMP07"    , {0},  4095}, /* BB7          */
    {"Engineering BB data"           , "TP_BB_TEMP08"    , {0},  4095}, /* BB8          */
    {"Engineering BB data"           , "TP_BB_TEMP09"    , {0},  4095}, /* BB9          */
    {"Engineering BB data"           , "TP_BB_TEMP10"    , {0},  4095}, /* BB10         */
    {"Engineering BB data"           , "TP_BB_TEMP11"    , {0},  4095}, /* BB11         */
    {"Engineering BB data"           , "TP_BB_TEMP12"    , {0},  4095}, /* BB12         */
    {"Telemetry Major Cycle 4A of 7" , "CR_PVLW_A_ON"    , {0},  1   }, /* LWIR_PRI_ADC A on */
    {"Telemetry Major Cycle 4A of 7" , "CR_PVLW_B_ON"    , {0},  1   }, /* LWIR_RED_ADC B on */
    {"Telemetry Major Cycle 4B of 7" , "CR_PVSM_A_ON"    , {0},  1   }, /* MWIR_PRI_ADC A on */
    {"Telemetry Major Cycle 4B of 7" , "CR_PVSM_B_ON"    , {0},  1   }, /* MWIR_RED_ADC B on */
    {"Telemetry Major Cycle 4B of 7" , "CR_PVNIR_A_ON"   , {0},  1   }, /* NIR_PRI_ADC  A on */
    {"Telemetry Major Cycle 4B of 7" , "CR_PVNIR_B_ON"   , {0},  1   }, /* NIR_RED_ADC  B on */
    {"Telemetry Major Cycle 4B of 7" , "CR_PVVIS_A_ON"   , {0},  1   }, /* VIS_PRI_ADC  A on */
    {"Telemetry Major Cycle 4B of 7" , "CR_PVVIS_B_ON"   , {0},  1   }, /* VIS_RED_ADC  B on */
    {"Telemetry Major Cycle 4A of 7" , "CR_PCLW_A_ON"    , {0},  1   }, /* PCLW_PRI_ADC A on */
    {"Telemetry Major Cycle 4A of 7" , "CR_PCLW_B_ON"    , {0},  1   }, /* PCLW_RED_ADC B on */
    {"Telemetry Major Cycle 18 of 63", "TA_RC_CS"        , {0},  4095}, /* RC1          */
    {"Telemetry Major Cycle 18 of 63", "TA_RC_CS_OG"     , {0},  4095}, /* RC2          */
    {"Telemetry Major Cycle 19 of 63", "TA_RC_IS"        , {0},  4095}, /* RC3          */
    {"Telemetry Major Cycle 19 of 63", "TA_RC_IS_OG"     , {0},  4095}, /* RC4          */
    {"Telemetry Major Cycle 19 of 63", "TA_RC_OS_OG"     , {0},  4095}, /* RC5          */
    {"Telemetry Major Cycle 30 of 63", "VR_RC_LW_FPA_HTR", {0},  1023}  /* VR_RC_LW_FPA */
  };

  int16 FP_SPS = 0; /*FP Set Point State; temporary value*/
  int16 thermistor;
  uint16 *BB_temp[NUM_BB_THERMISTORS];
  int16 S = 0; /*scan_index*/
  int16 XS = 0;

  float32 tcavsum;
  int32 j, ncav;
  uint16           CP_index[MAX_NUM_SCANS];
  float32          *C_CAV;
  float32          *C_INS;
  float32          *C_MIR;
  Engineering_Coefficients_t Eng_Coeffs;

  int32 vd_id;

/********* Macro: CONVERT_TEMP_POLY(tdk,tdn,c,dnsat,toffset) *************/
/* This macro converts the temperature DN using the 5th order polynomial
 * coefficients appropriate to the DN being converted.  It also checks for
 * invalid DN (zero or saturated).  If invalid, nothing is assigned to the
 * variable being computed (remains at its initialized state).
 *    tdk = temperature in degrees (Kelvin)
 *    tdn = temperature digital number
 *    c   = polynomial coefficients
 *    dnsat = saturation DN for this temperature DN
 *    toffset = offset in degrees to add to converted temperature.
 *              (this may be required to convert from Celsius to Kelvin)
 */
#define CONVERT_TEMP_POLY(tdk,tdn,c,dnsat,toffset) \
 if (tdn > 0 && tdn < dnsat) tdk = toffset +   \
 c[0] + tdn * (c[1] + tdn * (c[2] + tdn * (c[3] + tdn * (c[4] + tdn * c[5]))));

/***************************************************************************/


  /* Use different sets of engineering coefficients for
   * Terra (PFM) and Aqua (FM1)
   */

  if (satellite_ID == TERRA)        /* Satellite is Terra */
    Eng_Coeffs = Engineering_Coefficients_PFM;
  else                              /* Satellite is Aqua */
    Eng_Coeffs = Engineering_Coefficients_FM1;

     /* Initialize structure members */

  for (S = 0; S < MAX_NUM_SCANS; S++)
  {
    for (j = 0; j < NUM_FOCAL_PLANES; j++)temps->fp[j][S] = 0;
    temps->mir1[S] = 0;
    temps->mir2[S] = 0;
    temps->mir_avg[S] = 0;
    temps->scn_mtr[S] = 0;
    temps->cav[S] = 0;
    temps->ins[S] = temps->bb_avg[S] = 0;
    for (j = 0; j < NUM_BB_THERMISTORS; j++)temps->bb[j][S] = 0;
    temps->fp_set_point_state[S] = 0;
    temps->lwir_adc_index[S] = 0;
    temps->mwir_adc_index[S] = 0;
    for (j = 0; j < NUM_T_CAV_THERMISTORS; j++) temps->cav_temp[j][S] = 0;
    for (j = 0; j < NUM_T_INS_THERMISTORS; j++) temps->ins_temp[j][S] = 0;
    for (j = 0; j < NUM_T_RC_VALUES; j++) temps->rc_temp[j][S] = 0;
    temps->vr_lwir[S] = VOLTAGE_BAD_VALUE;
  }

  /* read in data */

  for (i = FP_T1SET; i < NUM_TEMPS; i++)
  {
    /* Must read entire vdata for INS1 through INS4 for extracts. */
    if (i >= INS1 && i <= INS4) {
    returnStatus = read_vdata (v_id, start_scan, -1,
                               temp[i].vname, temp[i].fname,
                               (VOIDP)temp[i].buffer);
    } else {
    returnStatus = read_vdata (v_id, start_scan, num_scans,
                               temp[i].vname, temp[i].fname,
                               (VOIDP)temp[i].buffer);
    }

    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "read_vdata", 0, NULL, False);
      return returnStatus;
    }

    /* Check values against valid range. */

    for (S = 0; S < num_scans; S++)
    {
      if (temp[i].buffer[S] > temp[i].max_valid_value)
      {
        char msgbuf[256];
        sprintf(msgbuf, "Value[%d] of \"%s\" is out of the range [0,%d].",
                (int) S, temp[i].fname, (int) temp[i].max_valid_value);
        returnStatus = MODIS_F_OUT_OF_RANGE;
        L1BErrorMsg(location, returnStatus, msgbuf, NULL, 0,
            Invalid_MOD01_Msg, False);
        return returnStatus;
      }
    }
  }

  /* fill in values for instrument temperatures */

  for (i = INS1; i <= INS4; i++)
  {
    vd_id = VSattach(v_id, VSfind(v_id, "Telemetry Major Cycle 0 of 63"), "r");
    returnStatus = Fill_Invalid_Temp_DNs(VSelts(vd_id), 0,
        DN_SAT_12BITS, temp[i].buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
          "Fill_Invalid_Temp_DNs", 0, NULL, False);
      return returnStatus;
    }
    VSdetach(vd_id);
    if (extract_line_offset > 0)
      memcpy(&temp[i].buffer[0], &temp[i].buffer[extract_line_offset],
             sizeof(int16)*extract_line_count);
  }


  /* Determine Primary or Redundant CP indicator for each scan. */

    returnStatus = Get_Electronics_index(num_scans, temp[CPA].buffer,
                   temp[CPB].buffer, CP_index);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not get CPA or CPB indexes",
                "Get_Electronics_index", 0, NULL, False);
    return returnStatus;
  }


    for (i = BB1; i <= BB12; i++)
  {
    thermistor = i - BB1;
    BB_temp[thermistor] = temp[i].buffer;
  }

  /*
   * Compute black_body temperature with outlier algorithm implemented
   */

   returnStatus = Compute_BB_Temperature(tables,
                                        start_scan,
                                        num_scans,
                                        CP_index,
                                        BB_temp,
                                        satellite_ID,
                                        temps);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
        "Compute_BB_Temperature", 0, NULL, False);
    return returnStatus;
  }

  /*
   * Determine primary or redundant mwir adc board value for each scan.
   */

  returnStatus = Get_Electronics_index(num_scans,
                                       temp[MWIR_ADC_PRI].buffer,
                                       temp[MWIR_ADC_RED].buffer,
                                       &temps->mwir_adc_index[start_scan]);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not get MWIR ADC indexes",
                "Get_Electronics_index", 0, NULL, False);
    return returnStatus;
  }

  /*
   * Determine primary or redundant lwir adc board value for each scan.
   */

  returnStatus = Get_Electronics_index(num_scans,
                                       temp[LWIR_ADC_PRI].buffer,
                                       temp[LWIR_ADC_RED].buffer,
                                       &temps->lwir_adc_index[start_scan]);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not get LWIR ADC indexes",
                "Get_Electronics_index", 0, NULL, False);
    return returnStatus;
  }

  /*
   * Determine primary or redundant nir adc board value for each scan.
   */

  returnStatus = Get_Electronics_index(num_scans,
                                       temp[NIR_ADC_PRI].buffer,
                                       temp[NIR_ADC_RED].buffer,
                                       &temps->nir_adc_index[start_scan]);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not get NIR ADC indexes",
                "Get_Electronics_index", 0, NULL, False);
    return returnStatus;
  }

  /*
   * Determine primary or redundant vis adc board value for each scan.
   */

  returnStatus = Get_Electronics_index(num_scans,
                                       temp[VIS_ADC_PRI].buffer,
                                       temp[VIS_ADC_RED].buffer,
                                       &temps->vis_adc_index[start_scan]);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not get VIS ADC indexes",
                "Get_Electronics_index", 0, NULL, False);
    return returnStatus;
  }

  /*
   * Determine primary or redundant pclw adc board value for each scan.
   */

  returnStatus = Get_Electronics_index_special(num_scans,
                                       temp[PCLW_ADC_PRI].buffer,
                                       temp[PCLW_ADC_RED].buffer,
                                       &temps->pclw_adc_index[start_scan]);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not get PCLW ADC indexes",
                "Get_Electronics_index", 0, NULL, False);
    return returnStatus;
  }

  for (S = 0; S < num_scans; S++)
  {
    XS = S + start_scan;

         /* Compute focal plane assembly 1 and 2 temperatures.
          * The polynomial result is in degrees Celsius, so we
          * use the CTOK macro.
          */

    CONVERT_TEMP_POLY(temps->fp[0][XS],temp[FP1].buffer[S],Eng_Coeffs.C_FP1,
                      DN_SAT_12BITS,Celsius_To_Kelvin_Offset);
    CONVERT_TEMP_POLY(temps->fp[1][XS],temp[FP2].buffer[S],Eng_Coeffs.C_FP2,
                      DN_SAT_12BITS,Celsius_To_Kelvin_Offset);

        /* Determine FP_Set_Point_State, needed for FP_3 & _4 */

    if (temp[LWHTR_ON].buffer[S] == 0 && temp[SMHTR_ON].buffer[S] == 0)
      FP_SPS = 0; /*FP heaters off, running open loop*/
    else if (temp[FP_T1SET].buffer[S] == 1 && temp[FP_T3SET].buffer[S] == 0)
      FP_SPS = 0;
    else if (temp[FP_T1SET].buffer[S] == 0 && temp[FP_T3SET].buffer[S] == 0)
      FP_SPS = 1;
    else if (temp[FP_T1SET].buffer[S] == 0 && temp[FP_T3SET].buffer[S] == 1)
      FP_SPS = 2;
    else   /*temp[FP_T1SET].buffer[S] == 1 && temp[FP_T3SET].buffer[S] ==1,
             not expected!*/
      return (MODIS_F_NOK);

    /* Keep FP_SPS
    */
    if (temp[LWHTR_ON].buffer[S] == 0 && temp[SMHTR_ON].buffer[S] == 0)
      temps->fp_set_point_state[XS] = FP_SPS;
    else
      temps->fp_set_point_state[XS] = FP_SPS + 1;

  /* Note that for focal planes 3 and 4, the polynomial evaluates
   * the temperature in degrees Kelvin and no conversion from
   * Celsius to Kelvin is needed.
   */

    CONVERT_TEMP_POLY(temps->fp[2][XS], temp[FP3].buffer[S],
        Eng_Coeffs.C_FP3[FP_SPS], DN_SAT_12BITS,0.);
    CONVERT_TEMP_POLY(temps->fp[3][XS], temp[FP4].buffer[S],
        Eng_Coeffs.C_FP4[FP_SPS], DN_SAT_12BITS,0.);

    /*   Coefficients which may differ for A/B side */
    if (CP_index[S] == ELECTRONICS_PRIMARY)
    {
      C_CAV = Eng_Coeffs.C_CAV_A;
      C_INS = Eng_Coeffs.C_INS_A;
      C_MIR = Eng_Coeffs.C_MIR_A;
    }
    else
    {
      C_CAV = Eng_Coeffs.C_CAV_B;
      C_INS = Eng_Coeffs.C_INS_B;
      C_MIR = Eng_Coeffs.C_MIR_B;
    }
        /* compute scan mirror temperature */

    CONVERT_TEMP_POLY(temps->mir1[XS], temp[MIR1].buffer[S],
        C_MIR, DN_SAT_12BITS, Celsius_To_Kelvin_Offset);
    CONVERT_TEMP_POLY(temps->mir2[XS], temp[MIR2].buffer[S],
        C_MIR, DN_SAT_12BITS, Celsius_To_Kelvin_Offset);

    /* Compute average mirror temperature. The algorithm similar
     * to that for the cavity below. Note, we use the local variables
     * ncav and tcavsum as a convenience only.
     */

    ncav = 0;
    tcavsum = 0.0;
    if (tables->T_mir_function_flag[0] == ON && temp[MIR1].buffer[S] > 0 &&
                                      temp[MIR1].buffer[S] < DN_SAT_12BITS)
    {
      ncav++;
      tcavsum += temps->mir1[XS];
    }

    if (tables->T_mir_function_flag[1] == ON && temp[MIR2].buffer[S] > 0 &&
                                       temp[MIR2].buffer[S] < DN_SAT_12BITS)
    {
      ncav++;
      tcavsum += temps->mir2[XS];
    }
    if (ncav > 0)
      temps->mir_avg[XS] = tcavsum/(float32)ncav;
    else
      temps->mir_avg[XS] = tables->T_mir_default;

    /* Compute the Scan Motor Temperature. */

    CONVERT_TEMP_POLY(temps->scn_mtr[XS], temp[SCN_MTR].buffer[S], C_CAV,
                      DN_SAT_9BITS, Celsius_To_Kelvin_Offset);

    /* Compute cavity temperature.  The algorithm is to average the valid
     * values (checking function flag and DN).  If there are no valid values,
     * use the default from the LUT.
     */

    ncav = 0;
    tcavsum = 0.;
    for (j = 0, i = CAV1; j < NUM_T_CAV_THERMISTORS && i <= CAV4; j++, i++)
    {
      if (temp[i].buffer[S] > 0 && temp[i].buffer[S] < DN_SAT_9BITS)
      {
        CONVERT_TEMP_POLY(temps->cav_temp[j][XS],
                          temp[i].buffer[S],
                          C_CAV,
                          DN_SAT_9BITS,
                          Celsius_To_Kelvin_Offset);
        if (tables->T_cav_function_flag[j] == ON)
        {
          ncav++;
          tcavsum += temps->cav_temp[j][XS];
        }
      }
    }
    if (ncav > 0)
      temps->cav[XS] = tcavsum / (float32) ncav;
    else
      temps->cav[XS] = tables->T_cav_default;


        /* Compute instrument temperature.  The algorithm is to use the first
         * valid value of the available thermistors.  If none are valid,
         * then use the default.
         */

    temps->ins[XS] = -1;
    for (j = 0, i = INS1; j < NUM_T_INS_THERMISTORS && i <= INS4; j++, i++)
    {
      if (temp[i].buffer[S] > 0 && temp[i].buffer[S] < DN_SAT_12BITS)
      {
        CONVERT_TEMP_POLY(temps->ins_temp[j][XS],
                          temp[i].buffer[S],
                          C_INS,
                          DN_SAT_12BITS,
                          Celsius_To_Kelvin_Offset + tables->T_ins_offset[j]);
        if (temps->ins[XS] <= 0 && tables->T_ins_function_flag[j] == ON)
          temps->ins[XS] = temps->ins_temp[j][XS];
      }
    }
    if (temps->ins[XS] <= 0)
      temps->ins[XS] = tables->T_ins_default;

      /* Compute the RC temperatures and single voltage */

    for (j = 0, i = RC1; j < NUM_T_RC_VALUES && i <= RC5; j++, i++)
    {
      CONVERT_TEMP_POLY(temps->rc_temp[j][XS], temp[i].buffer[S],
          Eng_Coeffs.C_RC[j], DN_SAT_12BITS, 0.);
    }
    CONVERT_TEMP_POLY(temps->vr_lwir[XS],
                      temp[VR_RC_LW_FPA].buffer[S],
                      Eng_Coeffs.C_RC_VR_FPA,
                      DN_SAT_10BITS,
                      0.);


  }/*Loop over scan*/


  return (returnStatus);
}


PGSt_SMF_status Compute_BB_Temperature(emiss_tables_t *tables,
                                       int32          start_scan,
                                       int32          num_scans,
                                       uint16         *CP_index,
                                       uint16         **raw_BB_temp,
                                       int32          satellite_ID,
                                       Temperatures_t *temps)
/*
!C**************************************************************************
!Description:   Calculate temperatures for each BB thermistor and the
                thermistor averaged temperature with outlier rejection.

!Input Parameters:
     emiss_tables_t * tables         Contains BB_Weight, used in here
     int32            start_scan     Start scan index fo filling arrays in temps
     int32            num_scans      Number of scans to fill in arrays in temps
     uint16           CP_index       CP A or B on (0 = "A", 1 = "B")
     uint16 **        raw_BB_temp    2D Array of BB temperature DNs
                                     [thermistor][scan]
     int32            satellite_ID   TERRA = 0, AQUA = 1,
                                         INVALID_SATELLITE_ID = -1

!Output Parameters:
     Temperatures_t *temps

!Revision History:
 Revision 02.01 March 2003, Razor Issue #173
 Changed the type of variable "thermistor" from "int8" to "int16" for ANSI-C 
 compliance.
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 02.00 January 4, 2002  Razor Issue #154
 Added Engineering Temperature structure to accommodate Terra/Aqua
 temperature coefficient differences and changed code to use the structure
 variables.
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.05, January 29, 2001
 Improved efficiency of calculation of TBB.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.04 November 22, 1999
 Added calls to L1BErrorMsg.  Changed logic of immediately error-out
 to instead simply write error message and return to parent (so that
 parent can identify LUN of file).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 September 30, 1999
 Added counting the number of BB thermistor outliers on each scan.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.02 August 1999
 Used Get_Electronics_index() to determine CPA or CPB and added a check for denom > 1
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 01.01 May 1999
 Added checks on BB DN.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 Aug 1997
 (This function is part of the function Read_Convert_Temperature,developed
  by Zhiong Hao, modified by Shi-Yue Qiu and David Catozzi)
 Initial development
 Zhenying Gu (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
  The calculations for standard deviation are accomplished in double
  precision due to the particular formula being used.

  If an error occurs in this function, write the error message and then
  return to the parent (so that the parent can identify the LUN of
  the file).

!END********************************************************************
*/

{
  PGSt_SMF_status  returnStatus  = MODIS_S_OK;
  int16            thermistor    = 0;
  int16            S             = 0;
  int16            XS            = 0;
  uint16           DN            = 0;
  float32          V0            = 0;
  float32          I             = 0;
  double           denom         = 0;
  double           sum1          = 0;
  double           sum2          = 0;
  float32          Rc            = 0;
  float32          R             = 0;
  float32          lnR           = 0;
  float32          TBB           = 0;
  float32          temp1         = 0;
  float32          mean          = 0;
  float32          sdev          = 0;
  float32          Tmin          = 0;
  float32          Tmax          = 0;
  Engineering_Coefficients_t Eng_Coeffs;
  char *location = "Compute_BB_Temperature";

  /* Use different sets of engineering coefficients for
   * Terra (PFM) and Aqua (FM1)
   */

  if (satellite_ID == TERRA)        /* Satellite is Terra */
    Eng_Coeffs = Engineering_Coefficients_PFM;
  else                              /* Satellite is Aqua */
    Eng_Coeffs = Engineering_Coefficients_FM1;


  for (S = 0; S < num_scans; S++)
  {
    XS = S + start_scan;
    if (CP_index[S] == ELECTRONICS_PRIMARY)
    {
      V0 = Eng_Coeffs.V0_CPA;
      I  = Eng_Coeffs.I_CPA;
    }
    else
    {
      V0 = Eng_Coeffs.V0_CPB;
      I  = Eng_Coeffs.I_CPB;
    }

    denom = 0;
    sum1 = 0;
    sum2 = 0;
    for (thermistor = 0; thermistor < NUM_BB_THERMISTORS; thermistor++)
    {
      if (raw_BB_temp[thermistor][S] > 0 &&
          raw_BB_temp[thermistor][S] < DN_SAT_12BITS)
      {
        DN = raw_BB_temp[thermistor][S];
        Rc = (5.0*DN/4096.0-(V0+2.5))/I;
        if(Eng_Coeffs.Rp - Rc == 0)
        {
          returnStatus = MODIS_F_NOK;
          L1BErrorMsg(location, returnStatus,
                      "Algorithm failure computing BB temperature.\n"
                      "Will get a divide by zero because of Rp-Rc.",
                      NULL, 0, NULL, False);
          return returnStatus;
        }
        R = Rc*(Eng_Coeffs.Rp)/(Eng_Coeffs.Rp - Rc);

        lnR=log((double) R);
        if(lnR == 0)
        {
          returnStatus = MODIS_F_NOK;
          L1BErrorMsg(location, returnStatus,
             "Algorithm failure computing BB temperature. Will get a log(0).",
             NULL, 0, NULL, False);
          return returnStatus;
        }

        TBB = Eng_Coeffs.BB_a[thermistor][0]
              + lnR * (Eng_Coeffs.BB_a[thermistor][1]
              + lnR * (Eng_Coeffs.BB_a[thermistor][2]
              + lnR * Eng_Coeffs.BB_a[thermistor][3]));
        temps->bb[thermistor][XS] = 1.0/TBB;

        temp1 = tables->BB_Weight[thermistor] * temps->bb[thermistor][XS];
        sum1 += (double) temp1;
        sum2 += (double) temp1 * (double) temp1;
        denom += (double) tables->BB_Weight[thermistor];
      }
      else
      {
        temps->bb[thermistor][XS] = 0;
      }
    }
    if (denom > 1)                        /* multiple values contribute */
    {
      mean = (float32) (sum1 / denom);
      sdev = (float32) sqrt(fabs(sum2 - sum1*sum1/denom) / (denom - 1.0));

      if ( sdev < 1.0E-6 )
      {
        temps->bb_avg[XS] = mean;
        temps->num_thermistor_outliers[XS] = 0;
      }
      else
      {
        Tmin = mean - 3.0*sdev;   /* lower bound on BB temperature to include */
        Tmax = mean + 3.0*sdev;   /* upper bound on BB temperature to include */

        temps->bb_avg[XS] = 0;
        temps->num_thermistor_outliers[XS] = 0;
        denom = 0;
        for (thermistor = 0; thermistor < NUM_BB_THERMISTORS; thermistor++)
        {
          if ( temps->bb[thermistor][XS] < Tmin ||
              temps->bb[thermistor][XS] > Tmax )
          {
            if (tables->BB_Weight[thermistor] > TOLERANCE)
              temps->num_thermistor_outliers[XS] =
                  temps->num_thermistor_outliers[XS] + 1;
            continue;
          }

          temps->bb_avg[XS] +=
              tables->BB_Weight[thermistor] * temps->bb[thermistor][XS];
          denom += tables->BB_Weight[thermistor];
        }
        if ( denom < TOLERANCE )
        {
          temps->bb_avg[XS] = mean;
          temps->num_thermistor_outliers[XS] = 0;
        }
        else
          temps->bb_avg[XS] /= denom;
      }
    }
    else if (denom > TOLERANCE)           /* only one value contributes */
    {
      mean = (float32) (sum1 / denom);
      temps->bb_avg[XS] = mean;
      temps->num_thermistor_outliers[XS] = 0;
    }
    else                                   /* no values */
    {
      temps->bb_avg[XS] = 0.;
      temps->num_thermistor_outliers[XS] = 0;
    }
  }   /* end loop for scan */
  return (returnStatus);

}

PGSt_SMF_status Write_Geo_OBC_SDS(int32 obc_sd_id, int32 num_scans)
/*
!C************************************************************************
!Description:  Write SD Sun Azimuth and SD Sun Zenigh into L1B_OBC file.

!Input Parameters:
   int32    obc_sd_id     file id of obc file
   int32    num_scans     number of scans in L1A file
!Output Parameters:

  Output file: L1B OBC/Eng file

!Revision History:
 Revision 01.01 November 22, 1999
 Added use of L1BErrorMsg
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 May 4 1999
 initial development
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

 This function writes to the OBC file, so that the LUN used for write
 error messages is assumed to be L1B_OBC_FILE.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32           num_scans_geo;
  char            fname_geo[512];
  int32           geo_sd_id;
  PGSt_integer    Version = 1;
  float32         sun_buffer[MAX_NUM_SCANS];
  intn            hdf_return;
  char *location = "Write_Geo_OBC_SDS";

  if (PGS_PC_GetReference(GEOLOCATION_FILE, &Version, fname_geo) !=
      PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus,
                "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  geo_sd_id = SDstart(fname_geo, DFACC_RDONLY);
  if (geo_sd_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus,
                "Could not open file for SD read access.",
                "SDstart", GEOLOCATION_FILE,
                "The file may be missing, corrupted or not an HDF-4 file.",
                 True);
    return returnStatus;
  }

  returnStatus = read_attribute(geo_sd_id, GEO_NUM_SCANS_NAME, DFNT_INT32,
                                &num_scans_geo);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "read_attribute", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  if (num_scans_geo != num_scans)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "Number of scans in GEO file does not match the current L1A file.",
                NULL, GEOLOCATION_FILE, "The PCF may have the wrong names.", True);
    return returnStatus;
  }

  returnStatus = read_sds_rank1(geo_sd_id,
                                GEO_SD_SUN_AZ_SDS_NAME,
                                num_scans_geo,
                                sun_buffer);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "read_sds_rank1", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  returnStatus = write_sds_rank1(obc_sd_id,
                                 OBC_SD_SUN_AZ_SDS_NAME,
                                 NUM_SCANS_DIM_NAME,
                                 num_scans,
                                 "float32",
                                 sun_buffer);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "write_sds_rank1", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  returnStatus = read_sds_rank1(geo_sd_id,
                                GEO_SD_SUN_ZEN_SDS_NAME,
                                num_scans_geo,
                                sun_buffer);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "read_sds_rank1", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  returnStatus = write_sds_rank1(obc_sd_id,
                                 OBC_SD_SUN_ZEN_SDS_NAME,
                                 NUM_SCANS_DIM_NAME,
                                 num_scans,
                                 "float32",
                                 sun_buffer);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "write_sds_rank1", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  hdf_return = SDend(geo_sd_id);
  if (hdf_return == FAIL)
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus, NULL, "SDend", GEOLOCATION_FILE,
                "Memory or the disk file must have become corrupted.", True);
    return returnStatus;
  }

  return (returnStatus);
}

PGSt_SMF_status Write_L1B_OBCEng (L1A_granule_t         *L1A_Gran,
                                  lookup_tables_t       *tables,
                                  Moon_arrays_t         *moon_arrays,
                                  DN_OBC_Avg_t          *DN_OBC_Avg)
/*
!C************************************************************************
!Description:  Write OBC target sdss and preprocess_data sdss, and
               copy scan_level_metadata sdss, pixel_quality sdss,
               engineering_memory sdss, engineering vdata, SD Sun
               Azimuth, SD Sun Zenith, and DN OBC averages into
               L1B_OBC file.

!Input Parameters:
    L1A_granule_t         *L1A_Gran
    lookup_tables_t       *tables
    Moon_arrays_t         *moon_arrays
    DN_OBC_Avg_t    *DN_OBC_Avg
!Output Parameters:

  Output file: L1B OBC/Eng file

!Revision History:
 Revision 02.16 May 4 1999
 Deleted function call to Write_Preprocess_Data() because SV_avg
 and SV_uncert are deleted. Added writing the moon vector and moon
 in KOB arrays. Added writing DN OBC averages. Added function
 Write_Geo_OBC_SDS()
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.10 July 1 1998
 Added function call to Write_Preprocess_Data() back
 (had been removed since V2.0)
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 June 1998
 added calculation of DN to DN_star for SD
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)

 Revision 01.00 Feb. 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  char *location = "Write_L1B_OBCEng";
  PGSt_SMF_status    returnStatus = MODIS_S_OK;
  intn               hdf_return   = 0;
  PGSt_integer       Version      = 1;
  char               file_name[PGSd_PC_FILE_PATH_MAX];
  int32              file_id      = 0; /* sd_id */
  int32              v_id         = 0;
  int16              T            = 0;/*Target index*/
  int16              R            = 0;/*Resolution index*/
  int32              band_dim     = 0;
  int32              track_dim    = 0;
  int32              frame_dim    = 0;
  int32              subsamp_dim  = 0;
  int16              num_scans    = 0;

  int16 *L1B_OBC;
  L1B_OBC = (int16*) malloc(NUM_250M_BANDS * MAX_250M_TRACK_DIM * MAX_250M_OBC_FRAME_DIM * sizeof(int16));

  char               *avg_names[NUM_L1A_RESOLUTIONS] =
                                   {"DN_obc_avg_250m",
                                    "DN_obc_avg_500m",
                                    "DN_obc_avg_1km_day",
                                    "DN_obc_avg_1km_night"};
  char               *var_names[NUM_L1A_RESOLUTIONS] =
                                   {"DN_obc_var_250m",
                                    "DN_obc_var_500m",
                                    "DN_obc_var_1km_day",
                                    "DN_obc_var_1km_night"};
  char               *omask_names[NUM_L1A_RESOLUTIONS] =
                                   {"DN_obc_outlier_mask_250m",
                                    "DN_obc_outlier_mask_500m",
                                    "DN_obc_outlier_mask_1km_day",
                                    "DN_obc_outlier_mask_1km_night"};

  char               *subsamp_dim_names[NUM_L1A_RESOLUTIONS] =
                                   {"250m_subsamples",
                                    "500m_subsamples",
                                    "1km_subsamples",
                                    "1km_subsamples"};

#define WRITE_DN_BG_AT_RES(r, avg, var, omask)                     \
  band_dim    = L1A_BANDS_AT_RES[r];                               \
  track_dim   = num_scans * DETECT_PER_BAND_AT_RES[r];             \
  subsamp_dim = BAND_RATIO_AT_RES[r];                              \
  returnStatus = write_sds_rank3(file_id, avg_names[r],            \
                                 L1A_TRACK_DIM_NAME[r],            \
                                 L1A_BAND_DIM_NAME[r],             \
                                 subsamp_dim_names[r],             \
                                 track_dim, band_dim, subsamp_dim, \
                                 "float32", avg);                  \
  if (returnStatus != MODIS_S_OK)                                  \
  {                                                                \
    L1BErrorMsg(location, returnStatus, NULL,                      \
                "write_sds_rank3", L1B_OBC_FILE, NULL, True);      \
    return returnStatus;                                           \
  }                                                                \
  returnStatus = write_sds_rank3(file_id, var_names[r],            \
                                 L1A_TRACK_DIM_NAME[r],            \
                                 L1A_BAND_DIM_NAME[r],             \
                                 subsamp_dim_names[r],             \
                                 track_dim, band_dim, subsamp_dim, \
                                 "float32", var);                  \
  if (returnStatus != MODIS_S_OK)                                  \
  {                                                                \
    L1BErrorMsg(location, returnStatus, NULL,                      \
                "write_sds_rank3", L1B_OBC_FILE, NULL, True);      \
    return returnStatus;                                           \
  }                                                                \
  returnStatus = write_sds_rank4(file_id, omask_names[r],          \
                                 L1A_TRACK_DIM_NAME[r],            \
                                 L1A_BAND_DIM_NAME[r],             \
                                 subsamp_dim_names[r], "2",        \
                                 track_dim, band_dim, subsamp_dim, \
                                 2, "uint32", omask);              \
  if (returnStatus != MODIS_S_OK)                                  \
  {                                                                \
    L1BErrorMsg(location, returnStatus, NULL,                      \
                "write_sds_rank4", L1B_OBC_FILE, NULL, True);      \
    return returnStatus;                                           \
  }

  /*
   * Prepare for opening file: convert logical ID to physical name
   */

  returnStatus = PGS_PC_GetReference (L1B_OBC_FILE, &Version, file_name);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus,
                "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /*
   * Create/open file for HDF science data (SD) access.
   */

  file_id = SDstart (file_name, DFACC_CREATE);
  if (file_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_CREATED;
    L1BErrorMsg(location, returnStatus,
                "Could not create file for SD write access.",
                "SDstart", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /*
   * Open file for HDF Vdata access and initialize Vdata interface.
   * Call Hopen() and Vstart() in the same place.
   */

  v_id = Hopen(file_name, DFACC_RDWR, 0);
  if (v_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus,
                "Could not open file for Vdata read/write access.",
                "Hopen", L1B_OBC_FILE,
                "The file may be corrupted or not an HDF-4 file.", True);
    return returnStatus;
  }

  hdf_return = Vstart (v_id);
  if (hdf_return == FAIL)
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not initialize Vdata interface.",
                "Vstart", L1B_OBC_FILE,
                "The file may be corrupted or not an HDF-4 file.", True);
    return returnStatus;
  }

  /*
   * Write OBC sector DN SDSs.  Loop through targets and resolutions and
   * read data from the L1A file and then write data to the OBC file.
   */

  num_scans = L1A_Gran->num_scans;

  for (T = 0; T < NUM_TARGETS - 1; T++) /*Excluding EV target*/
    for (R = 0; R < NUM_L1A_RESOLUTIONS; R++)
    {
      if ((RFLAG & (1 << R)) != 0) continue;

      band_dim  = L1A_BANDS_AT_RES[R];
      track_dim = L1A_Gran->num_scans * DETECT_PER_BAND_AT_RES[R];
      frame_dim = TARGET_1km_FRAMES[T] * BAND_RATIO_AT_RES[R];

      returnStatus = read_sds_rank3(L1A_Gran->sd_id, L1A_TARGET_SDS_NAME[T][R],
                                    track_dim, band_dim, frame_dim, L1B_OBC);
      if (returnStatus != MODIS_S_OK)
      {
        L1BErrorMsg(location, returnStatus, NULL,
                    "read_sds_rank3", FIRST_L1A_GRANULE, NULL, True);
        return returnStatus;
      }

      returnStatus = write_sds_rank3(file_id, L1A_TARGET_SDS_NAME[T][R],
                                     L1A_TRACK_DIM_NAME[R],
                                     L1A_BAND_DIM_NAME[R],
                                     L1A_FRAME_DIM_NAME[T][R],
                                     track_dim, band_dim, frame_dim,
                                     "int16", L1B_OBC);
      if (returnStatus != MODIS_S_OK)
      {
        L1BErrorMsg(location, returnStatus, NULL,
                    "write_sds_rank3", L1B_OBC_FILE, NULL, True);
        return returnStatus;
      }

    }

  /* Write attribute first SV or BB frame to use */

  hdf_return = SDsetattr(file_id, DN_BG_FIRST_FRAME_SDS_NAME, DFNT_INT16, 1,
                         &tables->refl.DN_obc_avg_first_frame_to_use);
  if (hdf_return == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Failed to write DN_obc_avg_first_frame_to_use.",
                "SDsetattr", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /* Write attribute number of SV or BB frames to use */

  hdf_return = SDsetattr(file_id, DN_BG_NUM_FRAMES_SDS_NAME, DFNT_INT16, 1,
                         &tables->refl.DN_obc_avg_number_of_frames_to_use);
  if (hdf_return == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Failed to write DN_obc_avg_number_of_frames_to_use.",
                "SDsetattr", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /* Write moon vector */

  returnStatus = write_sds_rank2(file_id,
                                 MOON_VECTOR_SDS_NAME,
                                 NUM_SCANS_DIM_NAME,
                                 VEC_DIM_NAME,
                                 L1A_Gran->num_scans,
                                 VEC_DIM,
                                 "float32",
                                 moon_arrays->moon_vector);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "write_sds_rank2", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /* Write moon in KOB */

  returnStatus = write_sds_rank2(file_id,
                                 MOON_IN_KOB_SDS_NAME,
                                 NUM_SCANS_DIM_NAME,
                                 NUM_BANDS_DIM_NAME,
                                 L1A_Gran->num_scans,
                                 NUM_BANDS,
                                 "int8",
                                 moon_arrays->moon_in_SV_KOB);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "write_sds_rank2", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /* Write SD sun azimuth and SD sun zenith */

  returnStatus = Write_Geo_OBC_SDS(file_id, L1A_Gran->num_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Write_Geo_OBC_SDS", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  WRITE_DN_BG_AT_RES(INDEX_250M,
                     DN_OBC_Avg->DN_obc_250m_avg,
                     DN_OBC_Avg->DN_obc_250m_var,
                     DN_OBC_Avg->DN_obc_250m_outlier_mask);

  WRITE_DN_BG_AT_RES(INDEX_500M,
                     DN_OBC_Avg->DN_obc_500m_avg,
                     DN_OBC_Avg->DN_obc_500m_var,
                     DN_OBC_Avg->DN_obc_500m_outlier_mask);

  WRITE_DN_BG_AT_RES(INDEX_1000M_DAY,
                     DN_OBC_Avg->DN_obc_1km_day_avg,
                     DN_OBC_Avg->DN_obc_1km_day_var,
                     DN_OBC_Avg->DN_obc_1km_day_outlier_mask);

  WRITE_DN_BG_AT_RES(INDEX_1000M_NIGHT,
                     DN_OBC_Avg->DN_obc_1km_night_avg,
                     DN_OBC_Avg->DN_obc_1km_night_var,
                     DN_OBC_Avg->DN_obc_1km_night_outlier_mask);


  /* Copy other L1A data......*/
  returnStatus = Copy_EngMemData (L1A_Gran->sd_id, file_id,
                     L1A_Gran->num_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Copy_EngMemData", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  returnStatus = Copy_ScanMetadata (L1A_Gran->sd_id, file_id,
                     L1A_Gran->num_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Copy_ScanMetadata", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  returnStatus = Copy_PixelQualityData (L1A_Gran->sd_id, file_id,
                     L1A_Gran->num_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Copy_PixelQualityData", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  returnStatus = Copy_EngVdata (L1A_Gran->v_id, v_id);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Copy_EngVdata", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /*
   * Close file and all accesses.
   */

  if (v_id != FAIL)
  {
    hdf_return = Vend(v_id);
    if (hdf_return == FAIL)
    {
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, "Vend",
                   L1B_OBC_FILE,
                   "Memory or the disk file must have become corrupted.",
                   True);
    }
    hdf_return = Hclose(v_id);
    if (hdf_return == FAIL)
    {
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, "Hclose",
                   L1B_OBC_FILE,
                   "Memory or the disk file must have become corrupted.",
                   True);
    }
  }
  if (file_id != FAIL)
  {
    if (SDend(file_id) == FAIL)
    {
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, "SDend",
                   L1B_OBC_FILE,
                   "Memory or the disk file must have become corrupted.",
                   True);
    }
  }

  free(L1B_OBC);
  return (returnStatus);
}

PGSt_SMF_status Copy_EngMemData (int32 in_sd_id,
                                 int32 out_sd_id,
                                 int32 num_scans)
/*
!C************************************************************************
!Description:  Copy all Engineering/Memory Data from L1A to L1B_OBC and
               calculate the QA "Moon within defined limits of SVP".

!Input Parameters:
    int32               in_sd_id       L1A file sd_id
    int32               out_sd_id      L1B file sd_id
    int32               num_scans

!Output Parameters:

!Revision History:
 Revision 01.01 March 2003, Razor Issue #173
 Deleted the value "NUM_ENG_MEM_SDS" from the enum "i", replaced it with a 
 new defined macro "NUM_ENG_MEM_SDS", and modified the corresponding "for" 
 loop for ANSI-C compliance.
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 01.00 June 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

 For writing error messages, it is assumed that this function reads from
 the current L1A granule and writes to the L1B OBC file (for setting the
 LUN).

!END********************************************************************
*/
{
#define NUM_ENG_MEM_SDS 12
  char *location = "Copy_EngMemData";
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char     *sds_name = NULL;
  char     *dim_name[2];
  int32    rank           = 2; /*max rank is 2*/
  int32    data_type      = 0;
  int32    start[2]       = {0, 0};
  int32    edge[2]        = {0, 0};
  char     buffer[MAX_NUM_SCANS * MAX_ENG_MEM_DATA_DIM];
  enum {
    FPA_AEM_CONFIG,
    SCIENCE_STATE,
    SCIENCE_ABNORMAL,
    FPA_DCR_OFFSET,
    RAW_MIR_ENC,
    RAW_VS_DEF,
    RAW_VS_ACT,
    RAW_SCI_ENG,
    RAW_HK_TELEM,
    RAW_SC_ANCIL,
    RAW_PARAM,
    RAW_PV_GAINS
  } i;

  edge[0]     = num_scans;
  dim_name[0] = "nscans";

  for (i = FPA_AEM_CONFIG; i < NUM_ENG_MEM_SDS; i++)
  {
    switch(i)
    {
      case FPA_AEM_CONFIG:
        sds_name = "fpa_aem_config";
        data_type = DFNT_INT8;
        rank = 2;
        edge[1] = 10;
        dim_name[1] = "10";
        break;
      case SCIENCE_STATE:
        sds_name = "science_state";
        data_type = DFNT_INT8;
        rank = 1;
        break;
      case SCIENCE_ABNORMAL:
        sds_name = "science_abnormal";
        data_type = DFNT_INT8;
        rank = 1;
        break;
      case FPA_DCR_OFFSET:
        sds_name = "fpa_dcr_offset";
        data_type = DFNT_INT8;
        rank = 2;
        edge[1] = 550;
        dim_name[1] = "550";
        break;
      case RAW_MIR_ENC:
        sds_name = "raw_mir_enc";
        data_type = DFNT_INT16;
        rank = 2;
        edge[1] = 78;
        dim_name[1] = "78";
        break;
      case RAW_VS_DEF:
        sds_name = "raw_vs_def";
        data_type = DFNT_INT16;
        rank = 2;
        edge[1] = 40;
        dim_name[1] = "40";
        break;
      case RAW_VS_ACT:
        sds_name = "raw_vs_act";
        data_type = DFNT_INT16;
        rank = 2;
        edge[1] = 24;
        dim_name[1] = "24";
        break;
      case RAW_SCI_ENG:
        sds_name = "raw_sci_eng";
        data_type = DFNT_INT8;
        rank = 2;
        edge[1] = 212;
        dim_name[1] = "212";
        break;
      case RAW_HK_TELEM:
        sds_name = "raw_hk_telem";
        data_type = DFNT_INT8;
        rank = 2;
        edge[1] = 128;
        dim_name[1] = "128";
        break;
      case RAW_SC_ANCIL:
        sds_name = "raw_sc_ancil";
        data_type = DFNT_INT16;
        rank = 2;
        edge[1] = 64;
        dim_name[1] = "64";
        break;
      case RAW_PARAM:
        sds_name = "raw_param";
        data_type = DFNT_INT8;
        rank = 2;
        edge[1] = 40;
        dim_name[1] = "40";
        break;
      case RAW_PV_GAINS:
        sds_name = "raw_pv_gains";
        data_type = DFNT_INT8;
        rank = 2;
        edge[1] = 550;
        dim_name[1] = "550";
        break;
    }

    /*
     * Read the SDS data from the L1A file.
     */

    returnStatus = read_sds_rankn(in_sd_id, sds_name, data_type, rank,
                             start, edge, buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "read_sds_rankn",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    /*
     * Write the data to the OBC file.
     */

    returnStatus = write_sds_rankn(out_sd_id, sds_name, data_type, rank,
                                   edge, dim_name, buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "write_sds_rankn",
                  L1B_OBC_FILE, NULL, True);
      return returnStatus;
    }

  }

  return returnStatus;
}

PGSt_SMF_status Copy_PixelQualityData (int32 in_sd_id,
                                       int32 out_sd_id,
                                       int32 num_scans)
/*
!C************************************************************************
!Description:  Copy L1A Pixel Quality Data into L1B_OBCEng file.

!Input Parameters:
    int32               in_sd_id     L1A file sd_id
    int32               out_sd_id    L1B_OBC file sd_id
    int32               num_scans    number of scans in current granule

!Output Parameters:

!Revision History:
 Revision 01.01 March 2003, Razor Issue #173
 Deleted the value "NUM_PIXEL_QUALITY_SDS" from the enum "i", replaced it 
 with a new defined macro "NUM_PIXEL_QUALITY_SDS", and modified the 
 corresponding "for" loops for ANSI-C compliance. 
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)
 
 Revision 01.00 June 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
#define NUM_PIXEL_QUALITY_SDS 5
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char *location = "Copy_PixelQualityData";
  char     *sds_name = NULL; 
  char     *dim_name[3];
  int32    j;
  int32    rank           = 3; /*max rank is 3*/
  int32    data_type      = 0;
  int32    start[3]       = {0, 0, 0};
  int32    edge[3]        = {0, 0, 0};
  int16    *buffer;
  buffer = (int16*) malloc(MAX_NUM_SCANS * 1400 * 2 * sizeof(int16));
  enum {
    SD,
    SRCA,
    BB,
    SV,
    EARTH
} i;

  data_type = DFNT_INT16;

  edge[0]     = num_scans;
  edge[2]     = 2;
  dim_name[0] = "nscans";
  dim_name[2] = "2";

  for (i = SD; i < NUM_PIXEL_QUALITY_SDS; i++)
  {
    switch(i)
    {
      case SD:
        sds_name    = "SD sector Pixel quality";
        edge[1]     = 64;
        dim_name[1] = "64";
        break;
      case SRCA:
        sds_name    = "SRCA sector Pixel quality";
        edge[1]     = 64;
        dim_name[1] = "64";
        break;
      case BB:
        sds_name    = "BB sector Pixel quality";
        edge[1]     = 64;
        dim_name[1] = "64";
        break;
      case SV:
        sds_name    = "SV sector Pixel quality";
        edge[1]     = 64;
        dim_name[1] = "64";
        break;
      case EARTH:
        sds_name    = "Earth sector Pixel quality";
        edge[1]     = 1400;
        dim_name[1] = "1400";
        if (extract_pixel_count > 0) edge[1] = extract_pixel_count;
        break;
    }

    /*
     * Read the SDS data from the L1A file.
     */

    returnStatus = read_sds_rankn(in_sd_id, sds_name, data_type, rank,
                             start, edge, buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "read_sds_rankn",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    if ((extract_pixel_offset > 0 ||
         extract_pixel_count > 0) &&
        i == EARTH) {
      edge[1] = 1400;  // reset to full width for output
      for (j=edge[0]-1; j>=0; j--) {
        memmove(&buffer[(j*edge[1]+extract_pixel_offset)*edge[2]],
               &buffer[j*extract_pixel_count*edge[2]],
               extract_pixel_count*edge[2]*sizeof(int16));
      }
      for (j=edge[0]-1; j>=0; j--) {
        memset(&buffer[j*edge[1]*edge[2]],
               0,
               extract_pixel_offset*edge[2]*sizeof(int16));
        memset(&buffer[(j*edge[1]+extract_pixel_offset+extract_pixel_count)*edge[2]],
               0,
               (edge[1]-extract_pixel_offset-extract_pixel_count)*edge[2]*sizeof(int16));
      }
    }

    /*
     * Write the data to the OBC file.
     */

    returnStatus = write_sds_rankn(out_sd_id, sds_name, data_type, rank,
                                   edge, dim_name, buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "write_sds_rankn",
                  L1B_OBC_FILE, NULL, True);
      return returnStatus;
    }

  }

  free(buffer);
  
  return returnStatus;
}

PGSt_SMF_status Copy_ScanMetadata (int32 in_sd_id,
                                   int32 out_sd_id,
                                   int32 num_scans)
/*
!C************************************************************************
!Description:  Copy L1A scan level metadata (SDSs) into L1B_OBCEng file.

!Input Parameters:
    int32               in_sd_id    L1A file sd_id
    int32               out_sd_id   L1B_OBC file sd_id
    int32               num_scans   number of scans in current granule

!Output Parameters:

!Revision History:
 Revision 01.01 March 2003, Razor Issue #173
 Deleted the value "NUM_SCANMETA_SDS" from the enum "i", replaced it 
 with a new defined macro "NUM_SCANMETA_SDS", and  modified the 
 corresponding "for" loops for ANSI-C compliance.
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 01.00 June 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
#define NUM_SCANMETA_SDS 14
  char *location = "Copy_ScanMetadata";
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char     *sds_name = NULL;
  char     *dim_name[2];
  int32    rank           = 2; /*max rank is 2*/
  int32    data_type      = 0;
  int32    start[2]       = {0, 0};
  int32    edge[2]        = {0, 0};
  int32    buffer[MAX_NUM_SCANS * 4];
  enum {
    SCAN_NUMBER,
    FRAME_COUNT_ARRAY,
    SCAN_TYPE,
    SD_START_TIME,
    SRCA_START_TIME,
    BB_START_TIME,
    SV_START_TIME,
    EV_START_TIME,
    SRCA_CALIBRATION_MODE,
    PACKET_SCAN_COUNT,
    CCSDS_APPLICATION_IDENTIFIERS,
    PACKET_EXPEDITED_DATA_FLAG,
    MIRROR_SIDE,
    SCAN_QUALITY_ARRAY
} i;

  edge[0]     = num_scans;
  dim_name[0] = "nscans";

  for (i = SCAN_NUMBER; i < NUM_SCANMETA_SDS; i++)
  {
    switch(i)
    {
      case SCAN_NUMBER:
        sds_name = "Scan number";
        data_type = DFNT_INT16;
        rank = 1;
        break;
      case FRAME_COUNT_ARRAY:
        sds_name = "Frame count array";
        data_type = DFNT_INT16;
        rank = 2;
        edge[1] = 6;
        dim_name[1] = "6";
        break;
      case SCAN_TYPE:
        sds_name = "Scan Type";
        data_type = DFNT_CHAR8;
        rank = 2;
        edge[1] = 10;
        dim_name[1] = "10";
        break;
      case SD_START_TIME:
        sds_name = "SD start time";
        data_type = DFNT_FLOAT64;
        rank = 1;
        break;
      case SRCA_START_TIME:
        sds_name = "SRCA start time";
        data_type = DFNT_FLOAT64;
        rank = 1;
        break;
      case BB_START_TIME:
        sds_name = "BB start time";
        data_type = DFNT_FLOAT64;
        rank = 1;
        break;
      case SV_START_TIME:
        sds_name = "SV start time";
        data_type = DFNT_FLOAT64;
        rank = 1;
        break;
      case EV_START_TIME:
        sds_name = "EV start time";
        data_type = DFNT_FLOAT64;
        rank = 1;
        break;
      case SRCA_CALIBRATION_MODE:
        sds_name = "SRCA calibration mode";
        data_type = DFNT_INT16;
        rank = 1;
        break;
      case PACKET_SCAN_COUNT:
        sds_name = "Packet scan count";
        data_type = DFNT_INT16;
        rank = 1;
        break;
      case CCSDS_APPLICATION_IDENTIFIERS:
        sds_name = "CCSDS Application Identifiers";
        data_type = DFNT_INT16;
        rank = 2;
        edge[1] = 3;
        dim_name[1] = "3";
        break;
      case PACKET_EXPEDITED_DATA_FLAG:
        sds_name = "Packet expedited data flag";
        data_type = DFNT_INT16;
        rank = 1;
        break;
      case MIRROR_SIDE:
        sds_name = "Mirror side";
        data_type = DFNT_INT16;
        rank = 1;
        break;
      case SCAN_QUALITY_ARRAY:
        sds_name = "Scan quality array";
        data_type = DFNT_INT32;
        rank = 2;
        edge[1] = 4;
        dim_name[1] = "4";
        break;
    }

    /*
     * Read the SDS data from the L1A file.
     */

    returnStatus = read_sds_rankn(in_sd_id, sds_name, data_type, rank,
                             start, edge, buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "read_sds_rankn",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    /*
     * Write the data to the OBC file.
     */

    returnStatus = write_sds_rankn(out_sd_id, sds_name, data_type, rank,
                                   edge, dim_name, buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "write_sds_rankn",
                  L1B_OBC_FILE, NULL, True);
      return returnStatus;
    }

  }

  return returnStatus;
}

PGSt_SMF_status Copy_EngVdata (int32 in_v_id,
                               int32 out_v_id)
/*
!C************************************************************************
!Description:  Copy L1A Engineering Vdata into L1B_OBCEng file.

!Input Parameters:
    int32    in_v_id      L1A file v_id
    int32    out_v_id     L1B_OBCEng file v_id

!Output Parameters:
    N/A

!Revision History:
 Revision 01.00 June 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char *location = "Copy_EngVdata";
  intn    hdf_return  = FAIL;
  int32   in_vdata_id = 0;
  int32   out_vdata_id= 0;
  int32   vdata       = 0; /*vdata index*/
  int32   field       = 0; /*field index*/
  int32   nfields     = 0;
  int32   nvdata      = 0;
  int32   nrecords    = 0;
  int32   interlace   = 0;
  int32   vd_size     = 0;
  int32   vdata_ref[MAX_NUM_VDATA];
  int32   f_type      = 0;
  int32   f_order     = 0;
  char    *f_name     = NULL; /* str pointer */
  char    fields[2048], vd_name[80];  /*max fields string less than 1024 bytes*/
  char    buffer[MAX_NUM_SCANS * MAX_VDATA_SIZE];

  /*Get nvdata in L1A file and the reference numbers*/
  nvdata = VSlone (in_v_id, NULL, 0);
  if (nvdata == FAIL)
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus, NULL, "VSlone(1st)",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  hdf_return = VSlone (in_v_id, vdata_ref, nvdata);
  if (hdf_return == FAIL)
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus, NULL, "VSlone(2nd)",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*Loop over all vdata*/
  for (vdata = 0; vdata < nvdata; vdata++)
  {
    /*attach a vdata to read*/
    in_vdata_id = VSattach (in_v_id, vdata_ref[vdata], "r");
    if (in_vdata_id == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSattach",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    /*inquire info*/
    hdf_return = VSinquire (in_vdata_id, &nrecords, &interlace,
                   fields, &vd_size, vd_name);
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSinquire",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }
    if (vd_size > MAX_VDATA_SIZE)
    {
      returnStatus = MODIS_F_NOK;
      L1BErrorMsg(location, returnStatus,
          "macro MAX_VDATA_SIZE is too small", NULL,
         FIRST_L1A_GRANULE,
         "Either there is something wrong with the file or there is a code bug.",
         True);
      return returnStatus;
    }

    if (strcmp(vd_name, "Discarded Packets") == 0) /* skip "Discarded Packets" */
    {
      hdf_return = VSdetach (in_vdata_id);
      if (hdf_return == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, NULL,
            "VSdetach(\"Discarded Packets\")",
                    FIRST_L1A_GRANULE, NULL, True);
        return returnStatus;
      }

      continue;
    }

    if (nrecords == 0) /*empty vdata: skip it*/
    {
      hdf_return = VSdetach (in_vdata_id);
      if (hdf_return == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, NULL, "VSdetach(empty vdata)",
                    FIRST_L1A_GRANULE, NULL, True);
        return returnStatus;
      }

      continue;
    }

    /*continue processing non_empty vdata*/
    nfields = VFnfields (in_vdata_id);
    if (nfields == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VFnfields",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    /*attach a vdata to write*/
    out_vdata_id = VSattach (out_v_id, -1, "w");
    if (out_vdata_id  == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSattach",
                  L1B_OBC_FILE, NULL, True);
      return returnStatus;
    }

    /*define fields*/
    for (field = 0; field < nfields; field++)
    {
      /*get info*/
      f_name = VFfieldname (in_vdata_id, field);
      if (f_name == NULL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, NULL, "VFfieldname",
                    FIRST_L1A_GRANULE, NULL, True);
        return returnStatus;
      }

      f_type = VFfieldtype (in_vdata_id, field);
      if (f_type == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, NULL, "VFfieldtype",
                    FIRST_L1A_GRANULE, NULL, True);
        return returnStatus;
      }

      f_order = VFfieldorder (in_vdata_id, field);
      if (f_order == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, NULL, "VFfieldorder",
                    FIRST_L1A_GRANULE, NULL, True);
        return returnStatus;
      }

      /*define field*/
      hdf_return = VSfdefine (out_vdata_id, f_name, f_type, f_order);
      if (hdf_return == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, NULL, "VSfdefine",
                    L1B_OBC_FILE, NULL, True);
        return returnStatus;
      }
    }

    /*read*/
    hdf_return = VSsetfields (in_vdata_id, fields);
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSsetfields",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    hdf_return = VSread (in_vdata_id, (void *)buffer, nrecords, interlace);
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_READ_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSread",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    /*write*/
    hdf_return = VSsetname (out_vdata_id, vd_name);
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSsetname",
                  L1B_OBC_FILE, NULL, True);
      return returnStatus;
    }

    hdf_return = VSsetfields (out_vdata_id, fields);
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSsetfields",
                  L1B_OBC_FILE, NULL, True);
      return returnStatus;
    }

    hdf_return = VSwrite (out_vdata_id, (void *)buffer, nrecords, interlace);
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_WRITE_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSwrite",
                  L1B_OBC_FILE, NULL, True);
      return returnStatus;
    }

    /*done with this vdata*/
    hdf_return = VSdetach (in_vdata_id);
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_READ_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSdetach",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    hdf_return = VSdetach (out_vdata_id);
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_READ_ERROR;
      L1BErrorMsg(location, returnStatus, NULL, "VSdetach",
                  L1B_OBC_FILE, NULL, True);
      return returnStatus;
    }

  } /*End of loop over all vdata*/

  return returnStatus;
}

PGSt_SMF_status Check_For_Moon_in_SV_KOB
                      (int32              num_scans,
                       common_QA_tables_t *common_QA_tables,
                       QA_Data_t          *QA,
                       Moon_arrays_t      *moon_arrays)
/*
!C**********************************************************************
!Description:
   For each scan of the middle granule, determine if the center of the
   moon lies within the space-view (SV) keep-out box (KOB).  This sets
   values to be used later in determining the DN offset and for writing
   one of the "bit QA flags" in the scan metadata.

!Input Parameters:
  int32  num_scans                        Number of scans in the middle L1A
                                             granule.
  common_QA_tables_t  *common_QA_tables
                                Contains: float32  moon_offset_limits
                                          Defines the offsets from the SV
                                          port center in 1km units.
                                          Dimensions: [NUM_BANDS][4]
!Output Parameters:
  QA_Data_t           *QA                 Used to set the appropriate
                                             "Bit QA Flags" value in L1B_Setup.
  Moon_arrays_t       *moon_arrays        keep the moon vector values and
                                             flags denoting if center of
                                             moon is within the keep-out
                                             box for every scan and band.
!Revision History:
 Revision 01.00 April 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
 According to James Kuyper (email, Jun 3, 1999) the magnitude of the
 moon vector should be 1 if successfully computed or 0 if not. There also
 could be "FloatInf" values (fill values) in some simulated data sets.
 The logic of this function is that a value of UNDETERMINED is assigned
 for either a zero magnitude moon vector or values of the components
 out of the range -1 through 1, inclusive.

!END********************************************************************
*/
{
  PGSt_SMF_status   returnStatus = MODIS_S_OK;
  PGSt_integer      Version       = 1;
  int32             geo_sd_id;
  char              fname[512];
  int32             B;
  int32             S;
  int32             num_scans_geo;
  float32           distance_square;
  float32           moon_scan_offset;
  float32           moon_track_offset;
  intn              hdf_return;
  int16             band_is_RSB[NUM_BANDS] = {1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 0, 0, 0,
                                              0, 0, 0, 1, 0, 0,
                                              0, 0,0, 0, 0, 0, 0, 0};

  float32           x, y, z;
  char *location = "Check_For_Moon_in_SV_KOB";

  if (PGS_PC_GetReference(GEOLOCATION_FILE, &Version, fname) != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus,
                "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  geo_sd_id = SDstart(fname, DFACC_RDONLY);
  if (geo_sd_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus,
                "Could not open file for SD read access.",
                "SDstart", GEOLOCATION_FILE,
                "The file may be missing, corrupted or not an HDF-4 file.",
                True);
    return returnStatus;
  }

  returnStatus = read_attribute(geo_sd_id,
                                GEO_NUM_SCANS_NAME,
                                DFNT_INT32,
                                &num_scans_geo);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "read_attribute", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  if (num_scans_geo != num_scans)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
      "Number of scans in GEO file does not match the current L1A file.",
      NULL, GEOLOCATION_FILE, "The PCF may have the wrong names.", True);
    return returnStatus;
  }

  returnStatus = read_sds_rank2(geo_sd_id,
                                MOON_VECTOR_SDS_NAME,
                                num_scans,
                                VEC_DIM,
                                moon_arrays->moon_vector);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "read_sds_rank2", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  hdf_return = SDend(geo_sd_id);
  if (hdf_return == FAIL)
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus, NULL, "SDend", GEOLOCATION_FILE,
                "Memory or the disk file must have become corrupted.", True);
    return returnStatus;
  }

  for (S = 0; S < num_scans; S++)
  {
    QA->QA_refl.moon_in_SV_KOB_RSB[S] = MOON_OUTSIDE_SV_KOB;
    QA->QA_emiss.moon_in_SV_KOB_TEB[S] = MOON_OUTSIDE_SV_KOB;

    x = moon_arrays->moon_vector[S][0];
    y = moon_arrays->moon_vector[S][1];
    z = moon_arrays->moon_vector[S][2];

      /*
       * Check that vector components are within range.
       * If so, then cannot compute.
       */

    if (x > 1.01 || x < -1.01 || y > 1.01 ||
           y < -1.01 || z > 1.01 || z < -1.01)
    {
      for (B = 0; B < NUM_BANDS; B++)
        moon_arrays->moon_in_SV_KOB[S][B] = UNDETERMINED_MOON_IN_SV_KOB;
      continue;
    }

    distance_square = x * x + y * y + z * z;

       /* Check that magnitude is not zero.  IF so, then cannot compute */

    if (distance_square < TOLERANCE)
    {
      for (B = 0; B < NUM_BANDS; B++)
        moon_arrays->moon_in_SV_KOB[S][B] = UNDETERMINED_MOON_IN_SV_KOB;
      continue;
    }
    moon_track_offset = - asin((double)x / sqrt((double)distance_square)) /
        PIXEL_WIDTH_RADIANS;
    moon_scan_offset = (atan2((double)z, (double)y) - SV_CNTR_ANG_FROM_Y) /
        PIXEL_WIDTH_RADIANS;
    for (B = 0; B < NUM_BANDS; B++)
    {
      if (moon_track_offset <=
             common_QA_tables->moon_offset_limits[B][TRK_UPPER] &&
          moon_track_offset >=
             common_QA_tables->moon_offset_limits[B][TRK_LOWER] &&
          moon_scan_offset  <=
             common_QA_tables->moon_offset_limits[B][SCN_UPPER] &&
          moon_scan_offset  >=
             common_QA_tables->moon_offset_limits[B][SCN_LOWER])
      {
        moon_arrays->moon_in_SV_KOB[S][B] = MOON_INSIDE_SV_KOB;
        if (band_is_RSB[B] == 1)
          QA->QA_refl.moon_in_SV_KOB_RSB[S] = MOON_INSIDE_SV_KOB;
        else
          QA->QA_emiss.moon_in_SV_KOB_TEB[S] = MOON_INSIDE_SV_KOB;
      }
      else
        moon_arrays->moon_in_SV_KOB[S][B] = MOON_OUTSIDE_SV_KOB;
    }/* end loop through bands */
  }/* end loop throught scans */

  return (returnStatus);
}

PGSt_SMF_status Fill_250m_DN_OBC_Avg(int32             sd_id,
                                     int32             num_scans,
                                     int16             *MirrorSide,
                                     refl_tables_t     *refl_tables,
                                     Moon_arrays_t     *moon_arrays,
                                     boolean           Ecal_On[][NUM_BANDS],
                                     DN_OBC_Avg_t      *DN_OBC_Avg,
                                     QA_Refl_t         *QA_refl)
/*
!C**********************************************************************
!Description:
   Compute the DN OBC averages for the 250m resolution bands.  Store in
   DN_OBC_Avg for writing to the OBC file and for use in the reflective
   algorithm.

!Input Parameters:
   int32           sd_id           SD file access ID for the L1A middle
                                      granule.
   int32           num_scans       Number of scans in the middle granule.
   int16           *MirrorSide     Mirror side indexes per scan
   refl_tables_t   *refl_tables    Structure which contains all reflective
                                      LUTs
   Moon_arrays_t   *moon_arrays    Structure containing the moon_in_SV_KOB
                                      array.
   boolean  Ecal_On[][NUM_BANDS]   Ecal on information for each band for
                                      each scan.
   Preprocess_Refl_t *PP_Refl      Structure containing DN saturation value.

!Output Parameters:
   DN_OBC_Avg_t   DN_OBC_Avg       Structure containing arrays that will
                                   ultimately be written to the OBC file.
   QA_Refl_t      QA_refl          Depending on processing values, some
                                   elements of these may be tripped.

!Revision History:


 Revision  01.10 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.05, March 26, 2001, Razor Issue #159
 Changed DN_upper_valid_limit to SATURATED_DN - 1.
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.04, Dec 7, 2000, Razor issue 146
 Removed DN_sat_prime.
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 Oct. 22, 1999
 Added checking Ecal on
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.02 August 28, 1999
 Added checking if the data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 August, 1999
 Implement changes to meet the requirements of the new ADC algorithm. See SDF.
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 April, 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;

  typedef int16 SV_250m_t[NUM_250M_BANDS][MAX_250M_OBC_FRAME_DIM];
  SV_250m_t *SV_250m;
  SV_250m = (SV_250m_t*) malloc(MAX_250M_TRACK_DIM * sizeof(SV_250m_t));
  typedef int16 BB_250m_t[NUM_250M_BANDS][MAX_250M_OBC_FRAME_DIM];
  BB_250m_t *BB_250m;
  BB_250m = (BB_250m_t*) malloc(MAX_250M_TRACK_DIM * sizeof(BB_250m_t));

  int32 S;
  int32 D;
  int32 T;
  int32 B;
  int32 subsamp;
  int32 increment;
  int16   rejects[MAX_1KM_OBC_FRAME_DIM];
  int16   DN_upper_valid_limit;
  int16   MS;
  float32 mean;
  float32 sdev;
  int32   track_dim;
  int32   band_dim;
  int32   frame_dim;
  int32   F;
  char *location = "Fill_250m_DN_OBC_Avg";

  track_dim = num_scans * DETECT_PER_BAND_AT_RES[INDEX_250M];
  band_dim  = NUM_250M_BANDS;
  frame_dim = SV_1km_FRAMES * BAND_RATIO_AT_RES[INDEX_250M];
  increment = BAND_RATIO_AT_RES[INDEX_250M];

  /*
   * Initialize arrays (in case we skip scans due to Mirror Side).
   */

  for (T = 0; T < MAX_250M_TRACK_DIM; T++) {
    for (B = 0; B < NUM_250M_BANDS; B++) {
      for (subsamp = 0; subsamp < NUM_250M_SUBSAMP; subsamp++) {
        DN_OBC_Avg->DN_obc_250m_avg[T][B][subsamp] = -1;
        DN_OBC_Avg->DN_obc_250m_var[T][B][subsamp] = 0;
        DN_OBC_Avg->DN_obc_250m_outlier_mask[T][B][subsamp][0] = 0;
        DN_OBC_Avg->DN_obc_250m_outlier_mask[T][B][subsamp][1] = 0;
      }
    }
  }

  /*
   * Read SV DNs from L1A middle granule and check for valid range.
   */

  returnStatus = read_sds_rank3(sd_id,
                                L1A_TARGET_SDS_NAME[SV_INDEX][INDEX_250M],
                                track_dim,
                                band_dim,
                                frame_dim,
                                SV_250m);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank3",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  returnStatus = Check_Valid_Range("SV_250m",
                                   DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   (track_dim*band_dim*frame_dim),
                                   (void *) SV_250m);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*
   * Read BB DNs from L1A middle granule and check for valid range.
   */

  returnStatus = read_sds_rank3(sd_id,
                                L1A_TARGET_SDS_NAME[BB_INDEX][INDEX_250M],
                                track_dim,
                                band_dim,
                                frame_dim,
                                BB_250m);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank3",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  returnStatus = Check_Valid_Range("BB_250m", DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   (track_dim*band_dim*frame_dim),
                                   (void *) BB_250m);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  for (S = 0; S < num_scans; S++)
  {
    MS = MirrorSide[S];
    if (MS != 0 && MS != 1) continue;              /* Skip this scan */

    for (D = 0; D < DETECTORS_PER_250M_BAND; D++)
    {
      T = S * DETECTORS_PER_250M_BAND + D;
      for (B = 0; B < NUM_250M_BANDS; B++)
      {
        for (subsamp = 0; subsamp < BAND_RATIO_AT_RES[INDEX_250M]; subsamp++)
        {
          for (F = 0; F < MAX_1KM_OBC_FRAME_DIM; F++)
            rejects[F] = 0;
          DN_upper_valid_limit = SATURATED_DN - 1;

          if (Ecal_On[S][B] == True)
          {
            mean = -1;
            sdev = 0;
          }
          else if (moon_arrays->moon_in_SV_KOB[S][B] == MOON_INSIDE_SV_KOB)
          {
            returnStatus =
              Get_DN_Avg_SDev_Rejects_LowN(
                               refl_tables->RSB_SV_DN_moon_include_frames,
                               subsamp,
                               SV_1km_FRAMES,
                               increment,
                               SV_250m[T][B],
                               DN_upper_valid_limit,
                               0,
                               &mean,
                               &sdev,
                               &rejects[0]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                "Get_DN_Avg_SDev_Rejects_LowN() in Fill_250m_DN_OBC_Avg()");
          }
          else
          {
            returnStatus =
                Get_DN_Avg_SDev_Rejects(
                     refl_tables->DN_obc_avg_first_frame_to_use
                          * BAND_RATIO_AT_RES[INDEX_250M] + subsamp,
                     refl_tables->DN_obc_avg_number_of_frames_to_use,
                     increment,
                     SV_250m[T][B],
                     DN_upper_valid_limit,
                     0,
                     &mean,
                     &sdev,
                     &rejects[refl_tables->DN_obc_avg_first_frame_to_use]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                        "Get_DN_Avg_SDev_Rejects() in Fill_250m_DN_OBC_Avg()");
          }
          if (mean < 0)
          {
            QA_refl->all_SV_DN_bad[S] = 1;
            returnStatus =
                Get_DN_Avg_SDev_Rejects(
                    refl_tables->DN_obc_avg_first_frame_to_use
                        * BAND_RATIO_AT_RES[INDEX_250M] + subsamp,
                    refl_tables->DN_obc_avg_number_of_frames_to_use,
                    increment,
                    BB_250m[T][B],
                    DN_upper_valid_limit,
                    0,
                    &mean,
                    &sdev,
                    &rejects[refl_tables->DN_obc_avg_first_frame_to_use]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                  "Get_DN_Avg_SDev_Rejects() in Fill_250m_DN_OBC_Avg()");
          }

          if (mean < 0)
            QA_refl->all_BB_DN_bad[S] = 1;

          DN_OBC_Avg->DN_obc_250m_avg[T][B][subsamp] = mean;
          DN_OBC_Avg->DN_obc_250m_var[T][B][subsamp] = sdev * sdev;
          DN_OBC_Avg->DN_obc_250m_outlier_mask[T][B][subsamp][0] = 0;
          DN_OBC_Avg->DN_obc_250m_outlier_mask[T][B][subsamp][1] = 0;
          returnStatus =
              Pack_Rejects_In_Outlier_Mask(SV_1km_FRAMES, rejects,
                  DN_OBC_Avg->DN_obc_250m_outlier_mask[T][B][subsamp]);
        }/* end loop through subsamples */
      }/* end loop through bands */
    }/* end loop through detector */
  }/* end loop through scans */

  free(SV_250m);
  free(BB_250m);

  return (returnStatus);
}

PGSt_SMF_status Fill_500m_DN_OBC_Avg(int32             sd_id,
                                     int32             num_scans,
                                     int16             *MirrorSide,
                                     refl_tables_t     *refl_tables,
                                     Moon_arrays_t     *moon_arrays,
                                     boolean           Ecal_On[][NUM_BANDS],
                                     DN_OBC_Avg_t      *DN_OBC_Avg,
                                     QA_Refl_t         *QA_refl)
/*
!C**********************************************************************
!Description:
   Compute the DN OBC averages for the 500m resolution bands.  Store in
   DN_OBC_Avg for writing to the OBC file and use in the reflective algorithm.

!Input Parameters:
   int32           sd_id           SD file access ID for the L1A middle
                                      granule.
   int32           num_scans       Number of scans in the middle granule.
   int16           *MirrorSide     Mirror side indexes per scan
   refl_tables_t   *refl_tables    Structure which contains all reflective
                                      LUTs
   Moon_arrays_t   *moon_arrays    Structure containing the moon_in_SV_KOB
                                      array.
   boolean  Ecal_On[][NUM_BANDS]   Ecal on information for each band for
                                      each scan.
   Preprocess_Refl_t *PP_Refl      Structure containing DN saturation value.

!Output Parameters:
   DN_OBC_Avg_t   DN_OBC_Avg       Structure containing arrays that will
                                   ultimately be written to the OBC file.
   QA_Refl_t      QA_refl          Depending on processing values, some
                                   elements of these may be tripped.

!Revision History:

 Revision 01.10 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.06, March 26, 2001, Razor Issue 159
 Changed DN_upper_valid_limit to SATURATED_DN - 1.
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.05, Dec 7, 2000, Razor issue 146
 Removed DN_sat_prime.
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.04 Oct. 22, 1999
 Added checking for ecal on
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.03 Sept. 2, 1999
 Implemented changes to meet the requirements of the new SWIR algorithm.
    See SDF.
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.02 August 28, 1999
 Added checking if the data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 August, 1999
 Implement changes to meet the requirements of the new ADC algorithm. See SDF.
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 April 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  
  typedef int16 SV_500m_t[NUM_500M_BANDS][MAX_500M_OBC_FRAME_DIM];
  SV_500m_t *SV_500m;
  SV_500m = (SV_500m_t*) malloc(MAX_500M_TRACK_DIM * sizeof(SV_500m_t));
  typedef int16 BB_500m_t[NUM_500M_BANDS][MAX_500M_OBC_FRAME_DIM];
  BB_500m_t *BB_500m;
  BB_500m = (BB_500m_t*) malloc(MAX_500M_TRACK_DIM * sizeof(BB_500m_t));
  
  int32 S;
  int32 D;
  int32 T;
  int32 B;
  int32 B_38;
  int32 subsamp;
  int32 increment;
  int16   rejects[MAX_1KM_OBC_FRAME_DIM];
  int16   DN_upper_valid_limit;
  int16   MS;
  float32 mean;
  float32 sdev;
  int32   track_dim;
  int32   band_dim;
  int32   frame_dim;
  int32   F;
  char *location = "Fill_500m_DN_OBC_Avg";

  track_dim = num_scans * DETECT_PER_BAND_AT_RES[INDEX_500M];
  band_dim  = NUM_500M_BANDS;
  frame_dim = SV_1km_FRAMES * BAND_RATIO_AT_RES[INDEX_500M];
  increment = BAND_RATIO_AT_RES[INDEX_500M];

  /*
   * Initialize arrays (in case we skip scans to due Mirror Side).
   */

  for (T = 0; T < MAX_500M_TRACK_DIM; T++) {
    for (B = 0; B < NUM_500M_BANDS; B++) {
      for (subsamp = 0; subsamp < NUM_500M_SUBSAMP; subsamp++) {
        DN_OBC_Avg->DN_obc_500m_avg[T][B][subsamp] = -1;
        DN_OBC_Avg->DN_obc_500m_var[T][B][subsamp] = 0;
        DN_OBC_Avg->DN_obc_500m_outlier_mask[T][B][subsamp][0] = 0;
        DN_OBC_Avg->DN_obc_500m_outlier_mask[T][B][subsamp][1] = 0;
      }
    }
  }

  /*
   * Read SV DNs from L1A middle granule and check for valid range.
   */

  returnStatus = read_sds_rank3(sd_id,
                                L1A_TARGET_SDS_NAME[SV_INDEX][INDEX_500M],
                                track_dim,
                                band_dim,
                                frame_dim,
                                SV_500m);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank3",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  returnStatus = Check_Valid_Range("SV_500m",
                                   DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   (track_dim*band_dim*frame_dim),
                                   (void *) SV_500m);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*
   * Read BB DNs from L1A middle granule and check for valid range.
   */

  returnStatus = read_sds_rank3(sd_id,
                                L1A_TARGET_SDS_NAME[BB_INDEX][INDEX_500M],
                                track_dim,
                                band_dim,
                                frame_dim,
                                BB_500m);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank3",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  returnStatus = Check_Valid_Range("BB_500m",
                                   DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   (track_dim*band_dim*frame_dim),
                                   (void *) BB_500m);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  for (S = 0; S < num_scans; S++)
  {
    MS = MirrorSide[S];
    if (MS != 0 && MS != 1) continue;              /* Skip this scan */

    for (D = 0; D < DETECTORS_PER_500M_BAND; D++)
    {
      T = S * DETECTORS_PER_500M_BAND + D;
      for (B = 0; B < NUM_500M_BANDS; B++)
      {
        B_38   = B + NUM_250M_BANDS;
        for (subsamp = 0;
               subsamp < BAND_RATIO_AT_RES[INDEX_500M]; subsamp++)
        {
          for (F = 0; F < MAX_1KM_OBC_FRAME_DIM; F++)
            rejects[F] = 0;
          DN_upper_valid_limit = SATURATED_DN - 1;

          if (Ecal_On[S][B_38] == True)
          {
            mean = -1;
            sdev = 0;
          }
          else if (moon_arrays->moon_in_SV_KOB[S][B_38] ==
                     MOON_INSIDE_SV_KOB)
          {
            returnStatus = Get_DN_Avg_SDev_Rejects_LowN(
                               refl_tables->RSB_SV_DN_moon_include_frames,
                               subsamp,
                               SV_1km_FRAMES,
                               increment,
                               SV_500m[T][B],
                               DN_upper_valid_limit,
                               0,
                               &mean,
                               &sdev,
                               &rejects[0]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                "Get_DN_Avg_SDev_Rejects_LowN() in Fill_500m_DN_OBC_Avg()");
          }
          else
          {
            returnStatus =
                Get_DN_Avg_SDev_Rejects(
                  refl_tables->DN_obc_avg_first_frame_to_use
                       * BAND_RATIO_AT_RES[INDEX_500M] + subsamp,
                  refl_tables->DN_obc_avg_number_of_frames_to_use,
                  increment,
                  SV_500m[T][B],
                  DN_upper_valid_limit,
                  0,
                  &mean,
                  &sdev,
                  &rejects[refl_tables->DN_obc_avg_first_frame_to_use]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                  "Get_DN_Avg_SDev_Rejects() in Fill_500m_DN_OBC_Avg()");
          }
          if (mean < 0)
            QA_refl->all_SV_DN_bad[S] = 1;
          if (mean < 0 && B < 2)  /* only use BB DNs for bands 3 and 4 */
          {
            returnStatus =
                Get_DN_Avg_SDev_Rejects(
                  refl_tables->DN_obc_avg_first_frame_to_use
                       * BAND_RATIO_AT_RES[INDEX_500M] + subsamp,
                  refl_tables->DN_obc_avg_number_of_frames_to_use,
                  increment,
                  BB_500m[T][B],
                  DN_upper_valid_limit,
                  0,
                  &mean,
                  &sdev,
                  &rejects[refl_tables->DN_obc_avg_first_frame_to_use]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                  "Get_DN_Avg_SDev_Rejects() in Fill_500m_DN_OBC_Avg()");
            if (mean < 0)
              QA_refl->all_BB_DN_bad[S] = 1;
          }

          DN_OBC_Avg->DN_obc_500m_avg[T][B][subsamp] = mean;
          DN_OBC_Avg->DN_obc_500m_var[T][B][subsamp] = sdev * sdev;
          DN_OBC_Avg->DN_obc_500m_outlier_mask[T][B][subsamp][0] = 0;
          DN_OBC_Avg->DN_obc_500m_outlier_mask[T][B][subsamp][1] = 0;
          returnStatus =
              Pack_Rejects_In_Outlier_Mask(SV_1km_FRAMES,
                  rejects,
                  DN_OBC_Avg->DN_obc_500m_outlier_mask[T][B][subsamp]);
        }/* end loop through subsamples */
      }/* end loop through bands */
    }/* end loop through detector */
  }/* end loop through scans */

  free(SV_500m);
  free(BB_500m);

  return (returnStatus);
}

PGSt_SMF_status Fill_1km_day_DN_OBC_Avg
                   (int32             sd_id,
                    int32             num_scans,
                    int16             *MirrorSide,
                    refl_tables_t     *refl_tables,
                    Moon_arrays_t     *moon_arrays,
                    boolean           Ecal_On[][NUM_BANDS],
                    DN_OBC_Avg_t      *DN_OBC_Avg,
                    QA_Refl_t         *QA_refl)
/*
!C**********************************************************************
!Description:
   Compute the DN OBC averages for the 1km day resolution bands.  Store in
   DN_OBC_Avg for writing to the OBC file for use in the reflective algorithm.

!Input Parameters:
   int32           sd_id           SD file access ID for the L1A middle
                                      granule.
   int32           num_scans       Number of scans in the middle granule.
   int16           *MirrorSide     Mirror side indexes per scan
   refl_tables_t   *refl_tables    Structure which contains all reflective
                                      LUTs
   Moon_arrays_t   *moon_arrays    Structure containing the moon_in_SV_KOB
                                      array.
   boolean  Ecal_On[][NUM_BANDS]   Ecal on information for each band for
                                      each scan.
   Preprocess_Refl_t *PP_Refl      Structure containing DN saturation value.

!Output Parameters:
   DN_OBC_Avg_t   DN_OBC_Avg       Structure containing arrays that will
                                   ultimately be written to the OBC file.
   QA_Refl_t            QA_refl    Depending on processing values, some
                                   elements of these may be tripped.

!Revision History:

 Revision 01.10 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.05, March 26, 2001, Razor Issue #159
 Changed DN_upper_valid_limit to SATURATED_DN - 1.
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.04, Dec 7, 2000, Razor issue 146
 Removed DN_sat_prime.
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 Oct. 22, 1999
 Added checking for Ecal on
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.02 August 28, 1999
 Added checking if the data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 August, 1999
 Implement changes to meet the requirements of the new ADC algorithm. See SDF.
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 April 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;

  typedef int16 SV_1km_day_t[NUM_1000M_DAY_BANDS][MAX_1KM_OBC_FRAME_DIM];
  SV_1km_day_t *SV_1km_day;
  SV_1km_day = (SV_1km_day_t*) malloc(MAX_1KM_TRACK_DIM * sizeof(SV_1km_day_t));
  typedef int16 BB_1km_day_t[NUM_1000M_DAY_BANDS][MAX_1KM_OBC_FRAME_DIM];
  BB_1km_day_t *BB_1km_day;
  BB_1km_day = (BB_1km_day_t*) malloc(MAX_1KM_TRACK_DIM * sizeof(BB_1km_day_t));
    
  int32 S;
  int32 D;
  int32 T;
  int32 B;
  int32 B_38;
  int32 subsamp;
  int32 increment;
  int16   rejects[MAX_1KM_OBC_FRAME_DIM];
  int16   DN_upper_valid_limit;
  int16   MS;
  float32 mean;
  float32 sdev;
  int32   track_dim;
  int32   band_dim;
  int32   frame_dim;
  int32   F;
  char *location = "Fill_1km_day_DN_OBC_Avg";

  track_dim = num_scans * DETECTORS_PER_1KM_BAND;
  band_dim  = NUM_1000M_DAY_BANDS;
  frame_dim = SV_1km_FRAMES * BAND_RATIO_AT_RES[INDEX_1000M_DAY];
  increment = BAND_RATIO_AT_RES[INDEX_1000M_DAY];

  /*
   * Initialize arrays (in case we skip scans to due Mirror Side).
   */

  for (T = 0; T < MAX_1KM_TRACK_DIM; T++) {
    for (B = 0; B < NUM_1000M_DAY_BANDS; B++) {
      for (subsamp = 0; subsamp < NUM_1KM_SUBSAMP; subsamp++) {
        DN_OBC_Avg->DN_obc_1km_day_avg[T][B][subsamp] = -1;
        DN_OBC_Avg->DN_obc_1km_day_var[T][B][subsamp] = 0;
        DN_OBC_Avg->DN_obc_1km_day_outlier_mask[T][B][subsamp][0] = 0;
        DN_OBC_Avg->DN_obc_1km_day_outlier_mask[T][B][subsamp][1] = 0;
      }
    }
  }

  /*
   * Read SV DNs from L1A middle granule and check for valid range.
   */

  returnStatus = read_sds_rank3(sd_id,
                                L1A_TARGET_SDS_NAME[SV_INDEX][INDEX_1000M_DAY],
                                track_dim,
                                band_dim,
                                frame_dim,
                                SV_1km_day);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank3",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  returnStatus = Check_Valid_Range("SV_1km_day", DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   (track_dim*band_dim*frame_dim),
                                   (void *) SV_1km_day);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*
   * Read BB DNs from L1A middle granule and check for valid range.
   */

  returnStatus = read_sds_rank3(sd_id,
                                L1A_TARGET_SDS_NAME[BB_INDEX][INDEX_1000M_DAY], track_dim,
                                band_dim,
                                frame_dim,
                                BB_1km_day);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank3",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  returnStatus = Check_Valid_Range("BB_1km_day",
                                   DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   (track_dim*band_dim*frame_dim),
                                   (void *) BB_1km_day);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  for (S = 0; S < num_scans; S++)
  {
    MS = MirrorSide[S];
    if (MS != 0 && MS != 1) continue;              /* Skip this scan */

    for (D = 0; D < DETECTORS_PER_1KM_BAND; D++)
    {
      T = S * DETECTORS_PER_1KM_BAND + D;
      for (B = 0; B < NUM_1000M_DAY_BANDS; B++)
      {
        B_38   = B + NUM_250M_BANDS + NUM_500M_BANDS;

        for (subsamp = 0;
               subsamp < BAND_RATIO_AT_RES[INDEX_1000M_DAY]; subsamp++)
        {
          for (F = 0; F < MAX_1KM_OBC_FRAME_DIM; F++)
            rejects[F] = 0;
          DN_upper_valid_limit = SATURATED_DN - 1;

          if (Ecal_On[S][B_38] == True)
          {
            mean = -1;
            sdev = 0;
          }
          else if (moon_arrays->moon_in_SV_KOB[S][B_38] ==
                     MOON_INSIDE_SV_KOB)
          {
            returnStatus = Get_DN_Avg_SDev_Rejects_LowN(
                               refl_tables->RSB_SV_DN_moon_include_frames,
                               subsamp,
                               SV_1km_FRAMES,
                               increment,
                               SV_1km_day[T][B],
                               DN_upper_valid_limit,
                               0,
                               &mean,
                               &sdev,
                               &rejects[0]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                "Get_DN_Avg_SDev_Rejects_LowN() in Fill_1km_day_DN_OBC_Avg()");
          }
          else
          {
            returnStatus =
               Get_DN_Avg_SDev_Rejects(
                   refl_tables->DN_obc_avg_first_frame_to_use
                      + subsamp,
                   refl_tables->DN_obc_avg_number_of_frames_to_use,
                   increment,
                   SV_1km_day[T][B],
                   DN_upper_valid_limit,
                   0,
                   &mean,
                   &sdev,
                   &rejects[refl_tables->DN_obc_avg_first_frame_to_use]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                  "Get_DN_Avg_SDev_Rejects() in Fill_1km_day_DN_OBC_Avg()");
          }
          if (mean < 0)
          {
            QA_refl->all_SV_DN_bad[S] = 1;
            returnStatus =
               Get_DN_Avg_SDev_Rejects(
                     refl_tables->DN_obc_avg_first_frame_to_use + subsamp,
                     refl_tables->DN_obc_avg_number_of_frames_to_use,
                     increment,
                     BB_1km_day[T][B],
                     DN_upper_valid_limit,
                     0,
                     &mean,
                     &sdev,
                     &rejects[refl_tables->DN_obc_avg_first_frame_to_use]);
            if (returnStatus != MODIS_S_OK)
              SMF_ERROR(returnStatus,
                  "Get_DN_Avg_SDev_Rejects() in Fill_1km_day_DN_OBC_Avg()");
          }

          if (mean < 0)
            QA_refl->all_BB_DN_bad[S] = 1;

          DN_OBC_Avg->DN_obc_1km_day_avg[T][B][subsamp] = mean;
          DN_OBC_Avg->DN_obc_1km_day_var[T][B][subsamp] = sdev * sdev;
          DN_OBC_Avg->DN_obc_1km_day_outlier_mask[T][B][subsamp][0] = 0;
          DN_OBC_Avg->DN_obc_1km_day_outlier_mask[T][B][subsamp][1] = 0;
          returnStatus =
              Pack_Rejects_In_Outlier_Mask(SV_1km_FRAMES,
                  rejects,
                  DN_OBC_Avg->DN_obc_1km_day_outlier_mask[T][B][subsamp]);
        }/* end loop through subsamples */
      }/* end loop through bands */
    }/* end loop through detector */
  }/* end loop through scans */

  free(SV_1km_day);
  free(BB_1km_day);

  return (returnStatus);
}

PGSt_SMF_status Fill_Band_26_DN_OBC_Avg
                     (int32             sd_id,
                      int32             num_scans,
                      int16             *MirrorSide,
                      refl_tables_t     *refl_tables,
                      Moon_arrays_t     *moon_arrays,
                      boolean           Ecal_On[][NUM_BANDS],
                      DN_OBC_Avg_t      *DN_OBC_Avg,
                      QA_Refl_t         *QA_refl)
/*
!C**********************************************************************
!Description:
   Compute the DN OBC averages for band 26. Store in DN_OBC_Avg structure
   for writing to the OBC file and for use in the reflective algorithm.

!Input Parameters:
   int32           sd_id           SD file access ID for the L1A middle
                                     granule.
   int32           num_scans       Number of scans in the middle granule.
   int16           *MirrorSide     Mirror side indexes per scan
   refl_tables_t   *refl_tables    Structure which contains all reflective
                                     LUTs
   Moon_arrays_t   *moon_arrays    Structure containing the moon_in_SV_KOB
                                     array.
   boolean  Ecal_On[][NUM_BANDS]   Ecal on information for each band for
                                     each scan.
   Preprocess_Refl_t *PP_Refl      Structure containing DN saturation value.

!Output Parameters:
   DN_OBC_Avg_t   DN_OBC_Avg       Structure containing arrays that will
                                   ultimately be written to the OBC file.
   QA_Refl_t      QA_refl          Depending on processing values, some
                                   elements of these may be tripped.

!Revision History:

 Revision 01.10 March 25, 2002   Razor Issue #178
 Strip out ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.06, March 26, 2001, Razor Issue #159
 Changed DN_upper_valid_limit to SATURATED_DN - 1.
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.05, Dec 7, 2000, Razor issue 146
 Removed DN_sat_prime.
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.04 Oct. 22, 1999
 Added checking for Ecal on
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.03 Sept. 2, 1999
 Implemented changes to meet the requirements of the new SWIR algorithm. See SDF.
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.02 August 28, 1999
 Added checking if the data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 August, 1999
 Implement changes to meet the requirements of the new ADC algorithm. See SDF.
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 April 1999
 Initial development
 Jim Rogers(rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  
  typedef int16 SV_1km_night_t[NUM_1000M_NIGHT_BANDS][MAX_1KM_OBC_FRAME_DIM];
  SV_1km_night_t *SV_1km_night;
  SV_1km_night = (SV_1km_night_t*) malloc(MAX_1KM_TRACK_DIM * sizeof(SV_1km_night_t));

  int32 S;
  int32 D;
  int32 T;
  int32 B;
  int32 B_38;
  int32 subsamp;
  int32 increment;
  int16   rejects[MAX_1KM_OBC_FRAME_DIM];
  int16   DN_upper_valid_limit;
  int16   MS;
  float32 mean;
  float32 sdev;
  int32   track_dim;
  int32   band_dim;
  int32   frame_dim;
  int32   F;
  char *location = "Fill_Band_26_DN_OBC_Avg";

  track_dim = num_scans * DETECTORS_PER_1KM_BAND;
  band_dim  = NUM_1000M_NIGHT_BANDS;
  frame_dim = SV_1km_FRAMES * BAND_RATIO_AT_RES[INDEX_1000M_NIGHT];
  increment = BAND_RATIO_AT_RES[INDEX_1000M_NIGHT];

  /*
   * Initialize arrays (in case we skip scans to due Mirror Side).
   */

  B = BAND26;  /* Index in the 1km night band group */

  for (T = 0; T < MAX_1KM_TRACK_DIM; T++) {
    for (subsamp = 0; subsamp < NUM_1KM_SUBSAMP; subsamp++) {
      DN_OBC_Avg->DN_obc_1km_night_avg[T][B][subsamp] = -1;
      DN_OBC_Avg->DN_obc_1km_night_var[T][B][subsamp] = 0;
      DN_OBC_Avg->DN_obc_1km_night_outlier_mask[T][B][subsamp][0] = 0;
      DN_OBC_Avg->DN_obc_1km_night_outlier_mask[T][B][subsamp][1] = 0;
    }
  }

  /*
   * Read SV DNs from L1A middle granule and check for valid range.
   */

  returnStatus = read_sds_rank3(sd_id,
                          L1A_TARGET_SDS_NAME[SV_INDEX][INDEX_1000M_NIGHT],
                          track_dim,
                          band_dim,
                          frame_dim,
                          SV_1km_night);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank3",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  returnStatus = Check_Valid_Range("SV_1km_night",
                                   DFNT_INT16,
                                   L1A_DN_SDS_LB,
                                   L1A_DN_SDS_UB,
                                   L1A_DN_SDS_FV,
                                   (track_dim*band_dim*frame_dim),
                                   (void *) SV_1km_night);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  B = BAND26;
  B_38   = MODIS_BAND26_INDEX;

  for (S = 0; S < num_scans; S++)
  {
    MS = MirrorSide[S];
    if (MS != 0 && MS != 1) continue;              /* Skip this scan */

    for (D = 0; D < DETECTORS_PER_1KM_BAND; D++)
    {
      T = S * DETECTORS_PER_1KM_BAND + D;

      for (subsamp = 0;
             subsamp < BAND_RATIO_AT_RES[INDEX_1000M_NIGHT]; subsamp++)
      {
        for (F = 0; F < MAX_1KM_OBC_FRAME_DIM; F++)
          rejects[F] = 0;
        DN_upper_valid_limit = SATURATED_DN - 1;

        if (Ecal_On[S][B_38] == True)
        {
          mean = -1;
          sdev = 0;
        }
        else if (moon_arrays->moon_in_SV_KOB[S][B_38] == MOON_INSIDE_SV_KOB)
        {
          returnStatus = Get_DN_Avg_SDev_Rejects_LowN(
                             refl_tables->RSB_SV_DN_moon_include_frames,
                             subsamp,
                             SV_1km_FRAMES,
                             increment,
                             SV_1km_night[T][B],
                             DN_upper_valid_limit,
                             0,
                             &mean,
                             &sdev,
                             &rejects[0]);
          if (returnStatus != MODIS_S_OK)
            SMF_ERROR(returnStatus,
                      "Get_DN_Avg_SDev_Rejects_LowN() in Fill_1km_day_DN_OBC_Avg()");
        }
        else
        {
          returnStatus =
              Get_DN_Avg_SDev_Rejects(
                refl_tables->DN_obc_avg_first_frame_to_use + subsamp,
                refl_tables->DN_obc_avg_number_of_frames_to_use,
                increment,
                SV_1km_night[T][B],
                DN_upper_valid_limit,
                0,
                &mean,
                &sdev,
                &rejects[refl_tables->DN_obc_avg_first_frame_to_use]);
          if (returnStatus != MODIS_S_OK)
            SMF_ERROR(returnStatus,
                "Get_DN_Avg_SDev_Rejects() in Fill_Band_26_DN_OBC_Avg()");
        }
        if (mean < 0)
          QA_refl->all_SV_DN_bad[S] = 1;

        DN_OBC_Avg->DN_obc_1km_night_avg[T][B][subsamp] = mean;
        DN_OBC_Avg->DN_obc_1km_night_var[T][B][subsamp] = sdev * sdev;
        DN_OBC_Avg->DN_obc_1km_night_outlier_mask[T][B][subsamp][0] = 0;
        DN_OBC_Avg->DN_obc_1km_night_outlier_mask[T][B][subsamp][1] = 0;
        returnStatus =
            Pack_Rejects_In_Outlier_Mask(SV_1km_FRAMES,
                 rejects,
                 DN_OBC_Avg->DN_obc_1km_night_outlier_mask[T][B][subsamp]);
      }/* end loop through subsamples */
    }/* end loop through detector */
  }/* end loop through scans */

  free(SV_1km_night);
  
  return (returnStatus);
}

PGSt_SMF_status Adjust_dn_star_Min(float32 *dn_star_Min,
                                   int32   L1A_v_id,
                                   int32   num_scans)
/*
!C**************************************************************************
!Description:  This routine check if the NAD door is closed. If the door is
               closed during any scan within the granule, set the dn_star_Min
               value to be -40 for those bands where dn_star_Min was set to 0.
               The reason is that if NAD door is closed, DN values should be
               very small and dn_star values for EV could be small negative
               values due to noise. So, set the minimum valid dn_star value to
               be -40.

!Input Parameters:
      float32  *    dn_star_Min      Mininum dn_star values for reflective
                                        bands
      int32         L1A_v_id         L1A file id for vdata interface
      int32         num_scans        number of scans in the L1A granule

!Output Parameters:
      float32  *    dn_star_Min      Mininum dn_star values for reflective
                                        bands

!Revision History:
   Revision 01.01  June 23, 2003, Razor Issue #192
    Changed the telemetry mnemonic which determines whether the nadir aperture
    door (NAD) is open from "CR_DR_NAD_OPEN" to the  equivalent "CR_DR_NAD_CLSD"
    mnemonic and altered the L1B code logic to correctly use the new value.  It
    was discovered that at two different times in the history of the MODIS/Terra
    instrument the mnemonic "CR_DR_NAD_OPEN" did not correctly reflect the state
    of the NAD, whereas the mnemonic "CR_DR_NAD_CLSD" has been consistently
    reliable.  
   Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.00 Oct. 12, 1999
   Initial development
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/

{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  uint16         NAD_door_clsd[MAX_NUM_SCANS];
  int16          S;
  int16          B;

  /* Check the input parameters are valid values */

  if (dn_star_Min == NULL)
    SMF_ERROR(MODIS_F_NOK,
        "dn_star_Min is NULL in Adjust_dn_star_Min(), Preprocess.c");

  if (num_scans < 0 || num_scans > MAX_NUM_SCANS)
    SMF_ERROR(MODIS_F_NOK,
        "invalid value for num_scans in Adjust_dn_star_Min(), Preprocess.c");

  /* Read the telemetry field for NAD door closed */

  returnStatus = read_vdata (L1A_v_id,
                             0,
                             num_scans,
                             "Telemetry Major Cycle 3A of 7",
                             "CR_DR_NAD_CLSD",
                             NAD_door_clsd);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg("Adjust_dn_star_Min", returnStatus, NULL, "read_vdata",
                FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  for (S = 0; S < num_scans; S++)
  {
    if (NAD_door_clsd[S] == 1)
    {
      for (B = 0; B < NUM_REFLECTIVE_BANDS; B++)
      {
        /*
         * If the NAD door is closed in any scan within the granule,
         * set the minimum value for dn_star to -40.
         */

        dn_star_Min[B] = -40;
      }

      break;
    }
  }

  return (returnStatus);
}

PGSt_SMF_status Check_For_Ecal_On(int32    lun,
                                  int32    num_scans,
                                  int32    v_id,
                                  boolean  Ecal_On[][NUM_BANDS])
/*
!C**************************************************************************
!Description:
   This function reads telemetry fields for Ecal-on for different
   kinds of bands, i.e. PV VIS, PV NIR, PV SM, PV LW and PC LW bands.
   There are two Ecal devices for each kind of bands, i.e. A and B.
   If one of them is on, set Ecal_On values for the corresponding bands
   for current scan to be true. If Ecal is on, the corresponding SV
   data can not be used for calibration.

!Input Parameters:
   int32     lun           The LUN of the L1A file being read (this is input
                           to allow identification in error messages).
   int32     num_scans     Number of scans in the L1A granule (could be the
                           leading, middle or trailing granule)
   int32     v_id          Vdata interface ID for the appropriate L1A file
                           (should be open for reading)

!Output Parameters:
   boolean   Ecal_On[][NUM_BANDS]  Ecal on or off information for every band
                                   for each scan
!Revision History:
   Revision 01.03  October 16, 2004  Razor Issue #200
   Casted Int32 variables in sprintf calls to "long" with the
   format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.02 March 2003, Razor Issue #173
   In the initializers of array-of-struct "temp", enclosed its rows with braces
   for ANSI-C compliance.
   Liqin Tan, SAIC GSO (ltan@saicmodis.com)

   Revision 01.01  Sep. 5, 2000
   Added "lun" to argument list to allow identification of L1A granule.
   Corrected the upper bound on Ecal telemetries (should be 1).
   Added error messages identifying telemetry point where failed.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.00  Oct. 21, 1999
   Initial development
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   When Ecal is ON, the SV DNs are invalid (these are use to telemeter Ecal
   results down).  However, the DNs from other sectors should be valid.

   a. For RSBs, the calculation of DN_OBC_avg should fail-over to use the
      BB DNs.
   b. For TEBs, a valid DN_OBC_avg (SV avg) cannot be calculated.  This should
      result in -1 in the DN_OBC_avg arrays for those affected bands and then
      the bad data pixel value will be set appropriately in Emissive_Cal.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32   S;
  int32   B_38;
  int32   i;

  int32   number_bands[NUM_BAND_TYPES] = {NUM_PV_VIS_BANDS,
                                          NUM_PV_NIR_BANDS,
                                          NUM_SM_BANDS,
                                          NUM_PV_LW_BANDS,
                                          NUM_PC_LW_BANDS};

  enum {
    PV_VIS,
    PV_NIR,
    PV_SM,
    PV_LW,
    PC_LW
  } band_type;

  int32   bands[NUM_BAND_TYPES][MAX_NUM_BANDS_PER_TYPE] = {
    /* VIS bands */
    {2, 3, 7, 8, 9, 10, 11},
    /* NIR bands */
    {0, 1, 12, 13, 14, 15, 16, 17, 18, 19, 20},
    /* SM bands */
    {4, 5, 6, 21, 22, 23, 24, 25, 26, 27},
    /* PV LW bands */
    {28, 29, 30, 31},
    /* PC LW bands */
    {32, 33, 34, 35, 36, 37}
  };

  /*
   * The number of the fields below, which equals the value of NUM_FIELDS
   * should be exactly twice as many as NUM_BAND_TYPES. Every two consecutive
   * fields with the same name except the last character, i.e. A or B, is
   * used for one kind of bands.
   * If this is changed, the logic below is not right.
   */

  enum {
    VIS_ECAL_A,
    VIS_ECAL_B,
    NIR_ECAL_A,
    NIR_ECAL_B,
    SM_ECAL_A,
    SM_ECAL_B,
    PV_LW_ECAL_A,
    PV_LW_ECAL_B,
    PC_LW_ECAL_A,
    PC_LW_ECAL_B,
    NUM_FIELDS
  } field;

  struct {
    char *vname;
    char *fname;
    uint16 buffer[MAX_NUM_SCANS];
  }temp[NUM_FIELDS] = {
    {"Telemetry Major Cycle 4B of 7", "CR_PVVISA_ECALON", {0}},
    {"Telemetry Major Cycle 4B of 7", "CR_PVVISB_ECALON", {0}},
    {"Telemetry Major Cycle 4B of 7", "CR_PVNIRA_ECALON", {0}},
    {"Telemetry Major Cycle 4B of 7", "CR_PVNIRB_ECALON", {0}},
    {"Telemetry Major Cycle 4B of 7", "CR_PVSMA_ECAL_ON", {0}},
    {"Telemetry Major Cycle 4B of 7", "CR_PVSMB_ECAL_ON", {0}},
    {"Telemetry Major Cycle 4A of 7", "CR_PVLWA_ECAL_ON", {0}},
    {"Telemetry Major Cycle 4A of 7", "CR_PVLWB_ECAL_ON", {0}},
    {"Telemetry Major Cycle 4A of 7", "CR_PCLWA_ECAL_ON", {0}},
    {"Telemetry Major Cycle 4A of 7", "CR_PCLWB_ECAL_ON", {0}}
  };
  char *location = "Check_For_Ecal_On";

  /* Check if the value of num_scans is valid */

  if (num_scans < 0 || num_scans > MAX_NUM_SCANS)
  {
    char errmsg[256];
    sprintf (errmsg, "The input value of \"num_scans\" (%ld) is invalid.",
             (long)num_scans);
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, errmsg, NULL, lun, NULL, True);
    return returnStatus;
  }

  /* Read the vdata and check the value */

  for (field = VIS_ECAL_A; field < NUM_FIELDS; field++)
  {
    returnStatus = read_vdata (v_id, 0, num_scans, temp[field].vname,
                               temp[field].fname,
                               (VOIDP)temp[field].buffer);
    if (returnStatus != MODIS_S_OK)
    {
      char errmsg[256];
      sprintf (errmsg, "Failed to read field \"%s\" from telemetry \"%s\".",
                       temp[field].fname, temp[field].vname);
      L1BErrorMsg(location, returnStatus, errmsg, "read_vdata",
                  lun, NULL, True);
      return returnStatus;
    }
    returnStatus = Check_Valid_Range(temp[field].fname, DFNT_UINT16,
                                     NULL, "1", NULL,
                                     num_scans, (void *) temp[field].buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                  lun, NULL, True);
      return returnStatus;
    }
  }

  for (S = 0; S < num_scans; S++)
  {
    for (band_type = PV_VIS, field = VIS_ECAL_A; band_type < NUM_BAND_TYPES;
          band_type++, field += 2)
    {
      /*
       * If Ecal A on or/and Ecal B on, it is considered to be Ecal on.
       * Set Ecal_On to be True for the corresponding bands for this scan.
       * Otherwise, set that to be False.
       */

      if (temp[field].buffer[S] == 1 || temp[field+1].buffer[S] == 1)
        for (i = 0; i < number_bands[band_type]; i++)
        {
          B_38 = bands[band_type][i];
          Ecal_On[S][B_38] = True;
        }

      else
        for (i = 0; i < number_bands[band_type]; i++)
        {
          B_38 = bands[band_type][i];
          Ecal_On[S][B_38] = False;
        }
    }
  }

  return (returnStatus);
}

PGSt_SMF_status Check_If_Sector_Rotated(int32    lun,
                                        int32    num_scans,
                                        int32    v_id,
                                        boolean  Sector_Rotated[])
/*
!C****************************************************************************
!Description:
   This function determines if the MODIS data collection sectors were
   rotated from their normal angular positions on any scan.  A sector rotation
   event, commanded from the ground, is used occasionally for calibration
   studies.

!Input Parameters:
   int32     lun           The LUN of the L1A file being read (this is input
                           to allow identification in error messages).
   int32     num_scans     Number of scans in the L1A granule (could be the
                           leading, middle or trailing granule)
   int32     v_id          Vdata interface ID for the appropriate L1A file
                           (should be open for reading)

!Output Parameters:
   boolean   Sector_Rotated[]  Array describing if there is any degree
                               of sector rotation for each scan.
                               Values of each element are "True" or "False".
                               (Array assumed allocated to MAX_NUM_SCANS).

!Revision History:
   Revision 01.01  October 16, 2004  Razor Issue #200
   Casted Int32 variables in sprintf calls to "long" with the
   format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.00  Sep. 5, 2000
   Initial development
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   The telemetry name is "Telemetry Major Cycle 6 of 7" and the field
   name is "CS_FR_ENC_DELTA". If the value is non-zero on any scan, then
   it is assumed that the sector is rotated for that scan.

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32  S;                              /* scan index */
  uint16 sect_rot_telem[MAX_NUM_SCANS];  /* holds telemetry data */
  char *location = "Check_If_Sector_Rotated";

      /* Check if the value of num_scans is valid */

  if (num_scans < 0 || num_scans > MAX_NUM_SCANS)
  {
    char errmsg[256];
    sprintf (errmsg, "The input value of \"num_scans\" (%ld) is invalid.",
             (long)num_scans);
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, errmsg, NULL, lun, NULL, True);
    return returnStatus;
  }


 /*--------------------------------------------------------------------- 
 * 03/01/2015 Xg 
 * The telemetry data CS_FR_ENC_DELTA may not represent correctly    
 * the sector rotation status. Our observations indicate that the 
 * SET_FR_ENC_DELTA in the Command Parameters field correlates well
 * with the sector rotation encoder position. 
 *
  returnStatus = read_vdata (v_id,
                             0,
                             num_scans,
                             "Telemetry Major Cycle 6 of 7",
                             "CS_FR_ENC_DELTA",
                             (VOIDP) sect_rot_telem);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus,
                "Failed to read telemetry field \"CS_FR_ENC_DELTA\".",
                "read_vdata", lun, NULL, True);
    return returnStatus;
  }
  */

      /* Check the upper bound.  If it is exceeded, then the L1A file
       * must be corrupted.
       */

  /*
  returnStatus = Check_Valid_Range("CS_FR_ENC_DELTA", DFNT_UINT16,
                                   NULL, "16383", NULL,
                                   num_scans, (void *) sect_rot_telem);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                lun, NULL, True);
    return returnStatus;
  }
  --------------------------------------------------------------------*/

  /* Read the telemetry data. */
  returnStatus = read_vdata (v_id,
                             0,
                             num_scans,
                             "Command Parameters",
                             "SET_FR_ENC_DELTA",
                             (VOIDP) sect_rot_telem);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus,
                "Failed to read telemetry field \"SET_FR_ENC_DELTA\".",
                "read_vdata", lun, NULL, True);
    return returnStatus;
  }


  /* Check the upper bound.  If it is exceeded, then the L1A file
   * must be corrupted.
   */

  returnStatus = Check_Valid_Range("SET_FR_ENC_DELTA", DFNT_UINT16,
                                    NULL, "16383", NULL,
                                    num_scans, (void *) sect_rot_telem);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                lun, NULL, True);
    return returnStatus;
  }
 
  /* Loop through scans, check each value, assign True or False. */
  for (S = 0; S < num_scans; S++)
  {
    if (sect_rot_telem[S] > 0) 
      Sector_Rotated[S] = True;
    else 
      Sector_Rotated[S] = False;
  }

  returnStatus = MODIS_S_OK;
  return returnStatus;
}


PGSt_SMF_status sort_int16_array (int32 n,
                                  int16 *a,
                                  int32 *indx)

/*
!C****************************************************************************
!Description:
    This function sorts an array of "int16" values from the lowest to
    highest.  The sorted values are returned in the original array (the
    original is changed upon exit).

!Input Parameters:
   int32 n       number of values in the array "a".
   int16 *a      array of n values to be sorted, from lowest to highest.
   int32 *indx   allocated space, contents unimportant

!Output Parameters:
   int16 *a      sorted array.
   int32 *indx   the indexes of the original array "a", after sorting

!Revision History:
   Revision 01.00  Sep. 11, 2000
   Initial development
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16 atemp;   /* Temporary variable */
  int32 itemp;   /* Temporary variable */
  int32 i, j;    /* Loop indices */

  for (i = 0; i < n; i++)
    indx[i] = i;
  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      if (a[j] < a[i]) {
        atemp   = a[i];
        a[i]    = a[j];
        a[j]    = atemp;
        itemp   = indx[i];
        indx[i] = indx[j];
        indx[j] = itemp;
      }
    }
  }
  return returnStatus;
}


PGSt_SMF_status Get_DN_Avg_SDev_Rejects_LowN(int32   N_include,
                                             int32   start_index,
                                             int32   N,
                                             int32   index_increment,
                                             int16   *DN_array,
                                             int16   DN_upper_valid_limit,
                                             int16   DN_lower_valid_limit,
                                             float32 *mean_DN,
                                             float32 *sdev_DN,
                                             int16   *rejects)
/*
!C****************************************************************************
!Description:
   This function has the same purpose as "Get_DN_Avg_SDev_Rejects" -- to
   compute the average, standard deviation and the rejected values for a set
   of DNs.  In fact, this function calls "Get_DN_Avg_SDev_Rejects".  However,
   prior to calling that function, the DNs are sorted from lowest to highest
   -- and only the lowest "N_include" values are passed into
   Get_DN_Avg_SDev_Rejects to compute the output values.

   This function is intended for use in calculating the space view (SV)
   average DN when the moon is in the space view port.  In that situation,
   several 1km frames will be corrupted by the moon, having abnormally high
   values.  By sorting the DNs from lowest to highest, we can eliminate
   those abnormally high values.

!Input Parameters:
   int32     N_include        Number of values to use to form the set of DNs
                              which are passed into "Get_DN_Avg_SDev_Rejects".
                              After sorting the values from lowest to highest,
                              the lowest N_include values will be used.
   int32     start_index      Starting index in the array "DN_array".
                              This value is normally the subsample index of
                              the DNs being averaged.
   int32     N                Number of DN values to form an array before
                              sorting.  This number should be greater than or
                              equal to N_include.  This number is normally
                              the number of 1km frames.
   int32     index_increment  Amount to increment when forming the array
                              of DNs to sort.  This value is normally the
                              number of subsamples appropriate to the
                              resolution of the data (1, 2 or 4).
   int16 *   DN_array         Array of digital numbers for all frames and
                              subsamples.  There should be
                              N * index_increment values in this array.
   int16     DN_upper_valid_limit   Upper limit of the valid range of
                                    DNs (typically 4094).
   int16     DN_lower_valid_limit   Lower limit of the valid range of
                                    DNs (typically 0).

!Output Parameters:
   float32 * mean_DN          Mean value of final set of data after all
                              invalid and outlier values are removed.  A
                              negative return value means that all values
                              were rejected for some reason.
   float32 * sdev_DN          standard deviation of the final set of values.
                              A negative return value means that all values
                              were rejected for some reason.
   int16 *   rejects          An array of dimension [N], where 1 denotes a
                              value was removed and 0 denotes that a value
                              was used in the final set.

!Revision History:
   Revision 01.11  October 16, 2004  Razor Issue #200
   Casted Int32 variables in sprintf calls to "long" with the
   format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.10 March 25, 2002   Razor Issue #178
   Strip out ADC Correction (delta_DN)
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.00  Sep. 13, 2000
   Initial development
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 i;                                   /* index */
  int16 DN_temp[MAX_1KM_OBC_FRAME_DIM];      /* Temp array for DNs */
  int32 DN_temp_indx[MAX_1KM_OBC_FRAME_DIM]; /* Temp array for indexes */
  int16 rejects_temp[MAX_1KM_OBC_FRAME_DIM]; /* Temp array for rejects */
  char *location = "Get_DN_Avg_SDev_Rejects_LowN";

  /*
   * Check dimensions and input values.
   */

  if (N > MAX_1KM_OBC_FRAME_DIM) {
    char errmsg[256];
    sprintf (errmsg, "Input N (%ld) exceeds MAX_1KM_OBC_FRAME_DIM (%ld)\n",
             (long)N, (long)MAX_1KM_OBC_FRAME_DIM);
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, errmsg, NULL, 0, NULL, True);
    return returnStatus;
  }
  if (N_include > N) {
    char errmsg[256];
    sprintf (errmsg, "Input N_include (%ld) exceeds N (%ld)\n",
             (long)N_include, (long)N);
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, errmsg, NULL, 0, NULL, True);
    return returnStatus;
  }
  if (N_include < 1) {
    char errmsg[256];
    sprintf (errmsg, "Input N_include (%ld) is invalid (should be > 0)\n",
             (long)N_include);
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, errmsg, NULL, 0, NULL, True);
    return returnStatus;
  }

  /*
   * Form the set of DNs and sort from lowest to highest.
   */

  for (i = 0; i < N; i++)
    DN_temp[i] = DN_array[start_index + i * index_increment];

  returnStatus = sort_int16_array (N, DN_temp, DN_temp_indx);
  if (returnStatus != MODIS_S_OK) {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus, NULL, "sort_int16_array",
                0, NULL, True);
    return returnStatus;
  }

  /*
   * Calculate <DN>, sdev and rejects.  Put rejects into temp array.
   */

  returnStatus = Get_DN_Avg_SDev_Rejects(0,
                                         N_include,
                                         1,
                                         DN_temp,
                                         DN_upper_valid_limit,
                                         DN_lower_valid_limit,
                                         mean_DN,
                                         sdev_DN,
                                         rejects_temp);
  if (returnStatus != MODIS_S_OK) {
    L1BErrorMsg(location, returnStatus, NULL, "Get_DN_Avg_SDev_Rejects",
                0, NULL, True);
    return returnStatus;
  }

  /*
   * Fill the rejects return variable.
   */

  for (i = 0; i < N_include; i++)
    rejects[DN_temp_indx[i]] = rejects_temp[i];
  for (  ; i < N; i++)
    rejects[DN_temp_indx[i]] = 1;

  return returnStatus;
}


PGSt_SMF_status Calculate_RVS_Correction(lookup_tables_t *tables,
                                         L1B_granule_t   *L1B_Gran)
/*
!C****************************************************************************
!Description:
 This function calculates the frame-by-frame RVS correction terms for
 the reflective solar bands and thermal emissive bands.  The calculated
 values are stored in the RSB_Cal_Coeff structure for use in
 Reflective_Cal and in the Emiss_Cal_Coeff structure for use
 in Emissive_Cal respectively.  The RVS correction values for the SV and
 BB are also used in Preprocess.c.
 The frame dependent uncertainty terms are also calculated in this function.

!Input Parameters:
 lookup_tables_t *  tables    Contains structure members for m1 and R*
                              time-dependent LUTs and dead detector list.
 L1B_granule_t      *L1B_Gran Will contain L1B data and also L1B Calibration
                              information.  Of interest are the RVS
                              correction terms which need to be calculated
                              so that they can be used in preprocessing.

!Output Parameters:
 L1B_granule_t      *L1B_Gran with the RVS correction inserted.

!Revision History:
 Revision 02.13 March 2003, Razor Issue #173
 Initialized variable "rvs_corr" to 0.0 for ANSI-C compliance.
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 1.01  May 14, 2002
 Replaced incorrect comparison of RVS coefficients for TEBs to RSB
 coefficents with a comparison to TEB coefficients.  This error never
 resulted in an incorrect answer but did lead to recalculating the
 coefficients for each detector in a thermal band even if the coefficients
 were the same for all 10 detectors.
 Alice Isaacman, SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Initial development  March 10, 2002      Razor Issue #174
 Alice Isaacman, SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END**************************************************************************
*/

{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16 band        = 0;     /* Band index */
  int16 ref_band    = 0;     /* Reference band index */         
  int16 det         = 0;     /* Detector */     
  int16 det_160     = 0;     /* Emissive detector index */         
  int16 sample      = 0;     /* Sample number */        
  int16 mirr_side   = 0;     /* Mirror side */           
  int16 frame       = 0;     /* Frame number */       
  int16 ref_frame   = 0;     /* Reference frame number */ 
  int16 BB_frame_no = 0;     /* Frame number to use for BB RVS Correction */
  int16 SV_frame_no = 0;     /* Frame number to use for SV RVS Correction */
  int16 i           = 0;     /* temporal variable */
  boolean coeffs_same = FALSE;
  boolean det_0_or_coeffs_diff = FALSE;
  float32 *rvs_coeffs;     /* Temporary pointer to current RVS coefficients */
  float32 *sigma_rvs_coeffs;     /* Temporary pointer to RVS uncertainty coefficients */
  float32 rvs_corr = 0.0;    /* Calculated RVS correction */
  float32 sigma_rvs = 0.0;       /* Calculated Sigma RVS */
  float32 *u2_samples;     /* Temporary pointer to RSB u2 uncertainty */
  int16   *u2_frames;       /* Temporary pointer to the specified frames loaded with RSB u2 uncertainty */
  float32 y1,y2;            /* temporary variables for u2 interpolation */
  int16   f1, f2, f_end;    /* temporary variables for frame index */
  char *location = "Calculate_RVS_Correction";

  int32 frame_no_squared[EV_1km_FRAMES];
                           /*  Array of squares of frames */

  /* Initialize the array of squares of frames */      
  for ( frame = 0; frame < EV_1km_FRAMES; frame++ ) 
      frame_no_squared[frame] = frame*frame;
      
  /* Frame numbers to use for BB and SV calculations */
    
  BB_frame_no = tables->emiss.RVS_BB_SV_Frame_No[0];
  SV_frame_no = tables->emiss.RVS_BB_SV_Frame_No[1];
            
  /* Calculate the frame-by-frame RVS Correction Terms for 250M Bands. */

  /* Loop through bands, detectors, and mirror sides: */
  
  for (band = 0; band < NUM_250M_BANDS; band++ ) 
  {
    for ( det = 0; det < DETECTORS_PER_250M_BAND; det++) 
    {
      for ( mirr_side = 0; mirr_side < NUM_MIRROR_SIDES; mirr_side++ ) 
      {
        
        /* Read the set of three rvs coefficients. */
        
        rvs_coeffs = tables->refl.RVS_RefSB[band][det][mirr_side];

      /*
       * If this is not the first detector, check to see whether these
       *     coefficients are the same as those for the detector before it.	      
       */
        
        if (det > 0)
        {
          i = 0;
          coeffs_same = TRUE;
          while (coeffs_same && i < NUM_RSB_RVS_COEFFS)
          {
            coeffs_same = (tables->refl.RVS_RefSB[band][det - 1][mirr_side][i] == rvs_coeffs[i]);
            i++;
          }
        }
        else
          coeffs_same = FALSE;
        
        det_0_or_coeffs_diff = (det == 0 || (coeffs_same == FALSE));    

        /* Run through each EV frame and calculate the correction. */            
        
        for ( ref_frame = 0; ref_frame < EV_1km_FRAMES; ref_frame++ )
        {
          
        /*
         * If this is the first detector or these coefficients are not the
         * same as those for the previous detector, calculate the new correction.
         */
          
          if (det_0_or_coeffs_diff)
          {    
            /* calculate RVS correction */
            rvs_corr = 0.;
            for (i = NUM_RSB_RVS_COEFFS-1; i >= 0; i -- )
            {
              rvs_corr = rvs_corr*ref_frame + rvs_coeffs[i];
            }

            if (rvs_corr < RVS_CORRECTION_LOWER_LIMIT ||
                rvs_corr > RVS_CORRECTION_UPPER_LIMIT)

            {
              returnStatus = MODIS_F_OUT_OF_RANGE;
              L1BErrorMsg(location, returnStatus, NULL, NULL, 
                  0, NULL, True);
              return returnStatus;
            }
          } /* Coefficients not the same */

          /* Set the actual frame number for the 250m resolution: */
          
          frame = NUM_250M_SUBSAMP*ref_frame;
          
          /* Loop through the subframes located here: */
          
          for ( sample = 0; sample < NUM_250M_SUBSAMP; sample++ ) 
          {
            /* If the correction was recalculated, use it now: */
            
            if (det_0_or_coeffs_diff)
              
                L1B_Gran->RSB_Cal_Coeff.RVS_250m
	                [band][det][frame][mirr_side] = rvs_corr;

			      else
              
            /* Otherwise, use the previously recorded correction term. */
              
		          L1B_Gran->RSB_Cal_Coeff.RVS_250m
                [band][det][frame][mirr_side] =
              L1B_Gran->RSB_Cal_Coeff.RVS_250m
		            [band][det - 1][frame][mirr_side];
			      
			      frame++;
            
          }  /* end loop through samples   */
        }    /* end loop through ref_frame */
      }      /* end loop through mirr_side */
    }        /* end loop through det       */
  }          /* end loop through bands     */

  /* Calculate the frame-by-frame RVS Correction Terms for 500M Bands. */

  /* Loop through bands, detectors, and mirror sides: */
  
  for (band = 0; band < NUM_500M_BANDS; band++ ) 
  {
    /* Find this band's index in the set of reflective solar bands. */
    
    ref_band = band + NUM_250M_BANDS; 
      
    for ( det = 0; det < DETECTORS_PER_500M_BAND; det++) 
    {
      for ( mirr_side = 0; mirr_side < NUM_MIRROR_SIDES; mirr_side++ ) 
      {
        
        /* Read the set of three rvs coefficients. */
        
        rvs_coeffs = tables->refl.RVS_RefSB[ref_band][det][mirr_side];

      /*
       * If this is not the first detector, check to see whether these
       *     coefficients are the same as those for the detector before it.	      
       */
        
        if (det > 0)
        {
          i = 0;
          coeffs_same = TRUE;
          while (coeffs_same && i < NUM_RSB_RVS_COEFFS)
          {
            coeffs_same = (tables->refl.RVS_RefSB[band][det - 1][mirr_side][i] == rvs_coeffs[i]);
            i++;
          }
        }
        else
          coeffs_same = FALSE;
        
        det_0_or_coeffs_diff = (det == 0 || (coeffs_same == FALSE));    
	
        /* Run through each EV frame and calculate the correction. */            
        
        for ( ref_frame = 0; ref_frame < EV_1km_FRAMES; ref_frame++ )
        {
          
        /*
         * If this is the first detector or these coefficients are not the
         * same as those for the previous detector, calculate the new correction.
         */
          
          if (det_0_or_coeffs_diff)
          {
            /* calculate RVS correction */
            rvs_corr = 0.;
            for (i = NUM_RSB_RVS_COEFFS-1; i >= 0; i -- )
            {
              rvs_corr = rvs_corr*ref_frame + rvs_coeffs[i];
            }
	      
            if (rvs_corr < RVS_CORRECTION_LOWER_LIMIT ||
                rvs_corr > RVS_CORRECTION_UPPER_LIMIT)

            {
              returnStatus = MODIS_F_OUT_OF_RANGE;
              L1BErrorMsg(location, returnStatus, NULL, NULL, 
                  0, NULL, True);
              return returnStatus;
            }
          }  /* Coefficients not same */
	    
          /* Set the actual frame number for the 500m resolution: */
          
          frame = NUM_500M_SUBSAMP*ref_frame;
          
          /* Loop through the subframes located here: */
          
          for ( sample = 0; sample < NUM_500M_SUBSAMP; sample++ ) 
          {
            /* If the correction was recalculated, use it now: */
            
            if (det_0_or_coeffs_diff)
              
                L1B_Gran->RSB_Cal_Coeff.RVS_500m
                  [band][det][frame][mirr_side] = rvs_corr;
            
            else
              
            /* Otherwise, use the previously recorded correction term. */
              
              L1B_Gran->RSB_Cal_Coeff.RVS_500m
                [band][det][frame][mirr_side] =
              L1B_Gran->RSB_Cal_Coeff.RVS_500m
                [band][det - 1][frame][mirr_side];
                               
            frame++;
            
          }  /* end loop through samples   */
        }    /* end loop through ref_frame */
      }      /* end loop through mirr_side */
    }        /* end loop through det       */
  }          /* end loop through bands     */


  /* Calculate the frame-by-frame RVS Correction Terms for 1KM RSBs. */

  /* Loop through bands, detectors, and mirror sides: */
  
  for (band = 0; band < NUM_1000M_REFL_BANDS; band++ ) 
  {
    /* Find this band's index in the set of reflective solar bands. */
    
    ref_band = band + (NUM_250M_BANDS + NUM_500M_BANDS);
      
    for ( det = 0; det < DETECTORS_PER_1KM_BAND; det++) 
    {
      for ( mirr_side = 0; mirr_side < NUM_MIRROR_SIDES; mirr_side++ ) 
      {
        /* Read the set of three rvs coefficients. */
        
        rvs_coeffs = tables->refl.RVS_RefSB[ref_band][det][mirr_side];
          
      /*
       * If this is not the first detector, check to see whether these
       *     coefficients are the same as those for the detector before it.	      
       */
        
        if (det > 0)
        {
          i = 0;
          coeffs_same = TRUE;
          while (coeffs_same && i < NUM_RSB_RVS_COEFFS)
          {
            coeffs_same = (tables->refl.RVS_RefSB[band][det - 1][mirr_side][i] == rvs_coeffs[i]);
            i++;
          }
        }
        else
          coeffs_same = FALSE;
        
        det_0_or_coeffs_diff = (det == 0 || (coeffs_same == FALSE));    
		
        /* Run through each EV frame and calculate the correction. */            
        
        for ( frame = 0; frame < EV_1km_FRAMES; frame++ )
        {
          
        /*
         * If this is the first detector or these coefficients are not the
         * same as those for the previous detector, calculate the new correction.
         */
          
          if (det_0_or_coeffs_diff)
          {
            /* calculate RVS correction */
            rvs_corr = 0.;
            for (i = NUM_RSB_RVS_COEFFS-1; i >= 0; i -- )
            {
              rvs_corr = rvs_corr*frame + rvs_coeffs[i];
            }

            if (rvs_corr < RVS_CORRECTION_LOWER_LIMIT ||
                rvs_corr > RVS_CORRECTION_UPPER_LIMIT)

            {
              returnStatus = MODIS_F_OUT_OF_RANGE;
              L1BErrorMsg(location, returnStatus, NULL, NULL, 
                  0, NULL, True);
              return returnStatus;
            }

	          /* Fill in the calculated correction: */
            
            L1B_Gran->RSB_Cal_Coeff.RVS_1km_RefSB
              [band][det][frame][mirr_side] = rvs_corr;
            
          }  /* Coefficients not same */
          
          else
            
            /* Otherwise, use the previously recorded correction term. */
              
            L1B_Gran->RSB_Cal_Coeff.RVS_1km_RefSB
              [band][det][frame][mirr_side] =
            L1B_Gran->RSB_Cal_Coeff.RVS_1km_RefSB
              [band][det - 1][frame][mirr_side] ;

        }    /* end loop through frame */
      }      /* end loop through mirr_side */
    }        /* end loop through det       */
  }          /* end loop through bands     */

  /* calculate the frame-dependent uncertainty */
  u2_frames = tables->refl.u2_frames;
  /* interpolation to get u2 on each frame */
  for ( det = 0; det < NUM_REFLECTIVE_DETECTORS; det++)
  {
    for ( mirr_side = 0; mirr_side < NUM_MIRROR_SIDES; mirr_side++ )
    {
      /* Read the set of three rvs coefficients. */
      u2_samples = tables->refl.u2_samples[det][mirr_side];
      for ( frame = 0; frame < EV_1km_FRAMES; frame++ )
      {
        if (frame <= u2_frames[0]) L1B_Gran->RSB_Cal_Coeff.u2[det][frame][mirr_side] = u2_samples[0];
        else if (frame >= u2_frames[NUM_U2_FRAME-1])
          L1B_Gran->RSB_Cal_Coeff.u2[det][frame][mirr_side] = u2_samples[NUM_U2_FRAME-1];
        else
        {
          /*
           * implement interpolation
           */
           f_end = 1;
           while ( frame > u2_frames[f_end] && f_end < (NUM_U2_FRAME-1) )
             f_end++;
           f1 = u2_frames[f_end-1];
           f2 = u2_frames[f_end];
           y1 = u2_samples[f_end-1];
           y2 = u2_samples[f_end];
           L1B_Gran->RSB_Cal_Coeff.u2[det][frame][mirr_side] = y1*((float32)(f2-frame))/((float32)(f2-f1))
                                        + y2*((float32)(frame-f1))/((float32)(f2-f1));
        }
      } /* end loop through frame */
    }   /* end loop through mirr_side */
  }     /* end loop through det       */


  /* Calculate the frame-by-frame RVS Correction Terms for Emissive Bands. */
  
  /* Set the emissive detector index to 0. */
  
  det_160 = 0;

  /* Loop through bands, detectors, and mirror sides: */
  
  for (band = 0; band < NUM_EMISSIVE_BANDS; band++ ) 
  {
    for ( det = 0; det < DETECTORS_PER_1KM_BAND; det++, det_160++ ) 
    {
      for ( mirr_side = 0; mirr_side < NUM_MIRROR_SIDES; mirr_side++ ) 
      {
        /* Read the set of three rvs coefficients. */
        
        rvs_coeffs = tables->emiss.RVS_TEB[band][det][mirr_side];

      /*
       * If this is not the first detector, check to see whether these
       *     coefficients are the same as those for the detector before it.	      
       */
        
        if (det > 0)
          coeffs_same = (
            tables->emiss.RVS_TEB[band][det - 1][mirr_side][0]
              == rvs_coeffs[0] &&
            tables->emiss.RVS_TEB[band][det - 1][mirr_side][1]
              == rvs_coeffs[1] &&
            tables->emiss.RVS_TEB[band][det - 1][mirr_side][2]
              == rvs_coeffs[2]);
        else
          coeffs_same = FALSE;
        
        det_0_or_coeffs_diff = (det == 0 || (coeffs_same == FALSE));    
	
        /* Run through each EV frame and calculate the correction. */            
        
        for ( frame = 0; frame < EV_1km_FRAMES; frame++ )
        {
          
        /*
         * If this is the first detector or these coefficients are not the
         * same as those for the previous detector, calculate the new correction.
         */
          
          if (det_0_or_coeffs_diff)
          {
            EVAL_2ND_ORDER_POLYNOMIAL(rvs_corr, rvs_coeffs, 
                frame, frame_no_squared[frame]);
          
            if (rvs_corr < RVS_CORRECTION_LOWER_LIMIT ||
                rvs_corr > RVS_CORRECTION_UPPER_LIMIT)

            {
              returnStatus = MODIS_F_OUT_OF_RANGE;
              L1BErrorMsg(location, returnStatus, NULL, NULL, 
                  0, NULL, True);
              return returnStatus;
            }

	          /* Fill in the calculated correction: */
            
            L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV
              [det_160][frame][mirr_side] = rvs_corr;
                
          }   /* Coefficients not same */
          
          else
            
            /* Otherwise, use the previously recorded correction term. */
              
            L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV
              [det_160][frame][mirr_side] =
            L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV
              [det_160 - 1][frame][mirr_side];
        }

        /* Calculate the sigma RVS for TEB */
        /* Read the set of three sigma rvs coefficients. */
        sigma_rvs_coeffs = tables->emiss.sigma_RVS_EV[band][det][mirr_side];

      /*
       * If this is not the first detector, check to see whether these
       *     coefficients are the same as those for the detector before it.
       */

        if (det > 0)
          coeffs_same = (
            tables->emiss.sigma_RVS_EV[band][det - 1][mirr_side][0]
              == sigma_rvs_coeffs[0] &&
            tables->emiss.sigma_RVS_EV[band][det - 1][mirr_side][1]
              == sigma_rvs_coeffs[1] &&
            tables->emiss.sigma_RVS_EV[band][det - 1][mirr_side][2]
              == sigma_rvs_coeffs[2]);
        else
          coeffs_same = FALSE;

        det_0_or_coeffs_diff = (det == 0 || (coeffs_same == FALSE));

        /* Run through each EV frame and calculate the correction. */
        for ( frame = 0; frame < EV_1km_FRAMES; frame++ )
        {

        /*
         * If this is the first detector or these coefficients are not the
         * same as those for the previous detector, calculate the new correction.
         */

          if (det_0_or_coeffs_diff)
          {
            EVAL_2ND_ORDER_POLYNOMIAL(sigma_rvs, sigma_rvs_coeffs,
                frame, frame_no_squared[frame]);

                  /* Fill in the calculated correction: */

            L1B_Gran->Emiss_Cal_Coeff.sigma_RVS_Emiss_EV
              [det_160][frame][mirr_side] = sigma_rvs;

          }   /* Coefficients not same */

          else

            /* Otherwise, use the previously recorded correction term. */

            L1B_Gran->Emiss_Cal_Coeff.sigma_RVS_Emiss_EV
              [det_160][frame][mirr_side] =
            L1B_Gran->Emiss_Cal_Coeff.sigma_RVS_Emiss_EV
              [det_160 - 1][frame][mirr_side];


        }    /* end loop through frame */
                   
      /* Fill in the BlackBody correction for this detector and mirror side: */
        
          L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_BB[det_160][mirr_side] = 
            L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV
              [det_160][BB_frame_no][mirr_side] ; 
        
      /* Fill in the Space View correction for this detector and mirror side: */
        
          L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_SV[det_160][mirr_side] = 
            L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV
              [det_160][SV_frame_no][mirr_side];  

      }      /* end loop through mirr_side */
    }        /* end loop through det       */
  }          /* end loop through bands     */

  return returnStatus;
}
              
