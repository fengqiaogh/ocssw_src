#include    "Preprocess.h"
#include    "Emissive_Cal.h"
#include    "HDF_Lib.h" 
#include    "PGS_Error_Codes.h"
#include    "FNames.h"
#include    <math.h>

#define FLOAT_TO_PIXEL(L, L_max, L_min, scale, offset, Scale_Integer)   \
if ((L) <= (L_min)) Scale_Integer = DN_MIN;                             \
else if ((L) >= (L_max)) Scale_Integer = DN15_SAT;                       \
else Scale_Integer = (uint16)( ( (L)/(scale) + (offset) ) + 0.5 );

/*
  Revision 05/10/06
  Skip "Bad_L1A_Error_Out" errors rather than exiting.
  J. Gales
*/

PGSt_SMF_status  Emissive_Cal (int16                S,
                               int16                MS,
                               L1B_granule_t        *L1B_Gran,
                               L1A_Scan_t           *L1A_Scan,
                               L1B_Scan_t           *L1B_Scan,
                               Preprocess_Emiss_t   *PP_emiss,
                               emiss_tables_t       *emiss_tables,
                               QA_tables_t          *QA_tables,
                               QA_Data_t            *QA,
                               L1B_Scan_Metadata_t  *L1B_Scan_Meta)
/*
!C**********************************************************************
!Description:
    This routine accomplishes the calibration of one scan of EV pixels for
    Emissive bands. Corrections to the digital numbers are applied for 
    electronic background,  for the angular dependence of the response of the
    scan mirror, for non-linearities in the Analog to Digital Converters. 
    Correction may also be applied for cross-talk for PC bands (the leak from
    band 31 to band 32, 33, 34, 35, 36) depending on the switch value in the
    lookup table. The radiance is computed after the corrections. If the 
    radiance can be computed and is in a valid range, it is converted
    to a scaled integer in the range of [0-32767]. If a valid value cannot be
    computed, the scaled integer is set to a value in the range of [32768-65535].
    Specific values in the range of [32768-65535] are used to denote why a valid
    value could not be obtained (a list of these is in the L1B file
    specifications).  This routine also computes the uncertainty in the
    radiance product for every pixel, and converts the uncertainty to a
    4-bit uncertainty index, stored in the 4 least significant bits of an 8-bit
    unsigned integer.

!Input Parameters:
      int16               S                 Current scan index 
      int16               MS                The corresponding mirror side
      L1B_granule_t       *L1B_Gran         contains scale and offset factors 
      L1A_Scan_t          *L1A_Scan         contains earth view DN values
      Preprocess_Emiss_t  *PP_emiss         contains ADC index, calibration coefficients a0,
                                            b1, a2 preprocessed data, radiance at scan mirror
                                            temperature, NEdL and electronic background DNs
      emiss_tables_t      *emiss_tables     contains emissive lookup tables values used in
                                            emissive calibration
      QA_tables_t         *QA_tables        contains the dead detector look up table value and
                                            NEdL look up table value             
      QA_Data_t           *QA               contains values to determine if the instrument is 
                                            in a sector rotation and if the moon is in SV KOB.                                  
!Output Parameters:
      L1B_Scan_t          *L1B_Scan         contains scaled integers and uncertainty indices         
      L1B_Scan_Metadata_t *L1B_Scan_Meta    contains the bit QA value for negative radiance
                                            beyond noise level. 
!Revision History:
 $Log: Emissive_Cal.c,v $
 Revision 1.13  2011-04-19 12:48:31-04  ltan
 Typo erro correction

 Revision 1.12  2011-04-07 14:38:15-04  xgeng
 TEB uncertainty algorithm update.

 Revision 1.11  2010-11-15 11:36:45-05  xgeng
 No calibration is performed for the scan if both sides of PCLW electronics are on at the same time.

 Revision 1.10  2005/01/18 21:57:33  ltan
 Added new file attributes prefixed with HDFEOS_FractionalOffset


 gevision 02.26, March 26, 2003  Razor Issue #191
 Change hard-coded band (28) whose EV data are saved for SWIR OOB correction to be 
   determined by a LUT instead.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.25, April 25, 2002  Razor Issue #180
 Moved checks on Dead detector, sector rotation, negative SV/BB average, and 
   negative b1 outside of Frame loop.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.24, March 25, 2001  Razor Issue #178
 Removed ADC Correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.23, February 2002  Razor Issue #180
 Inserted check on B1 to be sure it is not negative
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.22, November 9, 2001 Razor Issue #168
 Cleanup of spelling/grammar; changed TOLLERANCE to TOLERANCE
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.21, January 16, 2001
 Eliminated the "pow" functions in calculations of uncertainty.
 Also made other efficiency improvements and made variables a little
 more like the algorithm document.
 Changed frame index (F) to frame number in uncertainty calculation
 (minor bug in previous version)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.20 August 26, 1999
 Added checking if the EV data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.19 August 1999
 Removed interpolation correction.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.18 August 1999
 Implemented changes to meet the requirement of new ADC algorithms. See SDF.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov) 

 Revision 02.17 Aug 10, 1999
 Used L_Max and L_Min LUTs
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.16 May 23, 1999
 Many changes leading up to the May 30 delivery.  See SDF.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov) and Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.15  April 12, 1999
 Removed cubic term a3
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 02.12  Feb 9, 1999
 Added use of variables ADC_correction_switch, PCX_correction_switch and
 INT_correction_switch set to macro values to remove compiler warnings.
 Macros in upper case.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.11 July 1998
 modified EV, SI, UI to match L1B_Scan_t data structure change.
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Removed D_inv, which has already been taken care of in L1A
 Removed all the V vs L algorithm implementation.
 Added the DN vs. L algorithm implementation
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 02.00 Jan. 1997 
 Plugged in lookup tables and added Preprocess_data_t *PP as an input parameter;
 eliminated DED as an input parameter;
 expanded processing scope from one band to one scan.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.01 1996/04/05
 Update to match Version 1 Design Document
 Joan Baden(baden@highwire.gsfc.nasa.gov)
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Neal Devine(neal.devine@gsfc.nasa.gov)  -- debugging 1996/5/9-13

 Revision 01.00 1993
 Initial development
 Geir Kvaran(geir@highwire.gsfc.nasa.gov)

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
  int16     D            = 0; /* detector index within band (0,9) */
  int16     D_emiss      = 0; /* Emiss detector index (0,159) */
  int16     D_490        = 0; /* modis detector index (0,489) */
  int16     B            = 0; /* band index within L1A resolution (1km_Emiss res only here) */
  int16     B_emiss      = 0; /* Emiss band index (0,15) */
  int16     B_xt         = 0; /* PCX band index (32,36) */
  int16     B_38         = 0; /* MODIS band index (0,37) */
  int16     F            = 0; /* frame index (0 -> 1353) */
  int16     flag         = 0; /* for L < -NEdL */
  int16     F_shift      = 0; /* frame index in PCX */
  float32   L_ev         = 0; /* EV radiance */ 
  int16     DN_ev        = 0; /* must use signed integer (since MISSING FLAG is -1)*/
  float32   dn_ev        = 0; /* corrected Earth view DN */
  float32   dn_ev_31[10][EV_1km_FRAMES]; /* for PCX correction */
  float32   Fn           = 0;
  float32   dL_ev        = 0;

  /*
   * obsolete due to TEB UI algorithm update, 3/22/2011, Xu Geng
   */
  /*float32   Sigma_TEB_PV_resid_elec_sq; */ /* holds residual electronic crosstalk
                                                contribution to sigma^2 */
  /*float32   Sigma_TEB_ADC_sq;  */         /* holds RSB ADC contribution to sigma^2 */
  /*float32   Calibr_resid;      */         /* holds calibration polynomial fit 
                                               residual */
  /*float32   *calibr_resid_coeffs; */      /* temporary pointer to table for one 
                                               detector */
  /*float32   *si0, *si1;        */         /* temporary pointer to Ucoeff tables */
  /*float32   ri0, ri1;          */         /* summation variables for uncert coeff */

  float32   dn_PCX_corr = 0;            /* PC cross talk correction. Must be 
                                             initialized to 0. */
  float32   dtemp;                      /* temporary variable */

  uint16    bad_SI = 0;                 /* flag value for SI */
  boolean   bad_value = False;   
      
  int8      PCX_correction_switch;
  int16     uncert;

  int16 SV_frame_no = 0;     /* Frame number to use for SV RVS Correction */

  /* 
   * If the mirror side value is not 0 or 1, it is not a valid value.
   * It is most likely this scan is a missing scan and the following
   * calculation should not be proceeded.
   */

  if (MS != 0 && MS != 1)
    SMF_ERROR(MODIS_F_NOK, "Invalid mirror side value in Emissive_Cal()");

  /*
   * The Band chosen for the SWIR OOB correction may be any emissive band. 
   * Let "X" denote that band number in the following discussion. The values
   * of dn for band "X" are saved for use in "SWIR_out_of_band_correction".
   * It is assumed in that function that a negative value of dn_"X" is not
   * valid for making the SWIR OOB correction.  Thus, initialize the dn_"X"
   * to a negative value before the loop through individual pixels.  dn_"X"
   * will only be re-assigned if a value of dn is actually calculated.
   */

  for (D = 0; D < DETECTORS_PER_1KM_BAND; D++)
    for (F = 0; F < EV_1km_FRAMES; F++)
      L1B_Scan->dn_X[D][F] = -1.;
  
  /*Initialize extended indices, to be incremented within loops*/
  /*Note: band_26 is included in L1A_BANDS_AT_RES[INDEX_1000M_EMISS]*/
  B_38    = NUM_BANDS - NUM_EMISSIVE_BANDS - 1; 
  B_emiss = 0;
  D_emiss = 0;
  D_490   = NUM_REFLECTIVE_DETECTORS - 10; 
  SV_frame_no = emiss_tables->RVS_BB_SV_Frame_No[1];
 
  PCX_correction_switch = emiss_tables->PCX_correction_switch;

  for (B = 0; B < L1A_BANDS_AT_RES[INDEX_1000M_EMISS];  B++, B_38++)
  {
    /*Skip band_26*/
    if (B == MODIS_BAND26_INDEX_AT_RES)
    {
      D_490 += 10;   /* one band contains 10 detectors */
      continue;
    }
    /* Initialize the dn_ev_31 */
    if ( PCX_correction_switch == ON && B == BAND31)
    {
      for (D = 0; D < DETECT_PER_BAND_AT_RES[INDEX_1000M_EMISS]; D++)
      for (F = 0; F < EV_1km_FRAMES; F++)
        dn_ev_31[D][F] = 0;
    }
    for (D = 0; D < DETECT_PER_BAND_AT_RES[INDEX_1000M_EMISS]; 
                D++, D_emiss++, D_490++)
    {
      /*
       * obsolete due to TEB UI algorithm update, 3/22/2011, Xu Geng
       */
      /*
      Sigma_TEB_PV_resid_elec_sq = 
                         emiss_tables->Sigma_TEB_PV_resid_elec[D_emiss] *
                         emiss_tables->Sigma_TEB_PV_resid_elec[D_emiss];

      Sigma_TEB_ADC_sq = emiss_tables->Sigma_TEB_ADC[D_emiss] *
                         emiss_tables->Sigma_TEB_ADC[D_emiss];
      */

      /* Initialize the logical "bad value" signal: */
      bad_value = False;
            
        /* Check if the detector is dead */
              
        if (QA_tables->common_QA_tables.dead_detector[D_490] == 1)
        {
          bad_SI                                     = DEAD_DETECTOR_SI;    
          L1B_Gran->valid_pixels[B_38]              -= EV_1km_FRAMES;       
          L1B_Gran->dead_detector_pixels[B_38]      += EV_1km_FRAMES;       
          QA->QA_common.num_dead_detector_EV_data[D_490] += EV_1km_FRAMES;          
          bad_value = True;
        }

        /* Check if the instrument is in a sector rotation */

        else if (QA->QA_common.Sector_Rotation[S] == True)
        {
          bad_SI                                     = SECTOR_ROTATION_SI;
          L1B_Gran->valid_pixels[B_38]              -= EV_1km_FRAMES;
          L1B_Gran->bad_data_flag[B_38]              = 1;
          QA->QA_common.num_sector_rotation_EV_data[D_490] += EV_1km_FRAMES;
          bad_value = True;
        }

        /* Check if the instrument has both sides of PCLW electronics on at the same time */

        else if (QA->QA_common.Electronic_Anomaly[S] == True)
        {
          bad_SI                                     = UNABLE_CALIBRATE_SI;
          L1B_Gran->valid_pixels[B_38]              -= EV_1km_FRAMES;
          L1B_Gran->bad_data_flag[B_38]              = 1;
          bad_value = True;
        }

        /* Check if the zero point DN value is valid */

        else if ( PP_emiss->DN_sv[D_emiss][S] < 0)
        {
          bad_SI                                     = ZERO_POINT_DN_SI;
          L1B_Gran->valid_pixels[B_38]              -= EV_1km_FRAMES;
          L1B_Gran->bad_data_flag[B_38]              = 1;
          QA->QA_common.num_no_bg_DN_EV_data[D_490] += EV_1km_FRAMES;  
          bad_value = True;
        }
        
        /* Check to see that the b1 value is not negative. */        
        
        else if (PP_emiss->b1[D_emiss][S] < 0)
        {
          bad_SI                                     = TEB_B1_NOT_CALCULATED;
          L1B_Gran->valid_pixels[B_38]              -= EV_1km_FRAMES;
          L1B_Gran->bad_data_flag[B_38]              = 1;
          QA->QA_common.num_negative_b1[D_490]      += EV_1km_FRAMES;
          bad_value = True;
        }

        if (bad_value)
        {
          for (F = 0; F < EV_1km_FRAMES; F++)
          {
            L1B_Scan->SI.EV_1km_Emissive[B_emiss][D][F]    = (uint16) bad_SI;
            L1B_Scan->UI.EV_1km_Emissive_UI[B_emiss][D][F] = BAD_DATA_UI;
          }
          continue;
        }
                    
      for (F = 0; F < EV_1km_FRAMES; F++)
      {
	
        DN_ev = L1A_Scan->EV_1km_night[D][B][F]; 

        if (DN_ev < MISSING_L1A_FLAG || DN_ev > SATURATED_DN) {
          continue;
          /*
          Bad_L1A_Error_Out("EV_1km_night",
             " out of range in middle L1A file, Emissive_Cal(), Emissive_Cal.c");
          */
        }

        /* Skip missing data */

        if (DN_ev == MISSING_L1A_FLAG)
        {
          L1B_Scan->SI.EV_1km_Emissive[B_emiss][D][F] = L1A_DN_MISSING_SI;
          L1B_Scan->UI.EV_1km_Emissive_UI[B_emiss][D][F] = BAD_DATA_UI;
          L1B_Gran->missing_pixels[B_38]++;
          L1B_Gran->valid_pixels[B_38]--;
          L1B_Gran->bad_data_flag[B_38] = 1;
          QA->QA_common.num_missing_data_in_scans[D_490]++;
          continue;
        }
 
        /* Check if the EV DN read from L1A file is saturated */

        if (DN_ev == SATURATED_DN)
        {
          L1B_Scan->SI.EV_1km_Emissive[B_emiss][D][F] = SATURATED_DETECTOR_SI;
          L1B_Scan->UI.EV_1km_Emissive_UI[B_emiss][D][F] = BAD_DATA_UI;
          L1B_Gran->saturated_pixels[B_38]++;
          L1B_Gran->valid_pixels[B_38]--;
          L1B_Gran->bad_data_flag[B_38] = 1;
          QA->QA_common.num_saturated_EV_data[D_490]++;
          continue;
        }

        dn_ev = DN_ev - PP_emiss->DN_sv[D_emiss][S];
       
        /*
         * Save the value of dn for band X (any TEB band) to use for SWIR OOB
         * correction
         */

        if (B == L1B_Scan->band_X)
          L1B_Scan->dn_X[D][F] = dn_ev;

        /* check PCX Switch */
        if ( PCX_correction_switch == ON )
        {
          dn_PCX_corr = 0;
          if ( B == BAND31 )  dn_ev_31[D][F] = dn_ev;
          else if( B >= BAND32 && B <= BAND36 )
          {
            B_xt = B - BAND32;
            F_shift = F + emiss_tables->PC_XT[B_xt][D][1];
            if ( F_shift >= 0 && F_shift < EV_1km_FRAMES && 
                dn_ev_31[D][F_shift] > 0 )
            {
              dn_PCX_corr = -emiss_tables->PC_XT[B_xt][D][0]*0.01
                            *dn_ev_31[D][F_shift]
                            *emiss_tables->PC_XT[B_xt][D][2] + 
                            emiss_tables->PC_XT[B_xt][D][3];
              dn_ev += dn_PCX_corr;
            }
          }
        }


      	/* Compute L_ev */ 
        Fn = PP_emiss->a0[D_emiss][S] + 
             PP_emiss->b1[D_emiss][S] * dn_ev + 
             PP_emiss->a2[D_emiss][S] * dn_ev * dn_ev ;

        L_ev = (Fn - 
               (L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_SV[D_emiss][MS] - 
                L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS])
                * PP_emiss->Planck_mir[D_emiss][S] ) / 
                L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS];

        /* Check if the calculated radiance L exceeds the LUT L_Max.
         * If so, then set the unusable data value for SI and UI and
         * keep track of the number of instances of these pixels.
         */

        if (L_ev > emiss_tables->L_Max[B_emiss])
        {
          L1B_Scan->SI.EV_1km_Emissive[B_emiss][D][F] = TEB_OR_RSB_GT_MAX_SI;
          L1B_Scan->UI.EV_1km_Emissive_UI[B_emiss][D][F] = BAD_DATA_UI;
          L1B_Gran->valid_pixels[B_38]--;
          L1B_Gran->bad_data_flag[B_38] = 1;
          QA->QA_common.num_exceed_max_for_scaling[D_490]++;
          continue;
        }
 
        /* determine if the radiance is negative and below noise level */
        if (L_ev < -QA_tables->emiss_QA_tables.NEdL[D_emiss])
        {
          flag = 1; 
          L1B_Gran->negative_value_below_noise_pixels[B_38]++;
          L1B_Gran->valid_pixels[B_38]--;
          L1B_Gran->bad_data_flag[B_38] = 1;
        }

        /* Convert the product to scaled integer and store in the L1B scan */
       
        FLOAT_TO_PIXEL(L_ev,
                       emiss_tables->L_Max[B_emiss],
                       emiss_tables->L_Min[B_emiss],
                       L1B_Gran->SO.rad_scale_Emiss[B_emiss],
                       L1B_Gran->SO.rad_offset_Emiss[B_emiss],
                       L1B_Scan->SI.EV_1km_Emissive[B_emiss][D][F]);
        
        /* 
         * Check if NAD door is closed. If it is closed, set the most significant
         * bit to 1. If SI is larger than the upper limit for the NAD closed SI 
         * value, set it to that value.
         */
      
        if (QA->QA_common.NAD_Door_Open[S] == False)
        {
           L1B_Scan->SI.EV_1km_Emissive[B_emiss][D][F] |= 0x8000;
           QA->QA_common.num_nadir_door_closed_EV_data[D_490]++;         
  
           if (L1B_Scan->SI.EV_1km_Emissive[B_emiss][D][F] > NAD_CLOSED_UPPER_SI)
             L1B_Scan->SI.EV_1km_Emissive[B_emiss][D][F] = NAD_CLOSED_UPPER_SI;
        }

        /*
         * Set the uncertainty index.
         */

        if (L_ev <= TOLERANCE || QA->QA_common.NAD_Door_Open[S] == False
            || dn_ev <= TOLERANCE || PP_emiss->dn_bb[D_emiss][S] <= TOLERANCE )
        {
          L1B_Scan->UI.EV_1km_Emissive_UI[B_emiss][D][F] = BAD_DATA_UI;
        }
        else
        {

          /*
           * Compute uncertainty based on equation (3) in memo:
           * "PFM TEB Radiometric Uncertainty and LUT Format - Update 2",
           * Kwo-Fu (Vincent) Chiang, G. Godden and X. (Jack) Xiong, 
           * December 6, 1999.
           *********************************************************
           * The above algorithm is obsolete due to the recent update.
           * Please refer to the new algorithm in memo:
           * "TEB Radiometric Uncertainty and LUT Format",
           * by Kwo-Fu (Vincent) Chiang etc., March, 2011
           *********************************************************
           */

          dL_ev = 0.0;      /* use this variable to accumulate (dL_ev/L_ev)^2 */

          /***************************************************************
           * obsolete due to TEB UI algorithm update, 3/22/2011, Xu Geng
           */

          /*
           * Add in terms from Eq. (4).
           */

          /*
          frame_number = (float32) (F+1);
          for ( i = 0; i < NUM_UI_PARAMETERS; i++ )
          {
            si0 = emiss_tables->Ucoeff[D_emiss][i][0];
            si1 = emiss_tables->Ucoeff[D_emiss][i][1];
            EVAL_4TH_ORDER_POLYNOMIAL(ri0, si0, frame_number);
            EVAL_4TH_ORDER_POLYNOMIAL(ri1, si1, frame_number);
            dL_i_over_L_ev = (ri0 + ri1 * L_ev) / L_ev;
            dL_ev += (dL_i_over_L_ev * dL_i_over_L_ev);
          }
          */
 
          /*
           * Compute residual uncertainty to calibration fit (Eq. 9) and store
           * in variable "Calibr_resid".
           */
          /*
          calibr_resid_coeffs = emiss_tables->Ucoeff_Calibr_resid[D_emiss];
          EVAL_4TH_ORDER_POLYNOMIAL(Calibr_resid, calibr_resid_coeffs, L_ev);
          Calibr_resid /= L_ev;
          */

          /*
           * Add in remaining terms to (dL_ev/L_ev)^2
           */
          /*
          dtemp = PP_emiss->NEdL[D_emiss][S]/L_ev;
          dL_ev += (dtemp * dtemp) + Sigma_TEB_PV_resid_elec_sq +
                   (0.25 * dn_PCX_corr * dn_PCX_corr + Sigma_TEB_ADC_sq) / 
                   (dn_ev * dn_ev) + Calibr_resid * Calibr_resid;
          if (B_emiss == BAND21)
          {
            dtemp = emiss_tables->Band_21_Uncert_Lsat * L_ev;
            dL_ev += (dtemp * dtemp);
          }
          *******************************************************************/

          /*
           * Uncertainty due to a0
           */
          dtemp = PP_emiss->sigma_a0[D_emiss][S]*(1. - dn_ev/PP_emiss->dn_bb[D_emiss][S])
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to a2
           */
          dtemp = PP_emiss->sigma_a2[D_emiss][S]*(dn_ev-PP_emiss->dn_bb[D_emiss][S])*dn_ev
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to RVS_SV
           */
          dtemp = L1B_Gran->Emiss_Cal_Coeff.sigma_RVS_Emiss_EV[D_emiss][SV_frame_no][MS]*
                  (dn_ev/PP_emiss->dn_bb[D_emiss][S]-1.)*PP_emiss->Planck_mir[D_emiss][S]
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to RVS_EV
           */
          dtemp = L1B_Gran->Emiss_Cal_Coeff.sigma_RVS_Emiss_EV[D_emiss][F][MS]*
                  (PP_emiss->Planck_mir[D_emiss][S]/L_ev - 1.)
                  /L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS];
          dL_ev += (dtemp * dtemp);


          /*
           * Uncertainty due to epsilon_bb
           */
          dtemp = emiss_tables->sigma_epsilon_BB[B_emiss]*(PP_emiss->L_bb[D_emiss][S] -
                  emiss_tables->epsilon_cav[D_emiss]*PP_emiss->L_cav[D_emiss][S])*dn_ev
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev*
                  PP_emiss->dn_bb[D_emiss][S]);
          dL_ev += (dtemp * dtemp);
          /*
           * Uncertainty due to epsilon_cav
           */
          dtemp = emiss_tables->sigma_epsilon_CAV[B_emiss]*(1.-emiss_tables->epsilon_bb[D_emiss])
                  *PP_emiss->L_cav[D_emiss][S]*dn_ev
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev*
                  PP_emiss->dn_bb[D_emiss][S]);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to lambda
           */
          dtemp = (emiss_tables->sigma_L_lambda[B_emiss][0]+emiss_tables->sigma_L_lambda[B_emiss][1]*L_ev)
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to T_bb
           */
          dtemp = emiss_tables->sigma_L_Tbb[B_emiss]*emiss_tables->epsilon_bb[D_emiss]*dn_ev
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev*
                  PP_emiss->dn_bb[D_emiss][S]);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to T_sm
           */
          dtemp = emiss_tables->sigma_L_Tsm[B_emiss]*
                  (L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]-
                   L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_SV[D_emiss][MS]+
                   (L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_SV[D_emiss][MS] - 1)
                   *dn_ev/PP_emiss->dn_bb[D_emiss][S])
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to T_cav
           */
          dtemp = emiss_tables->sigma_L_Tcav[B_emiss]*(1.-emiss_tables->epsilon_bb[D_emiss])
                  *emiss_tables->epsilon_cav[D_emiss]*dn_ev/PP_emiss->dn_bb[D_emiss][S]
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to dn_EV
           */
          dtemp = PP_emiss->dn_bb_sdev[D_emiss][S]*(PP_emiss->b1[D_emiss][S]+
                  2.*PP_emiss->a2[D_emiss][S]*dn_ev)
                  /(L1B_Gran->Emiss_Cal_Coeff.RVS_1km_Emiss_EV[D_emiss][F][MS]*L_ev);
          dL_ev += (dtemp * dtemp);

          /*
           * Uncertainty due to band_21_b1
           */
          if (B_emiss == BAND21)
          {
            dtemp = emiss_tables->sigma_b1_B21[D][MS];
            dL_ev += (dtemp * dtemp);
          }

          /*
           * Uncertainty due to PC cross talk correction
           */
          if ( PCX_correction_switch == ON &&  B >= BAND32 && B <= BAND36 )
          {
            dtemp = emiss_tables->pcx_ui_factor[B_xt]*(dn_PCX_corr / dn_ev);
            dL_ev += (dtemp * dtemp);
          }

          /*
           * Compute fractional uncertainty (take square root of (dL_ev/L_ev)^2)
           */

          dL_ev = sqrt((double) dL_ev);

          /*
           * Calculate the uncertainty index using the formula from the
           * L1B Product User's Guide:
           *
           *            UI = n * log(percent_uncertainty / m)
           *
           * The value of n and m are band dependent and come from lookup tables.
           * The factor of 100 in the formula below converts the fractional
           * uncertainty to % uncertainty.
           */

          if (dL_ev > 0)
            uncert = emiss_tables->TEB_UI_scaling_factor[B_emiss] *
                     log((double) (100.0*dL_ev/emiss_tables->TEB_specified_uncertainty[B_emiss]));
          else
            uncert = BAD_DATA_UI;

          /* 
           * Ensure the uncertainty is in the range of 0-15. 
           */

          if (uncert > BAD_DATA_UI)
            L1B_Scan->UI.EV_1km_Emissive_UI[B_emiss][D][F] = BAD_DATA_UI;
          else if(uncert < 0)
            L1B_Scan->UI.EV_1km_Emissive_UI[B_emiss][D][F] = 0;
          else
            L1B_Scan->UI.EV_1km_Emissive_UI[B_emiss][D][F] = (uint8) uncert;
        }

      }/*End of loop over frames*/ 

    }/*End of loop over detectors within band*/ 

    B_emiss++;
  }/*End of loop over bands*/
  
  if (flag == 1)
    L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000008;

  return(MODIS_S_OK);
}


