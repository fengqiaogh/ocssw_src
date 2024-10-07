#ifndef	PREPROCESS_H
#define	PREPROCESS_H

#include    "L1B_Tables.h" /* contained variable definitions used here */

/*
!C-INC********************************************************************** 
!Description:      Header file Preprocess.h to be included in files 
                   where the preprocess routines and data structures 
                   are used. 
!Revision History:
 $Log: Preprocess.h,v $
 Revision 1.17  2016-07-22 13:23:32-04  xgeng
 Applied xtalk correction to Terra LWIR PV bands (27-30)

 Revision 1.16  2011-04-07 14:41:30-04  xgeng
 1. RSB &TEB uncertainty algorithm update; 2. The quadratic RSB RVS changed to 4th order.

 Revision 1.15  2006-10-30 10:00:13-05  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.

 Revision 02.20 June 28, 2002    Razor Issue #161
 Moved NUM_*_SUBSAMP definitions to Granule.h
 Gwyn Fireman, SAIC-GSO <Gwyn.Fireman@gsfc.nasa.gov>

 Revision 02.19  April 15, 2002  Razor Issue #182
 Applicable to MODIS/TERRA (PFM) processing only.
 Removed macro BAND_26_REFL_INDEX for Band 26 crosstalk correction
   and moved to L1B_SetupP.h
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.18, March 25, 2001  Razor Issue #178
 Removed ADC Correction passed parameters
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 02.17  March 20, 2002  Razor Issue #182
 Applicable to MODIS/TERRA (PFM) processing only.
 Added macro BAND_26_REFL_INDEX for Band 26 crosstalk correction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.16  March 11, 2002  Razor Issue #174
 Added L1B_Gran to Preprocess.c call
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.15, Dec. 7, 2000, Razor issue 146
 Removed DN_sat_prime arrays.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.14 August 23, 1999
 Removed DELTA_bb and dn_bb from Preprocess_Emiss_t.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.13 May 7, 1999
 Moved DN_to_DN_star to Reflective_Cal, because it is not used in Preporcess
 anymore
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 02.12 Feb 17, 1999
 Moved external function DN_to_DN_star into this module from Reflective_Cal to
 avoid circular include dependencies.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.11 Oct, 1998
 Added the prototype for Get_Leading_Gran_Temperatures(), Get_Middle_Gran_Temperatures(),
 Get_Trailing_and_Average_Temperatures(), Calculate_PP_Planck_Mir(), Get_Emissive_Bands_ADC_Board(),
 Get_Leading_Gran_Emiss_Coeff(), Get_Middle_Gran_Emiss_Coeff(), Get_Trailing_Gran_Emiss_Coeff(),
 Get_All_Emiss_Coeff()
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)
 
 Revision 02.10 June 2 , 1998
 Modified the prototype for Preprocess_L1A_Data(), Write_L1B_OBCEng() 
 ADDed the Calculate_SD_DN_star prototype.
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)
  
 Revision 02.10 Apr. 14, 1998
 Removed all old unused v2.0 prototype interfaces (e.g. cleanup.., malloc.. )
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Apr. 7, 1998
 Added MirrorSide to the Overlap_OBCEng_t
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Apr. 7, 1998
 Added the Get_Temps_From_Eng_Data() prototype.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Rewrite the L1B TEB algorithm using DN vs L
 Split the Preprocess_Data_t into Preprocess_Refl_t and Preprocess_Emiss_t
 Add the QA_Emiss_t to handle the new TEB QA data
 Rewrite the Temperatures_t to include more temperatures and lwir(mwir)_adc
 Add macro VAR_Calculation(..)
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)
 
 Revision 01.10 Aug. 1997
 Inlined DN_to_Voltage().
 Added Temperatures_t and the interface prototypes. 
 Zhidong Hao (hao@acrobat.gsfc.nasa.gov)

 Revision 01.00 Jan. 1997
 Initial development
 Zhidong Hao (hao@acrobat.gsfc.nasa.gov)

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END***********************************************************************
*/

/* If moon vector is bad value (0, 0, 0), set moon_in_SV_KOB to be -1 */
#define MOON_INSIDE_SV_KOB           1
#define MOON_OUTSIDE_SV_KOB          0
#define UNDETERMINED_MOON_IN_SV_KOB -1

typedef struct {
  float32 DN_obc_250m_avg[MAX_250M_TRACK_DIM][NUM_250M_BANDS][NUM_250M_SUBSAMP];
  float32 DN_obc_500m_avg[MAX_500M_TRACK_DIM][NUM_500M_BANDS][NUM_500M_SUBSAMP];
  float32 DN_obc_1km_day_avg[MAX_1KM_TRACK_DIM][NUM_1000M_DAY_BANDS][NUM_1KM_SUBSAMP];
  float32 DN_obc_1km_night_avg[MAX_1KM_TRACK_DIM][NUM_1000M_NIGHT_BANDS][NUM_1KM_SUBSAMP];
  float32 DN_obc_250m_var[MAX_250M_TRACK_DIM][NUM_250M_BANDS][NUM_250M_SUBSAMP];
  float32 DN_obc_500m_var[MAX_500M_TRACK_DIM][NUM_500M_BANDS][NUM_500M_SUBSAMP];
  float32 DN_obc_1km_day_var[MAX_1KM_TRACK_DIM][NUM_1000M_DAY_BANDS][NUM_1KM_SUBSAMP];
  float32 DN_obc_1km_night_var[MAX_1KM_TRACK_DIM][NUM_1000M_NIGHT_BANDS][NUM_1KM_SUBSAMP];
  uint32 DN_obc_250m_outlier_mask[MAX_250M_TRACK_DIM][NUM_250M_BANDS][NUM_250M_SUBSAMP][2];
  uint32 DN_obc_500m_outlier_mask[MAX_500M_TRACK_DIM][NUM_500M_BANDS][NUM_500M_SUBSAMP][2];
  uint32 DN_obc_1km_day_outlier_mask[MAX_1KM_TRACK_DIM][NUM_1000M_DAY_BANDS][NUM_1KM_SUBSAMP][2];
  uint32 DN_obc_1km_night_outlier_mask[MAX_1KM_TRACK_DIM][NUM_1000M_NIGHT_BANDS][NUM_1KM_SUBSAMP][2]; 
}DN_OBC_Avg_t;

typedef struct {
  float32 var_T_ins; /* holds variance of T_ins over the granule */
} Preprocess_Refl_t;

typedef struct {

  float32   Planck_mir[NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS]; 
  float32   a0        [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS];
  float32   sigma_a0  [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS];
  float32   b1        [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS];
  float32   a2        [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS];
  float32   sigma_a2  [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS];
  float32   DN_sv     [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS]; 
  float32   dn_bb     [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS]; 
  float32   dn_bb_sdev[NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS]; 
  float32   L_bb      [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS]; 
  float32   L_cav     [NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS]; 
  float32   NEdL[NUM_EMISSIVE_DETECTORS][MAX_NUM_SCANS];
  /* Temperatures running averaged XGranule
  */
  float32  T_fp       [NUM_FOCAL_PLANES][MAX_NUM_SCANS];    
  float32  T_bb                        [MAX_NUM_SCANS];
  float32  T_mir                       [MAX_NUM_SCANS];
  float32  T_ins                       [MAX_NUM_SCANS];
 

  /* Temperatures averaged over a granule, used in metadata
  */
  float32 T_bb_gran;
  float32 T_mir_gran;
  float32 T_fp1;
  float32 T_fp2;
  float32 T_fp3;
  float32 T_fp4;
  
  /* Focal_Plane_Set_Point_State
  */
  int8   fp_set_point_state[MAX_NUM_SCANS];

} Preprocess_Emiss_t;

typedef struct {
  DN_OBC_Avg_t        DN_OBC_Avg;
  Preprocess_Refl_t   PP_Refl;
  Preprocess_Emiss_t  PP_Emiss;
} Preprocess_Data_t;

/* Bands to apply the pvlw xtalk correction */
extern int16 PVLW_XT_BANDS[NUM_PVLW_XT_BANDS];

/* The relative position to reference optical axis for Bands 27 to 30 on the LWIR FPA */
extern int16 PVLW_XT_BANDS_FRAME_REL_POS[NUM_PVLW_XT_BANDS];

/*
 * Preprocess_L1A_Data -- called from main().
 */

PGSt_SMF_status  Preprocess_L1A_Data (lookup_tables_t     *tables,
                                      L1A_granule_t       *L1A_Gran, 
                                      L1B_granule_t       *L1B_Gran,
                                      QA_Data_t           *QA,
                                      Preprocess_Data_t   *PP); 


#endif














