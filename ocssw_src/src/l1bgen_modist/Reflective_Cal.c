#include    "Reflective_Cal.h"
#include    "Reflective_CalP.h"
#include    "L1B_Tables.h"
#include    "HDF_Lib.h"
#include    "PGS_Error_Codes.h"
#include    "FNames.h"
#include    <math.h>

/* The following value is currently the largest uncertainty index value,
 * for the reflective bands. If an index is calculated to be
 * larger than this value, then it is set to the bad data value.
 */
#define MAX_UNCERTAINTY_UI          15
static char *invalidlutfile =
               "This is most likely due to an invalid LUT file.";

extern int16 RFLAG;
extern int16 RSCL_FLAG;

PGSt_SMF_status Reflective_Cal (int16               S, 
                                L1A_granule_t       *L1A_Gran,
                                L1B_granule_t       *L1B_Gran,
                                L1A_Scan_t          *L1A_Scan,
                                L1B_Scan_t          *L1B_Scan,
                                Preprocess_Data_t   *PP,
                                refl_tables_t       *refl_tables,
                                common_QA_tables_t  *QA_tables,
                                QA_Common_t         *QA_Common)

/*
!C***************************************************************************** 
!Description:

 This module corrects the raw digital signals, DN, for known instrumental
 effects to produce corrected digital signals, dn* for every scan, frame,
 subframe, band and detector.  Corrections are applied for the effect of
 instrument and focal plane temperature on detector responsivity,  for the
 electronic background,  for the angular dependence of the response of the
 scan mirror, for non-linearities in the Analog to Digital Converters,  and
 for the effect of an out-of-band spectral leak in the SWIR bands 5, 6 and 7.
 Detectors within each spectral band are placed on a common scale by
 scaling dn* by relative calibration coefficients of the detectors in each
 band, to produce dn**.  Each final dn** value, if valid, is scaled to an
 integer in the range of [0-32767] and placed in an unsigned, 16-bit integer
 variable (which has a full range of [0-65535]).  If a valid value cannot be
 computed, the scaled integer is set to a value in the range of [32768-65535]. 
 Specific values in the range of [32768-65535] are used to denote why a valid
 value could not be obtained (a list of these is in the L1B file
 specifications).  This routine also computes the uncertainty in the
 reflectance product for every pixel, and converts the uncertainty to a
 4-bit uncertainty index, stored in the 4 least significant bits of an 8-bit
 unsigned integer.

!Input Parameters:

 int16               S               current scan index
 L1A_granule_t       *L1A_Gran       contains mirror side index of current scan
 L1B_granule_t       *L1B_Gran       contains calibration coefficients and scale, offset
                                     factors 
 L1A_Scan_t          *L1A_Scan       contains earth view DN values
 Preprocess_Data_t   *PP             contains focal plane and instrument temperatures,
                                     ADC indexes, average electronic background DNs,
                                     saturation DN values and scan-dependent parts of
                                     swir correction algorithms
 refl_tables_t       *refl_tables    contains the reflective lookup table values used 
                                     in reflective corrections
 common_QA_tables_t  *QA_tables      contains the dead detector lookup table value
 QA_Common_t         *QA_Common      contains NAD door open and sector rotation information

!Output Parameters:

 L1B_Scan_t          *L1B_Scan       contains scaled integers and uncertainty indices
 L1B_granule_t       *L1B_Gran       contains some granule level QA values

!Revision History:

 Revision 02.33  October 15, 2004  Razor Issue #199
 Changed the logic of setting the SWIR correction band detector index 
 for the SWIR band detector based on a new "swir_oob_sending_detector_table".
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 02.32 March 27, 2003, Razor Issue #173
 Initialized some variables for ANSI-C compliance.
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 02.31  March 26, 2003   Razor Issue #191
 Change use of Band 28 to band designated by Reflective LUT
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.30, Oct. 3, 2002  Razor Issue #187
 Remove R_star from code and replace it by 1/(m1*E_sun_over_pi) 
 Liqin Tan(ltan@saicmodis.com) 

 Revision 02.29, March 25, 2001  Razor Issue #178
 Removed ADC Correction
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 02.28, March 16, 2001, Razor issue #159
 Removed check on saturation of DN_ev using DN_sat LUT and replaced it with a
 check on saturation of dn_ev using a new LUT dn_sat_ev.  See MCST Memo M0991,
 "Update to Level 1B Handling of Reflective Solar Band Detector Electronics
 Saturation", G. Fireman and F. Adimi
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov)
 
 Revision 02.27, Dec 7, 2000, Razor issue 146, 147
 New SWIR algorithm, per Jack Xiong's memo.
 In-lined the SWIR correction algorithm since it is now only a couple of lines.
 Removed obsolete DN_sat_prime.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
  
 Revision 02.26 Oct. 6, 1999
 Corrected uncertainty calculation (missing the 100 for %) and restructured
 uncertainty calculations to make more understandable.  Added additional
 comments and refined comments for clarity.  Incorporated dn** nomenclature.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.25 Sept. 29, 1999
 Corrected the implementation of SWIR band correction to be on dn, not dn*.
 Changed the name of "SWIR_correction" to "SWIR_out_of_band_correction"
 as per Bruce B's request. Removed UI_ptr from argument list since we
 expect new SWIR uncertainty algorithm to require new module.
 Added and clarified documentation of the function.
 Rearranged some calculations to better match algorithm.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.24 Sept. 2, 1999
 Implemented changes to meet the requirements of the new SWIR algorithm. See SDF.
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov) 

 Revision 02.23 August 26, 1999
 Added checking if ev data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.22 August 1999
 Implemented changes to meet the requirement of new ADC algorithms. See SDF.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov) 

 Revision 02.21 May 23, 1999
 Many changes leading up to the May 30 delivery.  See SDF.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov) 
 
 Revision 02.11 Feb 22, 1999
 Changed the call to Scaled_Int_to_Radiance to a simple statement, using the
 new logic of how to calculate radiance from scaled integer.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 April, 15 1998
 Changed Preprocess_Refl_t *PP_refl paramter to Preprocess_Data_t *PP
 (since some temperatures in the PP->Preprocess_Emiss structure are needed).
 Changed the DN_to_DN_star() parameters for the new algorithm.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 April, 13 1998
 Removed the detector number inversion stuff (D_inv). The L1A
 data already is in inverted det order - no need to do it again.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 April, 8 1998
 Added bounds checking logic to the uncertainty (S_L) to 
 uncertainty_index conversion algorithm. 
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 April, 6 1998
 Got rid of the FLOAT_TO_PIXEL Macro call, replacing it with a 
 Float_To_Pixel() function call.
 Changed the rad_scale and rad_offset to counts_scale and counts_offset,
 and L_ev to DN_star in the Float_To_Pixel() parameter list since we want
 a scaled int of DN_star and not of L_ev.
 Added the get_uncertainties() call to calculate uncertainties.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Added the SWIR_correction(). Including adding:
 momo parameter, SWIR_correction_table to the refl_tables_t
 Changed Preprocess_Data_t *PP to Preprocess_Refl_t *PP_refl
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.00 Dec. 1996
 Plugged in lookup tables, trending_predictions, and 
 Preprocess_data_t *PP as an input parameter;
 eliminated DED as an input parameter;
 expanded processing scope from one band to one scan.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.01 1996/04/05
 Update to match Version 1 Design Document
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Joan Baden (baden@highwire.gsfc.nasa.gov)

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

 Summary of important symbology in this function:

 DN   = raw digital number, obtained directly from the Level 1A file.
 dn   = DN corrected for electronic zero point value (which usually
        comes from the averaged SV DN). 
 dn*  = dn corrected for instrument effects such as RVS and instrument 
        temperature.
 dn** = dn* placed on a common scale (one per band) by accounting for the
        relative calibration coefficients of the detectors in each band.
 SI   = scaled integer
 UI   = uncertainty index

 The calculations needed to calculate dn** and the uncertainty in the
 reflectance product are performed at four major levels of indexing: 
 (1) by L1A band groups, (2) by band within a group, (3) by detector
 within a band and (4) by frame within the resolution of a band group
 (an index that includes both 1km-frame and sub-frame).  Some
 calculations related to dn* but not dependent on the 1km-frame index
 are performed at level 3 to maintain efficiency.

 The following indicates the general flow of calculations:

 LOOP through L1A band groups
    LOOP through bands of this group
       - if thermal emissive band, skip to end of this loop
       - set various band-dependent indexes and quantities
       LOOP through detectors for this band
          - calculate and store some quantities independent of 1km-frame
          LOOP through frames at the resolution of this band group
             - extract values of pixel quantities from arrays
               that have band-group-dependent variable names
               (includes extracting addresses of where to assign
               scaled integer and uncertainty index)
             - check quantities to ensure that pixel can be calibrated
             - calculate dn
             - if SWIR band, apply a correction to dn to account for the
               out-of-band (OOB) leak
             - calculate dn* by applying RVS and temperature corrections
             - calculate dn** by applying relative calibration coefficient
               factor.
             - scale dn** to scaled integer
             - calculate reflectance uncertainty index
          END LOOP through frames at this resolution
       END LOOP through detectors for this band
    END LOOP through bands of this group
 END LOOP through L1A band groups

 Note that the detector order convention for all variables and indexes
 within variables is the "product" convention.

 The 38 MODIS bands are:
         1,2,3,4,5,6,7,8,9,10,11,12,13lo,13hi,14lo,14hi,
         15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,
         30,31,32,33,34,35,36
 The set of 490 MODIS detectors follows the band ordering indicated above
 with detectors for each band being inserted into the set (i.e., the
 first 40 detectors of the set are from band 1, the next 40 are from band 2,
 the next 20 are from band 3, and so on).

 The 22 reflective solar bands are MODIS bands:
         1,2,3,4,5,6,7,8,9,10,11,12,13lo,13hi,14lo,14hi,
         15,16,17,18,19,26
 The set of 330 reflective-band detectors follows a similar logic as that
 described above for the MODIS-detectors.

 Another band group used in this function consists of the SWIR bands,
 which are a subset of the reflective solar bands.  The SWIR bands are:
         5,6,7,26

 For an explanation of the "L1A" band groups and the "L1B" band groups,
 see the design notes in the file "Granule.h"

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char *location = "Reflective_Cal";

  /* Major loop indices */
  int16     R;                /* L1A band group index (0,3) */
  int16     B;                /* band index within one L1A band group */
  int16     D;                /* detector index within band */ 
  int16     F;                /* frame index within resolution of the band group */

  /* Other indices (see comments in design notes above) */
  int16     L1B_B;            /* band index within the L1B band group */
  int16     B_38;             /* band index within set of all MODIS bands (0-37) */
  int16     B_refl;           /* band index within set of reflective bands (0-21) */
  int16     D_refl;           /* detector index within set of refl. bands (0-329) */
  int16     D_490;            /* detector index within set of MODIS bands (0-489) */
  int16     B_swir;           /* SWIR band index (0-3) */
  int32     T;                /* track index (combined scan and detector) */
  int16     MS;               /* index of mirror side for this scan (0,1) */
  int16     subsamp;          /* sub-sample (or sub-frame) index */
  int16     focal_plane;      /* focal plane assembly index for a given band */

  /* Variables used in calculation of scaled integer and uncertainty index */
  int16     DN_ev;            /* earth view DN */
  float32   dn_ev;            /* DN with background subtracted. */
  float32   dn_star;          /* dn corrected for known instrument effects */
  float32   dn_star_star;     /* dn* adjusted for detector variablility within a
                                 band of the relative calibration coefficients. */
  float32   RVS;              /* Response vs. Scan Angle correction coefficient */
  float32   DN_obc_avg;             /* averaged electronics background DN to be
                                       subtract from DN_ev */ 
  float32   dT_inst;                /* difference in T_inst from reference */
  float32   dT_fpa;                 /* difference in T_fp from reference */
  float32   temp_correction[4];     /* holds correction to dn_ev */
  boolean   temp_correction_valid;  /* flags invalid T_ins or T_fp */
  /**************************************************************/
  /* obsolete due to UI algorithm update, 2/19/2011 by Xu Geng */   
  /*float32   sigma_RVS;           */    /* Uncertainty in RVS (from tables) */
  /*float32   sigma_m1_norm_sq[4]; */    /* holds m1 contribution to sigma^2 */
  /*float32   sigma_RVS_over_RVS;  */    /* temporary variable used in uncertainty */
  /*float32   sigma_Kinst_norm_sq[4]; */ /* holds Kinst contribution to sigma^2 */
  /*float32   sigma_PV_resid_elec_sq[4];*/  /* holds residual electronic crosstalk 
                                               contribution to sigma^2 */
  /*float32   RSB_NEdL_R_star_sq[4];*/   /* holds (NEdL * R_star)^2 */
  /*float32   sigma_R_EV;           */   /* deviation of R* linear fit residual contribution
                                            to sigma^2 */ 
  /*float32   sigma_Tinst_norm_sq[4]; */ /* holds Tinst contribution to sigma^2 */
  /*float32   sigma_RSB_ADC_sq = 0.0; */ /* holds RSB ADC contribution to sigma^2 */
  /**************************************************************/
  uint16    SI;                     /* scaled integer */
  uint16    *SI_ptr;                /* address of SI variable to set */
  uint8     *UI_ptr;                /* address of UI variable to set */
  float32   uncertainty;            /* fractional uncertainty in reflectance */
  float32   u2, u4;                 /* uncertainty terms */          
  float32   *u4_coeffs;             /* coefficients read from LUTs to calculate u4 */ 
  int32     ui;                     /* temporary uncertainty index value */
  float32   dn_swir_corr;           /* swir correction for dn_ev */
  boolean   dn_swir_corr_error;     /* flag for error in swir correction for dn_ev */
  float32   sigma_dn_swir;          /* holds part of uncertainty term */
  float32   reflectance;

  /* Variables requiring initialization (these remain unchanged in function) */

            /* mapping of reflective band index (B_refl) to FPA index */

  int16     fpa_index[NUM_REFLECTIVE_BANDS] = {
   /* band = 1   2   3   4   5    6    7   8   9  10  11  */
            NIR,NIR,VIS,VIS,SWIR,SWIR,SWIR,VIS,VIS,VIS,VIS,
   /* band = 12 13L 13H 14L 14H  15  16  17  18  19  26  */
            VIS,NIR,NIR,NIR,NIR,NIR,NIR,NIR,NIR,NIR,SWIR };

            /* Definition of number of subsamples at each resolution */

  int16     samples_per_frame_at_res[NUM_L1A_RESOLUTIONS] = {4,2,1,1};

  int16   bX_det_index = 0;   /* detector index for SWIR OOB Correction band */
  int16   bX_frame_index = 0; /* 1km frame index for SWIR OOB Correction band */
  float32 dn_X = -1.0;        /* dn of SWIR OOB Correction band pixel. */

    /* # subsamples in a swir band */
  int16   num_swir_subsamples[NUM_SWIR_BANDS] = {2, 2, 2, 1}; 
    /* ratio of # detectors to 1km num */
  int16   swir_detector_ratio[NUM_SWIR_BANDS] = {2, 2, 2, 1}; 
  SWIR_correction_tables_t *swir_tables = &refl_tables->SWIR_correction_tables;

  int16   saturated_dn[NUM_REFLECTIVE_BANDS] =
    {4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 4095, 32767,
     4095, 32767, 32767, 4095, 4095, 4095, 4095, 32767, 4095, 4095,
     4095, 4095};

  /************* END DECLARATIONS, BEGIN CODE STATEMENTS ********************/

  /*Initialize extended indices, to be incremented within loops*/

  D_refl = 0;
  D_490  = 0;
  B_refl = 0;
  B_38   = 0;

  /* Initialize the mirror side index for this scan. */

  MS = L1A_Gran->MirrorSide[S];

  /* 
   * If the mirror side is not 0 or 1, it is not a valid value. It is most
   * likely this is a missing scan. The calculation below should not be
   * proceeded.
   */
  
  if (MS != 0 && MS != 1)
    SMF_ERROR(MODIS_F_NOK, "Invalid mirror side value in Reflective_Cal()");

  /* If no rescaling then adjust saturated limits */
  /* -------------------------------------------- */
  if ((RSCL_FLAG & 1) == 0) saturated_dn[12] = 4095;
  if ((RSCL_FLAG & 2) == 0) saturated_dn[17] = 4095;
  if ((RSCL_FLAG & 4) == 0) saturated_dn[9]  = 4095;
  if ((RSCL_FLAG & 8) == 0) saturated_dn[11] = 4095;

  /*
   * Loop through L1A band groups.
   */

  for (R = 0; R < NUM_L1A_RESOLUTIONS; R++)
  {
    if ((RFLAG & (1 << R)) != 0) {
      B_refl += L1A_BANDS_AT_RES[R];
      B_38 += L1A_BANDS_AT_RES[R];
      D_490 += DETECT_PER_BAND_AT_RES[R]*L1A_BANDS_AT_RES[R];
      continue;
    }

    /*
     * Loop through bands within this L1A band group.  Note that 
     * we also increment the MODIS band index at the end of each loop.
     */

    for (B = 0; B < L1A_BANDS_AT_RES[R]; B++, B_38++)
    {

      /*
       * Skip thermal emissive bands (but, don't forget to increment D_490)
       */
   
      if (R == INDEX_1000M_EMISS && B != MODIS_BAND26_INDEX_AT_RES)
      {
        D_490 += DETECT_PER_BAND_AT_RES[R];
        continue;
      }

      /* 
       * Set various band-dependent (not detector and frame dependent) indexes 
       * and quantities. To be efficient, don't put these into inner loop. 
       */
  
      /* Determine band index within the L1B band group.  Only MODIS band 26
       * needs special treatment.  MODIS Band 26 is in the 4th and last L1A band
       * group but it gets placed as the last band in the 3rd L1B band group.
       */
   
      if (R == INDEX_1000M_EMISS)                   /* Last L1A band group */
      {
        L1B_B = L1B_BANDS_AT_RES[INDEX_1000M_REFL] - 1;
      }
      else                                          /* Any other L1A band group */
      {
        L1B_B = B;
      }

      /* obsolete, 2/19/2011 by Xu Geng 
       * Get the value of sigma RVS from tables.
       */
      /* sigma_RVS = refl_tables->Sigma_RVS_RSB[B_refl][MS]; */

      /*
       * Select index of focal plane array corresponding to this reflective
       * band.  Will be used later when calculating temperature corrections
       * for dn*.
       */

      focal_plane = fpa_index[B_refl];

      /*
       * Determine and set the SWIR band index.  If this band is not a SWIR
       * band, set -1 as a flag to know that the SWIR band correction is not
       * to be applied.
       */
   
      switch ( B_38 ) {
        case 4:             /* MODIS band 5 */
          B_swir = 0;
          break;

        case 5:             /* MODIS band 6 */
          B_swir = 1;
          break;

        case 6:             /* MODIS band 7 */
          B_swir = 2;
          break;

        case 27:            /* MODIS band 26 */
          B_swir = 3;
          break;

        default:            /* Any other reflective band that is not a SWIR band */
          B_swir = -1;
          break;

      } /* end switch */

      /*
       * Loop through Detectors for this band.  Note that we also increment
       * the reflective detector index and the MODIS detector index, D_490.
       */

      for (D = 0; D < DETECT_PER_BAND_AT_RES[R]; D++, D_refl++, D_490++)
      {


  /*
   * Calculate and store some quantities independent of 1km-frame.
   * (sub-frame quantities stored in arrays over sub-frame) 
   */

  /*
   * Compute temperature correction to dn_ev.   The correction factor is a
   * multiplicative one and is:
   * 1 + instrument temperature contribution + FPA temperature
   * contribution See "Level 1B Uncertainty Calculations for Reflective
   * Solar Bands" by Bruce Berriman and Gail Reichert, MCST Document #
   * M0703, May 27, 1999 for algorithm. 
   */

  /*
   * Use "temp_correction_valid" as a flag to denote that both instrument
   * temperature and FPA temperatures were valid (> 0).  If either of
   * these are not valid, the respective correction contribution is not
   * included in the correction factor (implying that there is no
   * correction to dn_ev). However, the uncertainty index for this pixel
   * is set to the bad data value.
   */

        temp_correction_valid = True;       /* Assume for now. */

        /*
         * First, set the temp_correction values to 1.0.
         */

        for (subsamp = 0; subsamp < SUBSAMPLES_AT_RES[R]; subsamp++) 
          temp_correction[subsamp] = 1.0;

  /*
   * Add in instrument temperature correction contribution if 
   * temperature is valid.
   */

        if (PP->PP_Emiss.T_ins[S] > 0)
        {
          dT_inst = PP->PP_Emiss.T_ins[S] - refl_tables->T_inst_ref;
          for (subsamp = 0; subsamp < SUBSAMPLES_AT_RES[R]; subsamp++) 
            temp_correction[subsamp] +=
                refl_tables->K_inst[B_refl][D][subsamp][MS] * dT_inst;
        }
        else
        {
          dT_inst = 0;
          temp_correction_valid = False;
        }

  /*
   * Add in FPA temperature correction contribution if 
   * temperature is valid.
   */

        if (PP->PP_Emiss.T_fp[focal_plane][S] > 0)
        {
          dT_fpa = PP->PP_Emiss.T_fp[focal_plane][S] - refl_tables->T_FPA_ref[focal_plane];
          for (subsamp = 0; subsamp < SUBSAMPLES_AT_RES[R]; subsamp++) 
            temp_correction[subsamp] +=
                refl_tables->K_FPA[B_refl][D][subsamp][MS] * dT_fpa;
        }
        else
        {
          dT_fpa = 0;
          temp_correction_valid = False;
        }

        
  /*
   * If the temperature correction is valid, compute quantities that are
   * used in the calculation of uncertainty.  See the memo: "PFM RSB
   * Reflectance Uncertainty Algorithm (ReflUncert) - update 1 and 2, K.
   * (Vincent) Chiang, G. Godden, etc. November 23 (24 for update 2), 1999. 
   */

  /***************************************************************** 
   * obsolete due to UI algorithm update, 2/19/2011 by Xu Geng 
   */
       /* 
        if (temp_correction_valid) {

          for (subsamp = 0; subsamp < SUBSAMPLES_AT_RES[R]; subsamp++) 
          { 
        */
            /*  
             * compute contribution of Kinst to square of uncertainty.
             */
       /* 
            t = refl_tables->Sigma_K_inst[B_refl][D][subsamp][MS] *
                dT_inst / temp_correction[subsamp];
            sigma_Kinst_norm_sq[subsamp] = t * t;
        */

            /*
             * compute contribution of Tinst to square of uncertainty.
             */

       /* 
            t = refl_tables->Sigma_T_inst * 
                refl_tables->K_inst[B_refl][D][subsamp][MS] / 
                temp_correction[subsamp];
            sigma_Tinst_norm_sq[subsamp] = t * t;
        */

            /*
             * compute PART contribution of NEdL to square of uncertainty.
             */

       /* 
            t = refl_tables->RSB_NEdL[B_refl][D][subsamp][MS] / 
                refl_tables->m1[B_refl][D][subsamp][MS] /
                refl_tables->E_sun_over_pi[D_refl];

            RSB_NEdL_R_star_sq[subsamp] = t * t;
        */

            /*
             * compute contribution of m1 to square of uncertainty.
             */
            
       /* 
            t = refl_tables->Sigma_m1[B_refl][D][subsamp][MS] /
                refl_tables->m1[B_refl][D][subsamp][MS];
            sigma_m1_norm_sq[subsamp] = t * t;
        */

            /*
             * compute variance of residual electronic crosstalk.
             */
          
       /*
            sigma_PV_resid_elec_sq[subsamp] = 
                refl_tables->Sigma_PV_Resid_Elec[B_refl][D][subsamp] *
                refl_tables->Sigma_PV_Resid_Elec[B_refl][D][subsamp];
           
          } 

          sigma_RSB_ADC_sq = refl_tables->Sigma_RSB_ADC[B_refl][D] *
                             refl_tables->Sigma_RSB_ADC[B_refl][D];
        }
        */
  /****************************************************************/

        /*
         * Set the track index for this scan and detector.
         */

        T = S * DETECT_PER_BAND_AT_RES[R] + D;

        /*
         * Set the SWIR Correction band detector index for this detector.
         * (this is only needed for SWIR bands).
         */

        if (B_swir >= 0)
          bX_det_index = D / swir_detector_ratio[B_swir];
        if ( swir_tables->SWIR_corr_sending_detector[bX_det_index] >= 0 &&
             swir_tables->SWIR_corr_sending_detector[bX_det_index] < DETECTORS_PER_1KM_BAND)
          bX_det_index = swir_tables->SWIR_corr_sending_detector[bX_det_index]; 
        else {
          returnStatus = MODIS_F_OUT_OF_RANGE;
          L1BErrorMsg(location, returnStatus,
                      "Detector to use for SWIR OOB Correction is out of range.",
                      NULL, REFLECTIVE_TABLES_FILE, invalidlutfile, True);
          return returnStatus;
        }

        
  /*
   * Loop through frames at this resolution.  This is equivalent to a loop
   * through 1km frames and a loop through subsamples, where the loop
   * through subsamples is the most rapidly varying one.
   */
        

        for (F = 0; F < EV_1km_FRAMES * BAND_RATIO_AT_RES[R]; F++) 
        {
/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS

          /*
           * Here, for computing band 26 only, we want to skip all other
           * reflective bands.
           */

          if (Reflective_Cal_Band_Flag == REFLECTIVE_BAND_26_ONLY
              && R != INDEX_1000M_EMISS)
            continue;

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

          /*
           * Set the subsample index for this frame
           */

          subsamp = (F % samples_per_frame_at_res[R]);

          /* 
           * Extract values of quantities from arrays that have 
           * band-group-dependent variable names (includes extracting 
           * addresses of where to assign scaled integer and uncertainty 
           * index). 
           */
     
          if (R == INDEX_250M)
          {
            DN_ev = L1A_Scan->EV_250m[D][B][F];
            DN_obc_avg = PP->DN_OBC_Avg.DN_obc_250m_avg[T][B][subsamp];
            RVS = L1B_Gran->RSB_Cal_Coeff.RVS_250m[B][D][F][MS];
            SI_ptr = &L1B_Scan->SI.EV_250m_RefSB[B][D][F];
            UI_ptr = &L1B_Scan->UI.EV_250m_RefSB_UI[B][D][F];
          }
          else if(R == INDEX_500M)
          {
            DN_ev = L1A_Scan->EV_500m[D][B][F];
            DN_obc_avg = PP->DN_OBC_Avg.DN_obc_500m_avg[T][B][subsamp];
            RVS = L1B_Gran->RSB_Cal_Coeff.RVS_500m[B][D][F][MS];
            SI_ptr = &L1B_Scan->SI.EV_500m_RefSB[B][D][F];
            UI_ptr = &L1B_Scan->UI.EV_500m_RefSB_UI[B][D][F];
          } 
          else if(R == INDEX_1000M_DAY)
          {
            DN_ev = L1A_Scan->EV_1km_day[D][B][F];
            DN_obc_avg = PP->DN_OBC_Avg.DN_obc_1km_day_avg[T][B][subsamp];
            RVS = L1B_Gran->RSB_Cal_Coeff.RVS_1km_RefSB[B][D][F][MS];
            SI_ptr = &L1B_Scan->SI.EV_1km_RefSB[B][D][F];
            UI_ptr = &L1B_Scan->UI.EV_1km_RefSB_UI[B][D][F];
          } 
          else
          {
            DN_ev = L1A_Scan->EV_1km_night[D][B][F];
            DN_obc_avg = PP->DN_OBC_Avg.DN_obc_1km_night_avg[T][BAND26][subsamp];
            RVS = L1B_Gran->RSB_Cal_Coeff.RVS_1km_RefSB[L1B_B][D][F][MS];
            SI_ptr = &L1B_Scan->SI.EV_1km_RefSB[L1B_B][D][F];
            UI_ptr = &L1B_Scan->UI.EV_1km_RefSB_UI[L1B_B][D][F];
          }

          /* Initialize *SI_ptr */
          *SI_ptr = 0;
         
  /*
   *  Check quantities to ensure that pixel can be calibrated. If the pixel
   *  cannot be calibrated, set a specific value of scaled integer to flag
   *  the reason why the pixel could not be calibrated.  Also, set the
   *  uncertainty index to the bad data value.  Specific values to be set
   *  are described in the L1B EV file specifications and come from macros
   *  in Granule.h.
   */

  /*
   * Check if the DN value read from the L1A granule is out of the range
   * [-1 to 4095].  If so, this indicates corrupted data or a code bug.
   */
        
          if (DN_ev < MISSING_L1A_FLAG || DN_ev > saturated_dn[B_refl])
            Bad_L1A_Error_Out("EV_250m or EV_500m or EV_1km_day or "
                "band 26 in EV_1km_night",
                " out of range in middle L1A files, Reflective_Cal(), "
                "Reflective_Cal.c");

          /*
           * Check for the DN value being flagged as missing by L1A.
           */

          if ( DN_ev == MISSING_L1A_FLAG )            
          { 
            *SI_ptr = (uint16)L1A_DN_MISSING_SI;
            *UI_ptr = (uint8)BAD_DATA_UI;
            L1B_Gran->missing_pixels[B_38]++;
            L1B_Gran->valid_pixels[B_38]--;
            L1B_Gran->bad_data_flag[B_38] = 1;
            QA_Common->num_missing_data_in_scans[D_490]++;
            continue;
          }
          
          /*
           * Check for a dead detector.
           */

          if ( QA_tables->dead_detector[D_490])
          {
            *SI_ptr = (uint16)DEAD_DETECTOR_SI;
            *UI_ptr = (uint8)BAD_DATA_UI;
            L1B_Gran->valid_pixels[B_38]--;
            L1B_Gran->dead_detector_pixels[B_38]++;
            QA_Common->num_dead_detector_EV_data[D_490]++;
            continue;
          }

          /*
           * Check for a dead subframe.
           * bits 0 to 3 for noisy subframes, bits 4 to 7 for dead subframes 
           */

          if(D_490 < NUM_HIGH_RESOLUTION_DETECTORS){ 
            if ( QA_tables->Detector_Quality_Flag2_Values[D_490][subsamp+4])
            {
              *SI_ptr = (uint16)DEAD_SUBFRAME_SI;
              *UI_ptr = (uint8)BAD_DATA_UI;
              L1B_Gran->valid_pixels[B_38]--;
              L1B_Gran->dead_subframe_pixels[B_38]++;
              QA_Common->num_dead_subframe_EV_data[D_490]++;
              continue;
            }
          }

          /* 
           * Check if the instrument is in a sector rotation. 
           */

          if (QA_Common->Sector_Rotation[S] == True)
          {
            *SI_ptr = (uint16)SECTOR_ROTATION_SI;
            *UI_ptr = (uint8)BAD_DATA_UI;
            L1B_Gran->valid_pixels[B_38]--;
            L1B_Gran->bad_data_flag[B_38] = 1;
            QA_Common->num_sector_rotation_EV_data[D_490]++;
            continue;
          }

          /*
           * Check if the instrument has both sides of PCLW electronics on at the same time.
           */

          if (QA_Common->Electronic_Anomaly[S] == True)
          {
            *SI_ptr = (uint16)UNABLE_CALIBRATE_SI;
            *UI_ptr = (uint8)BAD_DATA_UI;
            L1B_Gran->valid_pixels[B_38]--;
            L1B_Gran->bad_data_flag[B_38] = 1;
            continue;
          }
          
          /* Check that DN_ev is not saturated. */
          
          if (DN_ev == saturated_dn[B_refl])  {
            *SI_ptr = (uint16)SATURATED_DETECTOR_SI;
            *UI_ptr = (uint8)BAD_DATA_UI;
            L1B_Gran->saturated_pixels[B_38]++;
            L1B_Gran->valid_pixels[B_38]--;   
            L1B_Gran->bad_data_flag[B_38] = 1;      
            QA_Common->num_saturated_EV_data[D_490]++; 
            continue;
          }
          
  /*
   * Check that average electronic background DN is valid.
   * A negative value is assigned to average electronic background DN only if
   * a valid value could not be determined from SV or BB data.
   */
          
          if (DN_obc_avg < 0)
          {
            *SI_ptr = (uint16)ZERO_POINT_DN_SI;
            *UI_ptr = (uint8)BAD_DATA_UI;
            L1B_Gran->valid_pixels[B_38]--;
            L1B_Gran->bad_data_flag[B_38] = 1;
            QA_Common->num_no_bg_DN_EV_data[D_490]++;
            continue;
          }
         
          /*************** Checks complete ***********/


  /* 
   * Subtract average electronic background DN from earth view DN.
   * See memo "Level 1B Uncertainty Calculations for Reflective Solar 
   * Bands" by Bruce Berriman and Gail Reichert, MCST Document # M0703, 
   * May 27, 1999 for algorithm. 
   */

          dn_ev = DN_ev - DN_obc_avg;

          /* 
           * Determine if the dn_ev is saturated.
           */

          if (dn_ev >= refl_tables->dn_sat_ev[B_refl][D][subsamp][MS]) {
          if (B_refl != 9 && B_refl != 11 && B_refl != 12 && B_refl != 17) {
            *SI_ptr = (uint16)SATURATED_DETECTOR_SI;
            *UI_ptr = (uint8)BAD_DATA_UI;
            L1B_Gran->saturated_pixels[B_38]++;
            L1B_Gran->valid_pixels[B_38]--;   
            L1B_Gran->bad_data_flag[B_38] = 1;      
            QA_Common->num_saturated_EV_data[D_490]++; 
            continue;
          } else
             *SI_ptr = (uint16)UNRESCALED_HIGH_SI;
          }

          /* 
           * For SWIR band, if the switch is on, correct dn_ev.
           * See "SWIR Correction Algorithm Change in L1B", Draft,
           * Jack Xiong, Dec. 1, 2000.
           */

          if (B_swir >= 0 && swir_tables->SWIR_correction_switch == ON) {
            bX_frame_index = F / num_swir_subsamples[B_swir];
            dn_X = L1B_Scan->dn_X[bX_det_index][bX_frame_index];
            if (dn_X < 0)
            {
              dn_swir_corr = 0;
              dn_swir_corr_error = True;
              sigma_dn_swir = 0;
            }
            else {
              dn_swir_corr = swir_tables->X_OOB_0[B_swir][D][subsamp][MS] +
                             swir_tables->X_OOB_1[B_swir][D][subsamp][MS] * dn_X +
                             swir_tables->X_OOB_2[B_swir][D][subsamp][MS] * dn_X * dn_X;
              dn_swir_corr_error = False;
              sigma_dn_swir = refl_tables->swir_ui_factor[B_swir] * dn_swir_corr;
            }
            dn_ev -= dn_swir_corr;
          }
          else {
            dn_swir_corr = 0;
            dn_swir_corr_error = False;
            sigma_dn_swir = 0;
          }
        
  /* 
   * Convert dn to dn*.
   * See memo "Level 1B Uncertainty Calculations for Reflective Solar Bands",
   * Bruce Berriman and Gail Reichert, MCST Document # M0703, May 27, 1999.
   */
          
          dn_star = dn_ev * temp_correction[subsamp] / RVS;
   
  /* 
   * Convert dn* to dn**.
   * Use equation (5) in "Calculation of the Digital Signals Written to the 
   * Level 1B Data Products for the Reflective Solar Bands", Bruce Berriman
   * and Jim Rogers, MCST Document # M0825, October 27, 1999.
   */

          reflectance  = refl_tables->m0[B_refl][D][subsamp][MS] + dn_star *
              L1B_Gran->RSB_Cal_Coeff.m1_des_sq[B_refl][D][subsamp][MS];

          dn_star_star = reflectance / 
              L1B_Gran->RSB_Cal_Coeff.m1_des_sq_max[B_refl]; 

  /*
   * Check to see if dn** is below the bottom end of the range for scaling to
   * the scaled integer or is above the top end of the dynamic range for scaling.
   * Otherwise, convert dn** to the scaled integer.
   */

          if ( dn_star_star < L1B_Gran->SO.dn_star_Min[B_refl])
          {
            *SI_ptr = (uint16) RSB_DN_STAR_BELOW_MIN_SI;
            *UI_ptr = (uint8) BAD_DATA_UI;
            L1B_Gran->negative_value_below_noise_pixels[B_38]++;
            L1B_Gran->valid_pixels[B_38]--;
            L1B_Gran->bad_data_flag[B_38] = 1;
            QA_Common->num_bad_dn_star_star_RSB_EV_data[D_490]++;
            continue;
          }
          else if ( dn_star_star > L1B_Gran->SO.dn_star_Max[B_refl])
          {
            *SI_ptr = (uint16) TEB_OR_RSB_GT_MAX_SI;
            *UI_ptr = (uint8)BAD_DATA_UI;
            L1B_Gran->valid_pixels[B_38]--;
            L1B_Gran->bad_data_flag[B_38] = 1;
            QA_Common->num_exceed_max_for_scaling[D_490]++;

            if (B_refl == 9 && (RFLAG & 2) == 2 && (RSCL_FLAG & 4) == 4) {
              L1B_Scan->EV_500m_Aggr1km_RefSB[0][D][F] =
                (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[2]
                         + L1B_Gran->SO.counts_offset_RefSB[2] + 0.5);
              *SI_ptr = (uint16) RESCALED_L1B_SI;
            } else if (B_refl == 11 && (RFLAG & 2) == 2 && (RSCL_FLAG & 8) == 8) {
              L1B_Scan->EV_500m_Aggr1km_RefSB[1][D][F] =
                (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[3]
                         + L1B_Gran->SO.counts_offset_RefSB[3] + 0.5);
              *SI_ptr = (uint16) RESCALED_L1B_SI;
            } else if (B_refl == 12 && (RFLAG & 1) == 1 && (RSCL_FLAG & 1) == 1) {
              L1B_Scan->EV_250m_Aggr1km_RefSB[0][D][F] =
                (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[0]
                         + L1B_Gran->SO.counts_offset_RefSB[0] + 0.5);
              *SI_ptr = (uint16) RESCALED_L1B_SI;
            } else if (B_refl == 17 && (RFLAG & 1) == 1 && (RSCL_FLAG & 2) == 2) {
              L1B_Scan->EV_250m_Aggr1km_RefSB[1][D][F] =
                (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[1]
                         + L1B_Gran->SO.counts_offset_RefSB[1] + 0.5);
              *SI_ptr = (uint16) RESCALED_L1B_SI;
            }
            continue;
          }
          else if (*SI_ptr == (uint16)UNRESCALED_HIGH_SI)
            {
              if (B_refl == 9 && (RFLAG & 2) == 2 && (RSCL_FLAG & 4) == 4) {
                L1B_Scan->EV_500m_Aggr1km_RefSB[0][D][F] =
                  (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[2]
                           + L1B_Gran->SO.counts_offset_RefSB[2] + 0.5);
              } else if (B_refl == 11 && (RFLAG & 2) == 2 && (RSCL_FLAG & 8) == 8) {
                L1B_Scan->EV_500m_Aggr1km_RefSB[1][D][F] =
                  (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[3]
                           + L1B_Gran->SO.counts_offset_RefSB[3] + 0.5);
              } else if (B_refl == 12 && (RFLAG & 1) == 1 && (RSCL_FLAG & 1) == 1) {
                L1B_Scan->EV_250m_Aggr1km_RefSB[0][D][F] =
                  (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[0]
                           + L1B_Gran->SO.counts_offset_RefSB[0] + 0.5);
              } else if (B_refl == 17 && (RFLAG & 1) == 1 && (RSCL_FLAG & 2) == 2) {
                L1B_Scan->EV_250m_Aggr1km_RefSB[1][D][F] =
                  (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[1]
                           + L1B_Gran->SO.counts_offset_RefSB[1] + 0.5);
              }
              continue;
            }
          else
            SI = (uint16)(dn_star_star/L1B_Gran->SO.counts_scale_RefSB[B_refl]
                       + L1B_Gran->SO.counts_offset_RefSB[B_refl] + 0.5);
          
  /*
   * Check if the NAD door is closed. Set the most significant bit to 1 if
   * NAD door is closed.
   */

           if (QA_Common->NAD_Door_Open[S] == False)
           {
             SI = SI|0x8000;
             QA_Common->num_nadir_door_closed_EV_data[D_490]++;
 
             if (SI > NAD_CLOSED_UPPER_SI)
               SI = NAD_CLOSED_UPPER_SI;
           }

          *SI_ptr = (uint16)SI;

 /***********************************************************************
  * 2/19/2011, Xu Geng
  * The algorithm for uncertainty calculation has been updated.   
  * Please see Junqiang Sun, et. al., "RSB Uncertainty algorithm update", 
  * (MCST Internal Memo). 
  * Please ignore the old Memo mentioned below.  
  ***********************************************************************/

 /*
  * Calculate uncertainty in the reflectance product and convert to
  * uncertainty index.
  * See Kwo-Fu (Vincent) Chiang, et. al., "PFM Reflectance Uncertainty 
  * Algorithm (ReflUncert) - update 2", MCST Internal Memo, 
  * November 24, 1999, for algorithm.
  * Note that the last term of the uncertainty uses dn*, not dn**.
  */

 /*
  * The following are conditions under which the uncertainty index is set
  * to the bad data value (maximum percentage error):
  * - dn_ev <= TOLERANCE        (could get divide by zero, 
  *                                 negative or tiny number)
  * - dn_star <= TOLERANCE      (could get divide by zero, 
  *                                 negative or tiny number)
  * - temp. correction invalid  (unsure of how big correction 
  *                                 would have been)
  * - NAD is closed             (the scaled integers are set to 
  *                                 unusable values)
  * - SWIR band error           (could not use SWIR OOB Correction band
  *                                 to get correction)
  */

          if (dn_ev <= TOLERANCE || dn_star <= TOLERANCE
                      || !temp_correction_valid 
                      || QA_Common->NAD_Door_Open[S] == False
                      || dn_swir_corr_error)
          {
            *UI_ptr = (uint8)BAD_DATA_UI;
          } 
          else
          {

           /* 
            * Compute the residual uncertainty of the R* linear fit. 
            */
           /*
            *obsolete due to UI algorithm update, 2/19/2011, Xu Geng 
            */
            /*
            resid_coeff = 
                refl_tables->Sigma_R_Star_Lin_Resid_Ucoeff
                    [B_refl][D][subsamp][MS];
            EVAL_4TH_ORDER_POLYNOMIAL(sigma_R_EV,resid_coeff,reflectance);
            sigma_R_EV /= reflectance;
            */

            /*  
             * Compute RVS term
             */
            /* obsolete, 2/19/2011 by Xu Geng */
            /*sigma_RVS_over_RVS = sigma_RVS / RVS; */


            /*
             * Calculate U4 term 
             * u4 = (u4_coeffs[0]+u4_coeffs[1]*dn_ev+u4_coeffs[2]*dn_ev*dn_ev)/dn_ev 
             */
            u4_coeffs = refl_tables->u4_coeffs[B_refl][D][subsamp][MS]; 
            u4 = u4_coeffs[0]/dn_ev+u4_coeffs[1]+u4_coeffs[2]*dn_ev; 

            /*
             * Get the pre-calculated u2 
             */ 
            u2 = L1B_Gran->RSB_Cal_Coeff.u2[D_refl][F/samples_per_frame_at_res[R]][MS];

  /*
   * Compute the fractional uncertainty in the reflectance product.
   * This is the square root of the sum of the squares of the individual
   * contributions from the following 5 terms:
   * U1: The common term which is AOI and time independent. 
   * U2: The uncertainty due to calibrations using calibrators and the linear approximation in RVS AOI dependence. 
   * U3: The temperature impact.
   * U4: Scene dependent noise to signal ratio.
   * U5: The uncertainty due to cross talk correction for SWIR bands.   
   */

            uncertainty = sqrt( (double)
                                (refl_tables->u1[D_refl] * refl_tables->u1[D_refl] +  
                                u2 * u2 +
                                refl_tables->u3[D_refl][MS] * refl_tables->u3[D_refl][MS] + 
                                u4 * u4 +
                                (sigma_dn_swir * sigma_dn_swir) /
                                (dn_ev * dn_ev))
                              );
                                  
        
  /****************************************************************************
   * obsolete due to UI algorithm update, 2/19/2011, Xu Geng
    
            uncertainty = sqrt( (double)
                                (sigma_m1_norm_sq[subsamp] +
                                 sigma_RVS_over_RVS * sigma_RVS_over_RVS +
                                 sigma_Kinst_norm_sq[subsamp] +
                                 sigma_Tinst_norm_sq[subsamp] +
                                 sigma_PV_resid_elec_sq[subsamp] +
                                 (sigma_RSB_ADC_sq + sigma_dn_swir) / 
                                 (dn_ev * dn_ev) + 
                                 sigma_R_EV * sigma_R_EV +
                                 RSB_NEdL_R_star_sq[subsamp] / 
                                 (dn_star * dn_star))
                              );
  ****************************************************************************/

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

            if (uncertainty > TOLERANCE)
              ui = refl_tables->RSB_UI_scaling_factor[B_refl] * 
                   log((double) (100. * uncertainty / 
                       refl_tables->RSB_specified_uncertainty[B_refl]));
            else
              ui = 0;

            /*
             * Ensure that uncertainty index lies within the [0-15] range
             * and place value in the uncertainty index variable.
             * Note that MAX_UNCERTAINTY_UI is <= BAD_DATA_UI, which is 15.
             */

            if (ui > MAX_UNCERTAINTY_UI) ui = BAD_DATA_UI;
            else if(ui < 0) ui = 0;

            *UI_ptr = (uint8)ui; 
          }

        }                     /* End of loop over frames */
      }                       /* End of loop over detectors */
      B_refl++;               /* Increment B_refl*/
    }                         /* End of loop over bands */
  }                           /* End of loop over L1A band groups */


  return(MODIS_S_OK);
}

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS

int32 Reflective_Cal_Band_Flag;      /* external variable */

PGSt_SMF_status Copy_Band_26_Data (L1B_Scan_t *L1B_Scan)
/*
!C**********************************************************************
!Description:
    This function copies band 26 data from the EV_1km_RefSB structure member
    of L1B_Scan to the appropriate Band_26 structure member.

!Input Parameters:
    L1B_Scan_t * L1B_Scan   Contains all EV data for 1 scan, the structure
                            members "SI.EV_1kmRefSB" and "UI.EV_1km_RefSB_UI"
                            have been filled in function Reflective_Cal.
!Output Parameters:
    L1B_Scan_t * L1B_Scan   Structure members "Band26.SI" and
                            "Band26.UI" are now filled.

!Revision History:
   Revision 01.00 May 4, 1999
   Initial development, Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   Band 26 is the last band in the EV_1km_RefSB arrays in L1B_Scan.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 B = NUM_1000M_REFL_BANDS - 1;         /* index of band 26 data */
  int32 D;                                    /* detector index */
  int32 F;                                    /* frame index */

  for (D = 0; D < DETECTORS_PER_1KM_BAND; D++) {
    for (F = 0; F < EV_1km_FRAMES; F++) {
      L1B_Scan->Band26.SI[D][F] = L1B_Scan->SI.EV_1km_RefSB[B][D][F];
      L1B_Scan->Band26.UI[D][F] = L1B_Scan->UI.EV_1km_RefSB_UI[B][D][F];
    }
  }
  return returnStatus;
}
#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/



extern PGSt_SMF_status  Band_26_Crosstalk_Correction (
                   L1B_Scan_t            *L1B_Scan,
                   int16                 *b5_frame_offset,
#ifdef USE_B5_RAD_OFFSET
                   float32               b5_rad_offset,
#endif /* USE_B5_RAD_OFFSET */
                   float32               *b26_fr_b5_scaled_corr,
                   QA_Common_t           *QA_Common,
                   uint32                *valid_pixels,
                   uint32                *negative_value_below_noise_pixels,
                   int16                 *bad_data_flag,
                   boolean               isdaymode,
                   boolean               perform_B26_correction)
/*
!C**********************************************************************
!Description:
  For one scan of calibrated earth-view data, this routine performs
  a crosstalk correction to Band 26 based on the values of the 
  Band 5 Scaled Integers (SIs).  If no errors occur the corrected
  Band 26 value is output.  

!Input Parameters:
 L1B_Scan_t          *L1B_Scan       (->SI, ->UI) (Bands 5 & 26)
 L1A_Scan_t          *L1A_Scan       contains earth view DN values
 int16               *b5_frame_offset    Frame offset to use for Band 5 
 float32             b5_rad_offset   Band 5 radiance offset.  (It is not
                                       necessary to use this term unless
                                       the LUT value of dn* min for band
                                       5 is set to a value other than 0.)
 float32             *b26_fr_b5_scaled_corr
                                     already-scaled correction factors
 common_QA_tables_t  *QA_tables      contains the dead detector lookup table 
                                       values
 boolean             isdaymode       day mode indicator
 boolean             perform_B26_correction   correction switch

!Output Parameters:
 L1B_Scan_t          *L1B_Scan      (->SI, ->UI) (Band 26)
 common_QA_tables_t  *QA_tables      contains possibly updated error counts

!Revision History:
 $Log: Reflective_Cal.c,v $
 Revision 1.24  2011-04-07 14:41:47-04  xgeng
 RSB uncertainty algorithm update

 Revision 1.22  2010-11-15 10:40:39-05  xgeng
 No calibration is peformed for the scan if both sides of PCLW electronics are on at the same time.

 Revision 1.21  2008/11/18 19:49:15  xgeng
 merge branch for V6.0.0

 Revision 1.20.2.2  2008/06/05 14:06:46  xgeng
 Filled the dead subframe with value of DEAD_SUBFRAME_SI

 Revision 1.20  2005/01/18 19:34:49  ltan
 MOD_PR02_TERRA update to V5.0.4


 Revision 1.02  April 25, 2003   Razor Issue # 190
 Remove passed parameter L1A_Gran; no longer needed to differentiate which 
 platform is being processed because procedure is now applied to MODIS/Aqua
 as well.
 Alice Isaacman  SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 1.01    April 22, 2002
 Correct typo in original code; check on SI_26, not UI_26
 Alice Isaacman  SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Original Code     Razor Issue 182
 Applicable to MODIS/TERRA (PFM) processing only.
 Alice Isaacman  SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)
 
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
  int16   D               = 0;  /* detector index within band */
  int16   D_490           = 0;  /* detector index within all MODIS detectors*/
  int16   F               = 0;  /* Frame index within resolution of 1 km bands */
  int16   B5_F            = 0;  /* Calculated Band 5 frame to use */
  float32 det_corr_fac    = 0.; /* Scaled correction factor for a given detector */
  uint16  calc_corr_fac   = 0.; /* Calculated correction factor for a pixel */
  uint16  SI_5;                 /* Aggregated Band 5 scaled integer */
  uint8   UI_5;                 /* Band 5 uncertainty */
  uint16  *SI_26_ptr;           /* Pointer to Band 26 scaled integer */
  uint8   *UI_26_ptr;           /* Pointer to Band 26 uncertainty */

      
  /*
   *   Double check to be sure this is a day mode scan and that the Band 26
   *       correction switch is set.  
   *       Return without penalty if conditions are not met.
   */
  
  if ( (! isdaymode) || (! perform_B26_correction) )
    return MODIS_S_OK;
  
  /* Loop through Detectors */

  D_490 =  MODIS_BAND26_INDEX*DETECTORS_PER_1KM_BAND;
          
  for (D = 0; D < DETECTORS_PER_1KM_BAND; D++, D_490++)
  {
    
    det_corr_fac = b26_fr_b5_scaled_corr[D];
    for (F = 0; F < EV_1km_FRAMES; F++) 

    { 
      
  /*
   * Calculate the new frame for Band 5.  If out of bounds, leave the
   *           Band 26 SI and Band 26 UI alone.
   */ 
      B5_F = F + b5_frame_offset[D];
      if (B5_F < 0 || B5_F >= EV_1km_FRAMES) continue;
      
  /* Pull up the SI for the aggregated Band 5 and Band 26 radiance data. */  
      
      SI_5  = L1B_Scan->EV_500m_Aggr1km_RefSB[BAND_5_AGGR_INDEX][D][B5_F];
      UI_5  = L1B_Scan->EV_500m_Aggr1km_RefSB_UI[BAND_5_AGGR_INDEX][D][B5_F];
      SI_26_ptr = &L1B_Scan->SI.EV_1km_RefSB[ BAND_26_1KM_INDEX][D][F];
      UI_26_ptr = &L1B_Scan->UI.EV_1km_RefSB_UI[BAND_26_1KM_INDEX][D][F];
      
  /*
   *   If the Band 26 SI is invalid, do not correct the Band 26 value and
   *       do not change the Band 26 UI.
   */ 
      if ( *SI_26_ptr > (uint16) DN15_SAT) continue;

  /*
   *   If the Band 5 SI is invalid, set the Band 26 UI to 15 and do not
   *       correct the Band 26 value.
   */ 

      if ( SI_5 > (uint16) DN15_SAT)
      {
        *UI_26_ptr = (uint8) BAD_DATA_UI;
        continue;
      }
        
  /*   Calculate the correction = Band 5 SI * b26_fr_b5_scaled_corr[det] */
    
#ifdef USE_B5_RAD_OFFSET    
      
  /*
   * Use the Band 5 offset term in the calculation.  This will not 
   * become necessary unless the LUT value of dn* min is set
   * to a number other than 0.
   */
      
      calc_corr_fac = (uint16) (((float32) SI_5 - b5_rad_offset)
           * det_corr_fac + .5);
                      
#else
      
      calc_corr_fac = (uint16) ((float32) SI_5 * det_corr_fac + .5);
      
#endif /* USE_B5_RAD_OFFSET */
      
  /*
   *   If the correction which would be subtracted from the Band 26 SI value 
   *       is larger than the Band 26 SI value (i.e. if the resulting SI after
   *       "correction" were to be less than 0), set Band 26 UI to 15 and 
   *       set the Band 26 value to RSB_DN_STAR_BELOW_MIN_SI.
   */ 

      if ( calc_corr_fac > *SI_26_ptr)
      {
        *UI_26_ptr = (uint8) BAD_DATA_UI;
        *SI_26_ptr = (uint16) RSB_DN_STAR_BELOW_MIN_SI;
        negative_value_below_noise_pixels[MODIS_BAND26_INDEX]++;
        valid_pixels[MODIS_BAND26_INDEX]--;
        bad_data_flag[MODIS_BAND26_INDEX] = 1;
        QA_Common->num_bad_dn_star_star_RSB_EV_data[D_490]++;
        continue;
      }
        
       /* Perform the correction.   */
      
      *SI_26_ptr -= calc_corr_fac;

      /* If the Band 5 UI value was maximal set the Band 26 UI likewise. */
      
      if (UI_5 == (uint8) BAD_DATA_UI) *UI_26_ptr = (uint8) BAD_DATA_UI;
        
    }  /* Loop through frames    */
  }    /* Loop through detectors */
  
  return returnStatus;  
}

