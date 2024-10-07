#ifndef REFLECTIVE_CAL_H
#define REFLECTIVE_CAL_H

/*
 * UNCOMMENT THE FOLLOWING LINE TO USE THE BAND 5 RADIANCE OFFSET IN THE
 * CALCULATION OF THE BAND 26 CROSSTALK CORRECTION:
 */ 
/* #define USE_B5_RAD_OFFSET */

#include "L1B_Tables.h" /* contains declarations of look up table varables and
                           Granule.h */
#include "Preprocess.h" /* for Preprocess_Data_t definition */

/*
!C-INC***********************************************************************
!Description:   Contains external declarations for reflective calculation functions.

!Revision History:
 $Log: Reflective_Cal.h,v $
 Revision 1.10  2006-10-30 10:00:14-05  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.

 Revision 1.02  April 25, 2003   Razor Issue # 190
 Removed passed parameter L1A_Gran from Band_26_Crosstalk_Correction; no longer needed to
 differentiate which  platform is being processed because procedure is now applied to
 MODIS/Aqua as well.
 Alice Isaacman  SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.20  March 19, 2002    Razor Issue #182
 For MODIS/TERRA (PFM) processing only:
 Added function prototype for Band_26_Crosstalk_Correction
 Alice Isaacman   SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.10 May 13, 1999
 Added band 26 section
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 
 Revision 01.01 Feb 17, 1999
 Moved DN_to_DN_star to Preprocess.h as part of moving DN_to_DN_star
 to Preprocess.c to avoid circular include dependencies.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 
 Revision 01.00 Nov 1998
 initial development 
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST) for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS

/* Use the external variable Reflective_Cal_Band_Flag to signal to
 * Reflective_Cal whether we should calculate only band 26 or all
 * reflective bands.  The macros ALL_REFLECTIVE_BANDS and
 * REFLECTIVE_BAND_26_ONLY are the allowable values for
 * Reflective_Cal_Band_Flag.
 */
#define ALL_REFLECTIVE_BANDS    1
#define REFLECTIVE_BAND_26_ONLY 2
extern int32 Reflective_Cal_Band_Flag;
extern PGSt_SMF_status Copy_Band_26_Data (L1B_Scan_t *L1B_Scan);

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

PGSt_SMF_status  Reflective_Cal (int16               S, 
                                 L1A_granule_t       *L1A_Gran,
                                 L1B_granule_t       *L1B_Gran,
                                 L1A_Scan_t          *L1A_Scan,
                                 L1B_Scan_t          *L1B_Scan,
                                 Preprocess_Data_t   *PP,
                                 refl_tables_t       *refl_tables,
                                 common_QA_tables_t  *QA_tables,
                                 QA_Common_t         *QA_Common);

#define BAND_5_AGGR_INDEX  MODIS_BAND5_INDEX - NUM_250M_BANDS
#define BAND_26_1KM_INDEX  NUM_1000M_REFL_BANDS - 1

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
                   boolean               perform_correction);
#endif

