#ifndef GRANULE_H
#define GRANULE_H

#include    "PGS_SMF.h"    /* define PGSt_SMF_status */  
#include    "hdf.h"        /* contain prototype for hdf functions */
#include    "mfhdf.h"      /* contain prototype for hdf functions */
#include    "HdfEosDef.h"  /* contain prototype for eos hdf functions */

#define  MAX_NC_NAME   H4_MAX_NC_NAME
#define  MAX_VAR_DIMS  H4_MAX_VAR_DIMS

/*
!C-INC***********************************************************************
!Description:  Public header file containing macros, structure definitions 
               and external function declarations.

!Revision History:
   $Log: Granule.h,v $
   Revision 1.24  2011-04-07 14:40:12-04  xgeng
   1. RSB &TEB uncertainty algorithm update; 2. The quadratic RSB RVS changed to 4th order.

   Revision 1.23  2010-11-15 11:25:16-05  xgeng
   Defined a new fill value UNABLE_CALIBRATE_SI (65524); Added a new boolean member Electronic_Anomaly to structure QA_Common_t.

   Revision 1.22  2009/08/31 17:40:18  xgeng
   Change LOG to Log

   Revision 02.38, Janurary 7, 2008  Razor Issue #216
   Added a new constant "DEAD_SUBFRAME_SI" for dead subframe
   Added a new constant "NUM_HIGH_RESOLUTION DETECTORS"
   Added a new constant "NUM_HIGH_RESOLUTION SUBFRAMES"
   Added a new variable "dead_subframe_pixels" in L1B_granule_t
   Added a new variable "num_dead_subframe_EV_data" in QA_Common_t
   Xu Geng, SAIC GSO (xu.geng@saic.com)

   Revision 1.16  2006/10/27 15:02:54  ltan
   Added parameter PGE02_Version. Changed for ANSI-C compliance. Correction for the generation of code change log.

   Revision 02.37, October 31, 2003  Razor Issue #195
   Added ProcessingCenter to runtime parameters
   Liqin Tan, SAIC GSO (ltan@saicmodis.com)

   Revision 02.36 September 15, 2003  Razor Issue #181
   Added ProcessingEnvironment to runtime parameters
   Liqin Tan, SAIC GSO (ltan@saicmodis.com)

   Revision 02.35  March 26, 2003  Razor Issue #191
   Changed name of dn_28 to dn_X for use with any emissive band.
   Added band number to L1B_Scan structure type.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.32 June 28, 2002    Razor Issue #161
   Moved NUM_*_SUBSAMP definitions from Preprocess.h; added NUM_REFL_INDICES
   Gwyn Fireman, SAIC-GSO <Gwyn.Fireman@gsfc.nasa.gov>

   Revision 02.31   April 16, 2002  Razor Issue #166
   Added leading_granule_scan_gap and trailing_granule_scan_gap to 
      QA_Common_t
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.29 March 19, 2002  Razor Issue #182
   Applicable to MODIS/TERRA (PFM) processing only.
   Added b26_fr_b5_scaled_corr to RSB_Cal_Coeff_t.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.28 March 10, 2002  Razor Issue #174
   Added RVS tables to RSB_Cal_Coeff_t.  Defined TEB_Cal_Coeff_t to hold
      RVS tables.   
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.27 February 7, 2002  Razor Issue #180
   Added flag TEB_B1_NOT_CALCULATED and counter for it in QA_Common.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.26 December 14, 2001
   Changed MAX_*_TRACK_DIM parameters so that they are automatically determined
   by MAX_NUM_SCANS.  
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.25 November 12, 2001       (Razor Issue #169)
   Added boolean skip_night_hi_res to prototype for Close_L1B_Granule.
   Added Write_Night_Mode_HiRes_Data to Run_Time_Parameters_t
   Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

   Revision 02.24  November 6, 2001
   Added MCST_LUT_Version to Run_Time_Parameters_t (Razor issue #167)
   Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

   Revision 02.23  Sept. 28, 2001
   Eliminated use of defined PI value due to conflict with SDP Toolkit 5.2.7
   (PI was defined to be PGS_PI in any case) (Razor issue #168)
   Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

   Revision 02.22  Feb 22, 2001, Issue 155
   Added Run_Time_Parameters_t and Read_Run_Time_Parameters.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.21  Jan 29, 2001
   Corrected variable names and simplified L1A_Scan_t.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.21  Dec. 5, 2000
   Razor issue 147.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.20  Dec. 1, 2000
   Added dn_28 to L1B_Scan_t
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.19  October 29, 2000
   Added new unusable data value for L > L_Max as per issue 139.
   Added code for granule average QA values as per issue 140.
   Added structure element to L1A_granule_t, Razor issue 142.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.18  September 21, 2000
   Added counting num_rsb_at_night_scans as per issue 137.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   ... (many changes not logged) ...

 Revision 02.17 Nov 05, 1999
 Change momos_t to RSB_Cal_Coeff_t and add new variables
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.16 Sep 30, 1999
 Made num_thermistor_outliers dimensioned [MAX_NUM_SCANS] in QA_Emiss_t
 Jim Rogers (rogers@msct.gsfc.nasa.gov)

 Revision 02.15 May 13, 1999
 Added band 26 section.
 Jim Rogers (rogers@msct.gsfc.nasa.gov)

 Revision 02.14 Feb 22, 1999
 Added macros INDEX_1000M_DAY and INDEX_1000M_NIGHT to allow clarity
 uf usage elsewhere in the code.
 Jim Rogers (rogers@msct.gsfc.nasa.gov)

 Revision 02.13 Feb 8, 1999
 Moved definition of L1B_Scan_Metadata_t to L1B_SetupP.h since that structure
 is local in scope to L1B_Setup.c
 Jim Rogers (rogers@msct.gsfc.nasa.gov)

 Revision 02.12 Feb 2, 1999
 Separated out parts belonging to L1B_Setup as part of the
 process of separating L1B_Setup out into its own module.
 Added L1B_BANDS_AT_RES as an external declaration.
 Jim Rogers (rogers@msct.gsfc.nasa.gov)

 Revision 02.11 July 1998
 Add Write_L1B_SI_UI() prototype
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)
  
 Revision 02.10 April 13, 1998.
 Removed the SUN_EARTH_DIST_ADJUST() and DAYS_IN_YEAR macros.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 April 10, 1998.
 Changed the L1A_granule_t's EVStartTime_TAIsecond, MirrorSide, and
 ScanType from pointers to static arrays.
 Removed the L1B_granule_t's SunAdj, EVStartTime_TAIday, 
 and EVStartTime_TAIyear.
 Added momos_t variable to L1B_granule_t.
 Removed the following from L1B_ScaleOffset_t:
   rad_scale_RefSB
   rad_offset_RefSB
   refl_scale_RefSB
   refl_offset_RefSB
 (they're replaced by the momos_t)
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)
 
 Revision 02.10 April 1998
 Added L1B_Scan_metadat_t structure.
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)
 
 Revision 02.10 April 1998.
 Added the momos_t structure.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 March 1998.
 Add ADC, PCX and INT correction switch 
 Add emiss_band_index_t data structure
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 02.00 March 1997.
 Summarizes Geom_Param_Lib.h, Granule.h, L1A_to_L1B.h, ......
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.01 1996/04/09
 Updated to match Version 1 Design Document.
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Neal Devine(devine@ltpmail.gsfc.nasa.gov)
 Joan Baden (baden@highwire.gsfc.nasa.gov)

 Revision 01.00 1993
 Initial development
 Geir E. Kvaran(geir@highwire.gsfc.nasa.gov)

 Revision 01.00 Nov. 1998
 This is part of original Granule.h, L1B.h
 Initial development
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   
1. L1A Granule.  The L1A granule contains 4 HDF SDSs, each corresponding
   to a group of MODIS bands:

   L1A SDS name               MODIS bands
   ------------    ---------------------------------------------------------
   "XX_250m"       1,2
   "XX_500m"       3,4,5,6,7
   "XX_1km_day"    8,9,10,11,12,13lo,13hi,14lo,14hi,15,16,17,18,19
   "XX_1km_night"  20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36

   where:

   XX = SD        (Solar Diffusor)
   XX = SRCA      (Spectro-Radiometric Calibration Assembly)
   XX = BB        (Black Body)
   XX = SV        (Space View)
   XX = EV        (Earth View)

   "lo" and "hi" refer to two different gain levels.  Thus, 13lo and 13hi
   are treated as different bands even though they have the same center
   wavelegnth.  Similarly for 14lo and 14hi.  Thus, MODIS has 36 different
   center wavelengths but 38 different "bands" of data.

   Within the code, we refer to the above L1A SDSs as 4 "resolutions", even
   though the last two groups have the same resolution (1km).  All SDSs are
   present in one L1A granule HDF file.  The logic implemented in
   "Open_L1A_EV_SDS" is that each of the above SDSs must be present in the
   L1A granule.  If any is missing, then the function returns a
   MODIS_F_READ_ERROR.  The logic implemented in "Read_L1A_EV_Scan" (which
   reads one scan) is that the "EV_1km_night" data are always read but the
   other three are read only if the scan type is "DAY".

   Each L1A input granule covers a time period of approximately 5 minutes.
   The actual number of scans will be either 203 or 204, depending on how a
   scan fits into the granule (there are always whole scans in the granule).
   Thus, the maximum number of scans that can occur in an L1A granule is 204.

2. L1B Granules.  The L1B EV granules group the MODIS bands slightly
   differently than the above L1A SDSs.  Band 26 (in the L1A "night" group)
   is placed in the 1km reflective solar band (RefSB) group.

   L1B SDSs                   MODIS bands
   ------------       ------------------------------------------------------
   "EV_250_RefSB"     1,2
   "EV_500_RefSB"     3,4,5,6,7
   "EV_1KM_RefSB"     8,9,10,11,12,13lo,13hi,14lo,14hi,15,16,17,18,19,26
   "EV_1KM_Emissive"  20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36

   The L1B SDSs are placed in three L1B output HDF files which are based on
   resolution (see the file specifications).  This includes aggregation of
   the higher-resolution SDSs to a lower resolution.

   In the code, when we fill the L1B SDSs from the L1A data, there is logic
   which specifically deals with MODIS band 26. 

3. Detectors.  There are different numbers of detectors based on band
   resolution.  Also note that bands 13 and 14 also each have two sets
   of detectors.  In terms of the L1B SDSs:

   L1B SDSs           # of detectors      total in all bands of this res.
   ----------         --------------      ------------------------------
   "EV_250_RefSB"          40                        80
   "EV_500_RefSB"          20                       100
   "EV_1KM_RefSB"          10                       150
   "EV_1KM_Emissive"       10                       160

   Thus, the total number of detectors is 490.

4. Coordinate system:

      |---------------- ALONG SCAN ---------------->
      |
      |      --------------------------------------
      |      --------------------------------------
      |      --------------------------------------
    ALONG    -------------- SCAN LINES ------------
    TRACK    --------------------------------------
      |      --------------------------------------
      |      --------------------------------------
      |      --------------------------------------
     \|/     --------------------------------------

   For one sweep of the MODIS scan mirror, the number of scan lines generated
   equals the number of detectors (which is dependent on band). The number of
   samples generated along a scan line is also dependent on band (the 250m
   bands being sampled at 0.25km IFOV, the 500m bands being sampled at 0.5km
   IFOV and the 1km bands being sampled at 1km IFOV).

   Thus, for a given sweep of the mirror, 250m resolution bands have 4 times
   the number of scan lines as 1km resolution bands since they have 4 times
   the number of detectors.  Similarly, for a given aperture width (such as
   the earth view), 250m resolution bands have 4 times the number of samples
   ("FRAMES") as 1km resolution bands.

   In the macros below, the following terminology is used:
     "TRACK_DIM"    refers to scan lines in the along-track direction
     "FRAME_DIM"    refers to samples in the along-scan direction.
     "1km_FRAMES"   refers to samples taken at the rate of 1km IFOV
     "500m_FRAMES"  refers to samples taken at the rate of 0.5km IFOV
     "250m_FRAMES"  refers to samples taken at the rate of 0.25km IFOV

5. Aperture widths.  The following table relates the basic aperture widths 
   in terms of 1km frames:

      Aperture       acronym       # 1km_frames
     ----------      -------       ------------
     earth view        EV             1354
     space view        SV               50
     black body        BB               50
     SRCA              SRCA             10
     solar diffuser    SD               50      

!END***********************************************************************
*/

/***************************************************************************
Band/Detector/Frame Geometry constants
--------------------------------------
  MAX_NUM_SCANS             Maximum number of scans in one 5-minute granule.
  NUM_BANDS                 Number of bands (see note 1 above).
  NUM_DETECTORS             Total number of detectors for all bands.
  NUM_HIGH_RESOLUTION_DETECTORS   Total number of detectors in 250m and 500m bands.
  NUM_HIGH_RESOLUTION_SUBFRAMES   Total number of subframes in 250m and 500m bands.
  NUM_REFLECTIVE_DETECTORS  Total number of detectors for the L1B reflective
                            solar bands (RefSB).
  NUM_EMISSIVE_DETECTORS    Total number of detectors for the L1B emissive
                            bands (16 emissive bands, 10 detectors each).
  NUM_MIRROR_SIDES          Number of mirror sides (2).
  NUM_250M_BANDS            Number of 250m-resolution bands (2).
  NUM_500M_BANDS            Number of 500m-resolution bands (5).
  NUM_1000M_REFL_BANDS      Number of L1B reflective solar bands that have
                            1km resolution.
  NUM_REFLECTIVE_BANDS      Total number of L1B reflective solar bands.  This
                            is the sum of NUM_250M_BANDS, NUM_500M_BANDS and
                            NUM_1000M_REFL_BANDS.
  NUM_1000M_EMISS_BANDS     Number of L1B 1km-resolution emissive bands.
                            Since all emissive bands are 1km-resolution bands,
                            this macro means the same as the one below.
  NUM_EMISSIVE_BANDS        Number of L1B emissive bands.
  NUM_1000M_DAY_BANDS       Number of bands in the "EV_1km_day" L1A SDS.
  NUM_1000M_NIGHT_BANDS     Number of bands in the "EV_1km_night" L1A SDS.
  DETECTORS_PER_1KM_BAND    Number of detectors in any 1km-resolution band.
  DETECTORS_PER_500M_BAND   Number of detectors in any 500m-resolution band.
  DETECTORS_PER_250M_BAND   Number of detectors in any 250m-resolution band.
  NUM_1KM_SUBSAMP           Number of subsamples per frame in any 
                            1km-resolution band.
  NUM_500M_SUBSAMP          Number of subsamples per frame in any 
                            500m-resolution band.
  NUM_250M_SUBSAMP          Number of subsamples per frame in any 
                            250m-resolution band.
  NUM_REFL_INDICES          Total number of valid reflective 
                            band/detector/subframe/mirror side combinations.
  MODIS_BAND26_INDEX        For the set of 38 bands of this code (indexed 0
                            through 37), this is the index for MODIS band 26.
  MODIS_BAND26_INDEX_AT_RES For the set of 17 L1A 1km night bands (indexed
                            0 through 16), this is the index corresponding to
                            MODIS band 26.
  MODIS_BAND20_INDEX        For the set of 38 bands of this code (indexed 0
                            through 37), this is the index for MODIS band 20.
  MAX_250M_TRACK_DIM        Maximum number of scan lines spanning the along-
                            track dimension of the entire granule for any 250m
                            resolution band. This equals the maximum number of
                            scans times the number of detectors for 250m
                            resolution).
  MAX_500M_TRACK_DIM        Maximum number of scan lines for any 500m
                            resolution band. 
  MAX_1KM_TRACK_DIM         Maximum number of scan lines for any 1km
                            resolution band.
  MAX_SWIR_TRACK_DIM        Maximum number of scan lines for any SWIR band.
  MAX_250M_OBC_FRAME_DIM    Maximum number of 1km frames in any of the OBC
                            data sets (50) times the number of subsamples per
                            frame for a 250m-resolution band (4)
  MAX_500M_OBC_FRAME_DIM    Maximum number of 1km frames in any of the OBC
                            data sets (50) times the number of subsamples per
                            frame for a 500m-resolution band (2)
  MAX_1KM_OBC_FRAME_DIM     Maximum number of 1km frames in any of the OBC
                            data sets (50).
  SRCA_250M_FRAME_DIM       Number of 1km frames for the SRCA data sets (10)
                            times the number of subsamples per frame for a
                            250m-resolution band (4).
  SRCA_500M_FRAME_DIM       Number of 1km frames for the SRCA data sets (10)
                            times the number of subsamples per frame for a
                            500m-resolution band (2).
  SRCA_1KM_FRAME_DIM        Number of 1km frames for the SRCA data sets (10).
  TOTAL_1km_FRAMES          Total number of EV, SV, BB, SRCA and SD 1km frames.
  EV_1km_FRAMES             Number of earth view 1km frames.
  EV_250m_FRAMES            Number of earth view 250m frames
                            (EV_1km_FRAMES * 4).
  EV_500m_FRAMES            Number of earth view 500m frames
                            (EV_1km_FRAMES * 2).
  SD_1km_FRAMES             Number of solar diffusor 1km frames.
  SRCA_1km_FRAMES           Number of SRCA 1km frames.
  BB_1km_FRAMES             Number of black-body 1km frames.
  SV_1km_FRAMES             Number of space view 1km frames.
  NADIR_1km_FRAME_NUM       The 1km frame number nominally desginated as
                            the center frame number of the earth-view window.
  MAX_DETECTORS_PER_BAND    Maximum number of detectors in any band 
                            (occurs in the 250m bands).
  MAX_SAMPLES_PER_BAND      Maximum number of subsamples in a 1km frame for
                            any band (occurs in the 250m bands).
****************************************************************************/
#define    NUM_SCANS_DIM_NAME       "nscans"

/* Normally MAX_NUM_SCANS should be 204. However, the number of scans could
 * reach 208 due to a L1A abnormality. So, define it to be 208.
 */

#define    MAX_NUM_SCANS           1000
#define    NUM_BANDS                 38
#define    NUM_DETECTORS            490
#define    NUM_HIGH_RESOLUTION_DETECTORS            180
#define    NUM_HIGH_RESOLUTION_SUBFRAMES            520
#define    NUM_REFLECTIVE_DETECTORS 330
#define    NUM_EMISSIVE_DETECTORS   160
#define    NUM_MIRROR_SIDES           2
#define    NUM_250M_BANDS             2
#define    NUM_500M_BANDS             5
#define    NUM_1000M_REFL_BANDS      15
#define    NUM_REFLECTIVE_BANDS      22
#define    NUM_1000M_EMISS_BANDS     16
#define    NUM_EMISSIVE_BANDS        16
#define    NUM_1000M_DAY_BANDS       14
#define    NUM_1000M_NIGHT_BANDS     17
#define    DETECTORS_PER_1KM_BAND    10
#define    DETECTORS_PER_500M_BAND   20      
#define    DETECTORS_PER_250M_BAND   40 
#define    NUM_1KM_SUBSAMP            1
#define    NUM_500M_SUBSAMP           2
#define    NUM_250M_SUBSAMP           4
#define    NUM_REFL_INDICES        1340
#define    MODIS_BAND5_INDEX          4    
#define    MODIS_BAND26_INDEX        27
#define    MODIS_BAND26_INDEX_AT_RES  6
#define    MODIS_BAND20_INDEX        21

#define    L1B_1KM_EMISS_BAND28_INDEX  7    /* emissive band group */

/* The four macros below are related to the macro MAX_NUM_SCANS. */

#define    MAX_250M_TRACK_DIM      MAX_NUM_SCANS*DETECTORS_PER_250M_BAND  /* 8320 */  /* 208 * 40 */
#define    MAX_500M_TRACK_DIM      MAX_NUM_SCANS*DETECTORS_PER_500M_BAND  /* 4160 */
#define    MAX_1KM_TRACK_DIM       MAX_NUM_SCANS*DETECTORS_PER_1KM_BAND   /* 2080 */
#define    MAX_SWIR_TRACK_DIM      MAX_NUM_SCANS*DETECTORS_PER_500M_BAND  /* 4160 */

#define    MAX_250M_OBC_FRAME_DIM   200   /* 50 * 4 */
#define    MAX_500M_OBC_FRAME_DIM   100   /* 50 * 2 */
#define    MAX_1KM_OBC_FRAME_DIM     50 
#define    SRCA_250M_FRAME_DIM       40    /* 10 * 4 */
#define    SRCA_500M_FRAME_DIM       20
#define    SRCA_1KM_FRAME_DIM        10
#define    SD_1km_FRAMES             50
#define    SRCA_1km_FRAMES           10
#define    BB_1km_FRAMES             50
#define    SV_1km_FRAMES             50
#define    EV_1km_FRAMES           1354
#define    EV_250m_FRAMES          5416
#define    EV_500m_FRAMES          2708 
#define    TOTAL_1km_FRAMES        1514
#define    NADIR_1km_FRAME_NUM      677
#define    MAX_DETECTORS_PER_BAND    40
#define    MAX_SAMPLES_PER_BAND       4


/* Three correction switch in Emissive Band */
typedef enum { OFF, ON } switch_t ;

#define  SCAN_TYPE_TEXT_SIZE       10
#define  COMMON_TEXT_SIZE          10
#define  MAX_ERROR_MESSAGE_LENGTH  1024 

#define  NUM_BB_THERMISTORS        12
#define  NUM_FOCAL_PLANES           4

#define SCAN_QUALITY_ARRAY_NUM_ELEMENTS 4

#define PGS_PI          3.14159265358979323846

#define PGS_EXP    	2.7182818284590452354

#define SAME       0 /*Used with strcmp()*/

/*
 * Macros for L1B sector DNs, to be used with the Check_Valid_Range
 * function.  Note that the HDF default is currently used in the
 * L1A code (does not comply with file specs, but we must allow it).
 *   L1A_DN_SDS_LB -- lower bound
 *   L1A_DN_SDS_UB -- upper bound
 *   L1A_DN_SDS_FV -- fill value
 * These are strings since that is the interface to Check_Valid_Range.
 */

#define L1A_DN_SDS_LB        "-1"
#define L1A_DN_SDS_UB        "4095"
#define L1A_DN_SDS_FV        "-32767"  /* HDF default */

#define HDF_DEFAULT_FILLVALUE    -32767
#define SATURATED_DN      4095

#define MISSING_L1B_DATA  65535 
#define MISSING_L1B_UI      255
#define MISSING_L1A_FLAG    -1 

#define L1A_SCAN_DATA_MISSING_SI 65535
#define L1A_DN_MISSING_SI        65534
#define SATURATED_DETECTOR_SI    65533
#define ZERO_POINT_DN_SI         65532
#define DEAD_DETECTOR_SI         65531
#define RSB_DN_STAR_BELOW_MIN_SI 65530    
#define TEB_OR_RSB_GT_MAX_SI     65529
#define AGGREGATION_FAIL_SI      65528
#define SECTOR_ROTATION_SI       65527
#define TEB_B1_NOT_CALCULATED    65526
#define DEAD_SUBFRAME_SI         65525
#define UNABLE_CALIBRATE_SI      65524
#define UNRESCALED_HIGH_SI       65521
#define RESCALED_L1B_SI          65520
#define NAD_CLOSED_UPPER_SI      65500
#define L1A_SCAN_DATA_MISSING_UI 255
#define BAD_DATA_UI              15

#define TOLERANCE      1.0E-20  /* defined zero */

#define True                         1
#define False                        0

/* Focal Plane Indices */

#define VIS  0
#define NIR  1
#define SWIR 2
#define LWIR 3

/* satellite ID */

#define TERRA 0
#define AQUA  1
#define INVALID_SATELLITE_ID -1

/* MOD_PR02 exit code */

extern int32 MOD_PR02_Failure_Exit_Code;

/*****************************************************************************
 The following defines macros for the number of L1A or L1B "resolutions" (or
 groupings of bands) and the indices representing the ordering of the
 resolutions in arrays.  Even though INDEX_1000M_DAY and INDEX_1000M_REFL are
 both the same value, they represent different band groupings because of
 MODIS band 26 (see earlier comments).  Currently, the code uses these macros
 somewhat inconsistently, so beware.  Similar comments apply to
 INDEX_1000M_NIGHT and INDEX_1000M_EMISS.
 *****************************************************************************/
typedef  int8  boolean;

typedef enum {
  INDEX_250M,
  INDEX_500M,
  INDEX_1000M_DAY,
  INDEX_1000M_NIGHT,
  NUM_L1A_RESOLUTIONS
} resolution_index_t;
#define INDEX_1000M_REFL  INDEX_1000M_DAY
#define INDEX_1000M_EMISS INDEX_1000M_NIGHT

/*****************************************************************************
 The following defines macros for the number of earth-view products generated
 and indices for each product.
 *****************************************************************************/

typedef enum {
  INDEX_L1B_250m,
  INDEX_L1B_500m,
  INDEX_L1B_1km,
  NUM_L1B_EV_FILES
} L1B_EV_FILE_index_t;

/***************************************************************************
 L1B_BANDS_AT_RES         Number of bands in each L1B "resolution"
                          (group of bands -- see notes at top)
 L1A_BANDS_AT_RES         Number of bands in each L1A "resolution"
                          (group of bands -- see notes at top)
 DETECT_PER_BAND_AT_RES   Number of detectors in any band of the 4 different
                          "resolutions" (same for both L1A and L1B even though
                          the "groupings" are defined slightly differently).
 BAND_RATIO_AT_RES        This is the ratio of the #samples at a resolution to
                          the number of samples at a 1km resolution (also the
                          ratio of the number of detectors relative to the
                          number of detectors for a 1km resolution band).
                          Same for both L1A and L1B "resolutions".
 **************************************************************************/
extern int32 SUBSAMPLES_AT_RES[NUM_L1A_RESOLUTIONS];

extern int16 L1B_BANDS_AT_RES[NUM_L1A_RESOLUTIONS];

extern int16 L1A_BANDS_AT_RES[NUM_L1A_RESOLUTIONS]; 

extern int16 DETECT_PER_BAND_AT_RES [NUM_L1A_RESOLUTIONS]; 

extern int16 BAND_RATIO_AT_RES[NUM_L1A_RESOLUTIONS];

typedef enum {
  Orbital_Node,
  Number_Of_Scans,
  Number_Day_Scans,
  Number_Night_Scans,
  Total_Frames,
  Earth_Frames,
  SD_Frames,
  SRCA_Frames,
  BB_Frames,
  SV_Frames,
  Scan_Type,
  Incomplete_Scans,
  Missing_Packets,
  Packets_With_Bad_CRC,
  Discarded_Packets, 
  NUM_L1A_ATTRIBUTES
} L1A_attr_index_t;

typedef enum {
  SD_INDEX,
  SRCA_INDEX,
  BB_INDEX,
  SV_INDEX,
  EV_INDEX,
  NUM_TARGETS
} target_index_t;

typedef enum {
  BAND20, BAND21, BAND22, BAND23, BAND24, BAND25, BAND26, BAND27, BAND28,
  BAND29, BAND30, BAND31, BAND32, BAND33, BAND34, BAND35, BAND36, NUM_NIGHT_BANDS 
} night_band_index_t;

/******************************************************************************
 * Define macros for granule average quantities.
 *
 * NOTE 1: the code uses loops when assigning some of these quantities.
 * In particular, the number and order of BB, INS, CAV and RC temperatures
 * must match the arrays defined in Temperatures_t in PreprocessP.h
 * and the calculations in Read_Convert_Temperatures.
 * Consequently, these cannot be arbitrarily rearranged without inspecting
 * and possibly changing the code.
 *
 * NOTE 2: some elements of "granule_averages" (a structure member in
 * QA_Common_t) are unused.  This variable is currently
 * initialized in Calculate_Temp_QA, where all non-zero values
 * are assigned.  The initialization sets the unused elements to zero.
 * If other granule averages are added, then the location of the
 * initialization may need to be changed.
 */

typedef enum {
  GRAN_AVG_TP_BB_TEMP01,        /* BB Thermistor 1 */
  GRAN_AVG_TP_BB_TEMP02,        /* BB Thermistor 2 */
  GRAN_AVG_TP_BB_TEMP03,        /* BB Thermistor 3 */
  GRAN_AVG_TP_BB_TEMP04,        /* BB Thermistor 4 */
  GRAN_AVG_TP_BB_TEMP05,        /* BB Thermistor 5 */
  GRAN_AVG_TP_BB_TEMP06,        /* BB Thermistor 6 */
  GRAN_AVG_TP_BB_TEMP07,        /* BB Thermistor 7 */
  GRAN_AVG_TP_BB_TEMP08,        /* BB Thermistor 8 */
  GRAN_AVG_TP_BB_TEMP09,        /* BB Thermistor 9 */
  GRAN_AVG_TP_BB_TEMP10,        /* BB Thermistor 10 */
  GRAN_AVG_TP_BB_TEMP11,        /* BB Thermistor 11 */
  GRAN_AVG_TP_BB_TEMP12,        /* BB Thermistor 12 */
  GRAN_AVG_TA_AO_VIS_FPAE,      /* (FPA1) VIS Focal Plane */
  GRAN_AVG_TA_AO_NIR_FPAE,      /* (FPA2) NIR Focal Plane */
  GRAN_AVG_TA_RC_SMIR_CFPAE,    /* (FPA3) SMIR Focal plane Temp */
  GRAN_AVG_TA_RC_LWIR_CFPAE,    /* (FPA4) LWIR Focal plane Temp */
  GRAN_AVG_TP_SA_RCT1_MIRE,     /* Scan Mirror Temp 1 */
  GRAN_AVG_TP_SA_RCT2_MIRE,     /* Scan Mirror Temp 2 */
  GRAN_AVG_TP_SA_A_MTR,         /* Scan motor temperature A */
  GRAN_AVG_TP_MF_CALBKHD_SR,    /* (CAV1) Cal bulkhead below SRCA mount */
  GRAN_AVG_TP_SR_SNOUT,         /* (CAV2) SRCA Bulkhead Temp */
  GRAN_AVG_TP_MF_Z_BKHD_BB,     /* (CAV3) Mid zenith bulkhead near BB */
  GRAN_AVG_TP_MF_CVR_OP_SR,     /* (CAV4) cover opposite SRCA */
  GRAN_AVG_TP_AO_SMIR_OBJ,      /* (INS1) SMIR Objective Lens Temp */
  GRAN_AVG_TP_AO_LWIR_OBJ,      /* (INS2) LWIR Objective Lens Temp */
  GRAN_AVG_TP_AO_SMIR_LENS,     /* (INS3) SMIR Eye Assy Temp */
  GRAN_AVG_TP_AO_LWIR_LENS,     /* (INS4) LWIR Eye Assy Temp */
  GRAN_AVG_TA_RC_CS,            /* (RC1) Cold stage Temp */
  GRAN_AVG_TA_RC_CS_OG,         /* (RC2) Cold stage outgas Temp */
  GRAN_AVG_TA_RC_IS,            /* (RC3) Intermediate Stage Temp */
  GRAN_AVG_TA_RC_IS_OG,         /* (RC4) Intermediate Stage outgas Temp */
  GRAN_AVG_TA_RC_OS_OG,         /* (RC5) Outerstage outgas Temp */
  GRAN_AVG_VR_RC_LW_FPA_HTR,    /* LWIR heater voltage */
  NUM_GRAN_AVERAGES
} gran_average_def_t;

/* Define the dimension of the attribute.  Unused elements are set to
 * zero within the attribute.
 */

#define MAX_NUM_GRAN_AVERAGES   50

/* Define the number of RC values in the above definition.  This value
 * is also used in Preprocess when calculating and temporarily storing
 * the temperatures.
 */

#define NUM_T_RC_VALUES  5

/* Define a bad value for voltages
 */

#define VOLTAGE_BAD_VALUE  -1000.

/***************************************************************************/

/*--------------------------------------------------------------------------
  Data Structures
--------------------------------------------------------------------------*/

  /* Upper bound on length of ProcessingEnvironment archive
   * metadata item
   */

#define  MAX_PROCESSING_ENVIRONMENT_STRLEN 200

#define MAX_RUNTIME_PARAM_SIZE 256

typedef struct {
  char SatelliteInstrument[MAX_RUNTIME_PARAM_SIZE];
  char ReprocessingPlanned[MAX_RUNTIME_PARAM_SIZE];
  char ReprocessingActual[MAX_RUNTIME_PARAM_SIZE];
  char PGE02_Version[MAX_RUNTIME_PARAM_SIZE];
  char MCST_LUT_Version[MAX_RUNTIME_PARAM_SIZE];
  char Write_Night_Mode_HiRes_Data[MAX_RUNTIME_PARAM_SIZE];
  char ProcessingEnvironment[MAX_RUNTIME_PARAM_SIZE];
  char ProcessingCenter[MAX_RUNTIME_PARAM_SIZE];
  } Run_Time_Parameters_t;

typedef  struct {
  int32     v_id;
  int32     sd_id;
  int32     satellite_id;
  int32     num_scans;
  int32     num_day_scans;
  int16     MirrorSide[MAX_NUM_SCANS]; 
  char      ScanType[MAX_NUM_SCANS][SCAN_TYPE_TEXT_SIZE];
  float64   EVStartTime_TAIsecond[MAX_NUM_SCANS];
  int32     scan_quality[MAX_NUM_SCANS][SCAN_QUALITY_ARRAY_NUM_ELEMENTS];
  float64   data_collection_time;

  /* For L1B metadata */

  int32     num_night_scans;
  int32     incomplete_scans;
  int32     max_ev_frames;
  int32     Extract_Pixel_Offset;
  int32     Extract_Pixel_Count;
  int32     Extract_Line_Offset;
  int32     Extract_Line_Count;

  /* For identifying completely missing scans */

  boolean   missing_scan[MAX_NUM_SCANS];

  } L1A_granule_t;

typedef struct 
{
  int32 sds_id[NUM_L1A_RESOLUTIONS];

  int16 EV_250m     [DETECTORS_PER_250M_BAND]
                    [NUM_250M_BANDS]
                    [EV_250m_FRAMES];

  int16 EV_500m     [DETECTORS_PER_500M_BAND]
                    [NUM_500M_BANDS]
                    [EV_500m_FRAMES];

  int16 EV_1km_day  [DETECTORS_PER_1KM_BAND]
                    [NUM_1000M_DAY_BANDS]
                    [EV_1km_FRAMES];

  int16 EV_1km_night[DETECTORS_PER_1KM_BAND]
                    [NUM_1000M_NIGHT_BANDS]
                    [EV_1km_FRAMES];
} L1A_Scan_t;

/* Level 1B scale_offset data */

typedef struct {
  float32 dn_star_Max[NUM_REFLECTIVE_BANDS];
  float32 dn_star_Min[NUM_REFLECTIVE_BANDS]; 

  float32 rad_scale_RefSB[NUM_REFLECTIVE_BANDS]; 
  float32 rad_offset_RefSB[NUM_REFLECTIVE_BANDS]; 
  float32 refl_scale_RefSB[NUM_REFLECTIVE_BANDS];
  float32 refl_offset_RefSB[NUM_REFLECTIVE_BANDS];

  float32 counts_scale_RefSB[NUM_REFLECTIVE_BANDS];
  float32 counts_offset_RefSB[NUM_REFLECTIVE_BANDS];

  float32 rad_scale_Emiss[NUM_EMISSIVE_BANDS];
  float32 rad_offset_Emiss[NUM_EMISSIVE_BANDS];
} L1B_ScaleOffset_t;


typedef struct 
{
  float32 m1_des_sq[NUM_REFLECTIVE_BANDS]
                   [MAX_DETECTORS_PER_BAND]
                   [MAX_SAMPLES_PER_BAND]
                   [NUM_MIRROR_SIDES];

  float32 m1_des_sq_max[NUM_REFLECTIVE_BANDS];

  float32 RVS_250m[NUM_250M_BANDS]
                  [DETECTORS_PER_250M_BAND]
                  [EV_250m_FRAMES]
                  [NUM_MIRROR_SIDES];

  float32 RVS_500m[NUM_500M_BANDS]
                  [DETECTORS_PER_500M_BAND]
                  [EV_500m_FRAMES]
                  [NUM_MIRROR_SIDES];

  float32 RVS_1km_RefSB[NUM_1000M_REFL_BANDS]
                       [DETECTORS_PER_1KM_BAND]
                       [EV_1km_FRAMES]
                       [NUM_MIRROR_SIDES];

  float32 u2[NUM_REFLECTIVE_DETECTORS]
            [EV_1km_FRAMES]
            [NUM_MIRROR_SIDES]; 

} RSB_Cal_Coeff_t;

typedef  struct {
  float32 RVS_1km_Emiss_EV[NUM_EMISSIVE_DETECTORS]
                          [EV_1km_FRAMES]
                          [NUM_MIRROR_SIDES];
  float32 RVS_1km_Emiss_SV[NUM_EMISSIVE_DETECTORS]
                          [NUM_MIRROR_SIDES];
  float32 RVS_1km_Emiss_BB[NUM_EMISSIVE_DETECTORS]
                          [NUM_MIRROR_SIDES];
  float32 sigma_RVS_Emiss_EV[NUM_EMISSIVE_DETECTORS]
                            [EV_1km_FRAMES]
                            [NUM_MIRROR_SIDES];
} Emiss_Cal_Coeff_t;

typedef  struct {
  int32     sw_f_id[NUM_L1B_EV_FILES];
  int32     v_id[NUM_L1B_EV_FILES];
  int32     sd_id[NUM_L1B_EV_FILES];
  int32     num_scans;
  int32     num_day_scans;
  
  RSB_Cal_Coeff_t   RSB_Cal_Coeff;
  Emiss_Cal_Coeff_t Emiss_Cal_Coeff;
  L1B_ScaleOffset_t SO;

  float32 b26_fr_b5_scaled_corr[DETECTORS_PER_1KM_BAND];
  
  /* For L1B metadata */

  float32 validEVPercent[NUM_BANDS];
  float32 satEVPercent[NUM_BANDS];
  float32 missEVPercent[NUM_BANDS];
  uint32  elecRedVec[2];
  int8    FPSetPointState[NUM_FOCAL_PLANES];

  uint32 total_pixels[NUM_BANDS];
  uint32 valid_pixels[NUM_BANDS];
  uint32 saturated_pixels[NUM_BANDS];
  uint32 missing_pixels[NUM_BANDS]; 
  uint32 interpolated_pixels[NUM_BANDS];
  uint32 dead_detector_pixels[NUM_BANDS];
  uint32 dead_subframe_pixels[NUM_BANDS];
  uint32 negative_value_below_noise_pixels[NUM_BANDS];
  int16  bad_data_flag[NUM_BANDS];
  float32 Earth_Sun_Dist;
} L1B_granule_t;

typedef struct 
{
  uint16 EV_250m_RefSB  [NUM_250M_BANDS]
                        [DETECTORS_PER_250M_BAND]
                        [EV_250m_FRAMES];

  uint16 EV_500m_RefSB  [NUM_500M_BANDS]
                        [DETECTORS_PER_500M_BAND]
                        [EV_500m_FRAMES];

  uint16 EV_1km_RefSB   [NUM_1000M_REFL_BANDS]
                        [DETECTORS_PER_1KM_BAND]
                        [EV_1km_FRAMES];

  uint16 EV_1km_Emissive[NUM_1000M_EMISS_BANDS]
                        [DETECTORS_PER_1KM_BAND]
                        [EV_1km_FRAMES];
} L1B_Scan_SI_t;

typedef struct 
{
  uint8 EV_250m_RefSB_UI  [NUM_250M_BANDS]
                          [DETECTORS_PER_250M_BAND]
                          [EV_250m_FRAMES];

  uint8 EV_500m_RefSB_UI  [NUM_500M_BANDS]
                          [DETECTORS_PER_500M_BAND]
                          [EV_500m_FRAMES];

  uint8 EV_1km_RefSB_UI   [NUM_1000M_REFL_BANDS]
                          [DETECTORS_PER_1KM_BAND]
                          [EV_1km_FRAMES];

  uint8 EV_1km_Emissive_UI[NUM_1000M_EMISS_BANDS]
                          [DETECTORS_PER_1KM_BAND]
                          [EV_1km_FRAMES];
} L1B_Scan_UI_t;

/************************* Begin Band 26 Section **************************/
/*
 * If the macro WRITE_BAND_26_SDS is defined, then the MODIS band 26 SDSs will
 * be written to the EV 1km file.  Two SDSs are written: (1) scaled integer,
 * (2) uncertainty index.  These are written regardless of whether the scan
 * is day or night.
 */

/*** UNCOMMENT THE FOLLOW STATEMENT TO ENABLE WRITING BAND 26 SDS ***/

#define WRITE_BAND_26_SDS  

#ifdef WRITE_BAND_26_SDS

#define BAND_26_SI_SDS_NAME      "EV_Band26"
#define BAND_26_UI_SDS_NAME      "EV_Band26_Uncert_Indexes"
#define BAND_26_SI_SDS_LONG_NAME "Earth View Band 26 Scaled Integers"
#define BAND_26_UI_SDS_LONG_NAME "Earth View Band 26 Uncertainty Indexes"

typedef struct {
  int32 SI_sds_id;
  int32 UI_sds_id;
  uint16 SI[DETECTORS_PER_1KM_BAND][EV_1km_FRAMES];
  uint8  UI[DETECTORS_PER_1KM_BAND][EV_1km_FRAMES];
} Band_26_t;

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

typedef struct {
/*-------------------------
  SI == Scaled Integer
  UI == Uncert Index
  SU == Samples Used
-------------------------*/
  int32    SI_sds_id[NUM_L1A_RESOLUTIONS];
  int32    EV_250m_Aggr500m_RefSB_sds_id;
  int32    EV_250m_Aggr1km_RefSB_sds_id;
  int32    EV_500m_Aggr1km_RefSB_sds_id;
  int32    UI_sds_id[NUM_L1A_RESOLUTIONS];
  int32    EV_250m_Aggr500m_RefSB_UI_sds_id;
  int32    EV_250m_Aggr1km_RefSB_UI_sds_id;
  int32    EV_500m_Aggr1km_RefSB_UI_sds_id;
  int32    EV_250m_Aggr500m_RefSB_SU_sds_id;
  int32    EV_250m_Aggr1km_RefSB_SU_sds_id;
  int32    EV_500m_Aggr1km_RefSB_SU_sds_id;
  L1B_Scan_SI_t   SI;
  /* uint16   ***SI[NUM_L1A_RESOLUTIONS];*/    /*[band][det][frame]*/
  uint16   EV_250m_Aggr500m_RefSB   [NUM_250M_BANDS]
                                    [DETECTORS_PER_500M_BAND]
                                    [EV_500m_FRAMES];

  uint16   EV_250m_Aggr1km_RefSB    [NUM_250M_BANDS]
                                    [DETECTORS_PER_1KM_BAND]
                                    [EV_1km_FRAMES];

  uint16   EV_500m_Aggr1km_RefSB    [NUM_500M_BANDS]
                                    [DETECTORS_PER_1KM_BAND]
                                    [EV_1km_FRAMES];

  L1B_Scan_UI_t   UI;

  /* uint8    ***UI[NUM_L1A_RESOLUTIONS];
 */   /*
                                    [band][det][frame]*/
  uint8    EV_250m_Aggr500m_RefSB_UI[NUM_250M_BANDS]
                                    [DETECTORS_PER_500M_BAND]
                                    [EV_500m_FRAMES];
         
  uint8    EV_250m_Aggr1km_RefSB_UI [NUM_250M_BANDS]
                                    [DETECTORS_PER_1KM_BAND]
                                    [EV_1km_FRAMES];
           
  uint8    EV_500m_Aggr1km_RefSB_UI [NUM_500M_BANDS]
                                    [DETECTORS_PER_1KM_BAND]
                                    [EV_1km_FRAMES];
             
  int8     EV_250m_Aggr500m_RefSB_SU[NUM_250M_BANDS]
                                    [DETECTORS_PER_500M_BAND]
                                    [EV_500m_FRAMES];

  int8     EV_250m_Aggr1km_RefSB_SU [NUM_250M_BANDS]
                                    [DETECTORS_PER_1KM_BAND]
                                    [EV_1km_FRAMES];

  int8     EV_500m_Aggr1km_RefSB_SU [NUM_500M_BANDS]
                                    [DETECTORS_PER_1KM_BAND]
                                    [EV_1km_FRAMES];

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
/*
 * This structure member holds arrays dimensioned appropriately for one
 * scan of the Band 26 SDSs.
 */
  Band_26_t Band26;
#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

  int16   band_X;
  float32 dn_X[DETECTORS_PER_1KM_BAND][EV_1km_FRAMES]; 
} L1B_Scan_t;

typedef struct 
{
  uint8 noise_T_bb[12];
  uint8 noise_T_bb_avg;
  uint8 noise_T_lwir;
  uint8 noise_T_mwir;
  uint8 noise_T_mir1;
  uint8 noise_T_mir2;
  uint8 noise_T_mir_avg;
  uint8 noise_T_ins;
  uint8 noise_T_cav;
  uint8 NEdL[NUM_EMISSIVE_DETECTORS];        /*noise in system*/
  uint8 change_b1[NUM_EMISSIVE_DETECTORS];
  uint8 num_thermistor_outliers[MAX_NUM_SCANS];     /* range [0-12] */
  int8  change_dc_restore[MAX_NUM_SCANS]
                         [NUM_EMISSIVE_BANDS]
                         [DETECTORS_PER_1KM_BAND];
  int8  emiss_band_bad_data[NUM_EMISSIVE_BANDS];
  int8  moon_in_SV_KOB_TEB[MAX_NUM_SCANS]; 
} QA_Emiss_t;

typedef struct 
{
  uint8 noise_T_vis;
  uint8 noise_T_nir;
  int8  change_dc_restore_250m[MAX_NUM_SCANS][NUM_250M_BANDS]
                              [DETECTORS_PER_250M_BAND];
  int8  change_dc_restore_500m[MAX_NUM_SCANS][NUM_500M_BANDS]
                              [DETECTORS_PER_500M_BAND];
  int8  change_dc_restore_1km[MAX_NUM_SCANS][NUM_1000M_REFL_BANDS]
                             [DETECTORS_PER_1KM_BAND];
  int8  refl_band_bad_data[NUM_REFLECTIVE_BANDS];
  int8  moon_in_SV_KOB_RSB[MAX_NUM_SCANS];
  int8  all_SV_DN_bad[MAX_NUM_SCANS];
  int8  all_BB_DN_bad[MAX_NUM_SCANS]; 
} QA_Refl_t;

typedef struct 
{
  int8  all_l1b_error_flags_off;
  boolean missing_leading_granule;   /* True or False */
  boolean missing_trailing_granule;  /* True or False */
  boolean leading_granule_scan_gap;  /* True or False */
  boolean trailing_granule_scan_gap; /* True or False */
  boolean NAD_Door_Open[MAX_NUM_SCANS];
  boolean Sector_Rotation[MAX_NUM_SCANS];
  boolean Electronic_Anomaly[MAX_NUM_SCANS];
  int32 num_missing_scans;
  int32 num_rsb_at_night_scans;
  int32 num_missing_data_in_scans[NUM_DETECTORS];
  int32 num_dead_detector_EV_data[NUM_DETECTORS];
  int32 num_dead_subframe_EV_data[NUM_HIGH_RESOLUTION_DETECTORS];
  int32 num_sector_rotation_EV_data[NUM_DETECTORS];
  int32 num_saturated_EV_data[NUM_DETECTORS];
  int32 num_no_bg_DN_EV_data[NUM_DETECTORS];
  int32 num_moon_in_SVP_TEB_EV_data[NUM_DETECTORS];
  int32 num_bad_dn_star_star_RSB_EV_data[NUM_DETECTORS];
  int32 num_exceed_max_for_scaling[NUM_DETECTORS];
  int32 num_nadir_door_closed_EV_data[NUM_DETECTORS];
  int32 num_negative_b1[NUM_DETECTORS];
  uint32 bit_QA_flags_last_value;
  uint32 bit_QA_flags_change;
  float32 granule_averages[MAX_NUM_GRAN_AVERAGES];
} QA_Common_t;

typedef struct 
{
  QA_Common_t  QA_common;
  QA_Emiss_t   QA_emiss;
  QA_Refl_t    QA_refl;
} QA_Data_t;

/* 16 Bit Constants for Re-Scaling of Radiance */
#define DN15_SAT                32767
#define DN_MIN                  0

/************************ Error message utilities ************************/
/*
 * The following are global values to use for the "other_msg" input to
 * L1BErrorMsg. These errors may occur many times.
 */

extern char Invalid_MOD01_Msg[];  /* Invalid or out of bounds L1A data */

/*
 * Function protos
 */
int safe_strcat(
  char *buf,     /* A character buffer that we want to concatenate to */
  char *str,     /* The string to be concatenated */
  int  buflen    /* The known memory limitation of the character buffer */
);

void L1BErrorMsg(
       /* name of L1B function that error occurred in */    
  char          *L1B_location,    
    
       /* associated error code for this error */    
  PGSt_SMF_code code,             
    
       /* short (usually 1-line) description of error */    
  char          *input_message,   
    
       /* name of the associated function that failed */    
  char          *assoc_function,  
  
       /* associated LUN for the file being accessed */  
  int32         lun,            
  
       /* other message to add such as probable cause */  
  char          *other_msg,     
  
       /* flag to tell whether or not to call SMF_ERROR */  
  boolean       error_out       
  
);

void            SMF_ERROR  
                  (PGSt_SMF_code         code, 
                   char                  *messagestring);

void            Bad_L1A_Error_Out
                  (char                  *name,
                   char                  *message);

/********************** Other function Prototypes **************************/

PGSt_SMF_status Read_Run_Time_Parameters 
                  (Run_Time_Parameters_t *runtime_params);

PGSt_SMF_status Open_and_Read_L1A 
                  (Run_Time_Parameters_t *runtime_params,
                   L1A_granule_t         *L1A_Gran,
                   boolean               *skip_night_hi_res);


PGSt_SMF_status Get_Satellite_ID 
                  (PGSt_PC_Logical lun, 
                   int32                 *satellite_ID);


PGSt_SMF_status Read_L1A_EV_Scan 
                  (int16                 S,
                   L1A_granule_t         *L1A_Gran,
                   L1A_Scan_t            *L1A_Scan);


PGSt_SMF_status Aggregate_L1B 
                  (L1B_Scan_t            *L1B_Scan);


PGSt_SMF_status Fill_Dead_Detector_SI 
                  (boolean               isdaymode,
                   int8                  *dead_detector,
                   L1B_Scan_t            *L1B_Scan,
                   L1B_granule_t         *L1B_Gran,
                   QA_Common_t           *QA_Common);


PGSt_SMF_status Write_L1B_EV_Scan 
                  (int16                 S,
                   L1B_granule_t         *L1B_Gran,
                   L1B_Scan_t            *L1B_Scan,
                   boolean               isdaymode);


PGSt_SMF_status Close_L1A_Granule 
                  (L1A_granule_t         *L1A_Gran,
                   L1A_Scan_t            *L1A_Scan);


PGSt_SMF_status Close_L1B_Granule 
                  (L1B_granule_t         *L1B_Gran,
                   L1B_Scan_t            *L1B_Scan,
                   boolean               skip_night_hi_res);


#endif

