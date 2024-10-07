/*****************************************************************************

File: Granule.c

External functions:
   Aggregate_L1B
   Close_L1A_Granule
   Close_L1B_Granule
   Read_L1A_EV_Scan
   safe_strcat
   L1BErrorMsg
   SMF_ERROR
   Bad_L1A_Error_Out
   Write_L1B_EV_Scan
   Fill_Dead_Detector_SI
   Open_and_Read_L1A
   Get_Satellite_ID
   Read_Run_Time_Parameters

Other functions:
   Compute_Aggregates
   Write_L1B_SI_UI

Revision History:
   $Log: Granule.c,v $
   Revision 1.27  2012-06-27 15:58:24-04  xgeng
   correced an error message

   Revision 1.26  2008-01-31 16:45:17-05  ltan
   Relax the RVS correction limit range. Removed the ScanType of "Mixed" from the code.

   Revision 1.24  2006/10/30 14:53:59  ltan
   Changed to write PGEVersion value extracted from the PCF file (LUN 800500). Also check to make sure the PGE Version from the PCF matchs code macro PGE02_VERSION.

   Revision 1.22  2005/01/18 19:34:42  ltan
   Added new file attributes prefixed with HDFEOS_FractionalOffset

*****************************************************************************/

#include    "L1B_Tables.h"
#include    "GranuleP.h"
#include    "HDF_Lib.h"
#include    "PGS_PC.h"
#include    "PGS_TD.h"
#include    "PGS_MET.h"
#include    "PGS_Error_Codes.h"
#include    "FNames.h"
#include    <string.h>
#include    <math.h>
#include    <time.h>
#include    "granule_extract_metadata.h"

/* MOD_PR02 exit code */

int32 MOD_PR02_Failure_Exit_Code = 1;

/*----------------------------------------------------------------------------
     The following are external variables used by other modules
 ---------------------------------------------------------------------------*/

int16 L1B_BANDS_AT_RES[NUM_L1A_RESOLUTIONS] = {
  NUM_250M_BANDS,       NUM_500M_BANDS,
  NUM_1000M_REFL_BANDS, NUM_1000M_EMISS_BANDS
};

int16 L1A_BANDS_AT_RES[NUM_L1A_RESOLUTIONS] = {
  NUM_250M_BANDS,       NUM_500M_BANDS,
  NUM_1000M_DAY_BANDS,  NUM_1000M_NIGHT_BANDS
};

int16 DETECT_PER_BAND_AT_RES [NUM_L1A_RESOLUTIONS] = {
  DETECTORS_PER_1KM_BAND * 4, DETECTORS_PER_1KM_BAND * 2,
  DETECTORS_PER_1KM_BAND    , DETECTORS_PER_1KM_BAND
};

int16 BAND_RATIO_AT_RES[NUM_L1A_RESOLUTIONS] = {4,2,1,1};

int16 RFLAG;
int16 RSCL_FLAG;

PGSt_SMF_status  Aggregate_L1B (L1B_Scan_t  *L1B_Scan)
/*
!C**********************************************************************
!Description:
  For one scan of calibrated earth-view data, this routine performs
  spatial integration (or aggregation) of each of the higher
  resolution bands (250m or 500m) so that they appear at a lower
  resolution appropriate for the 500m or 1km L1B granule products.

!Input Parameters:
      L1B_Scan_t  *L1B_Scan   (->SI, ->UI)

!Output Parameters:
      L1B_Scan_t  *L1B_Scan   (->SI_aggr, ->UI_aggr, ->SU_aggr)

!Revision History:
 (continue at top of the file)

 Revision 02.11 July 1998
 Modified the routine to match the change of L1B_Scant_t change.
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)

 Revision 01.00 April 1997
 Initial development
 (based on the routine Compute_Aggregates, which 
 is a C translation of the FORTRAN code written 
 by Vicky Lin & Rich Hucek)
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
  int16   B         = 0;
  int16   line_dim  = 0;
  int16   pixel_dim = 0;

  /*------------------
    250_to_500m
  -------------------*/
  line_dim = DETECT_PER_BAND_AT_RES[INDEX_500M];
  pixel_dim = BAND_RATIO_AT_RES[INDEX_500M] * EV_1km_FRAMES;

  for (B = 0; B < L1B_BANDS_AT_RES[INDEX_250M]; B++)
  {
    returnStatus = Compute_Aggregates 
                     (2, line_dim, pixel_dim,
                      (VOIDP)L1B_Scan->SI.EV_250m_RefSB[B],
                      (VOIDP)L1B_Scan->UI.EV_250m_RefSB_UI[B],
                      (VOIDP)L1B_Scan->EV_250m_Aggr500m_RefSB[B],
                      (VOIDP)L1B_Scan->EV_250m_Aggr500m_RefSB_UI[B],
                      (VOIDP)L1B_Scan->EV_250m_Aggr500m_RefSB_SU[B]);
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR(returnStatus, 
                "Compute_Aggregates(250_to_500m) in "
                "Aggregate_L1B(), Granule.c");
  }


  /*---------------------------
    250_to_1km & 500_to_1km
  ---------------------------*/
  line_dim = DETECT_PER_BAND_AT_RES[INDEX_1000M_REFL];
  pixel_dim = EV_1km_FRAMES;
  
  if ((RFLAG & 1) == 0) {
  for (B = 0; B < L1B_BANDS_AT_RES[INDEX_250M]; B++)
  {
    returnStatus = Compute_Aggregates 
                    (4, line_dim, pixel_dim,
                     (VOIDP)L1B_Scan->SI.EV_250m_RefSB[B],
                     (VOIDP)L1B_Scan->UI.EV_250m_RefSB_UI[B],
                     (VOIDP)L1B_Scan->EV_250m_Aggr1km_RefSB[B],
                     (VOIDP)L1B_Scan->EV_250m_Aggr1km_RefSB_UI[B],
                     (VOIDP)L1B_Scan->EV_250m_Aggr1km_RefSB_SU[B]);
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR(returnStatus,
                "Compute_Aggregates(250_to_1km) in "
                "Aggregate_L1B(), Granule.c");
  }
  }

  if ((RFLAG & 2) == 0) {
  for (B = 0; B < L1B_BANDS_AT_RES[INDEX_500M]; B++)
  {
    returnStatus = Compute_Aggregates 
                    (2, line_dim, pixel_dim,
                     (VOIDP)L1B_Scan->SI.EV_500m_RefSB[B],
                     (VOIDP)L1B_Scan->UI.EV_500m_RefSB_UI[B],
                     (VOIDP)L1B_Scan->EV_500m_Aggr1km_RefSB[B],
                     (VOIDP)L1B_Scan->EV_500m_Aggr1km_RefSB_UI[B],
                     (VOIDP)L1B_Scan->EV_500m_Aggr1km_RefSB_SU[B]);
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR(returnStatus,
                "Compute_Aggregates(500_to_1km) in "
                "Aggregate_L1B(), Granule.c");
  }
  }

  return(MODIS_S_OK);
}


PGSt_SMF_status Close_L1A_Granule (L1A_granule_t *L1A_Gran,
                                   L1A_Scan_t    *L1A_Scan)
/*
!C**********************************************************************
!Description:
    This routine ends access to all L1A SDSs and ends the SD and Vdata
    interfaces to the L1A middle granule file.

!Input Parameters:
      L1A_granule_t *L1A_Gran
      L1A_Scan_t    *L1A_Scan

!Output Parameters:
      L1A_granule_t *L1A_Gran
      L1A_Scan_t    *L1A_Scan

!Revision History:
 (continue at top of the file)

 Revision 02.11 July 1998
 Removed Cleanup_L1A_Granule ( the arrays are changed to static allocated )
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)

 Revision 02.00 March 1997
 Modified to match changes in L1A_granule_t & L1A_Scan_t.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 (Neal's version between 01.01 & 02.00 ......???)

 Revision 01.01 1996/04/05
 Updated to match Version Design Document.
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Joan Baden (baden@highwire.gsfc.nasa.gov)

 Revision 01.00 1993
 Initial development
 Geir E. Kvaran(geir@highwire.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16     R   = 0;
 
  /*Close L1A EV SDS's*/
  for (R = 0; R < NUM_L1A_RESOLUTIONS; R++) {
    if ((RFLAG & (1 << R)) != 0) continue;
    if (SDendaccess(L1A_Scan->sds_id[R]) == FAIL)
      return(MODIS_F_NOK);
  }

  /*Clode L1A file*/
  if (SDend(L1A_Gran->sd_id) == FAIL) /* close the SD interface */
    return(MODIS_F_NOK);

  if (Vend(L1A_Gran->v_id) == FAIL) 
    return(MODIS_F_NOK);

  if (Hclose(L1A_Gran->v_id) == FAIL) /* close the Vdata interface */
    return(MODIS_F_NOK);

  return(returnStatus);
}

PGSt_SMF_status Close_L1B_Granule (L1B_granule_t         *L1B_Gran,
                                   L1B_Scan_t            *L1B_Scan,
                                   boolean               skip_night_hi_res)
/*
!C**********************************************************************
!Description:
    This routine ends SDS access to all open EV SDSs and ends swath interface
    to all L1B EV granule files.

!Input Parameters:
     L1B_granule_t         *L1B_Gran
     L1B_Scan_t            *L1B_Scan
     Run_Time_Parameters_t *runtime_params    Values read from PCF file.
     boolean               skip_night_hi_res  Logical; TRUE if and only 
                                              if all scans in the granule 
                                              are in NIGHT mode.                        

!Output Parameters:
     L1B_granule_t *L1B_Gran
     L1B_Scan_t    *L1B_Scan

!Revision History:
 (continue at top of the file)

 Revision 02.21  November 13, 2001 (Razor Issue #169)
  Added skip_night_hi_res to input variables. 
  Changed logic so that 250m and 500m data are not written when granule
  has no day mode scans and runtime parameter
  Write_Night_Mode_HiRes_Data is False.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 02.20 May   1999
 Added band 26 section.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.11 July  1998
 Removed Cleanup_L1B_Granule(), because all arrays are changed to 
    static allocated.
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)

 Revision 02.00 March 1997
 Modified to match changes in L1B_granule_t & L1B_Scan_t.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 (Neal's version between 01.01 & 02.00 ......???)

 Revision 01.01 1996/04/05
 Uppdated to match Version 1 Design Document.
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Joan Baden (baden@highwire.gsfc.nasa.gov)

 Revision 01.00 1993
 Initial development
 Geir E. Kvaran(geir@highwire.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   Code to make Vgroup has been removed.  It is not currently implemented 
   in SDS_Services, although it could be.
!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16     R                  = 0;
  int16     f                  = 0;
  int16     start_output_index = 0;
  
  if (skip_night_hi_res == True)
    start_output_index = INDEX_L1B_1km;
  else
    start_output_index = INDEX_L1B_250m;

  /*Close L1B EV SDS's*/
    
  for (R = start_output_index; R < NUM_L1A_RESOLUTIONS; R++)
  {
    if ((RFLAG & (1 << R)) != 0) continue;
    if (SDendaccess(L1B_Scan->SI_sds_id[R]) == FAIL) return(MODIS_F_NOK);
    if (SDendaccess(L1B_Scan->UI_sds_id[R]) == FAIL) return(MODIS_F_NOK);
  }
  if (skip_night_hi_res == False)
  {
    if (SDendaccess(L1B_Scan->EV_250m_Aggr500m_RefSB_sds_id) == FAIL) 
      return(MODIS_F_NOK);
    if (SDendaccess(L1B_Scan->EV_250m_Aggr500m_RefSB_UI_sds_id) == FAIL) 
      return(MODIS_F_NOK);
    if (SDendaccess(L1B_Scan->EV_250m_Aggr500m_RefSB_SU_sds_id) == FAIL) 
      return(MODIS_F_NOK);
  }

  if (SDendaccess(L1B_Scan->EV_250m_Aggr1km_RefSB_sds_id) == FAIL) 
    return(MODIS_F_NOK);
  if (SDendaccess(L1B_Scan->EV_250m_Aggr1km_RefSB_UI_sds_id) == FAIL) 
    return(MODIS_F_NOK);
  if (SDendaccess(L1B_Scan->EV_250m_Aggr1km_RefSB_SU_sds_id) == FAIL) 
    return(MODIS_F_NOK);
  if (SDendaccess(L1B_Scan->EV_500m_Aggr1km_RefSB_sds_id) == FAIL) 
    return(MODIS_F_NOK);
  if (SDendaccess(L1B_Scan->EV_500m_Aggr1km_RefSB_UI_sds_id) == FAIL) 
    return(MODIS_F_NOK);
  if (SDendaccess(L1B_Scan->EV_500m_Aggr1km_RefSB_SU_sds_id) == FAIL) 
    return(MODIS_F_NOK);

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
/*
 * Close SDS access to the band 26 SDSs.
 */
  if (SDendaccess(L1B_Scan->Band26.SI_sds_id) == FAIL)
    return(MODIS_F_NOK);
  if (SDendaccess(L1B_Scan->Band26.UI_sds_id) == FAIL)
    return(MODIS_F_NOK);
#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

  /*
   * If this granule has no day mode scans and the runtime parameter
   * Write_Night_Mode_HiRes_Data is False, the 250m and 500m
   * resolution files were not written, and so start with output 
   * file number 2, not 0.
   */
   
  if (skip_night_hi_res == True)
    start_output_index = INDEX_L1B_1km;
  else
    start_output_index = INDEX_L1B_250m;
  

  /*Close L1B EV files*/
  for (f = start_output_index; f < NUM_L1B_EV_FILES; f++)
    if (SWclose(L1B_Gran->sw_f_id[f]) == FAIL) return(MODIS_F_NOK);

  return(returnStatus);
}


PGSt_SMF_status  Compute_Aggregates (int16     scale,
                                     int16     line_dim_lower,
                                     int16     frame_dim_lower,
                                     uint16    *SI_in,
                                     uint8     *UI_in,
                                     uint16    *SI_out,
                                     uint8     *UI_out,
                                     int8      *SU_out)
/*
!C****************************************************************************
!Description:
   Perform spatial integration (aggregation) of higher resolution L1B data to
   appear as lower resolution data. This function computes the aggregated
   scaled integer, the aggregated uncertainty index and the number of samples
   used in aggregation.

!Input Parameters:
   int16  scale           spatial integration factor, in both the along track
                          and along scan directions of the data. (there are
                          limited valid values -- see the local variables)
   int16  line_dim_lower  line dimension (along track) of lower resolution data
   int16  frame_dim_lower frame dimension (along scan) of lower resolution data
   uint16 *SI_in          scaled integer data of higher resolution
                          (size = line_dim_lower * scale * frame_dim_lower * scale)
   uint8  *UI_in          uncertainty index data of higher resolution
                          (size = line_dim_lower * scale * frame_dim_lower * scale)

!Output Parameters:
   uint16 *SI_out         aggregated scaled integer data at lower resolution
                          (array size = line_dim_lower * frame_dim_lower).
   uint8  *UI_out         aggregated uncertainty index data at lower resolution
                          (array size = line_dim_lower * frame_dim_lower).
   int8   *SU_out         number of samples used in each aggregation
                          (array size = line_dim_lower * frame_dim_lower).

!Revision History:
 (continue at top of the file)

   Revision 01.01, January 9, 2001, Razor issue 149
   The assumption on how the first samples of different resolution data are
   registered was corrected. Also added design notes, made variable names more
   meaningful, improved efficiency slightly and added comments.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.00 April 1997
   Initial development
   (this is a translation of the FORTRAN code, 
   Aggregate_L1B, of Vicky Lin & Richard Hucek)
   Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   This code is translated from the FORTRAN code written by 
        Vicky Lin     (vlin@ltpmail.gsfc.nasa.gov)
        Richard Hucek (rhucek@ltpmail.gsfc.nasa.gov)
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   1. For the input and output arrays of data, each array represents a
      2D array, similar to:
         A[line_dim_lower][frame_dim_lower]               (lower resolution)
         A[line_dim_lower*scale][frame_dim_lower*scale]   (higher resolution)
      Thus, each higher resolution array has scale^2 times the values of
      a lower resolution array.
   2. In the line dimension (along track), each set of scale lines of the
      higher resolution are averaged (using uniform weighting) to form
      one aggregated line at the lower resolution.
   3. In the frame dimension (along scan), the effective response of a
      detector of either resolution is a triangle function, with the peak of
      the triangle at the effective location of the frame, and the base
      of the triangle having a total width of 2 times the frame width.  The
      triangular responses are overlapping by 1/2 the frame width
      (see the Geolocation ATBD for an explanation and diagram).
   4. The frames of different resolutions are registered such that the 1st
      frame of any line at one resolution is aligned with the first frame
      in a line of a different resolution.
   5. To aggregate in the frame dimension, (2*scale-1) consecutive values
      are chosen with triangular weighting to form the aggregated response.
      Because of the assumption in item 4 above, the first (scale-1) frames
      are missing from the first aggregate in any given line.

!END**************************************************************************
*/
{
  int16    lines_to_aggregate;   /* number of lines of higher resolution data
                                  * to aggregate to one line of lower resolution
                                  * data (same as input value of "scale").
                                  */
  int16    frames_to_aggregate;  /* number of frames (triangle functions) of
                                  * higher resolution data to aggregate to one
                                  * frame of lower resolution data (2*scale-1).
                                  */
  int16    line_offset;          /* essentially this is the index of the first
                                  * line of higher resolution data to use in
                                  * one aggregation across the set of lines.
                                  */
  int16    frame_offset;         /* essentially this is the index of the first
                                  * frame of higher resolution data to use in
                                  * one aggregation across the set of frames.
                                  */
  int16    frame_dim_higher;     /* number of higher resolution frames in each
                                  * line of higher resolution data (must be
                                  * frame_dim_lower * scale).
                                  */
  int16    L1;                   /* line index in lower resolution data */
  int16    F1;                   /* frame index in lower resolution data */
  int16    L2;                   /* line index in higher resolution data */
  int16    L;                    /* line index in one set of higher resolution
                                  * lines to aggregate.
                                  */
  int16    F;                    /* frame index in one set of higher resolution
                                  * frames to aggregate.
                                  */
  int16    F_start;              /* starting frame index within one set of
                                  * higher resolution frames to aggregate.
                                  * (not necessarily 0 because of design note 4).
                                  */
  int16    scale_index;          /* index in 1st dimension of W and SqW arrays */
  int32    index_higher;         /* array index within any higher resolution
                                  * array (must be int32).
                                  */
  int32    index_lower;          /* array index within any lower resolution
                                  * array (must be int32).
                                  */
  int8     samples;              /* used to count number of valid samples used
                                  * in one aggregation.
                                  */
  float32  sum_W;                /* accumulates the weight factors */
  float32  sum_WSI;              /* accumulates the weighted scaled integer */
  float32  sum_SqW;              /* accumulates the square of the weight factors */
  float32  sum_SqWUI;            /* accumulates the square of the weighted UI */

  /*
   * The following define how many valid scale factors are handled by this
   * function and the maximum number of frames for triangle weighting.
   */

#define NUM_SCALE_FACTORS  2
#define MAX_AGGR_FRAMES    7

  /*
   * The following arrays hold the valid scale factors and the associated
   * weighting factors, for aggregating along the frame (scan) dimension.
   *   W   = weight
   *   SqW = square of weight.
   * The number of frames to aggregate is (2 * scale - 1).
   */

  int16 valid_scale_values[NUM_SCALE_FACTORS] = {2, 4};

  float32  W[NUM_SCALE_FACTORS][MAX_AGGR_FRAMES] = {
    {1,  2,  1,  0,  0,  0,  0}, 
    {1,  2,  3,  4,  3,  2,  1}
  };

  float32  SqW[NUM_SCALE_FACTORS][MAX_AGGR_FRAMES] = {
    {1,  4,  1,  0,  0,  0,  0}, 
    {1,  4,  9, 16,  9,  4,  1}
  };

  /*************************************************************************/

  /*
   * Determine the "scale_index", the index in the W and SqW arrays.
   * This checks the value of "scale" in the process.
   */

  for (scale_index = 0; scale_index < NUM_SCALE_FACTORS; scale_index++)
  {
    if (scale == valid_scale_values[scale_index])
      break;
  }
  if (scale_index == NUM_SCALE_FACTORS)
  {
    char errmsg[512];
    sprintf (errmsg, "Input argument (scale = %d) is invalid.", scale);
    L1BErrorMsg("Compute_Aggregates", MODIS_F_INVALID_ARGUMENT, errmsg, NULL, 0,
                NULL, True);
  }


  lines_to_aggregate  = scale;
  frames_to_aggregate = 2 * scale - 1;
  frame_dim_higher    = frame_dim_lower * scale; 

  line_offset = 0;
  for (L1 = 0; L1 < line_dim_lower; L1++, line_offset += scale)
  {

    /*
     * frame_offset starts negative.  For the first aggregation of a line,
     * we will start at the peak of the triangle weight to compensate.
     */

    frame_offset = - (scale - 1);
    for (F1 = 0; F1 < frame_dim_lower; F1++, frame_offset += scale)
    {

      /*
       * Because the first frames of different resolutions align, for the
       * first aggregate of a line we are missing some frames on the left side
       * of the triangle weighting function.  Set F_start so that we will
       * have a valid frame index later.
       */

      if (F1 == 0)
        F_start = scale - 1;    /* start at the peak of triangle weighting */
      else
        F_start = 0;            /* start at beginning of triangle weighting */

      samples   = 0;
      sum_W     = 0;
      sum_WSI   = 0;
      sum_SqW   = 0;
      sum_SqWUI = 0;

      for (L = 0; L < lines_to_aggregate; L++)
      {
        L2 = line_offset + L;

        index_higher = L2 * frame_dim_higher + frame_offset + F_start;
        for (F = F_start; F < frames_to_aggregate; F++, index_higher++)
        {
          if (SI_in[index_higher] <= DN15_SAT)
          {
            samples++;
            sum_W     += W[scale_index][F];
            sum_WSI   += W[scale_index][F] * SI_in[index_higher];
            sum_SqW   += SqW[scale_index][F];
            sum_SqWUI += SqW[scale_index][F] * UI_in[index_higher] * UI_in[index_higher];
          }
        }
      }

      index_lower = L1 * frame_dim_lower + F1; 
      if (samples > 0)
      {
        SI_out[index_lower] = (uint16)(sum_WSI / sum_W + 0.5);
        UI_out[index_lower] = (uint8)(sqrt((double)(sum_SqWUI / sum_SqW)) + 0.5);
        SU_out[index_lower] = samples;
      }
      else
      {
        SI_out[index_lower] = AGGREGATION_FAIL_SI;
        UI_out[index_lower] = BAD_DATA_UI;
        SU_out[index_lower] = 0;
      }
    }
  }

  return(MODIS_S_OK);
}

PGSt_SMF_status Read_L1A_EV_Scan (int16            S,
                                  L1A_granule_t    *L1A_Gran,
                                  L1A_Scan_t       *L1A_Scan)
/*
!C**************************************************************************
!Description:
   For the input scan index, S, this routine reads earth-view (EV) scaled
   integer data for each of the four L1A EV SDSs: EV_250m, EV_500m,
   EV_1km_day, and EV_1km_night.  If the scan is not a "Day" mode scan,
   then only the data from the EV_1km_night SDS is read.

!Input Parameters:
   int16            S           current scan index
   L1A_granule_t    *L1A_Gran   contains day-night flag, number of scans
                                and valid SDS ids for each L1A EV SDS.

!Output Parameters:
   L1A_Scan_t       *L1A_Scan   filled SDS data for all resolutions of 1 scan

!Revision History:
 (continue at top of the file)

   Revision 02.1.2, January 22, 2001
   Corrected variable names and simplified the logic within this function
   to remove unnecessary lines. Added design notes.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.1.1 July 1997
   Modified to match change in L1A_Scan_t;
   Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)

   Revision 02.00 March 1997
   Modified to match change in L1A_Scan_t;
   limited processing to EV data.
   Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

   (Neal's version between 01.01 & 02.00......???)

   Revision 01.01 1996/04/05
   Update to match Version 1 Design Document
   John Hannon(hannon@highwire.gsfc.nasa.gov)
   Joan Baden (baden@highwire.gsfc.nasa.gov)

   Revision 01.00 1993
   Initial development
   Geir E. Kvaran(geir@highwire.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   1.  We only need to read the 250m, 500m and 1km_day SDS data if the scan
       is a "Day" scan, according to the day-night flag.  We always read the
       1km_night data.
   2.  The L1A SDSs are oversized in terms of the number of frames.  We
       read the full set of data for 1 scan to make the internal read
       more efficient (each call will read one stream of data from the file
       rather than having to read many disjointed streams of data).
!END************************************************************************
*/
{
  int16  R;                     /* L1A resolution index */
  int16  D;                     /* Detector index */
  int16  B;                     /* Band index */
  int32  start[3] = {0, 0, 0};  /* starting indices, each dimension */
                                /* (must be initialized to zero) */
  int32  edge[3];               /* # of elements to read, each dimension */
  intn   hdf_return;
  char   *location = "Read_L1A_EV_Scan";    /* This function name */
  int16  *ev_data_p[NUM_L1A_RESOLUTIONS];   /* pointers to each SDS data */
  int32 offset;

  /*
   * Check gross validity of inputs.
   */

  if (S < 0 || S >= L1A_Gran->num_scans || !L1A_Scan || !L1A_Gran) {
    L1BErrorMsg(location, MODIS_F_INVALID_ARGUMENT,
                "Scan index (S), L1A_Gran or L1A_Scan value is invalid.",
                NULL, FIRST_L1A_GRANULE, NULL, True);
    return MODIS_F_INVALID_ARGUMENT;
  }

  /*
   * Assign data pointers.
   */

  ev_data_p[INDEX_250M]        = (int16 *) L1A_Scan->EV_250m;
  ev_data_p[INDEX_500M]        = (int16 *) L1A_Scan->EV_500m;
  ev_data_p[INDEX_1000M_DAY]   = (int16 *) L1A_Scan->EV_1km_day;
  ev_data_p[INDEX_1000M_NIGHT] = (int16 *) L1A_Scan->EV_1km_night;

  /*
   * Loop through L1A resolutions and read one scan of L1A data.
   * (see design notes regarding which SDSs to read).
   */

  for (R = 0; R < NUM_L1A_RESOLUTIONS; R++)
  {
    if ((RFLAG & (1 << R)) != 0) continue;

    if (strcmp(L1A_Gran->ScanType[S],"Day") == SAME ||
        R == INDEX_1000M_EMISS)
    {
      if (L1A_Gran->Extract_Pixel_Offset == -1 &&
          L1A_Gran->Extract_Pixel_Count == -1) {

      start[0] = S * DETECT_PER_BAND_AT_RES[R];
      edge[0]  = DETECT_PER_BAND_AT_RES[R];
      edge[1]  = L1A_BANDS_AT_RES[R];
      edge[2]  = EV_1km_FRAMES * BAND_RATIO_AT_RES[R];

      hdf_return = SDreaddata (L1A_Scan->sds_id[R], start, NULL, edge,
                               ev_data_p[R]);
      if (hdf_return == FAIL) {
        char errmsg[256];
        char *resnames[NUM_L1A_RESOLUTIONS] = {
          "250m", "500m", "1km_day", "1km_night"
        };
        sprintf (errmsg, "Could not read L1A SDS: EV_%s.", resnames[R]);
        L1BErrorMsg(location, MODIS_F_READ_ERROR, errmsg,
                    "SDreaddata", FIRST_L1A_GRANULE, NULL, True);
        return MODIS_F_READ_ERROR;
      }
      } else {

        /* Clear buffer */
        memset(ev_data_p[R], 0, DETECT_PER_BAND_AT_RES[R] * L1A_BANDS_AT_RES[R] *
               EV_1km_FRAMES * BAND_RATIO_AT_RES[R]);

        for (D = 0; D < DETECT_PER_BAND_AT_RES[R]; D++)

          for (B = 0; B < L1A_BANDS_AT_RES[R]; B++) {
            start[0] = S*DETECT_PER_BAND_AT_RES[R] + D;
            start[1] = B;
            start[2] = 0;
            edge[0]  = 1;
            edge[1]  = 1;
            edge[2]  = L1A_Gran->Extract_Pixel_Count * BAND_RATIO_AT_RES[R];

            offset = (D*L1A_BANDS_AT_RES[R]*EV_1km_FRAMES +
                      B*EV_1km_FRAMES + L1A_Gran->Extract_Pixel_Offset) *
                      BAND_RATIO_AT_RES[R];

            hdf_return = SDreaddata (L1A_Scan->sds_id[R], start, NULL, edge,
                                     &ev_data_p[R][offset]);
          }
      }
    }
  }

  return MODIS_S_OK;
}


#include <string.h>  /* for the strcat function */

int safe_strcat(char *buf, char *str, int buflen)
/*
!C*****************************************************************************
!Description:
   Safely concatenate a string onto a buffer.  If neccessary, truncate the
   string to fit onto the end of the buffer, allowing a null to be at the end.
   Return values:
     0: the string was completely and successfully concatenated onto the end
        of the buffer.
     1: the buffer did not exist, the string did not exist or the buffer was
        to small to hold the entire string -- in which case as many characters
        as possible from the string are concatenated onto the end of the buffer.

!Input Parameters:
   char *buf       Buffer holding characters to which we are concatenating
                   a string onto the end of.
   char *str       The string to be placed at the end of the buffer.
   int  buflen     The memory length of the buffer (not the string length).

!Output Parameters:
   (none)

!Revision History:
 (continue at top of the file)

   Revision 1.0.1  January 17, 2002   Razor Issue #172
   Improve portability of code to 64-bit mode.
   Change variables j, nb, and ns to type size_t.
   Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)
   
   Revision 1.0 Nov. 1, 1999
   Initial development.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:


!END******************************************************************************
*/
{
  int i;           /* indexes */
  size_t j;
  size_t nb;       /* string length of buf */
  size_t ns;       /* string length of str */
  int mx;          /* maximum number of non-null chars allowed in buffer */

  if (!buf || !str || buflen <= 0) return 1;

    /* Compute total length of required memory. */

  nb = strlen(buf);
  ns = strlen(str);
  if ((nb + ns) < buflen) {
    strcat(buf, str);
    return 0;
  }
  else {
    mx = buflen - 1;
    for (i = 0, j = nb; j < mx; i++, j++) buf[j] = str[i];
    buf[mx] = '\0';
    return 1;
  }
}


/*
 * The following are global values to use for the "other_msg" input to
 * L1BErrorMsg. These errors may occur many times.
 */

char Invalid_MOD01_Msg[] = "\
The MOD01 granule appears to be invalid.  Data values are present\n\
that are not allowed by the MOD01 file specifications.";


void L1BErrorMsg(
  char *L1B_location,    /* name of L1B function that error occurred in */
  PGSt_SMF_code code,    /* associated error code for this error */
  char *input_message,   /* short (usually 1-line) description of error */
  char *assoc_function,  /* name of an associated function that failed */
  int32 lun,             /* associated LUN for the file being accessed */
  char *other_msg,       /* other message to add such as probable cause */
  boolean error_out      /* flag to tell whether or not to call SMF_ERROR */
)
/*
!C**********************************************************************
!Description:
   Format an error message and write it to the LogReport and LogStatus
   files. If the flag "error_out" is True, then also call SMF_ERROR.
   Otherwise, this function will simply return to continue processing.
   The general form of the message created is:

   ERROR in [L1B_location]():
   [input_message]
   Call to [assoc_function]() failed.
   File LUN: [lun] ([type of file for this LUN])
   [other_msg]

   If an input is NULL (for char *) or zero (for int32), then that part of
   the message is not included.

!Input Parameters:
   char *L1B_location    name of L1B function that error occurred in
   PGSt_SMF_code code    associated error code for this error
   char *input_message   short (usually 1-line) description of error
   char *assoc_function  name of an associated function that failed
   int32 LUN             associated LUN for the file being accessed
   char *other_msg       other message to add such as probable cause
   boolean error_out     flag to tell whether or not to call SMF_ERROR 

!Output Parameters:
   (none)

!Revision History:
 (continue at top of the file)

   Revision 1.02  October 16, 2004  Razor Issue #200
   Casted Int32 variables in sprintf calls to "long" with the
   format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 1.01  August 15, 2002
   Add messages regarding "ProcessingEnvironment", "MCST LUT Version", 
   and "Write_Night_Mode_HiRes_Data".
   Alice Isaacman, SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 1.0 Nov. 4, 1999
   Initial development.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   The messages going to the LogStatus and LogReport are formatted to 
   look similar to each other.

   The message size should be smaller than PGS_SMF_MAX_MSGBUF_SIZE
   due to internal limits in PGS_SMF_SetDynamicMsg.  This routine will
   fill up the message buffer until there is no more room and then
   call PGS_SMF_SetDynamicMsg.  Thus, a message may appear truncated.

   The same size message buffer is assumed for the LogReport message.
   This could result in the LogReport message being truncated from the
   message written to the log status file (since the L1B_location is added
   at the beginning of the LogReport file).

   This restriction may not be necessary.  The source code for
   PGS_SMF_GenerateStatusReport does not appear to have any string
   length restrictions.

   The lun (if one of the list of L1B LUNs) is used to determine the
   the category of the file (e.g., QA_Lookup_Tables_file).

   All "Concatenate" operations will use the "safe_strcat" function that
   ensures that the string can be concatenated onto the end and that flags
   when a message is truncated.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus;           /* Return value for some PGS functs */
  char msg[PGS_SMF_MAX_MSGBUF_SIZE];      /* Buffer for LogStatus message */
  char report[PGS_SMF_MAX_MSGBUF_SIZE];   /* Buffer for LogReport message */
  char lunstr[PGS_SMF_MAX_MSGBUF_SIZE];   /* Buffer for LUN string */
  char category[PGS_SMF_MAX_MSGBUF_SIZE]; /* Buffer for file category */
  PGSt_SMF_boolean internal_failure;      /* Flags internal function failure. */

  /*
   * First, format the message to be written to the LogStatus file.
   * The L1B function name is not needed for that message (it is an argument
   * to the PGS_SMF_SetDynamicMsg function).
   */

  /* initialize the message buffer */

  strcpy(msg, "");

    /* Add the specific input message */

  if (input_message) {
    safe_strcat(msg, input_message, PGS_SMF_MAX_MSGBUF_SIZE);
    if (msg[strlen(msg)-1] != '\n')
      safe_strcat(msg, "\n", PGS_SMF_MAX_MSGBUF_SIZE);
  }

    /* Add the associated function name */

  if (assoc_function) {
    safe_strcat(msg, "Call to ", PGS_SMF_MAX_MSGBUF_SIZE);
    safe_strcat(msg, assoc_function, PGS_SMF_MAX_MSGBUF_SIZE);
    safe_strcat(msg, "() failed.\n", PGS_SMF_MAX_MSGBUF_SIZE);
  }

    /* Add the LUN information */

  if (lun != 0) {
    sprintf(lunstr, "File LUN:  %ld\nFile type: ", (long)lun);
    safe_strcat(msg, lunstr, PGS_SMF_MAX_MSGBUF_SIZE);
    switch (lun) {
      case MCF_FILE_QKM:           
        strcpy(category, "MOD02QKM MCF\n");                         break;
      case MCF_FILE_HKM:           
        strcpy(category, "MOD02HKM MCF\n");                         break;
      case MCF_FILE_1KM:           
        strcpy(category, "MOD021KM MCF\n");                         break;
      case MCF_FILE_OBC:           
        strcpy(category, "MOD02OBC MCF\n");                         break;
      case L1B_EV_250M_FILE:       
        strcpy(category, "MOD02QKM output granule\n");              break;
      case L1B_EV_500M_FILE:       
        strcpy(category, "MOD02HKM output granule\n");              break;
      case L1B_EV_1000M_FILE:      
        strcpy(category, "MOD021KM output granule\n");              break;
      case L1B_OBC_FILE:           
        strcpy(category, "MOD02OBC output granule\n");              break;
      case REFLECTIVE_TABLES_FILE: 
        strcpy(category, "Reflective_Lookup_Tables_file\n");        break;
      case EMISSIVE_TABLES_FILE:   
        strcpy(category, "Emissive_Lookup_Tables_file\n");          break;
      case QA_TABLES_FILE:         
        strcpy(category, "QA_Lookup_Tables_file\n");                break;
      case LEADING_L1A_GRANULE:    
        strcpy(category, "previous MOD01 input granule\n");         break;
      case FIRST_L1A_GRANULE:      
        strcpy(category, "current MOD01 input granule\n");          break;
      case TRAILING_L1A_GRANULE:   
        strcpy(category, "following MOD01 input granule\n");        break;
      case GEOLOCATION_FILE:       
        strcpy(category, "MOD03 input granule\n");                  break;
      case PROCESSING_ENVIRONMENT_LUN: strcpy(category, 
          "PCF File; ProcessingEnvironment; not required.\n");      break;      
      case MCST_LUT_VERSION_LUN: strcpy(category, 
          "PCF File; invalid MCSTLUTVersion value\n");              break;
      case WRITE_NIGHT_HIRES_LUN: strcpy(category, 
          "PCF File; invalid Write_Night_Mode_HiRes_Data value\n"); break;
      default:                     
        strcpy(category, "invalid value, L1B defect\n");            break;
    }
    safe_strcat(msg, category, PGS_SMF_MAX_MSGBUF_SIZE);
    
  }

  /*
   * Add the other message
   */

  if (other_msg) {
    safe_strcat(msg, other_msg, PGS_SMF_MAX_MSGBUF_SIZE);
  }

  /*
   * Write the message to the LogStatus file.
   */

  internal_failure = PGS_FALSE;
  returnStatus = PGS_SMF_SetDynamicMsg(code, msg, L1B_location);
  if (returnStatus != PGS_S_SUCCESS)
    internal_failure = PGS_TRUE;

  /*
   * Now format the same message to be written to the LogReport file.
   * Make the appearance of the message similar in the LogReport file.
   */

  strcpy(report, "");

    /* Add the L1B location */

  if (L1B_location) 
  {
    safe_strcat(report, L1B_location, PGS_SMF_MAX_MSGBUF_SIZE);
    safe_strcat(report, "():", PGS_SMF_MAX_MSGBUF_SIZE);
  }
  else 
  {
    safe_strcat(report, ":", PGS_SMF_MAX_MSGBUF_SIZE);
  }

  /* Add the textual and numerical equivalents of the exit code. */

  switch (code) {
    case MODIS_F_OUT_OF_MEMORY:
      safe_strcat(report, "MODIS_F_OUT_OF_MEMORY:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_MEM_FREE_FAIL:
      safe_strcat(report, "MODIS_F_MEM_FREE_FAIL:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_FILE_NOT_FOUND:
      safe_strcat(report, "MODIS_F_FILE_NOT_FOUND:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_READ_ERROR:
      safe_strcat(report, "MODIS_F_READ_ERROR:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_WRITE_ERROR:
      safe_strcat(report, "MODIS_F_WRITE_ERROR:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_OUT_OF_RANGE:
      safe_strcat(report, "MODIS_F_OUT_OF_RANGE:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_W_OUT_OF_RANGE:
      safe_strcat(report, "MODIS_W_OUT_OF_RANGE:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_W_TIME_INCORRECT:
      safe_strcat(report, "MODIS_W_TIME_INCORRECT:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_NOK:
      safe_strcat(report, "MODIS_F_NOK:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_HDF_ERROR:
      safe_strcat(report, "MODIS_F_HDF_ERROR:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_NO_MORE:
      safe_strcat(report, "MODIS_F_NO_MORE:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_S_NO_MORE:
      safe_strcat(report, "MODIS_S_NO_MORE:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_FILE_NOT_OPENED:
      safe_strcat(report, "MODIS_F_FILE_NOT_OPENED:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_FILE_NOT_CREATED:
      safe_strcat(report, "MODIS_F_FILE_NOT_CREATED:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_F_INVALID_ARGUMENT:
      safe_strcat(report, "MODIS_F_INVALID_ARGUMENT:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_E_TESTING:
      safe_strcat(report, "MODIS_E_TESTING:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    case MODIS_S_OK:
      safe_strcat(report, "MODIS_S_OK:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
    default:
      safe_strcat(report, "UNKNOWN_CODE:", 
          PGS_SMF_MAX_MSGBUF_SIZE); break;
  }
  sprintf(lunstr, "%d\n", code);
  safe_strcat(report, lunstr, PGS_SMF_MAX_MSGBUF_SIZE);

  /*
   * Add the previously formatted message written to LogStatus file.
   * Write the total message to the LogReport file.
   */

  safe_strcat(report, msg, PGS_SMF_MAX_MSGBUF_SIZE);
  returnStatus = PGS_SMF_GenerateStatusReport(report);
  if (returnStatus != PGS_S_SUCCESS)
    internal_failure = PGS_TRUE;

  /*
   * If there is an internal failure, write fatal error message and
   * error out.  (The way that SMF_ERROR is designed, the error message
   * will probably not get written anyway.)
   */

  if (internal_failure == PGS_TRUE) 
  {
    SMF_ERROR(MODIS_F_NOK,
      "Fatal error in HDFFuncErrorMsg:  Failed to write to\n"
      "PGS_SMF_SetDynamicMsg or to PGS_SMF_GenerateStatusReport");
  }

  /*
   * Call SMF_ERROR to error out if the flag is set to True.
   */

  if (error_out == True)
    SMF_ERROR(code, NULL);
}


/*
 * The following variable contains error message information used in
 * SMF_ERROR.  The first entry in the variable is the
 * default (in case of input error) and the last entry in the table
 * must be NULL.
 */

  static struct 
  {
    char *name;         /* macro name of the error code */
    int32 value;        /* numerical value of the code */
    char *actions;      /* operator actions to take */
  } err_msg_info[] = 
  {

    { "MODIS_F_NOK",              MODIS_F_NOK, "\
  Operator Actions:\n\
  Contact MCST.\n"
    },

    { "MODIS_F_OUT_OF_MEMORY",    MODIS_F_OUT_OF_MEMORY, "\
  Operator Actions:\n\
  Check the system health to see if it is overloaded.\n\
  Try re-running the PGE with system loading reduced.\n\
  Contact MCST.\n"
    },

    { "MODIS_F_MEM_FREE_FAIL",    MODIS_F_MEM_FREE_FAIL, "\
  Operator Actions:\n\
  Contact MCST.\n"
    },

    { "MODIS_F_FILE_NOT_FOUND",   MODIS_F_FILE_NOT_FOUND, "\
  Operator Actions:\n\
  Check previous messages to identify LUN of the file and other possible causes.\n\
  Check the PCF file to ensure that UR exists and PCF format is correct.\n\
  Check the PCF file to ensure that file names and directories are correct.\n\
  Contact MCST.\n"
    },

    { "MODIS_F_READ_ERROR",       MODIS_F_READ_ERROR, "\
  Operator Actions:\n\
  Check previous messages to identify LUN of the file and other possible causes.\n\
  Check the PCF file to ensure that file names and directories are correct.\n\
  Check the file type and size to see if it is consistent with file category.\n\
  Contact MCST.\n"
    },

    { "MODIS_F_WRITE_ERROR",      MODIS_F_WRITE_ERROR, "\
  Operator Actions:\n\
  Check previous messages to identify LUN of the file and other possible causes.\n\
  Check the PCF file to ensure that file names and directories are correct.\n\
  Check available system disk space.\n\
  Check to see if file has been moved or corrupted during execution.\n\
  Contact MCST.\n"
    },

    { "MODIS_F_OUT_OF_RANGE",     MODIS_F_OUT_OF_RANGE, "\
  Operator Actions:\n\
  Check previous messages to identify LUN of the file and other possible causes.\n\
  Contact MCST.\n\
  Also contact SDST if the LUN corresponds to a MOD01 or MOD03 granule.\n"
    },

    { "MODIS_W_OUT_OF_RANGE",     MODIS_W_OUT_OF_RANGE, "\
  Operator Actions:\n\
  Contact MCST.\n"
    },

    { "MODIS_W_TIME_INCORRECT",   MODIS_W_TIME_INCORRECT, "\
  Operator Actions:\n\
  Check Leap Seconds and other toolkit files and environment variables.\n\
  Notify L1B product users that time in product may be incorrect.\n\
  Contact MCST.\n"
    },

    { "MODIS_F_HDF_ERROR",        MODIS_F_HDF_ERROR, "\
  Operator Actions:\n\
  Check previous messages to identify LUN of the file and other possible causes.\n\
  Contact MCST.\n"
    },

    { "MODIS_F_NO_MORE",          MODIS_F_NO_MORE, "\
  Operator Actions:\n\
  Contact MCST.\n"
    },

    { "MODIS_S_NO_MORE",          MODIS_S_NO_MORE, NULL },

    { "MODIS_F_FILE_NOT_OPENED",  MODIS_F_FILE_NOT_OPENED, "\
  Operator Actions:\n\
  Check previous messages to identify LUN of the file and other possible causes.\n\
  Check the PCF file to ensure that file names and directories are correct.\n\
  Contact MCST.\n"
    },

    { "MODIS_F_FILE_NOT_CREATED", MODIS_F_FILE_NOT_CREATED, "\
  Operator Actions:\n\
  Check previous messages to identify LUN of the file and other possible causes.\n\
  Check the PCF file to ensure that file names and directories are correct.\n\
  Check file permissions of the system.\n\
  Contact MCST.\n"
    },

    { "MODIS_F_INVALID_ARGUMENT", MODIS_F_INVALID_ARGUMENT, "\
  Operator Actions:\n\
  Contact MCST.\n"
    },

    { "MODIS_E_TESTING",          MODIS_E_TESTING, "\
  Operator Actions:\n\
  Contact MCST.\n"
    },

    { "MODIS_S_OK",               MODIS_S_OK, NULL },

    { NULL, 0, NULL }       /*** Must be NULL terminated ***/
  };


void SMF_ERROR (PGSt_SMF_code code, char *messagestring)
/*
!C**************************************************************************
!Description:
   Provides the interface to the Status Message Facility in the SDP Toolkit.
   A message is formed using the format string from the L1B "seed" file,
   the input message string and the local time.  The message is written to
   both the LogStatus file and the LogReport file.  Additionally, operator
   actions are written as the last messages to the files.

   IMPORTANT: This function will invoke "exit" for fatal error codes or
   internal failures that should not occur (such as the environment or
   seed messages do not function).  See other comments in design notes.

!Input Parameters:
   PGSt_SMF_code code          One of the L1B message codes defined in
                               PGS_Error_Codes.h" and in the seed file.
   char *       messagestring  A string which contains a specific error
                               message.  The function name should be already
                               part of the string, if appropriate.

!Output Parameters:
   (none)

!Revision History:
 (continue at top of the file)

   Revision 2.0.1  January 17, 2002   Razor Issue #172
   Improve portability of code to 64-bit mode.
   Change variable n to type size_t.  Use n instead of i in strlen operations.
   Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.00 Oct. 27, 1999
   Added checks on input validity (NULL, string length).  Reformatted
   message so that it does not appear all on one line.  Added "MOD_PR02"
   to the beginning to distinguish the source from other processes that
   may write to the same Log files.  Added additional error messages and
   warnings.  Function was completely restructured. Added design notes.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.00 199604/05
   Initial development, new function for Version 1.
   John Hannon(hannon@highwire.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   This routine is based on the CHECK_ERR function used by Geir Kvaran
   in the Beta version of the code.
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

 The final form of the message consists of three parts, concatenated together:
 1. leading part: string from the seed file associated with this code
 2. middle part:  the input "messagestring"
 3. last part:    the local time that the error occurred
 Note that the string from the seed file actually is the format string for the
 sprintf function (it has two "%s" fields for the middle and trailing parts).

 If the total message is too large to fit into the maximum size allowed by the
 Toolkit (PGS_SMF_MAX_MSGBUF_SIZE - 1), then the user supplied message is
 truncated so that it will fit.  In this case, the process will exit only if the
 original message code was fatal.

 Side effect:  The input messagestring may be modified by this function if the
 string length is too long.  Thus, the calling function should not plan on using
 the string afterward.

 For a NULL or zero-length messagestring, assume that the seed file message is
 the complete message.

 Note: The following conditions could cause an exit invocation:
 1.  The input "code" corresponds to a fatal error, or
 2.  There is an internal failure inside this function.
 An internal failure is a failure of one of the PGS toolkit functions.  Other
 errors that occur may generate warnings, but will not halt the process
 in and of themselves. 

 Depending on the internal failure, there may be no error messages written.
 If absolutely no message of any kind comes out, perhaps the LogStatus and
 LogReport files are not set properly in the PCF.

!END***************************************************************************
*/
{
  PGSt_SMF_status returnStatus;          /* Return value for some PGS functs */
  char msgfmt[PGS_SMF_MAX_MSGBUF_SIZE];  /* Holds the format retrieved from
                                          * the seed file (format associated with
                                          * the particular input error code).
                                          */
#define ASCTIMEBUF 26
  char timestr[ASCTIMEBUF];              /* Holds the local time in a string.
                                          * The time function "asctime_r" sets a
                                          * null-terminated string of length 25.
                                          */
  char msgstr[PGS_SMF_MAX_MSGBUF_SIZE];  /* Holds the modified input message (could
                                          * be truncated or have newlines added)
                                          */
  char buf[PGS_SMF_MAX_MSGBUF_SIZE];     /* Holds the final concatenated message. */
  struct tm *tptr;                       /* used in getting time string */
  time_t local;                          /* used in getting time string */
  PGSt_SMF_boolean internal_failure;     /* Flags any internal function failure. */
  size_t n;                               /* temporary variable */

  char *defaultmsgfmt = "Error running...%sTIME:%s";  /* backup value */
  PGSt_SMF_code internal_failure_code = MODIS_F_NOK;  /* value to use for internal
                                                       * failure.
                                                       */
  char *s;                               /* temporary variable */
  int32 i;                               /* loop index */
  int32 imsg;                            /* index in "err_mgs_info" for this code */

/*************************** MACRO WRITE_MSG(c,m) *********************************
 Writes a string to both the LogStatus and LogReport.  For writing to the
 LogStatus, the function name is not input (set to NULL).  Note that if either
 call fails, then we probably will not get any messages at all.  There is one
 final check on the string length, although this should be OK.
    c = SMF message code
    m = message string to be written
**********************************************************************************/
#define WRITE_MSG(c,m)                               \
  if (strlen(m) > (PGS_SMF_MAX_MSGBUF_SIZE - 1)) {   \
    m[strlen(m)-1] = '\0';                           \
  }                                                  \
  returnStatus = PGS_SMF_SetDynamicMsg(c, m, NULL);  \
  if (returnStatus != PGS_S_SUCCESS)                 \
    internal_failure = PGS_TRUE;                     \
  returnStatus = PGS_SMF_GenerateStatusReport(m);    \
  if (returnStatus != PGS_S_SUCCESS)                 \
    internal_failure = PGS_TRUE;

  /******************* Begin code ***************/

  /*
   * Initialize flag to false.
   */

  internal_failure = PGS_FALSE;

  /*
   * Get the seed file message format string for this error code.
   */

  returnStatus = PGS_SMF_GetMsgByCode(code, msgfmt);
  if (returnStatus != PGS_S_SUCCESS)
  {
    char errmsg[] = "MOD_PR02 Warning in SMF_ERROR: Call to PGS_SMF_GetMsgByCode failed.\n"
                    "Seed file may not be set up properly or code is invalid.\n"
                    "The input message will be printed below and process will exit.\n";
    WRITE_MSG(internal_failure_code,errmsg);
    internal_failure = PGS_TRUE;
    strcpy(msgfmt, defaultmsgfmt);
  }

  /*
   * Form the local time and copy the time string into the timestr.
   * If the time functions fail, copy a small error message into the timestr.
   */

  local = time(NULL);
  tptr = localtime(&local);
  s = asctime(tptr);
  if (s)
  {
    strcpy(timestr, s);
  }
  else
  {
    strcpy(timestr,"WARNING: asctime failed");  /* this string must be < 26 chars! */
  }
  if (strlen(timestr) > (ASCTIMEBUF-1))
  {
    strcpy(timestr,"WARNING: asctime failed");  /* this string must be < 26 chars! */
  }

  /*
   * Put the input message into the msgstr.  Truncate if necessary.
   */

  if (!messagestring)
  {
    strcpy(msgstr, "MOD_PR02. ");
  }
  else if (!strlen(messagestring))
  {
    strcpy(msgstr, "MOD_PR02. ");
  }
  else
  {
    strcpy(msgstr, "MOD_PR02:");
    if (messagestring[0] != '\n')
      strcat(msgstr, "\n");
    n = PGS_SMF_MAX_MSGBUF_SIZE - 1 - strlen(msgfmt) - 
        strlen(timestr) - strlen(msgstr);
    if (strlen(messagestring) > n)
      messagestring[n] = '\0';
    strcat(msgstr, messagestring);
    n = strlen(msgstr);
    if (msgstr[n-1] != '\n')
      strcat(msgstr, "\n");
  }

  /*
   * Form the final message string.
   */

  sprintf(buf, msgfmt, msgstr, timestr);

  /*
   * Determine index of operator actions, based on exit code.
   */

  i = 0; imsg = 0;
  while (err_msg_info[i].name)
  {
    if (code == err_msg_info[i].value)
    {
      imsg = i;
      break;
    }
    i++;
  }

  /*
   * Determine string length of operator actions.
   */

  if (err_msg_info[imsg].actions)
    n = strlen(err_msg_info[imsg].actions);
  else
    n = 0;

  /*
   * If there is room, append the operator actions to the end of the message
   * buffer to be written.  Otherwise, write the message buffer and the write
   * the operator actions separately.
   */

  if (n == 0)
  {
    WRITE_MSG(code,buf);
  }
  else if ((strlen(buf)+n) < PGS_SMF_MAX_MSGBUF_SIZE)
  {
    strcat(buf, err_msg_info[imsg].actions);
    WRITE_MSG(code,buf);
  }
  else
  {
    WRITE_MSG(code,buf);
    WRITE_MSG(code,err_msg_info[imsg].actions);
  }

  /*
   * Test to see if the code corresponds to a fatal error or if there was any
   * internal failure and then exit.
   */

  if ((PGS_SMF_TestFatalLevel(code)) == PGS_TRUE || internal_failure == PGS_TRUE)
  {

    exit (MOD_PR02_Failure_Exit_Code);

  }
}

void Bad_L1A_Error_Out(char *name, char *message)
/*
!C*******************************************************************
!Description:   This function posts a error message to log file including
                the name of the data and the description of the problem
                when bad L1A data are read.

!Input Parameters:
     char *    name      the name of the data being read
     char *    message   error message

!Output Parameters:

!Revision History:
 (continue at top of the file)

 Revision 01.00 August 27, 1999
 Initial development
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
 

!END********************************************************************
*/
{
  char buffer [MAX_ERROR_MESSAGE_LENGTH];
  char *default_message = "Bad L1A data.\n";
  char *fmt = "\
ERROR READING DATA FROM ONE OF THE MOD01 FILES INPUT TO PGE02.\n\
   Name of MOD01 data: %s\n\
   Error:              %s\n\
This error is due to:\n\
(1) the MOD01 file is corrupted,\n\
(2) the process MOD_PR01 (in PGE01) has a bug in it, or\n\
(3) there is an address bug in MOD_PR02.";

  if (strlen(fmt) + strlen(name) + strlen(message) < MAX_ERROR_MESSAGE_LENGTH)
    sprintf(buffer, fmt, name, message);
  else
    sprintf(buffer, "%s", default_message);

  SMF_ERROR(MODIS_F_NOK, buffer);
}


PGSt_SMF_status Write_L1B_EV_Scan (int16            S,
                                   L1B_granule_t    *L1B_Gran,
                                   L1B_Scan_t       *L1B_Scan,
                                   boolean          isdaymode)
/*
!C**********************************************************************
!Description:
    This routine writes one scan of L1B EV data, including scaled integers,
    uncertainty indices and samples used, as appropriate.  For a day-mode
    scan, all resolutions are written.  For a night mode scan, only data
    for the emissive bands are written.

!Input Parameters:
      int16            S              current scan index
      L1B_granule_t    *L1B_Gran
      L1B_Scan_t       *L1B_Scan
      boolean          isdaymode      Flag denoting that scan is a day mode scan
                                      (if True)
!Output Parameters:
      
!Revision History:
 (continue at top of the file)

 Revision 03.10 May 13, 1999
 Added band 26 section.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 03.00 March 1997
 Modified to match changes in L1B_Scan_t and in file spec; 
 limited processing to EV data. 
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 02.00 1996/06/28
 New version using SDS_Services           Currently corresponds to version 1 specification
 Neal Devine(neal.devine@gsfc.nasa.gov)

 Revision 01.01 1996/04/05
 Updated to match Version q Design Document
 Joan Baden(baden@highwire.gsfc.nasa.gov)
 John Hannon(hannon@highwire.gsfc.nasa.gov)

 Revision 01.00 1993
 Initial development
 Geir E. Kvaran(geir@highwire.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   New version subsumes version 1 functions Write_L1B_PixelQuality, 
   Write_L1B_EngMem, and Write_L1B_Slice, which are no longer needed.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus;
  int16    R            = 0;
  int32    start[3]     = {0, 0, 0};
  int32    edge[3]      = {0, 0, 0};

/*------------------------------------------------------------------------------
  This macro writes a set of three EV_Aggr SDS's (SI & UI & SU)
------------------------------------------------------------------------------*/
#define  WRITE_EV_Aggr_SDS(SI_id, SI, UI_id, UI, SU_id, SU)                    \
    if(SDwritedata(SI_id,start,NULL,edge,SI) == FAIL)         \
      return(MODIS_F_WRITE_ERROR);                                              \
    if(SDwritedata(UI_id,start,NULL,edge,UI) == FAIL)         \
      return(MODIS_F_WRITE_ERROR);                                              \
    if(SDwritedata(SU_id,start,NULL,edge,SU) == FAIL)         \
      return(MODIS_F_WRITE_ERROR)
/*----------------------------------------------------------------------------*/

  if (S >= L1B_Gran->num_scans)
    return(MODIS_F_NO_MORE);

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
      /*
       * Write this scan's worth of Band 26 data.  These are always written,
       * regardless of day or night.  Note that since these are 2D SDSs,
       * we pass in the address of start[1] and edge[1], since the last
       * two dimensions of the 3D EV SDSs are the same as the two 
       * dimensions for these SDSs.
       */
  start[1] = S * DETECT_PER_BAND_AT_RES[INDEX_1000M_REFL];
  edge[1]  = DETECT_PER_BAND_AT_RES[INDEX_1000M_REFL];
  edge[2]  = EV_1km_FRAMES * BAND_RATIO_AT_RES[INDEX_1000M_REFL];
  if (SDwritedata(L1B_Scan->Band26.SI_sds_id, &start[1], NULL, &edge[1],
                  (VOIDP)L1B_Scan->Band26.SI) == FAIL)
     return(MODIS_F_WRITE_ERROR);
  if (SDwritedata(L1B_Scan->Band26.UI_sds_id, &start[1], NULL, &edge[1],
                  (VOIDP)L1B_Scan->Band26.UI) == FAIL)
     return(MODIS_F_WRITE_ERROR);

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

  for (R = 0; R < NUM_L1A_RESOLUTIONS; R++)
  {
    start[1] = S * DETECT_PER_BAND_AT_RES[R];
    edge[0]  = L1B_BANDS_AT_RES[R];
    edge[1]  = DETECT_PER_BAND_AT_RES[R];
    edge[2]  = EV_1km_FRAMES * BAND_RATIO_AT_RES[R];
    if (!isdaymode && R != INDEX_1000M_EMISS)
      continue;

    returnStatus = Write_L1B_SI_UI(S, L1B_Scan, R);
    if (returnStatus != MODIS_S_OK)
       SMF_ERROR(returnStatus, "Write_L1B_SI_UI() in Write_L1B_EV_Scan(), Granule.c"); 

    switch(R)
    {
      case INDEX_250M:
        break;

      case INDEX_500M:
        edge[0] = L1B_BANDS_AT_RES[INDEX_250M];
        WRITE_EV_Aggr_SDS (L1B_Scan->EV_250m_Aggr500m_RefSB_sds_id, 
                             L1B_Scan->EV_250m_Aggr500m_RefSB,
                           L1B_Scan->EV_250m_Aggr500m_RefSB_UI_sds_id, 
                             L1B_Scan->EV_250m_Aggr500m_RefSB_UI,
                           L1B_Scan->EV_250m_Aggr500m_RefSB_SU_sds_id, 
                             L1B_Scan->EV_250m_Aggr500m_RefSB_SU);
        break;

      case INDEX_1000M_REFL:
        edge[0] = L1B_BANDS_AT_RES[INDEX_250M];
        WRITE_EV_Aggr_SDS (L1B_Scan->EV_250m_Aggr1km_RefSB_sds_id, 
                             L1B_Scan->EV_250m_Aggr1km_RefSB,
                           L1B_Scan->EV_250m_Aggr1km_RefSB_UI_sds_id, 
                             L1B_Scan->EV_250m_Aggr1km_RefSB_UI,
                           L1B_Scan->EV_250m_Aggr1km_RefSB_SU_sds_id, 
                             L1B_Scan->EV_250m_Aggr1km_RefSB_SU);

        edge[0] = L1B_BANDS_AT_RES[INDEX_500M];
        WRITE_EV_Aggr_SDS (L1B_Scan->EV_500m_Aggr1km_RefSB_sds_id, 
                             L1B_Scan->EV_500m_Aggr1km_RefSB,
                           L1B_Scan->EV_500m_Aggr1km_RefSB_UI_sds_id, 
                             L1B_Scan->EV_500m_Aggr1km_RefSB_UI,
                           L1B_Scan->EV_500m_Aggr1km_RefSB_SU_sds_id, 
                             L1B_Scan->EV_500m_Aggr1km_RefSB_SU);
        break;

      case INDEX_1000M_EMISS:
        ;
    }/*End of switch*/
  }/*End of loop over R: L1A_resolution*/

  return(MODIS_S_OK);
}


PGSt_SMF_status Write_L1B_SI_UI (int16 S, L1B_Scan_t *L1B_Scan, int16 R)
/*
!C**********************************************************************
!Description:
    For a given resolution, this routine writes one scan of L1B EV scaled
    integer and uncertainty index data to the appropriate L1B granule.

!Input Parameters:
      int16            S              current scan index
      L1B_Scan_t       *L1B_Scan
      int16            R
!Output Parameters:
      
!Revision History:
 (continue at top of the file)

 Revision 01.00 1998
 Initial development
 Zhenying Gu(zgu@ltpmail.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  int32    start[3]     = {0, 0, 0};
  int32    edge[3]      = {0, 0, 0};
  start[1] = S * DETECT_PER_BAND_AT_RES[R];
  edge[0]  = L1B_BANDS_AT_RES[R];
  edge[1]  = DETECT_PER_BAND_AT_RES[R];
  edge[2]  = EV_1km_FRAMES * BAND_RATIO_AT_RES[R];
  
  if ((RFLAG & (1 << R)) != 0) return(MODIS_S_OK);

  switch (R){
    case 0:
      if (SDwritedata(L1B_Scan->SI_sds_id[R], start, NULL, edge, 
                               (VOIDP)L1B_Scan->SI.EV_250m_RefSB) == FAIL)
         return(MODIS_F_WRITE_ERROR);
      if (SDwritedata(L1B_Scan->UI_sds_id[R], start, NULL, edge, 
                               (VOIDP)L1B_Scan->UI.EV_250m_RefSB_UI) == FAIL)
         return(MODIS_F_WRITE_ERROR);
      break;
    
    case 1:
      if (SDwritedata(L1B_Scan->SI_sds_id[R], start, NULL, edge, 
                               (VOIDP)L1B_Scan->SI.EV_500m_RefSB) == FAIL)
         return(MODIS_F_WRITE_ERROR);
      if (SDwritedata(L1B_Scan->UI_sds_id[R], start, NULL, edge, 
                               (VOIDP)L1B_Scan->UI.EV_500m_RefSB_UI) == FAIL)
         return(MODIS_F_WRITE_ERROR);
      break;
    
    case 2:
      if (SDwritedata(L1B_Scan->SI_sds_id[R], start, NULL, edge, 
                               (VOIDP)L1B_Scan->SI.EV_1km_RefSB) == FAIL)
         return(MODIS_F_WRITE_ERROR);
      if (SDwritedata(L1B_Scan->UI_sds_id[R], start, NULL, edge, 
                               (VOIDP)L1B_Scan->UI.EV_1km_RefSB_UI) == FAIL)
         return(MODIS_F_WRITE_ERROR);
      break;

    case 3: 
      if (SDwritedata(L1B_Scan->SI_sds_id[R], start, NULL, edge, 
                               (VOIDP)L1B_Scan->SI.EV_1km_Emissive) == FAIL)
         return(MODIS_F_WRITE_ERROR);
      if (SDwritedata(L1B_Scan->UI_sds_id[R], start, NULL, edge, 
                               (VOIDP)L1B_Scan->UI.EV_1km_Emissive_UI) == FAIL)
         return(MODIS_F_WRITE_ERROR);
      break;

    default: break;
  }
  
  return(MODIS_S_OK);
} 

PGSt_SMF_status Fill_Dead_Detector_SI (boolean isdaymode,
                                       int8 *dead_detector,
                                       L1B_Scan_t *L1B_Scan,
                                       L1B_granule_t *L1B_Gran,
                                       QA_Common_t *QA_Common)
/*
!C************************************************************************
!Description:
  This function fills in reasonable pixel values in one Level 1B EV
  product file SDSs, for pixels that correspond to dead detectors.
  Values from adjacent (live) detectors are used to determine the values
  to assign to the dead-detector pixels.  If possible, a linear average
  from adjacent pixels is calculated.  This operation is applied only to
  the native resolution L1B data sets, not aggregated data sets.

!Input Parameters:
  boolean isdaymode      A flag that denotes that the scan is a
                         day-mode scan.  Some RSB SDSs are computed only
                         if the scan is a day-mode scan.
  int8 *dead_detector    Array spanning 490 MODIS detectors which
                         identifies a non-functional detector (if = 1).
  L1B_Scan_t *L1B_Scan   Scan data to be written to the EV products.
                         Includes the scaled integers that may be
                         modified.

!Output Parameters:
  L1B_Scan_t *L1B_Scan   SI pixel values corresponding to dead detectors
                         may be changed.  Other values should be unchanged.
  L1B_granule_t *L1B_Gran  Contains interpolated pixel values for each band.
  QA_Common_t *QA_Common   Contains number of dead detector EV data for each
                           detector.

!Revision History:
 (continue at top of the file)

  Revision 01.01 March 2003, Razor Issue #173
  Initialized SI, SI_prev,and SI_subs to 0 for ANSI-C compliance.
  Liqin Tan, SAIC GSO (ltan@saicmodis.com)

  Revision 01.00 June 15, 2000
  Initial development
  Jim Rogers (rogers@mcst.gfsc.nasa.gov)

!Team-unique Header:
  This software is developed by the MODIS Science Data Support
  Team for the National Aeronautics and Space Administration,
  Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
  HDF portions developed at the National Center for Supercomputing
  Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

  The change is accomplished on a scan-by-scan basis on the native-resolution
  SDSs only.  No artificial values are used in forming aggregated SI values.
  There is no "cross-over" of information between scans.  The logic is as
  follows.  If a detector is flagged as "dead" by the detector quality flag
  LUT, then the SI values of nearest "live" detectors in the same frame and
  subsample are used to form a linear average value for the dead detector SI
  value.  This linear average is accomplished independently for each frame and
  subsample of the scan line of the dead detector.  There are conditions
  checked for:

  1)  if the dead detector is either the first or last detector in the set,
  then the second or next-to-last value is simply copied (no linear average).

  2)  if an adjacent value to be used in either averaging or copying is an
  "unusable" value (>32767), then it is not used to form the artificial value.
  For example, if we started off with two values to form a linear average but
  one of the values is unusable, then the remaining valid value is simply
  copied.

  3)  If there are no usable values available to form the artificial value,
  then the pixel value remains at 65533 (the "dead detector" SI value).

  In the product, users may read the "Dead Detector List" attribute to
  determine which detectors were classified as dead and thus which scan lines
  may have the artificial values filled in.

!END**********************************************************************
*/
{
        /* Define order and number of SDSs that may need to be modified.
         * These are in different multi-dimensional arrays.
         */

typedef enum {
#ifdef WRITE_BAND_26_SDS
  FDD_RSB_BAND26,
#endif /* WRITE_BAND_26_SDS */
  FDD_RSB_250M,
  FDD_RSB_500M,
  FDD_RSB_1KM_DAY,
  FDD_TEB_1KM,
  NUM_FDD_DATA_SETS
} data_set_enum_t;
#define MAX_BANDS_PER_FDD_SET 16

        /* Declare arrays which describe features (such as dimensions)
         * of the SDSs.
         */

  int32 num_bands_per_set[NUM_FDD_DATA_SETS] = 
  {
#ifdef WRITE_BAND_26_SDS
    1,
#endif /* WRITE_BAND_26_SDS */
    NUM_250M_BANDS,
    NUM_500M_BANDS,
    NUM_1000M_REFL_BANDS,
    NUM_1000M_EMISS_BANDS
  };

  int32 num_dets_per_set[NUM_FDD_DATA_SETS] = 
  {
#ifdef WRITE_BAND_26_SDS
    DETECTORS_PER_1KM_BAND,
#endif /* WRITE_BAND_26_SDS */
    DETECTORS_PER_250M_BAND,
    DETECTORS_PER_500M_BAND,
    DETECTORS_PER_1KM_BAND,
    DETECTORS_PER_1KM_BAND
  };

  int32 num_frames_per_set[NUM_FDD_DATA_SETS] = 
  {
#ifdef WRITE_BAND_26_SDS
    EV_1km_FRAMES,
#endif /* WRITE_BAND_26_SDS */
    EV_250m_FRAMES,
    EV_500m_FRAMES,
    EV_1km_FRAMES,
    EV_1km_FRAMES
  };

        /* The following denotes whether the data set is "Day" only */

  boolean isdayonly[NUM_FDD_DATA_SETS] = 
  {
#ifdef WRITE_BAND_26_SDS
    False,
#endif /* WRITE_BAND_26_SDS */
    True,
    True,
    True,
    False
  };
    
  boolean qa_flag; /* to check if the qa parameter values for band 26 is
                    * already adjusted.
                    */

        /* We will use dead_detector to fill the array below.
         * This is static because values are determined one time
         * and then saved for subsequent scans.
         */

  static int8 dead_det_array [NUM_FDD_DATA_SETS][MAX_BANDS_PER_FDD_SET]
                             [MAX_DETECTORS_PER_BAND];

        /* The following static variable is used to determine if we have
         * already defined the above array (we only need to fill it for
         * one scan -- it is the same for all scans).
         */

  static boolean dead_det_array_is_filled = False;

        /* The following static variables is used to map the index of
         * bands through all bands to the index of bands within each L1A
         * resolution.
         */

  static int8 band_38[NUM_FDD_DATA_SETS][MAX_BANDS_PER_FDD_SET];

        /* The following static variables is used to map the index of
         * detectors through all detectors to the index of detectors within
         * each band within each L1A resolution.
         */

  static int16 det_490[NUM_FDD_DATA_SETS][MAX_BANDS_PER_FDD_SET]
                       [MAX_DETECTORS_PER_BAND];

  int32 R;          /* Index through L1A resolutions */
  int32 B;          /* Index through bands within one L1A resolution */
  int32 B_38;       /* Index through all bands, in order */
  int32 D;          /* Index through detectors within a band */
  int32 D_490;      /* Index through all detectors, in order */
  int32 B_emiss;    /* Index through L1B emissive bands */
  int32 B26_RSB;    /* index of band 26 inside L1B RSB band set */
  int32 D_prev;     /* index of previous detector. */
  int32 D_subs;     /* index of subsequent detector. */
  int32 iset;       /* index of data set being operated on */
  int32 F;          /* Frame index (frame and sample) */

  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char *location = "Fill_Dead_Detector_SI";

        /* Utility pointers -- to an array of frames */

  uint16 *SI = 0;       /* SI for scan line in question */
  uint16 *SI_prev = 0;  /* SI for previous scan line */
  uint16 *SI_subs = 0;  /* SI for subsequent scan line */

  /************************* Declarations complete ***********************/


  if (dead_det_array_is_filled == False) 
    {

        /* Assemble new dead detector array having more convenient
         * dimensionality.
         */

        /* Initialize to a bad value so we can check later */

      for (iset = 0; iset < NUM_FDD_DATA_SETS; iset++) 
      {
        for (B = 0; B < MAX_BANDS_PER_FDD_SET; B++) 
        {
          for (D = 0; D < MAX_DETECTORS_PER_BAND; D++) 
          {
            dead_det_array [iset][B][D] = -1;
          }
        }
      }

          /* Fill the 3D dead detector array. */

      B_38 = 0;
      D_490 = 0;
      R = INDEX_250M;
      for (B = 0; B < L1A_BANDS_AT_RES[R]; B++, B_38++) 
      {
        band_38[FDD_RSB_250M][B] = B_38;
        for (D = 0; D < DETECT_PER_BAND_AT_RES[R]; D++, D_490++) 
        {
          dead_det_array[FDD_RSB_250M][B][D] = dead_detector[D_490];
          det_490[FDD_RSB_250M][B][D] = D_490;
        }
      }
      R = INDEX_500M;
      for (B = 0; B < L1A_BANDS_AT_RES[R]; B++, B_38++) 
      {
        band_38[FDD_RSB_500M][B] = B_38;
        for (D = 0; D < DETECT_PER_BAND_AT_RES[R]; D++, D_490++) 
        {
          dead_det_array[FDD_RSB_500M][B][D] = dead_detector[D_490];
          det_490[FDD_RSB_500M][B][D] = D_490;
        }
      }
      R = INDEX_1000M_DAY;
      for (B = 0; B < L1A_BANDS_AT_RES[R]; B++, B_38++) 
      {
        band_38[FDD_RSB_1KM_DAY][B] = B_38;
        for (D = 0; D < DETECT_PER_BAND_AT_RES[R]; D++, D_490++) 
        {
          dead_det_array[FDD_RSB_1KM_DAY][B][D] = dead_detector[D_490];
          det_490[FDD_RSB_1KM_DAY][B][D] = D_490;
        }
      }
      R = INDEX_1000M_NIGHT;
      B_emiss = 0;
      B26_RSB = L1B_BANDS_AT_RES[INDEX_1000M_DAY] - 1;
      for (B = 0; B < L1A_BANDS_AT_RES[R]; B++, B_38++) 
      {
        if (B == MODIS_BAND26_INDEX_AT_RES) 
        {
          band_38[FDD_RSB_1KM_DAY][B26_RSB] = B_38;
#ifdef WRITE_BAND_26_SDS
          band_38[FDD_RSB_BAND26][0] = B_38;
#endif /* WRITE_BAND_26_SDS */

          for (D = 0; D < DETECT_PER_BAND_AT_RES[R]; D++, D_490++) 
          {
            dead_det_array[FDD_RSB_1KM_DAY][B26_RSB][D] = dead_detector[D_490];
            det_490[FDD_RSB_1KM_DAY][B26_RSB][D] = D_490;
#ifdef WRITE_BAND_26_SDS
            dead_det_array[FDD_RSB_BAND26][0][D] = dead_detector[D_490];
            det_490[FDD_RSB_BAND26][0][D] = D_490;
#endif /* WRITE_BAND_26_SDS */
          }
        }
        else 
        {
          band_38[FDD_TEB_1KM][B_emiss] = B_38;
          for (D = 0; D < DETECT_PER_BAND_AT_RES[R]; D++, D_490++) 
          {
            dead_det_array[FDD_TEB_1KM][B_emiss][D] = dead_detector[D_490];
            det_490[FDD_TEB_1KM][B_emiss][D] = D_490;
          }
          B_emiss++;
        }
      }

          /* Check for a bad value -- make sure we assigned each one. */

      for (iset = 0; iset < NUM_FDD_DATA_SETS; iset++) 
      {
        for (B = 0; B < num_bands_per_set[iset]; B++) 
        {
          for (D = 0; D < num_dets_per_set[iset]; D++) 
          {
            if (dead_det_array [iset][B][D] == -1)
            {
              returnStatus = MODIS_F_NOK;
              L1BErrorMsg(location, returnStatus,
                          "Dead detector array values not properly assigned.",
                          NULL, 0,  NULL, True);
              return returnStatus;
            }
          }
        }
      }

      dead_det_array_is_filled = True;
    }

        /* Main loop through all data sets. */

    for (iset = 0; iset < NUM_FDD_DATA_SETS; iset++) 
      {

        /* Skip this set if it is day only and scan is not day */

        if (isdaymode == False && isdayonly[iset] == True)
          continue;

        /* Loop through bands in this data set */

        for (B = 0; B < num_bands_per_set[iset]; B++) 
          {
            B_38 = band_38[iset][B];

            qa_flag = True;

#ifdef WRITE_BAND_26_SDS
            if (iset != FDD_RSB_BAND26 && B_38 == MODIS_BAND26_INDEX)
              qa_flag = False;
#endif /* WRITE_BAND_26_SDS */

            /* Loop through detectors */

            for (D = 0; D < num_dets_per_set[iset]; D++) 
              {

                /* If this detector is not "dead", skip it.
                 * (leave pixels as they are).
                 */

                if (dead_det_array[iset][B][D] == 0)
                  continue;

                        /* Determine D_prev by going back
                         * and finding the 1st non-dead detector.
                         */

                D_prev = D - 1;
                while (D_prev >= 0) 
                  {
                    if (dead_det_array[iset][B][D_prev] == 0)
                      break;
                    D_prev--;
                  }

                        /* Determine D_subs in a similar way. */

                D_subs = D + 1;
                while (D_subs < num_dets_per_set[iset]) 
                  {
                    if (dead_det_array[iset][B][D_subs] == 0)
                      break;
                    D_subs++;
                  }

                if (D_subs == num_dets_per_set[iset])
                  D_subs = -1;

                        /* If both D_prev and D_subs are invalid, continue
                         * (omit this detector)
                         */

                if (D_prev == -1 && D_subs == -1)
                  continue;

                if (D_prev == D_subs)       /*** Should not happen ***/
                  continue;

              /* Assign the data to be used for filling in values */

#ifdef WRITE_BAND_26_SDS
                if (iset == FDD_RSB_BAND26) 
                  {
                    SI = L1B_Scan->Band26.SI[D];
                    if (D_prev >= 0)
                      SI_prev = L1B_Scan->Band26.SI[D_prev];
                    if (D_subs >= 0)
                      SI_subs = L1B_Scan->Band26.SI[D_subs];
                  } 

                  else
#endif /* WRITE_BAND_26_SDS */
                  if (iset == FDD_RSB_250M) 
                    {
                      SI = L1B_Scan->SI.EV_250m_RefSB[B][D];
                      if (D_prev >= 0)
                        SI_prev = L1B_Scan->SI.EV_250m_RefSB[B][D_prev];
                      if (D_subs >= 0)
                        SI_subs = L1B_Scan->SI.EV_250m_RefSB[B][D_subs];
                    }

                  else if (iset == FDD_RSB_500M) 
                    {
                      SI = L1B_Scan->SI.EV_500m_RefSB[B][D];
                      if (D_prev >= 0)
                        SI_prev = L1B_Scan->SI.EV_500m_RefSB[B][D_prev];
                      if (D_subs >= 0)
                        SI_subs = L1B_Scan->SI.EV_500m_RefSB[B][D_subs];
                    }

                  else if (iset == FDD_RSB_1KM_DAY) 
                    {
                      SI = L1B_Scan->SI.EV_1km_RefSB[B][D];
                      if (D_prev >= 0)
                        SI_prev = L1B_Scan->SI.EV_1km_RefSB[B][D_prev];
                      if (D_subs >= 0)
                        SI_subs = L1B_Scan->SI.EV_1km_RefSB[B][D_subs];
                    }

                  else if (iset == FDD_TEB_1KM) 
                    {
                      SI = L1B_Scan->SI.EV_1km_Emissive[B][D];
                      if (D_prev >= 0)
                        SI_prev = L1B_Scan->SI.EV_1km_Emissive[B][D_prev];
                      if (D_subs >= 0)
                        SI_subs = L1B_Scan->SI.EV_1km_Emissive[B][D_subs];
                    }

                /* At long last, fill the array of values to be written */

                D_490 = det_490[iset][B][D];

                if (D_prev == -1) 
                  { 
                    for (F = 0; F < num_frames_per_set[iset]; F++) 
                      {
                        if (SI_subs[F] <= DN15_SAT && SI[F] == DEAD_DETECTOR_SI) 
                          {
                            SI[F] = SI_subs[F];
                            if (qa_flag) 
                              {
                                L1B_Gran->interpolated_pixels[B_38]++;
                                L1B_Gran->valid_pixels[B_38]++;
                                QA_Common->num_dead_detector_EV_data[D_490]--;
                              }
                          }
                      }
                  }

                else if (D_subs == -1) 
                  {
                    for (F = 0; F < num_frames_per_set[iset]; F++) 
                      {
                        if (SI_prev[F] <= DN15_SAT && SI[F] == DEAD_DETECTOR_SI) 
                          {
                            SI[F] = SI_prev[F];
                            if (qa_flag) 
                              {
                                L1B_Gran->interpolated_pixels[B_38]++;
                                L1B_Gran->valid_pixels[B_38]++;
                                QA_Common->num_dead_detector_EV_data[D_490]--;
                              }
                          }
                      }
                  }

                else 
                  {
                    float32 fract = ((float32) (D - D_prev)) /
                                        ((float32) (D_subs - D_prev));
                    for (F = 0; F < num_frames_per_set[iset]; F++) 
                      {
                        if (SI_subs[F] <= DN15_SAT && SI_prev[F] <= DN15_SAT &&
                            SI[F] == DEAD_DETECTOR_SI) 
                          {
                            SI[F] = (uint16)((float32)SI_prev[F] +
                                    (fract * ((float32)SI_subs[F] - 
                                    (float32)SI_prev[F])));
                            if (qa_flag) 
                              {
                                L1B_Gran->interpolated_pixels[B_38]++;
                                L1B_Gran->valid_pixels[B_38]++;
                                QA_Common->num_dead_detector_EV_data[D_490]--;
                              }
                          }

                        else if (SI_subs[F] <= DN15_SAT && 
                                   SI[F] == DEAD_DETECTOR_SI) 
                          {
                            SI[F] = SI_subs[F];
                            if (qa_flag) 
                              {
                                L1B_Gran->interpolated_pixels[B_38]++;
                                L1B_Gran->valid_pixels[B_38]++;
                                QA_Common->num_dead_detector_EV_data[D_490]--;
                              }
                          }

                        else if (SI_prev[F] <= DN15_SAT && 
                             SI[F] == DEAD_DETECTOR_SI) 
                          {
                            SI[F] = SI_prev[F];
                            if (qa_flag) 
                              {
                                L1B_Gran->interpolated_pixels[B_38]++;
                                L1B_Gran->valid_pixels[B_38]++;
                                QA_Common->num_dead_detector_EV_data[D_490]--;
                              }
                          }
                      }
                  }

              }      /* loop through D */
          }          /* loop through B */
      }              /* loop through iset */


  return returnStatus;
}

PGSt_SMF_status Open_and_Read_L1A 
                         (Run_Time_Parameters_t *runtime_params,
                          L1A_granule_t         *L1A_Gran,
                          boolean               *skip_night_hi_res)
    
/*
!C**************************************************************************
!Description: This function opens the L1A granule to be processed and reads
              in data for all members of the L1A_granule_t structure. The
              file remains open (both SD and Vdata) upon function exit.

!Input Parameters:
   Run_Time_Parameters_t *runtime_params   Values read from PCF file.

!Output Parameters:
     L1A_granule_t *L1A_Gran         contains SD file interface ID, Vdata
                                     file interface ID, number of scans and
                                     other data relating to the middle
                                     granule

!Revision History:
 (continue at top of the file)

   Revision 01.02 June 4, 2002    Razor Issue #169
   Changed logic to check whether any scans are in day mode and if so, set 
   skip_night_hi_res to False.
   Gwyn Fireman & Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)
   
   Revision 01.01 November 12, 2001       (Razor Issue #169)
   Added boolean skip_night_hi_res as output variable.
   Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

   Revision 01.00 June 27, 2000
   Initial development.
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
!END********************************************************************
*/
{
  PGSt_SMF_status    returnStatus = MODIS_S_OK;
  intn               hdf_return   = 0;
  PGSt_integer       Version      = 1;
  char               file_name[PGSd_PC_FILE_PATH_MAX];
  int32              S;                  /* scan index */
  int32              middle_scan;
  int32              geo_satellite_id;
  char               *location = "Open_and_Read_L1A";
  char               *satelliteInstrument;

      /* 
       * Will be set to TRUE by PCF file only if writing 250m and 
       * 500m files when in night mode is desired. 
       */                        
  boolean            write_nightmode_hires = 0;  
  boolean            is_nightmode_data     = 0;
  
  /*
   * Prepare for opening file: convert logical ID to physical name
   */

  returnStatus = PGS_PC_GetReference (FIRST_L1A_GRANULE, &Version, file_name);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus, "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*
   * Open file for HDF science data (SD) access.
   */

  L1A_Gran->sd_id = SDstart(file_name,DFACC_RDONLY); /*for SD interface */
  if (L1A_Gran->sd_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus, "Could not open file for SD read access.",
                "SDstart", FIRST_L1A_GRANULE,
                "The file may be missing, corrupted or not an HDF-4 file.", True);
    return returnStatus;
  }

  /*
   * Open file for HDF Vdata access and initialize Vdata interface.
   * Call Hopen() and Vstart() in the same place.
   */

  L1A_Gran->v_id = Hopen(file_name,DFACC_RDONLY,0);
  if (L1A_Gran->v_id == FAIL)
  {
    SDend (L1A_Gran->sd_id);
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus, "Could not open file for Vdata read access.",
                "Hopen", FIRST_L1A_GRANULE,
                "The file may be corrupted or not an HDF-4 file.", True);
    return returnStatus;
  }

  hdf_return = Vstart (L1A_Gran->v_id); /*initialize internal structures*/
  if (hdf_return == FAIL)
  {
    SDend (L1A_Gran->sd_id); Hclose (L1A_Gran->v_id);
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus, "Could not initialize Vdata interface.",
                "Vstart", FIRST_L1A_GRANULE,
                "The file may be corrupted or not an HDF-4 file.", True);
    return returnStatus;
  }

  /*
   * Read num_scans
   */

  returnStatus = read_attribute (L1A_Gran->sd_id, "Number of Scans",
                                 DFNT_INT32, (void *)&L1A_Gran->num_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not read Number of Scans.",
                "read_attribute", FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  if (L1A_Gran->num_scans <= 0)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    MOD_PR02_Failure_Exit_Code = 233;
    L1BErrorMsg(location, returnStatus,
                "Number of Scans in the MOD01 granule is less than "
                "or equal to 0.\n"
                "The MOD01 granule is probably empty (no valid data)",
                NULL, FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  if (L1A_Gran->num_scans > MAX_NUM_SCANS)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "Number of Scans in the MOD01 granule is greater "
                "than MAX_NUM_SCANS.",
                NULL, FIRST_L1A_GRANULE,
                "If the file is valid, the code macro must be increased.", True);
    return returnStatus;
  }

  /*
   * Read num_day_scans
   */

  returnStatus = read_attribute (L1A_Gran->sd_id, "Number of Day mode scans",
                                 DFNT_INT32, (VOIDP)&L1A_Gran->num_day_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not read Number of Day "
                "mode scans.",
                "read_attribute", FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  /*
   * Set skip_night_hi_res flag
   */

  if (strcmp(runtime_params->Write_Night_Mode_HiRes_Data, "0") == 0)
    write_nightmode_hires = False;  
  else
    write_nightmode_hires = True;

  if (L1A_Gran->num_day_scans == 0)
    is_nightmode_data =  True;
  else
    is_nightmode_data =  False;

  if (is_nightmode_data == True && write_nightmode_hires == False) 
    *skip_night_hi_res = True;
  else
    *skip_night_hi_res = False;

  /*
   * Read L1A ScanType
   */

  returnStatus = read_sds_rank2 (L1A_Gran->sd_id, "Scan Type",
                                 L1A_Gran->num_scans, SCAN_TYPE_TEXT_SIZE,
                                 (void *)L1A_Gran->ScanType);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not read Scan Type.",
                "read_sds_rank2", FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  for (S = 0; S < L1A_Gran->num_scans; S++)
  {
    if (strcmp(L1A_Gran->ScanType[S], "Day") != 0   &&
        strcmp(L1A_Gran->ScanType[S], "Night") != 0 &&
        strcmp(L1A_Gran->ScanType[S], "Other") != 0)
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus,
                  "Scan Type is not \"Day\", \"Night\", or \"Other\".",
                  NULL, FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
      return returnStatus;
    }
  }

  /*
   * Read L1A EVStartTime_TAIsecond
   */

  returnStatus = read_sds_rank1 (L1A_Gran->sd_id, "EV start time",
                                 L1A_Gran->num_scans,
                                 (void *)L1A_Gran->EVStartTime_TAIsecond);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not read EV start time.",
                "read_sds_rank1", FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  /* 
   * Determine data collection time. The valid EV start time for the scan
   * closest to the middle scan is chosen to be the data collection time.
   * This time is used for time-dependent look up tables.
   */

  middle_scan = L1A_Gran->num_scans/2;
  L1A_Gran->data_collection_time = 0;
  for (S = 0; middle_scan-S >= 0 && 
      middle_scan+S < L1A_Gran->num_scans; S++) {
    if (L1A_Gran->EVStartTime_TAIsecond[middle_scan-S] > 0) {
      L1A_Gran->data_collection_time = 
          L1A_Gran->EVStartTime_TAIsecond[middle_scan-S];
      break;
    }
    else if (L1A_Gran->EVStartTime_TAIsecond[middle_scan+S] > 0) {
      L1A_Gran->data_collection_time = 
          L1A_Gran->EVStartTime_TAIsecond[middle_scan+S];
      break;
    }
  }

  if (L1A_Gran->data_collection_time == 0) {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "No valid or only one valid EV start time in the "
                "array \"EV start time\"",
                NULL, 0, NULL, True);
    return returnStatus;
  }
 
  /*
   * Read L1A MirrorSide
   */

  returnStatus = read_sds_rank1 (L1A_Gran->sd_id, "Mirror side",
                                 L1A_Gran->num_scans,
                                 (void *)L1A_Gran->MirrorSide);
  if ( returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not read Mirror side.",
                "read_sds_rank1", FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  for (S = 0; S < L1A_Gran->num_scans; S++)
  {
    if (L1A_Gran->MirrorSide[S] != 0 && L1A_Gran->MirrorSide[S] != 1
                                     && L1A_Gran->MirrorSide[S] != -1)
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus, "Mirror side is not [-1, 0, 1].",
                  NULL, FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
      return returnStatus;
    }
  }

  /*
   * Read num_night_scans
   */

  returnStatus = read_attribute (L1A_Gran->sd_id, 
                                 "Number of Night mode scans",
                                 DFNT_INT32, 
                                 (VOIDP)&L1A_Gran->num_night_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, 
                "Could not read Number of Night mode scans.",
                "read_attribute", FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  /*
   * Read incomplete_scans
   */

  returnStatus = read_attribute (L1A_Gran->sd_id, 
                                 "Incomplete Scans", DFNT_INT32,
                                 (VOIDP)&L1A_Gran->incomplete_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, 
                "Could not read number of incomplete scans.",
                "read_attribute", 
                FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  /*
   * Read max_ev_frames
   */

  returnStatus = read_attribute (L1A_Gran->sd_id, 
                                 "Max Earth Frames", DFNT_INT32,
                                 (VOIDP)&L1A_Gran->max_ev_frames);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, 
                "Could not read maximum number of EV frames.",
                "read_attribute", FIRST_L1A_GRANULE, 
                Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  /* Read Extract Metadata */
  read_extract_metadata(L1A_Gran->sd_id,
                        &L1A_Gran->Extract_Pixel_Offset,
                        &L1A_Gran->Extract_Pixel_Count,
                        &L1A_Gran->Extract_Line_Offset,
                        &L1A_Gran->Extract_Line_Count);

  /*
   * Read scan quality array
   */

  returnStatus = read_sds_rank2(L1A_Gran->sd_id, "Scan quality array",
                                L1A_Gran->num_scans, 
                                SCAN_QUALITY_ARRAY_NUM_ELEMENTS,
                                (VOIDP)L1A_Gran->scan_quality);
  if (returnStatus != MODIS_S_OK) {
    L1BErrorMsg(location, returnStatus, "Could not read scan quality array",
                "read_sds_rank2", FIRST_L1A_GRANULE, NULL, True);
    return returnStatus;
  }

  /*
   * Determine if a scan is completely missing based on the second element
   * of Scan quality array.
   */

  for (S = 0; S < L1A_Gran->num_scans; S++)
  {
    if (L1A_Gran->scan_quality[S][0] == 0)
      L1A_Gran->missing_scan[S] = True;
    else
      L1A_Gran->missing_scan[S] = False;
  }

  /*
   * Determine the satellite
   */

  returnStatus = Get_Satellite_ID(FIRST_L1A_GRANULE, &L1A_Gran->satellite_id);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not get satellite ID.",
                "Get_Satellite_ID", FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  /*
   * Determine the satellite instrument value in geo file.
   */

  returnStatus = Get_Satellite_ID(GEOLOCATION_FILE, &geo_satellite_id);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, "Could not get satellite ID.",
                "Get_Satellite_ID", GEOLOCATION_FILE,
                "The geolocation file is invalid.", True);
    return returnStatus;
  }
  
  /*
   * Check if the geo granule and L1A granule is for the same instrument.
   */

  if (geo_satellite_id != L1A_Gran->satellite_id)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "*** INCORRECT GEOLOCATION FILE ***\n" 
                "The satellite instrument in the geolocation file (lun 600000) is\n"
                "different from that in the middle L1A granule (lun 500001).", 
                NULL, 0, NULL, True);
    return returnStatus;
  } 

  /*
   * Check the run time parameter SatelliteInstrument from the PCF file against
   * the satellite inherent in the L1A data.
   */

  satelliteInstrument = runtime_params->SatelliteInstrument;
  if (!((strcmp(satelliteInstrument, "AM1M") == 0 && 
         L1A_Gran->satellite_id == TERRA) ||
        (strcmp(satelliteInstrument, "PM1M") == 0 &&
         L1A_Gran->satellite_id == AQUA)))
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus, 
                "The satellite instrument set in the PCF does not match\n"
                "the value inherent in the middle L1A granule.", 
                NULL, 0, NULL, True);
    return returnStatus;
  }

  RFLAG = 0;
  if (SDnametoindex(L1A_Gran->sd_id, "EV_250m") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "BB_250m") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SV_250m") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SD_250m") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SRCA_250m") == -1) {
    RFLAG = RFLAG | 1;
  }

  if (SDnametoindex(L1A_Gran->sd_id, "EV_500m") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "BB_500m") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SV_500m") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SD_500m") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SRCA_500m") == -1) {
    RFLAG = RFLAG | 2;
  }

  if (SDnametoindex(L1A_Gran->sd_id, "EV_1km_day") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "BB_1km_day") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SV_1km_day") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SD_1km_day") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SRCA_1km_day") == -1) {
    RFLAG = RFLAG | 4;
  }

  if (SDnametoindex(L1A_Gran->sd_id, "EV_1km_night") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "BB_1km_night") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SV_1km_night") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SD_1km_night") == -1 ||
      SDnametoindex(L1A_Gran->sd_id, "SRCA_1km_night") == -1) {
    RFLAG = RFLAG | 8;
  }

  if (SDfindattr(L1A_Gran->sd_id, "B_13_scale_parm") != -1)
    RSCL_FLAG = RSCL_FLAG | 1;
  if (SDfindattr(L1A_Gran->sd_id, "B_16_scale_parm") != -1)
    RSCL_FLAG = RSCL_FLAG | 2;
  if (SDfindattr(L1A_Gran->sd_id, "B_10_scale_parm") != -1)
    RSCL_FLAG = RSCL_FLAG | 4;
  if (SDfindattr(L1A_Gran->sd_id, "B_12_scale_parm") != -1)
    RSCL_FLAG = RSCL_FLAG | 8;

  returnStatus = MODIS_S_OK;  /* Successful completion of function */
  return returnStatus;        /* L1A granule remains open */
}

PGSt_SMF_status Get_Satellite_ID (PGSt_PC_Logical lun, int32 *satellite_ID)
/*
!C**************************************************************************
!Description: This function reads the metadata "SHORTNAME" from the file
              designated by the lun and determines the satelitte name.

!Input Parameters:
   int32 lun            logic identifier of the file.

!Output Parameters:
   int32 *satellite_ID  identifier of the satellite (TERRA = 0, AQUA = 1,
                        INVALID_SATELLITE_ID = -1)

!Revision History:
 (continue at top of the file)

   Revision 01.00 June 27, 2000
   Initial development.
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  PGSt_integer       Version      = 1;
  char shortname[10];
  char *shortnameptr;
  char *attrName = "SHORTNAME";
  char *location = "Get_Satellite_ID";

  /* To use the address of the shortname[], a pointer is needed. */

  shortnameptr = shortname;

  /* Get the shortname */

  returnStatus = PGS_MET_GetPCAttr(lun, Version, "CoreMetadata.0",
                                   attrName, &shortnameptr);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_READ_ERROR;
    L1BErrorMsg (location, returnStatus, "Could not get metadata \"SHORTNAME\".",
                 "PGS_MET_GetPCAttr", lun, NULL, True);
  }

  if (strncmp(shortname, "MOD", 3) == 0)
    *satellite_ID = TERRA;
  else if (strncmp(shortname, "MYD", 3) == 0)  
    *satellite_ID = AQUA;
  else
    *satellite_ID = INVALID_SATELLITE_ID; 

  return MODIS_S_OK;
}

PGSt_SMF_status Read_Run_Time_Parameters (Run_Time_Parameters_t *runtime_params)
/*
!C**************************************************************************
!Description:
   Read all run-time parameters from the PCF file and set the values in the
   output variable.

!Input Parameters:
   none

!Output Parameters:
   Run_Time_Parameters_t *runtime_params  
                 contains SatelliteInstrument,
                          ReprocessingPlanned,
                          ReprocessingActual,
                          PGE02_Version,
                          MCST_LUT_Version,
                          Write_Night_Mode_HiRes_Data,
                          ProcessingEnvironment
!Revision History:
 (continue at top of the file)

   Revision 01.07  October 31, 2003 Razor Issue #195
   Added the code of reading the ProcessingCenter value from the
   Run_Time_Parameters_t structure.
   Liqin Tan, SAIC GSO (ltan@saicmodis.com)

   Revision 01.06  October 2, 2002  Razor Issue #183
   Bug fix to setting ProcessingEnvironment.  The function used, while POSIX
     compliant, was on the list of functions forbidden by the SDP Toolkit.  
     "getenv" was subsituted.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)  
 
   Revision 01.05  August 19, 2002  Razor Issue #183
   Revised setting ProcessingEnvironment.  It is now set by calling
     the appropriate POSIX function from within the module, rather
     than read from the PCF file.  This is to avoid seeing error messages
     from the SDP Toolkit for each granule when the field is left blank
     in the PCF file, which is allowable.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)  

   Revision 01.04  April 25, 2002  Razor Issue #183
   Corrected syntax for setting ProcessingEnvironment to blank
   Gwyn Fireman, SAIC GSO (fireman@mcst.gsfc.nasa.gov)

   Revision 01.03  March 21, 2002  Razor Issue #181
   Added ProcessingEnvironment parameter as per SDST
   request for Collection 4
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)  
  
   Revision 01.02  November 12, 2001
   Added Write_Night_Mode_HiRes_Data (Razor issue #169)
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)
   
   Revision 01.01  November 6, 2001
   Added MCST_LUT_Version (Razor issue #167)
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.00 February 21, 2001
   Initial development.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus;

  char  *location = "Read_Run_Time_Parameters";
        /* Default ProcessingEnvironment identifier */
  char  thisname[MAX_PROCESSING_ENVIRONMENT_STRLEN] = "";     
  char  *ename;  /* Value of environment variable from "getenv" calls */
  char  *osnames[4] = {"OS", "OSTYPE", "BRAND", "MACHINE"};
  char  *hostnames[3] = {"HOST", "HOSTNAME", "HOSTTYPE"};
  char  *revision = "REVISION";
  
  int16 j = 0;      /* A counter */  

  if (!runtime_params)
  {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, "Input runtime_params is NULL.",
                NULL, 0, NULL, True);
    return returnStatus;
  }

  /*
   * Read and check the SatelliteInstrument value.
   */

  returnStatus = PGS_PC_GetConfigData(SATELLITE_INSTRUMENT_LUN,
                                      runtime_params->SatelliteInstrument);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "Failed to get parameter SatelliteInstrument from PCF file.",
                "PGS_PC_GetConfigData", SATELLITE_INSTRUMENT_LUN,
                "Invalid PCF file", True);
    return returnStatus;
  }
  if (!(strcmp(runtime_params->SatelliteInstrument, "AM1M") == 0 ||
        strcmp(runtime_params->SatelliteInstrument, "PM1M") == 0))
  {
    char errmsg[256];
    sprintf (errmsg, 
             "The SatelliteInstrument value (%s) in the PCF is invalid.",
             runtime_params->SatelliteInstrument);
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus, errmsg,
                NULL, SATELLITE_INSTRUMENT_LUN, NULL, True);
    return returnStatus;
  }

  /*
   * Set the ProcessingEnvironment value.  This 
   * is not a mandatory entry and it is acceptable that 
   * the value be missing or blank, but in practice it is easier to set it
   * within L1B and avoid seeing error messages when it is not in the PCF file.
   */

   strcpy((runtime_params->ProcessingEnvironment), "");

    /*
     *    Obtain the processing environment fields 
     *        by calling the POSIX "getenv" function
     */ 

   for (j = 0; j < 4; j++)
   {
     ename = getenv(osnames[j]);

     if (ename != NULL)

     {
      /* Concatenate the desired names onto the ProcessingEnvironment string. */
      safe_strcat(thisname, ename, MAX_PROCESSING_ENVIRONMENT_STRLEN); 
       safe_strcat(thisname, " ", MAX_PROCESSING_ENVIRONMENT_STRLEN);
       
       break;
     }
   }

   for (j = 0; j < 3; j++)
   {
     ename = getenv(hostnames[j]);

     if (ename != NULL)

     {
      /* Concatenate the desired names onto the ProcessingEnvironment string. */
       safe_strcat(thisname, ename, MAX_PROCESSING_ENVIRONMENT_STRLEN); 
       safe_strcat(thisname, " ", MAX_PROCESSING_ENVIRONMENT_STRLEN);
       
       break;
     }
   }

   ename = getenv(revision);

   if (ename != NULL)

   {
    /* Concatenate the desired names onto the ProcessingEnvironment string. */
     safe_strcat(thisname, ename, MAX_PROCESSING_ENVIRONMENT_STRLEN); 
     safe_strcat(thisname, " ", MAX_PROCESSING_ENVIRONMENT_STRLEN);
   }

    /*
     * If the function getenv fails to return any values it is OK.  
     */
   
    safe_strcat(runtime_params->ProcessingEnvironment, thisname,
       MAX_PROCESSING_ENVIRONMENT_STRLEN);
    
  /*
   * Read the ReprocessingPlanned value.  This is not checked for specific
   * values because it is not clear that there is a constant set of valids
   * for this field (some could be added in the future).
   */

  returnStatus = PGS_PC_GetConfigData(REPROCESSING_PLANNED_LUN,
                                      runtime_params->ReprocessingPlanned);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "Failed to get parameter ReprocessingPlanned from PCF file.",
                "PGS_PC_GetConfigData", REPROCESSING_PLANNED_LUN,
                "Invalid PCF file", True);
    return returnStatus;
  }
  if (strlen(runtime_params->ReprocessingPlanned) == 0)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "The ReprocessingPlanned value in the PCF is invalid (no string).",
                NULL, REPROCESSING_PLANNED_LUN, NULL, True);
    return returnStatus;
  }

  /*
   * Read the ReprocessingActual value.  This is not checked for specific
   * values because it is not clear that there is a constant set of valids
   * for this field (some could be added in the future).
   */

  returnStatus = PGS_PC_GetConfigData(REPROCESSING_ACTUAL_LUN,
                                      runtime_params->ReprocessingActual);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "Failed to get parameter ReprocessingActual from PCF file.",
                "PGS_PC_GetConfigData", REPROCESSING_ACTUAL_LUN,
                "Invalid PCF file", True);
    return returnStatus;
  }
  if (strlen(runtime_params->ReprocessingActual) == 0)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "The ReprocessingActual value in the PCF is invalid (no string).",
                NULL, REPROCESSING_ACTUAL_LUN, NULL, True);
    return returnStatus;
  }

  /*
   * Read the Write_Night_Mode_HiRes_Data value.  This value must be 
   * True (1) or False (0).
   */
  returnStatus = PGS_PC_GetConfigData(WRITE_NIGHT_HIRES_LUN,
                         runtime_params->Write_Night_Mode_HiRes_Data);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "Failed to get parameter Write_Night_Mode_HiRes_Data "
                "from PCF file.",
                "PGS_PC_GetConfigData", WRITE_NIGHT_HIRES_LUN,
                "Invalid PCF file", True);
    return returnStatus;
  }
  if (!(strcmp(runtime_params->Write_Night_Mode_HiRes_Data, "0") == 0 ||
        strcmp(runtime_params->Write_Night_Mode_HiRes_Data, "1") == 0))
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "The Write_Night_Mode_HiRes_Data value in the PCF file is "
                "invalid (must be 0 (False) or 1 (True)).",
                NULL, WRITE_NIGHT_HIRES_LUN, NULL, True);
    return returnStatus;
  }
  
  /*
   * Read the PGE02 Version value from PCF file.
   */

  returnStatus = PGS_PC_GetConfigData(PGE02_VERSION_LUN,
                                      runtime_params->PGE02_Version);

  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "Failed to get parameter PGE02_Version from PCF file.",
                "PGS_PC_GetConfigData", PGE02_VERSION_LUN,
                "Invalid PCF file", True);
    return returnStatus;
  }
  if (strlen(runtime_params->PGE02_Version) == 0)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "The PGE02_Version value in the PCF is invalid (no string).",
                NULL, PGE02_VERSION_LUN, NULL, True);
    return returnStatus;
  }

  /*
   * Check the PGE02 version given in the PCF file against the code macro.
   */

  if (strcmp(runtime_params->PGE02_Version, PGE02_VERSION))
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    char errmsg[512];
    sprintf(errmsg, "The PGE02 version in PCF file (%s) does not match the expected value (%s).",
                    runtime_params->PGE02_Version, PGE02_VERSION);
    L1BErrorMsg(location, returnStatus,
      errmsg, NULL, 0,
      "Code files must have the same PGE02 version as specified in the PCF file.",
      True);

    return returnStatus;
  }

  /*
   * Read the MCST_LUT_Version value.  This is not checked for specific
   * values because the values vary with the MCST LUT Version numbers.
   */

  returnStatus = PGS_PC_GetConfigData(MCST_LUT_VERSION_LUN,
                                      runtime_params->MCST_LUT_Version);

  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "Failed to get parameter MCST_LUT_Version from PCF file.",
                "PGS_PC_GetConfigData", MCST_LUT_VERSION_LUN,
                "Invalid PCF file", True);
    return returnStatus;
  }
  if (strlen(runtime_params->MCST_LUT_Version) == 0)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "The MCST_LUT_Version value in the PCF is invalid (no string).",
                NULL, MCST_LUT_VERSION_LUN, NULL, True);
    return returnStatus;
  }

  /*
   * Read the ProcessingCenter value from the Run_Time_Parameters_t
   * structure.
   */

  returnStatus = PGS_PC_GetConfigData(PROCESSING_CENTER_LUN,
                                      runtime_params->ProcessingCenter);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "Failed to get parameter ProcessingCenter from PCF file.",
                "PGS_PC_GetConfigData", PROCESSING_CENTER_LUN,
                "Invalid PCF file", True);
    return returnStatus;
  }
  if (strlen(runtime_params->ProcessingCenter) == 0)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus,
                "The ProcessingCenter value in the PCF is invalid (no string).",
                NULL, PROCESSING_CENTER_LUN, NULL, True);
    return returnStatus;
  }

    return MODIS_S_OK;
}

