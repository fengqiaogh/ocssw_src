#ifndef REFLECTIVE_CALP_H
#define REFLECTIVE_CALP_H

#include "L1B_Tables.h"  /* Contains declaration of refl_tables_t and 
                            common_QA_tables_t */

/*
!C-INC***********************************************************************
!Description:   Private header file to be used only in Reflective_Cal.c.
                Contains prototype of functions.

!Revision History:
 $Log: Reflective_CalP.h,v $
 Revision 1.3  2006-10-27 11:01:38-04  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.

 Revision 01.05, Dec 7, 2000, Razor issue 146, 147
 New SWIR algorithm, in-lined the SWIR algorithm (removed
 SWIR_out_of_band_correction).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.04 Sep 29, 1999
 Changed the name of "SWIR_correction" to "SWIR_out_of_band_correction"
 as per Bruce B's request.  Removed UI_ptr from argument list since we
 expect new SWIR uncertainty algorithm to require new module.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 May 7, 1999
 Moved DN_to_DN_star in because it is only used in Reflective_Cal now.
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.02 Feb 22, 1999
 Two function protos were deleted: Scaled_Int_to_Radiance (replaced by single
 statement) and Set_L1B_UI_Left_Bit (never used, obsolete).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 Feb 17, 1999
 Changed names of Get_L1B_SI and Get_L1B_UI to "Set...".  Deleted
 duplicate declarations that were present in Reflective_Cal.h.
 Moved ADC_correction to PreprocessP.h as part of moving DN_to_DN_star
 to Preprocess.c to avoid circular include dependencies.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 Nov 1998 
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

#endif /* REFLECTIVE_CALP_H */


