#include "L1A_prototype.h"
#include "SC_scan.h"
#include "MD_metadata.h"


void   initialize_scan (SC_SCAN_DATA_t           *L1A_scan,
                        SC_PIXEL_QUALITY_DATA_t  *scan_pixel,
                        MD_SCAN_MET_t            *scan_meta)

/*
!C*****************************************************************************

!Description:  
        This function calls other subroutines to fill in the Scan Level
        Metadata (section 2 of the MODIS Level 1A Data Product Format), Pixel
        Quality Data (section 3 of the MODIS Level 1A Data Product Format), and
        Scan Data (section 4 of the MODIS Level 1A Data Product Format) with
        the appropriate fill data.

!Input Parameters:
               None

!Output Parameters:
               SC_SCAN_DATA_t           *L1A_scan    ** Scan Data           **
               SC_PIXEL_QUALITY_DATA_t  *scan_pixel  ** Pixel Quality Data  **
               MD_SCAN_MET_t            *scan_meta   ** Scan Level Metadata **

Return Values: 
               None

Externally Defined:  
               SC_SCAN_DATA_t           (SC_scan.h)
               SC_PIXEL_QUALITY_DATA_t  (SC_scan.h)
               MD_SCAN_MET_t            (MD_metadata.h)

Called By:
               process_a_scan

Routines Called:
               initialize_scan_metadata
	       initialize_scan_data

!Revision History:
  $Log: initialize_scan.c,v $
  Revision 4.1  2003/03/10 16:05:54  vlin
  Touch up description in the prologue.

  Revision 4.0  2002/12/03 16:27:17  vlin
  called initialize_scan_data() to fill the Scan Data structure
  with fill value (-1).
  vlin@saicmodis.com

               Revision 2.2  2001/01/04
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Simplified code. Got rid of 2 function calls.

               Revision 2.1  1997/09/08  17:50 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code. Baseline from Version 1.

!Team-unique Header:

        This software is developed by the MODIS Science Data Support Team 
        for the National Aeronautics and Space Administration, 
        Goddard Space Flight Center, under contract NAS5-32373.

References and Credits: None

Design Notes: 

!END
***************************************************************************/

{

  initialize_scan_metadata (scan_meta);
  initialize_pixel_qual_data(scan_pixel);
  initialize_scan_data(L1A_scan);

}
