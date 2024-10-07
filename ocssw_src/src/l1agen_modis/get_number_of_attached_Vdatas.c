#include "L1A_prototype.h"
#include "hdfi.h"
#include "VU_vdata_utility.h"


int16   get_number_of_attached_Vdatas (void)

/*
!C************************************************************************

!Description:  This function is part of the Vdata Utility Package. It returns 
               the number of Vdatas currently attached to the HDF file.
               As an internal consistency check it is useful to check that 
               there are no attached Vdatas prior to closing the HDF file 
               (since detaching Vdatas is required by HDF and one may forget 
               to do it).

!Input Parameters:
               None

!Output Parameters:
               None

!Input/Output Parameters:
               None

Return Values: 
               int16    result           ** the number of attached Vdatas ** 
                                         ** (a non-negative integer)      **

Externally Defined:
               int16                        (hdfi.h)
               VU_REPORT                    (VU_vdata_utility.h)

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/10/02  10:45 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

               Revision 1.0  1997/07/14  15:58 EDT
               David Catozzi/SAIC/GSC (cato@ltpmail.gsfc.nasa.gov)
               Original design.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes:
               None

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/

  int16  result;



  /**************************************************************************/
  /*                                                                        */
  /* return  attached_Vdata_counter( VU_REPORT )                            */
  /*                                                                        */
  /**************************************************************************/

     result = attached_Vdata_counter(VU_REPORT);

     return result;

 } /* End of routine get_number_of_attached_Vdatas */
