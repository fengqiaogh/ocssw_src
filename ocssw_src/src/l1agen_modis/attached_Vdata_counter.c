#include "L1A_prototype.h"
#include "hdf.h"
#include "hdfi.h"
#include "VU_vdata_utility.h"


int16  attached_Vdata_counter ( int16 action )

/*
!C************************************************************************

!Description:  This function keeps track of how many Vdatas are currently 
               attached to the HDF file.  As an internal consistency check
               it is useful to check that there are no attached Vdatas prior
               to closing the HDF file (since detaching Vdatas is required 
               by HDF and one may forget to do it).

!Input Parameters:
               int16    action              ** either VU_INCREMENT, **
               				    **        VU_DECREMENT, **
                                            **     or VU_REPORT     **
!Output Parameters:
               None

Return Values: 
               FAIL                         (hdf.h)
               int16    result              ** the number of attached Vdatas **
                                            ** (a non-negative integer) or   **
                                            **  FAIL (-1)                    **
Externally Defined:
               int16                        (hdfi.h)
               VU_MAX_NUMBER_OF_VDATAS      (VU_vdata_utility.h)
               VU_INCREMENT                 (VU_vdata_utility.h)
               VU_DECREMENT                 (VU_vdata_utility.h)
               VU_REPORT                    (VU_vdata_utility.h)

Called By:
               remember
               forget
               end_vdata_access_to_file

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/10/02  11:35 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

               Revision 1.1  1997/09/03  10:55
               Tom Johnson/GSC     (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate walkthrough comments

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
               This PDL cannot be implemented in a language that doesn't 
               provide persistent local variables (i.e. local variables that
               maintain their values between calls).

               This is a private function - only to be used by functions 
               internal to the Vdata_id_table.

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /* declare total_number_of_attached_Vdatas to be a persistent local       */
  /* variable of type int16 with initial value of 0                         */
  /*                                                                        */
  /**************************************************************************/

  static int16 total_number_of_attached_Vdatas = 0;



  /**************************************************************************/
  /* SWITCH (action)                                                        */
  /*    case VU_INCREMENT:                                                  */
  /*       IF (total_number_of_attached_Vdatas >= VU_MAX_NUMBER_OF_VDATAS)  */
  /*       THEN                                                             */
  /*          set total_number_of_attached_Vdatas to FAIL                   */
  /*       ELSE                                                             */
  /*          increment total_number_of_attached_Vdatas                     */
  /*       ENDIF                                                            */
  /*                                                                        */
  /*       break                                                            */
  /*                                                                        */
  /**************************************************************************/

     switch (action) {
        case VU_INCREMENT:
           if (total_number_of_attached_Vdatas >= VU_MAX_NUMBER_OF_VDATAS)
              total_number_of_attached_Vdatas = FAIL;
           else
              total_number_of_attached_Vdatas++;
  
           break;



  /**************************************************************************/
  /*    case VU_DECREMENT:                                                  */
  /*       IF (total_number_of_attached_Vdatas <= 0)                        */
  /*       THEN                                                             */
  /*          set total_number_of_attached_Vdatas to FAIL                   */
  /*       ELSE                                                             */
  /*          decrement total_number_of_attached_Vdatas                     */
  /*       ENDIF                                                            */
  /*                                                                        */
  /*       break                                                            */
  /*                                                                        */
  /**************************************************************************/

        case VU_DECREMENT:
           if (total_number_of_attached_Vdatas <= 0)
              total_number_of_attached_Vdatas = FAIL;
           else
              total_number_of_attached_Vdatas--;
           break;


  /**************************************************************************/
  /*    default:                                                            */
  /*       {the default case is VU_REPORT}                                  */
  /*       {it returns the number of attached Vdatas}                       */
  /*                                                                        */
  /*       break                                                            */
  /*                                                                        */
  /* END_SWITCH                                                             */
  /*                                                                        */
  /* return total_number_of_attached_Vdatas                                 */
  /*                                                                        */
  /**************************************************************************/

        default:

        break;
     }

     return total_number_of_attached_Vdatas;

 } /* End of routine attached_Vdata_counter */
