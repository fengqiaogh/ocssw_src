#include "L1A_prototype.h"
#include "hdf.h"
#include "hdfi.h"


int16   equal_strings ( char  *a, 
                        char  *b )

/*
!C************************************************************************

!Description:  This functions takes two strings up to 79 characters long and
	       determines if they are equal.

!Input Parameters:
               char  *a                   ** the first of the two strings    **
               				  ** to compare                      **

               char  *b                   ** the second of the two strings   **
               				  ** to compare                      **

!Output Parameters:
               None

Return Values: 
               TRUE   (1)                     (hdf.h)
               FALSE  (0)                     (hdf.h)

Externally Defined:  
               int16                          (hdfi.h)

Called By:
               assign_data_type
               get_index
               recall_id

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/10/01  17:25 EDT
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
  /*                                                                        */
  /* Set value to FALSE                                                     */
  /*                                                                        */
  /**************************************************************************/
  
  int16  value = FALSE;



  /**************************************************************************/
  /*                                                                        */
  /* IF string a is equal to string b (for up to 79 characters)             */
  /* THEN                                                                   */
  /*    Set value to TRUE                                                   */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* return value                                                           */
  /*                                                                        */
  /**************************************************************************/

     if ((a && b) != FALSE)            /* if both a and b is non-NULL then  */
        value = (!strncmp(a, b, 79));  
     
     return value;

 }  /* End of routine equal_strings */
