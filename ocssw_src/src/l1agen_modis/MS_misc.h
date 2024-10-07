#ifndef MS_MISC_H    
#define MS_MISC_H    

/*
!C-INC************************************************************************

!Description:  This header file contains miscellaneous MACRO definitions
               for the L1A software.

!Input Parameters:
               N/A

!Output Parameters:
               N/A

Return Values: 
               N/A

Externally Defined:
               N/A

Called By:
               N/A

Routines Called:
               N/A

!Revision History:
               Revision 1.0  1997/09/24  14:45 EDT
               Timi Adelekan/GSC/SAIC (adelekan@ltpmail.gsfc.nasa.gov)
               Created include file from Version 1 include files:
               modis_init.h,

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The ".h" file below was specifically written for development
               in C. Any other language choice may require reworking of the
               ".h" file before coding can begin.

!END**********************************************************************
*/


#define	MS_HEADER_BUFF_SIZE		10000 
#define	MS_FOOTER_BUFF_SIZE 		    1

#define MS_MT_EXIT_FATAL_FAILURE            1
#define MS_MT_EXIT_SUCCESS                  0

#define MS_BLANK                          ' '
#define MS_COMMENT                        '-'
#define MS_NEW_LINE                      '\n'

#endif /* MS_MISC_H */
