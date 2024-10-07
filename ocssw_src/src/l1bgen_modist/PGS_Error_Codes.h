#ifndef PGS_ERROR_CODES_H
#define PGS_ERROR_CODES_H

/*
!C-INC***********************************************************************
!Description:   Header file PGS_Error_Codes.h - defines error codes.
                (see design notes).

!Revision History:
 $Log: PGS_Error_Codes.h,v $
 Revision 1.3  2006-10-30 10:07:47-05  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.

 Revision 02.11, Feb 2, 1999
 Made this the standard include for all L1B .c files.  Invoked the include
 for PGS_MODIS_36100.h within this file.  Defined fallbacks for all macros.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 April, 13 1998
 Added this standard prolog.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 01.00 March 1997
 Initial Development.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST) for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   All error codes should be defined in the PGS_MODIS_36100.h header file.
   On acrobat, and possibly other machines, this file had incorrect names.
   Thus, fallbacks were defined for all L1B error codes.  The fallbacks
   match the correct defines in PGS_MODIS_36100.h.

!END********************************************************************
*/

#include    "PGS_MODIS_36100.h"

/*
 * Fallbacks:
 */

#ifndef MODIS_F_OUT_OF_MEMORY
#define MODIS_F_OUT_OF_MEMORY                  295735296 /* 0x11a09000L */
#endif

#ifndef MODIS_F_MEM_FREE_FAIL
#define MODIS_F_MEM_FREE_FAIL                  295735297 /* 0x11a09001L */
#endif

#ifndef MODIS_F_FILE_NOT_FOUND
#define MODIS_F_FILE_NOT_FOUND                 295735298 /* 0x11a09002L */
#endif

#ifndef MODIS_F_READ_ERROR
#define MODIS_F_READ_ERROR                     295735299 /* 0x11a09003L */
#endif

#ifndef MODIS_F_WRITE_ERROR
#define MODIS_F_WRITE_ERROR                    295735300 /* 0x11a09004L */
#endif

#ifndef MODIS_F_OUT_OF_RANGE
#define MODIS_F_OUT_OF_RANGE                   295735301 /* 0x11a09005L */
#endif

#ifndef MODIS_W_OUT_OF_RANGE
#define MODIS_W_OUT_OF_RANGE                   295734278 /* 0x11a08c06L */
#endif

#ifndef MODIS_W_TIME_INCORRECT
#define MODIS_W_TIME_INCORRECT                 295734279 /* 0x11a08c07L */
#endif

#ifndef MODIS_F_NOK
#define MODIS_F_NOK                            295735304 /* 0x11a09008L */
#endif

#ifndef MODIS_F_HDF_ERROR
#define MODIS_F_HDF_ERROR                      295735305 /* 0x11a09009L */
#endif

#ifndef MODIS_F_NO_MORE
#define MODIS_F_NO_MORE                        295735306 /* 0x11a0900aL */
#endif

#ifndef MODIS_S_NO_MORE
#define MODIS_S_NO_MORE                        295731723 /* 0x11a0820bL */
#endif

#ifndef MODIS_F_FILE_NOT_OPENED
#define MODIS_F_FILE_NOT_OPENED                295735308 /* 0x11a0900cL */
#endif

#ifndef MODIS_F_FILE_NOT_CREATED
#define MODIS_F_FILE_NOT_CREATED               295735309 /* 0x11a0900dL */
#endif

#ifndef MODIS_F_INVALID_ARGUMENT
#define MODIS_F_INVALID_ARGUMENT               295735310 /* 0x11a0900eL */
#endif

#ifndef MODIS_E_TESTING
#define MODIS_E_TESTING                        295734799 /* 0x11a08e0fL */
#endif

#ifndef MODIS_S_OK
#define MODIS_S_OK                             295731728 /* 0x11a08210L */
#endif

#endif  /* _PGS_ERROR_CODES_H */ /* DO NOT ADD ANYTHING PAST THIS LINE */

