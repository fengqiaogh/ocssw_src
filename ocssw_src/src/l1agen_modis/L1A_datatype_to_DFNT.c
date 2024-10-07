#include "L1A_prototype.h"
#include <string.h>
#include "mfhdf.h"
#include "mapi.h"

#define NULLstring(c) (((c)==NULL) || (*(c) == '\0'))


int32   L1A_datatype_to_DFNT(char *datatype)

/*
!C************************************************************************

!Description:  L1A_datatype_to_DFNT takes a M-API data type string and
               returns the HDF number type integer used by HDF routines 
               to identify number types.
  
!Input Parameters:
               char *datatype      pointer to M-API data type string 
                                   permitted types: "int8"
                                                    "uint8"
                                                    "int16"
                                                    "uint16"
                                                    "int32"
                                                    "uint32"
                                                    "int64"
                                                    "uint64"
                                                    "float32"
                                                    "float64"
                                                    "char *"
  
!Output Parameters:
               None

Return Values:
               HDF number type or MFAIL if an error occurs

Externally Defined:
               DFNT_INT8          (hdf.h)
               DFNT_UINT8         (hdf.h)
   	       DFNT_INT16         (hdf.h)
               DFNT_UINT16        (hdf.h)
               DFNT_INT32         (hdf.h)
               DFNT_UINT32        (hdf.h)
               DFNT_INT64         (hdf.h)
               DFNT_UINT64        (hdf.h)
               DFNT_FLOAT32       (hdf.h)
               DFNT_FLOAT64       (hdf.h)
               DFNT_CHAR8         (hdf.h)
 	       I8	          (mapi.h)
               UI8                (mapi.h)
 	       I16                (mapi.h)
               UI16               (mapi.h)
               I32                (mapi.h)
               UI32               (mapi.h)
               I64                (mapi.h)
               UI64               (mapi.h)
               R32                (mapi.h)
               R64                (mapi.h)
               TXT                (mapi.h)
               MFAIL              (mapi.h)

Called By:
               create_Vdata_field

Routines Called:
               None

!Revision History:
               Revision 1.0  1997/10/01  11:25 EDT
               Tom Johnson/GSC     (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate walkthrough comments

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes:
               This routine is a exact copy of the MAPI routine
               datatype_to_DFNT.

!END************************************************************************
*/

{
  int32 output=0;   /* used to return HDF datatype */

  /* if datatype points to a nullstring return MFAIL */  
  if (NULLstring(datatype)) return(MFAIL);

  /* select correct HDF datatype */
  if ((strcmp(datatype,I8)==0)) output=DFNT_INT8;
  else if ((strcmp(datatype,UI8)==0)) output=DFNT_UINT8;
  else if ((strcmp(datatype,I16)==0)) output=DFNT_INT16;
  else if ((strcmp(datatype,UI16)==0)) output=DFNT_UINT16;
  else if ((strcmp(datatype,I32)==0)) output=DFNT_INT32;
  else if ((strcmp(datatype,UI32)==0)) output=DFNT_UINT32;
/************************************************************************ 
   64 bit integer support disabled until HDF decides to support it.
   To re-enable, remove comments around the following block of code.

  else if ((strcmp(datatype,I64)==0)) output=DFNT_INT64;
  else if ((strcmp(datatype,UI64)==0)) output=DFNT_UINT64;
************************************************************************/
  else if ((strcmp(datatype,R32)==0)) output=DFNT_FLOAT32;
  else if ((strcmp(datatype,R64)==0)) output=DFNT_FLOAT64;
  else if ((strcmp(datatype,TXT)==0)) output=DFNT_CHAR8;
  else output=MFAIL;

  return(output);

} /* end L1A_datatype_to_DFNT */

