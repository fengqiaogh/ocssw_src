#ifndef   EN_ENG_DATA_H
#define   EN_ENG_DATA_H

/*
!C-INC**********************************************************************

!Description:  This include file contains information about eng_data (vdata)
               structures.

!Input Parameters:  N/A

!Output Parameters: N/A

Externally Defined:
               M01LAST_VALID_SCAN        (mapiL1A.h)
               int8                      (hdfi.h)
               uint8                     (hdfi.h)
               int16                     (hdfi.h)
               uint16                    (hdfi.h)
               int32                     (hdfi.h)
               uint32                    (hdfi.h)

!Revision History:
  $Log: EN_eng_data.h,v $
  Revision 3.1  2002/09/20 21:36:02  vlin
  2 fields "revision" and "instrument" added to Vdata structure.

               Revision 2.1  2000/12/27
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Added #define for EN_NUM_MAJ_CYCLES 32 for use in
               module update_eng_data_for_maj_cycle_n to get rid of
               magic numbers.

               Revision 2.0 1998/10/26   11:29 EST
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Added #defines to mark new vdatas in the EN_DATA_LIST
               and to use #defines from mapiL1A.h where possible. Also
               added order and type to eng_data field structure to 
               handle data with different data sizes and orders. Also
               several #defines to assist in picking the needed Vdatas
               out of the eng_data structure.

               Revision 1.0  1997/07/21  15:58 EDT
               David Catozzi/SAIC/GSC (cato@ltpmail.gsfc.nasa.gov)
               Originated code development.

!Team-unique Header:

    This software is developed by the MODIS Science Data Support Team 
    for the National Aeronautics and Space Administration, 
    Goddard Space Flight Center, under contract NAS5-32373.

!END
*************************************************************************/

#include "hdfi.h"
#include "mapiL1A.h"
#include "PGS_PC.h"
#include "PGS_TYPES.h"


/********************** GLOBAL CONSTANTS (MACROS): *****************************/

#define  EN_NUM_VDATAS                   68
#define  EN_ENG_VDATA_START_INDEX        51  /* the vdata_array index of the first
                                             (eng pkt 1-2) engineering Vdata */

#define  EN_NUM_MAJ_CYCLES               32

#define  EN_MAX_FIELDS_PER_VDATA         70
#define  EN_MAX_VDATA_NAME_LENGTH        80
#define  EN_MAX_FIELD_NAME_LENGTH        50

/* LAST_VALID_SCAN flag values */
#define  EN_INITIAL_LAST_VALID_SCAN_VALUE      65535
#define  EN_ENG_DATA_IS_FROM_PREVIOUS_GRANULE      0

#define  EN_LAST_VALID_SCAN             M01LAST_VALID_SCAN

#define  EN_INITIAL_FIELD_VALUE                    0

#define  EN_CP_HK_TLMY_PRIOR_OFFSET     512 

#define  EN_EMPTY_SLOT                  "EMPTY_SLOT"

#define  EN_MAX_LINE_LEN                255

#define  EN_SC_ANCILLARY_VDATA_START    48
#define  EN_SC_ANCILLARY_VDATA_END      49
#define  EN_COMM_PROC_VDATA_NUMBER      50
#define  EN_ANCIL_VDATA_BUFFER_SIZE     71

#define  EN_MAX_VDATA_NUM_BITS          64
#define  EN_MIN_VDATA_NUM_BITS           1
#define  EN_NUM_SCAN_ELEMENTS            5
#define  EN_MAX_VDATA_ORDER              8
#define  EN_MIN_VDATA_ORDER              1
#define  EN_MAX_VDATA_TYPE              25
#define  EN_MIN_VDATA_TYPE              20
/****************************************************************************
 ***                           VDATA STRUCTURE:                           ***
 ****************************************************************************/

/* Vdata field structure: */
typedef struct{
  char    field_name[EN_MAX_FIELD_NAME_LENGTH];
  uint16  num_bits;
  uint16  start_bit_pos;
  uint16  order;
  uint16  type;
  union {
    int8     i8type;
    uint8    ui8type;
    int16    i16type;
    uint16   ui16type;
    int32    i32type;
    uint32   ui32type;
  } union_value[EN_MAX_VDATA_ORDER];
  uint16   value;
} EN_FIELD_TYPE_t;

/* Vdata structure */
typedef struct{
  char             vdata_name[EN_MAX_VDATA_NAME_LENGTH];
  uint8            num_fields;
  EN_FIELD_TYPE_t  field[EN_MAX_FIELDS_PER_VDATA];
  char revision[PGSd_PC_VALUE_LENGTH_MAX];
  PGSt_tag instrument;
} EN_VDATA_TYPE_t;


#endif  /* EN_ENG_DATA_H */
