#include "L1A_prototype.h"
#include "EN_eng_data.h"
#include "PGS_IO_L0.h"
#include "hntdefs.h"
#include "PGS_SMF.h"
#include "PGS_IO.h"
#include "PGS_MODIS_35005.h"

void process_group2_packet1_vdata ( PGSt_IO_L0_Packet *pkt, EN_VDATA_TYPE_t *eng_data)

/*
!C************************************************************************
!Description:  This function adds the data for the Prior and Current Ancillary
               Vdata to the Vdata array (eng_data). This procedure extracts the
               data from the packet, and stores the data in the eng_data structure.

!Input Parameters:
               PGSt_IO_L0_Packet    *pkt         ** The current eng     **
                                                 ** packet              **

!Output Parameters:
               EN_VDATA_TYPE_t      *eng_data  ** Stores Vdata values **

!Input/Output Parameters:
               None

Return Values:
               None

Externally Defined:
               PGSt_IO_L0_Packet               (PGS_IO_L0.h)
               EN_VDATA_TYPE_t                 (EN_eng_data.h)
               EN_SC_ANCILLARY_VDATA_START     (EN_eng_data.h)
               EN_SC_ANCILLARY_VDATA_END       (EN_eng_data.h)
               PD_NUM_BITS_IN_BYTE             (PD_pkt_data.h)
               PD_PKT_CONTENTS_BYTE_OFFSET     (PD_pkt_data.h)
               DFNT_INT8                       (hntdefs.h)
               DFNT_UINT8                      (hntdefs.h)
               DFNT_INT16                      (hntdefs.h)
               DFNT_UINT16                     (hntdefs.h)
               DFNT_INT32                      (hntdefs.h)
               DFNT_UINT32                     (hntdefs.h)
               EN_MAX_VDATA_ORDER              (EN_eng_data.h)
               MODIS_E_INVALID_VDATA_ORDER     (PGS_MODIS_35005.h)
               MODIS_E_INVALID_VDATA_TYPE      (PGS_MODIS_35005.h)

Called By:
               process_eng_packet

Routines Called:
               extr_bits
               log_fmt_msg

!Revision History:
               Revision 1.0  1998/10/26  11:45:00 EST
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)   
               Original design.
    
!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.
        
!References and Credits:
               None
                        
!Design Notes:
               This routine is written in ANSI C.
               This routine assumes that the vdatas with orders greater than 1
                 all have type of uint8. If this changes, code changes will need
                 to be made to the start_pos calculation to replace k with
                 the number of bytes in the field.
               This routine assumes that the data currently stored in the
                 int16 union value field is really 12 bit signed data. Sign
                 extension is performed on all data stored in the int16 union value.

!END************************************************************************
*/
{
  int i,j, k;                                     /* loop control variables */
  int start_pos;                        /* starting position of data in pkt */
  int num_bits;                       /* number of bits to extract from pkt */
  int start_byte;             /* starting byte for data extraction from pkt */
  int start_bit;               /* starting bit for data extraction from pkt */
  char msg[300];                           /* array to store error messages */
  char *routine = "process_group2_packet1_vdata";
  uint16 ui16_value;      /* unsigned int16 value for proper sign extension */

/****************************************************************************/
/* Check for NULL input paramaters.                                         */
/****************************************************************************/

if ((pkt == NULL) || (eng_data == NULL)) {
  log_fmt_msg(MODIS_E_NULL_POINTER, routine, " ");
  return;
}

/****************************************************************************/
/* Loop through the 2 S/C Ancillary Data Vdatas.                            */
/****************************************************************************/
  for (i=EN_SC_ANCILLARY_VDATA_START; i <= EN_SC_ANCILLARY_VDATA_END; i++) { 

/****************************************************************************/
/* Loop through all fields for each vdata.                                  */
/****************************************************************************/
    for (j=0; j< eng_data[i].num_fields; j++) {

/****************************************************************************/
/* Loop through all orders for each field.                                  */
/****************************************************************************/
      if ((eng_data[i].field[j].order > 0) && (eng_data[i].field[j].order <= 
                                              EN_MAX_VDATA_ORDER)) {

/****************************************************************************/
/* Calculate where in the packet the data for this vdata field exists.      */
/****************************************************************************/
        num_bits = eng_data[i].field[j].num_bits / eng_data[i].field[j].order;
        for (k=0; k < eng_data[i].field[j].order; k++) {
          start_pos = eng_data[i].field[j].start_bit_pos;
          start_pos += (PD_NUM_BITS_IN_BYTE * k);
          start_byte = (start_pos / PD_NUM_BITS_IN_BYTE) +
                                                  PD_PKT_CONTENTS_BYTE_OFFSET;
          start_bit = start_pos % PD_NUM_BITS_IN_BYTE;

/****************************************************************************/
/* Extract data from packet into proper union type value field.             */
/****************************************************************************/
          switch (eng_data[i].field[j].type) {
            case DFNT_INT8 : eng_data[i].field[j].union_value[k].i8type = 
                              (int8)  extr_bits(pkt,start_bit,start_byte,num_bits);
                     break;
            case DFNT_UINT8 :
                     eng_data[i].field[j].union_value[k].ui8type = 
                              (uint8) extr_bits(pkt,start_bit,start_byte,num_bits);
                     break;
            case DFNT_INT16 :
                     ui16_value = (uint16) extr_bits(pkt,start_bit,start_byte,num_bits);
                     if (ui16_value > 0x800)
                       eng_data[i].field[j].union_value[k].i16type = 
                             (int16) (ui16_value - 0x1000);
                     else
                       eng_data[i].field[j].union_value[k].i16type = (int16) ui16_value;
                     break;
            case DFNT_UINT16 :
                     eng_data[i].field[j].union_value[k].ui16type = 
                             (uint16) extr_bits(pkt,start_bit,start_byte,num_bits);
                     break;
            case DFNT_INT32 :
                     eng_data[i].field[j].union_value[k].i32type = 
                             (int32)  extr_bits(pkt,start_bit,start_byte,num_bits);
                     break;
            case DFNT_UINT32 :
                     eng_data[i].field[j].union_value[k].ui32type = 
                             (uint32) extr_bits(pkt,start_bit,start_byte,num_bits);
                     break;
            default          :
                     sprintf(msg, "Vdata Name = %s, Field Name = %s", eng_data[i].vdata_name,
                                  eng_data[i].field[j].field_name);
                     log_fmt_msg(MODIS_E_INVALID_VDATA_TYPE, routine, msg);
                     break;

          } /* end switch */

        } /* for k */

      } /* if order */

/**********************************************************************************/
/* Invalid order value error.                                                     */
/**********************************************************************************/
      else {
        sprintf(msg, "Vdata Name = %s, Field Name = %s", eng_data[i].vdata_name,
                                  eng_data[i].field[j].field_name);
        log_fmt_msg(MODIS_E_INVALID_VDATA_ORDER, routine, msg);
      } /* else */

    } /* end for j*/

  } /* end for i */

  return;

}/* eng process_group2_packet1_vdata.c */

