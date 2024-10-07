#include "L1A_prototype.h"
#include "hdf.h"
#include "PGS_IO_L0.h"
#include "EN_eng_data.h"
#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"

void   update_eng_data (uint16             index,
                        PGSt_IO_L0_Packet  *eng_packet,
                        uint16             scan_number,
                        EN_VDATA_TYPE_t    *eng_data, 
                        int                use_cp_prior_offset)

/*
!C*****************************************************************************

!Description:  This function updates the Vdata referenced by 'index' in
               the eng_data (array) structure with field values extracted
               from the eng_packet, and updates the Vdata's LAST_VALID_SCAN
               field with the scan_number. If use_cp_prior_offset is True
               then the bit offset EN_CP_HK_TLMY_PRIOR_OFFSET will be added
               to the start bit positions of the Vdata fields.

!Input Parameters:
               uint16              index        ** The eng_data (array) index 
               					   of the eng Vdata to 
               					   update **
               PGSt_IO_LO_Packet  *eng_packet   ** The eng packet from which 
               					   to extract values to update
               					   the eng data with **
               uint16              scan_number  ** The scan number (counting 
               					   from 1) of the current scan
               					   within the current 
               					   granule **

               boolean             use_cp_prior_offset  ** whether or not to use
                                                      the bit offset: 
                                                      EN_CP_HK_TLMY_PRIOR_OFFSET
                                                      (should only be set to True
                                                      when updating eng_data for 
                                                      the "prior" part of eng pkt
                                                      2-1 **

!Output Parameters:
               None

!Input/Output Parameters:
               EN_VDATA_TYPE_t     *eng_data    ** The Vdata array structure **

Return Values: 
               None

Externally Defined:
               EN_VDATA_START_INDEX         (EN_eng_data.h)
               EN_CP_HK_TLMY_PRIOR_OFFSET   (EN_eng_data.h)
               EN_VDATA_TYPE_t              (EN_eng_data.h)
               PGSt_IO_LO_Packet            (PGS_IO.h)
               PD_NUM_BITS_IN_BYTE          (PD_pkt_data.h)
               PD_PKT_CONTENTS_BYTE_OFFSET  (PD_pkt_data.h)

Called By:
               process_cp_hk_tlmy
               process_sci_eng_data
               update_eng_data_for_maj_cycle_n

Routines Called:
               extr_bits

!Revision History:
               revision 1.0 1997/09/12  17:30:00
               Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
               Original development

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes:    
               (c.f. CDRL Table 30-5D and the Vdata_list file)

!END************************************************************************
*/
{
  /******************************************************************************/
  /*                                                                            */
  /*               Define and Initialize Local Variables                        */
  /*                                                                            */
  /******************************************************************************/

  int                  bit_offset;
  int                  num_fields_this_vdata;
  int                  start_pos;
  int                  num_bits;
  int                  start_byte;
  int                  start_bit;
  int                  i;


  /******************************************************************************/
  /*                                                                            */
  /*   The sci eng Vdatas' start bit positions are relative to the start of the */
  /*   packet header, while those of the cp hk telemetry Vdatas' are relative   */
  /*   to the start of the packet's data field... so if the vdata to update is  */
  /*   one of the former (as shown by the index being >= ENG_VDATA_START_INDEX),*/
  /*   then subtract the number of header bits (144) to get the packet data     */
  /*   start bit position.                                                      */
  /*                                                                            */
  /*        -----------------------------------------------------------------   */
  /*                                                                            */
  /*  IF ( index < EN_ENG_VDATA_START_INDEX ) THEN                              */
  /*    set bit_offset to 0                                                     */
  /*  ELSE                                                                      */
  /*    set bit_offset to -144                                                  */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (index < EN_ENG_VDATA_START_INDEX)
    bit_offset = 0;
  else
    bit_offset = -144;


  /******************************************************************************/
  /*                                                                            */
  /*     if use_cp_prior_offset is True then add EN_CP_HK_TLMY_PRIOR_OFFSET to  */
  /*     the bit_offset                                                         */
  /*                                                                            */
  /*        -----------------------------------------------------------------   */
  /*                                                                            */
  /*  IF ( use_cp_prior_offset )                                                */
  /*  THEN                                                                      */
  /*    increment bit_offset by EN_CP_HK_TLMY_PRIOR_OFFSET                      */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (use_cp_prior_offset)
    bit_offset += EN_CP_HK_TLMY_PRIOR_OFFSET;


  /******************************************************************************/
  /*                                                                            */
  /*     set the LAST_VALID_SCAN field (the first field in                      */
  /*     every Vdata) to scan_number                                            */
  /*                                                                            */
  /*        -----------------------------------------------------------------   */
  /*                                                                            */
  /*  set eng_data[index].field[0] to scan_number                               */
  /*                                                                            */
  /*  set num_fields_this_vdata to eng_data[index].num_fields                   */
  /*                                                                            */
  /******************************************************************************/

  eng_data[index].field[0].value = scan_number;
  num_fields_this_vdata = eng_data[index].num_fields;


  /******************************************************************************/
  /*                                                                            */
  /*  DO FOR ( i = 1 to (num_fields_this_vdata - 1) )                           */
  /*                                                                            */
  /*    set start_pos to (eng_data[index].field[i].start_bit_pos + bit_offset)  */
  /*                                                                            */
  /*    set num_bits to eng_data[index].field[i].num_bits                       */
  /*                                                                            */
  /*    calculate start_byte and start_bit from start_pos                       */
  /*                                                                            */
  /*    CALL extr_bits to get the bits from the packet                          */
  /*      INPUTS:  eng_packet, start_byte, start_bit, num_bits                  */
  /*      OUTPUTS: None                                                         */
  /*      RETURN:  value                                                        */
  /*    set eng_data[index].field[i].value to the returned value                */
  /*                                                                            */
  /*  END DO                                                                    */
  /*                                                                            */
  /******************************************************************************/

  for (i=1; i<num_fields_this_vdata; i++)
    {
      start_pos = eng_data[index].field[i].start_bit_pos + bit_offset;
      num_bits = eng_data[index].field[i].num_bits;
      start_byte = (start_pos / PD_NUM_BITS_IN_BYTE) + PD_PKT_CONTENTS_BYTE_OFFSET;
      start_bit = start_pos % PD_NUM_BITS_IN_BYTE;
      eng_data[index].field[i].value = extr_bits(eng_packet,start_bit,start_byte,
                                                 num_bits);
    }

  return;
}
