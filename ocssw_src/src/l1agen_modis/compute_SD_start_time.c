#include "L1A_prototype.h"
#include "hdf.h"
#include "PH_pkt_hdr.h"
#include "PGS_TD.h"


void  compute_SD_start_time (PH_PACKET_HEADER_t *pkt_header,
                             PGSt_double *SD_start_time)

/*
!C**********************************************************************

!Description:  This function will compute the SD start time for the input 
               packet.

!Input Parameters:
               PH_PACKET_HEADER_t   pkt_header      **  Packet header  **

!Output Parameters:
               PGSt_double          SD_start_time   **  SD start time  **

Return Values:
               None

Externally Defined:
               PH_SEC_PKT_TYPE_ENG1_GROUP           (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG2_GROUP           (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH     (PH_pkt_hdr.h)
               PH_PACKET_HEADER_t                   (PH_pkt_hdr.h)
               PGSt_double                          (PGS_TD.h)
               global_time_offset_array             (level1a)

Called By:
               load_eng_data
               process_a_granule
               process_a_scan

Routines Called:
               None

!Revision History:
               revision 1.0 1997/08/29  17:30:00
               Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
               Original development

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:

!Design Notes: 

!END***************************************************************************
*/
{
  /*****************************************************************************/
  /*                                                                           */
  /*                   Declare Global Variables                                */
  /*                                                                           */
  /*****************************************************************************/

  extern PGSt_double  global_time_offset_array[5];  
					 /*  array containing the time         */
				         /*   offsets for each sector from     */
				         /*   the beginning of the scan        */
				         /*   (ie. SD sector)                  */


  /*****************************************************************************/
  /*                                                                           */
  /*          Declare the local variables and initialize them.                 */
  /*                                                                           */
  /*****************************************************************************/

  int8             index;


  /*****************************************************************************/
  /*                                                                           */
  /* IF (PH_PACKET_HEADER_t.pkt_type equals PH_SEC_PKT_TYPE_ENG1_GROUP or      */
  /*    PH_SEC_PKT_TYPE_ENG2_GROUP) or (PH_PACKET_HEADER_t.source_ID_type_flag */
  /*    equals PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH)                               */
  /* THEN                                                                      */
  /*    index = 4                                                              */
  /* ELSE                                                                      */
  /*    index = PH_PACKET_HEADER_t.cal_type                                    */
  /* ENDIF                                                                     */
  /*                                                                           */
  /*****************************************************************************/

  if ((pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG1_GROUP) ||
      (pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG2_GROUP) ||
      (pkt_header->source_ID_type_flag == PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH))
    index = 4;
  else
    index = pkt_header->cal_type;


  /*****************************************************************************/
  /*                                                                           */
  /* SD_start_time = PH_PACKET_HEADER_t.pkt_TAI_time -                         */
  /*    global_time_offset_array[index]                                        */
  /*                                                                           */
  /*****************************************************************************/

  *SD_start_time = pkt_header->pkt_TAI_time - global_time_offset_array[index];

  return;
}  
