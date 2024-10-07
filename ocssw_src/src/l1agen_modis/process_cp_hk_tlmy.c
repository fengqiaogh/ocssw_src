#include "L1A_prototype.h"
#include "hdf.h"
#include "PGS_IO_L0.h"
#include "EN_eng_data.h"


void   process_cp_hk_tlmy (EN_VDATA_TYPE_t   *eng_data,
                           PGSt_IO_L0_Packet *eng_pkt_2_1,
                           uint16             scan_number)
/*
!C*****************************************************************************

!Description:  This function processes eng group 2 pkt #1. The appropriate
               eng_data Vdatas are updated for this scan using the data in 
               the engineering packet's Current/Prior HK Tlmy section.
               (c.f. CDRL Table 30-6A and the Vdata_list file).

!Input Parameters:
               PGSt_IO_LO_Packet  *eng_pkt_2_1   ** eng grp 2 packet #1 
               					    for the current scan **

               uint16             scan_number    ** The scan number (counting 
               					    from 1) of the current scan
               					    within the current 
               					    granule **
!Output Parameters:
               None

!Input/Output Parameters:
               EN_VDATA_TYPE_t    *eng_data      ** The Vdata array 
                     				    structure  **

Return Values:
               None

Externally Defined:
               EN_VDATA_TYPE_t                 (EN_eng_data.h)
               PGSt_IO_LO_Packet               (PGS_IO.h)

Called By:
               process_eng_packet

Routines Called:
               extr_bits
               update_eng_data
               update_eng_data_for_maj_cycle_n

!Revision History:
               revision 1.0 1997/09/11  17:30:00
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

      Use the current Major Cycle (first 6 bits of eng_pkt_2_1's data) to
    decide which Vdatas to update:

   ------------------------------------------------------------------------
     If Major Cycle is:      | then update Vdatas at vdata_array indices:
   ------------------------------------------------------------------------
        all                  |    0, 1, 2
   ------------------------------------------------------------------------
     Major Cycle (MOD 8) =   |   
         0                   |    3
         1                   |    4
         2                   |    5  
         3                   |    6, 7, 8
         4                   |    9, 10
         5                   |    11, 12
         6                   |    13
         7                   |    14
   ------------------------------------------------------------------------
       0..32                 |    15..47  
   ------------------------------------------------------------------------  

    Then use the 'prior' Major Cycle (first 6 bits of the "prior" section (bits
  512..517) of eng_pkt_2_1->data (i.e. bits 656..661 from the start of the 
  packet header)) to decide which Vdatas to update using the above map, except
  this time do not update the "all" category (Vdatas 0, 1, or 2).

      (c.f. CDRL Table 30-6A, CDRL Table 20-4, and the Vdata_list) 

!END************************************************************************
*/
{
  /*****************************************************************************/
  /*                                                                           */
  /*               Define and Initialize Local Variables                       */
  /*                                                                           */
  /*****************************************************************************/

  int             is_cp_hk_prior_section = FALSE;
  int             start_byte;
  int             start_bit;
  int             num_bits;
  int             i;
  uint16          major_cycle;
  

  /*****************************************************************************/
  /*                                                                           */
  /*  set is_cp_hk_prior_section to FALSE   (done during declaration)          */
  /*                                                                           */
  /*****************************************************************************/

  /*****************************************************************************/
  /*                                                                           */
  /*  DO FOR ( i = 0 to 2 )    to update the three "Telemetry Major Cycle All" */ 
  /*                           Vdatas                                          */
  /*    CALL update_eng_data to update the eng_data (vdata_array)              */
  /*       structure with the eng data values in the eng pkt                   */
  /*      INPUTS:  i, eng_pkt_2_1, scan_number, eng_data,                      */
  /*               is_cp_hk_prior_section                                      */
  /*      OUTPUTS: eng_data                                                    */
  /*      RETURN:  None                                                        */
  /*  END DO                                                                   */
  /*                                                                           */
  /*****************************************************************************/

  for (i=0; i<3; i++)
    update_eng_data(i,eng_pkt_2_1,scan_number,eng_data,is_cp_hk_prior_section);


  /*****************************************************************************/
  /*                                                                           */
  /*  set start_byte to 18      <-- 144/8 = 18 = the start of the data section */
  /*  set start_bit to 0                                                       */
  /*  set num_bits to 6                                                        */
  /*                                                                           */
  /*****************************************************************************/

  start_byte = 18;
  start_bit = 0;
  num_bits = 6;


  /*****************************************************************************/
  /*                                                                           */
  /*  CALL extr_bits to get the current Major Cycle (the first 6 bits of       */
  /*                                          eng_pkt_2_1 's data section)     */
  /*    INPUTS:  num_bits, start_byte, start_bit, eng_pkt_2_1                  */
  /*    OUTPUTS: None                                                          */
  /*    RETURN:  major_cycle                                                   */
  /*                                                                           */
  /*****************************************************************************/

  major_cycle = extr_bits(eng_pkt_2_1,start_bit,start_byte,num_bits);


  /*****************************************************************************/
  /*                                                                           */
  /*  CALL update_eng_data_for_maj_cycle_n to update the appropriate "Telemetry*/
  /*     Major Cycle X of 63" and "Telemetry Major Cycle X of 7" Vdata(s) for  */
  /*     the "current" section of eng packet 2-1                               */
  /*    INPUTS:  major_cycle, eng_pkt_2_1, scan_number, eng_data,              */
  /*                                       is_cp_hk_prior_section              */
  /*    OUTPUTS: eng_data                                                      */
  /*    RETURN:  None                                                          */
  /*                                                                           */
  /*****************************************************************************/

  update_eng_data_for_maj_cycle_n( major_cycle,eng_pkt_2_1,scan_number,eng_data,
                                   is_cp_hk_prior_section);


  /*****************************************************************************/
  /*                                                                           */
  /*  set start_byte to 82   <-- 656/8 = 82 = the start of the "prior" part of */
  /*                            eng packet 2-1's data section                  */
  /*  set start_bit to 0                                                       */
  /*  set num_bits to 6                                                        */
  /*                                                                           */
  /*****************************************************************************/

  start_byte = 82;
  start_bit = 0;
  num_bits = 6;


  /*****************************************************************************/
  /*                                                                           */
  /*  CALL extr_bits to get the prior Major Cycle for the eng packet's         */
  /*     "prior" section                                                       */
  /*    INPUTS:  num_bits, start_byte, start_bit, eng_pkt_2_1                  */
  /*    OUTPUTS: None                                                          */
  /*    RETURN:  major_cycle                                                   */
  /*                                                                           */
  /*****************************************************************************/

  major_cycle = extr_bits(eng_pkt_2_1,start_bit,start_byte,num_bits);


  /*****************************************************************************/
  /*                                                                           */
  /*  set is_cp_hk_prior_section to True                                       */
  /*                                                                           */
  /*****************************************************************************/

  is_cp_hk_prior_section = TRUE;


  /*****************************************************************************/
  /*                                                                           */
  /*  CALL update_eng_data_for_maj_cycle_n to update the appropriate "Telemetry*/
  /*     Major Cycle X of 63" and "Telemetry Major Cycle X of 7" Vdata(s) for  */
  /*     the "prior" section of eng packet 2-1                                 */
  /*    INPUTS:  major_cycle, eng_pkt_2_1, scan_number, eng_data,              */
  /*                                       is_cp_hk_prior_section              */
  /*    OUTPUTS: eng_data                                                      */
  /*    RETURN:  None                                                          */
  /*                                                                           */
  /*****************************************************************************/

  update_eng_data_for_maj_cycle_n(major_cycle,eng_pkt_2_1,scan_number,eng_data,
                                   is_cp_hk_prior_section);

  return;
}
