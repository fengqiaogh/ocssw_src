#include "L1A_prototype.h"
#include "hdf.h"
#include "PGS_IO_L0.h"
#include "EN_eng_data.h"
#include "PH_pkt_hdr.h"

void  update_eng_data_for_maj_cycle_n ( uint16             major_cycle,
                                        PGSt_IO_L0_Packet *eng_pkt_2_1,
                                        uint16             scan_number,
                                        EN_VDATA_TYPE_t   *eng_data,
                                        int                is_cp_hk_prior_section )


/*
!C***********************************************************************************

!Description:  This function processes the major cycle dependent parts of group 2 
               pkt #1's Current/Prior HK Tlmy section. The appropriate eng_data Vdatas 
               are updated for this scan using the data in the engineering packet's 
               Current/Prior HK Tlmy section.
               (c.f. CDRL Table 30-6A and the Vdata_list file).

               The following table shows the Major Cycle to vdata_array index
               relationship.
 
               ----------------------------------------------------------------------
                If Major Cycle is:      | then update Vdatas at vdata_array indices:
               ----------------------------------------------------------------------
                Major Cycle (MOD 8) =   |   
                    0                   |    3
                    1                   |    4
                    2                   |    5  
                    3                   |    6, 7, 8
                    4                   |    9, 10
                    5                   |    11, 12
                    6                   |    13
                    7                   |    14
               ----------------------------------------------------------------------
                  0..32                 |    15..47  
               ----------------------------------------------------------------------  
                 (c.f. CDRL Table 30-6A, CDRL Table 20-4, and the Vdata_list) 



!Input Parameters:
               uint16              major_cycle   ** The major cycle for which 
                                                    to update the eng_data **

               PGSt_IO_LO_Packet  *eng_pkt_2_1   ** eng grp 2 packet #1 for the 
                                                    current scan **

               uint16              scan_number   ** The scan number (counting 
                                                    from 1) of the current scan
                                                    within the current granule **

               boolean             is_cp_hk_prior_section  ** whether or not 
                                                              this is the "prior"
                                                              part of eng pkt
                                                              2-1's current/prior
                                                              section **
!Output Parameters:
               None

!Input/Output Parameters:
               EN_VDATA_TYPE_t    *eng_data     ** The Vdata array structure **

Return Values:
               None

Externally Defined:
               EN_VDATA_TYPE_t                 (EN_eng_data.h)
               PGSt_IO_LO_Packet               (PGS_IO.h)
               EN_NUM_MAJ_CYCLES               (EN_eng_data.h)

Called By:
               process_cp_hk_tlmy

Routines Called:
               update_eng_data

!Revision History:
               Revision 2.0  2001/01/04
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Simplified code. Removed call to is_in_range

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
               None

!END************************************************************************
*/
{
  int     index;
  uint16  cycle_x_of_7;

  /*******************************************************************************/
  /*  IF major_cycle is less than or equal to EN_NUM_MAJ_CYCLES (32)             */
  /*  THEN                                                                       */
  /*    Set index to major_cycle + 15                                            */
  /*    CALL update_eng_data to update the appropriate                           */
  /*          "Telemetry Major Cycle X of 63" Vdata in the eng_data structure    */
  /*      INPUTS:  index,eng_pkt_2_1,scan_number,eng_data,is_cp_hk_prior_section */
  /*      OUTPUT:  eng_data                                                      */
  /*      RETURN:  None                                                          */
  /*  ENDIF                                                                      */
  /*******************************************************************************/

  if (major_cycle <= EN_NUM_MAJ_CYCLES)
    {
      index = major_cycle + 15;
      update_eng_data(index,eng_pkt_2_1,scan_number,eng_data,is_cp_hk_prior_section);
    }


  /****************************************************************************/
  /*  SWITCH ( major_cycle MODULO 8 )                                         */
  /*                                                                          */
  /*    CASE 0:                                                               */
  /*      Set index to 3                                                      */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 0                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*      BREAK                                                               */
  /*                                                                          */
  /*    CASE 1:                                                               */
  /*      Set index to 4                                                      */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 1                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*      BREAK                                                               */
  /*                                                                          */
  /*    CASE 2:                                                               */
  /*      Set index to 5                                                      */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 2                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*      BREAK                                                               */
  /*                                                                          */
  /*    CASE 3:                                                               */
  /*      FOR ( index = 6 to 8 )                                              */
  /*         CALL update_eng_data to update the the appropriate "Telemetry    */
  /*            Major Cycle X of 7" Vdata in the eng_data structure with the  */
  /*            eng data values in the eng pkt for major cycle 3              */
  /*           INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,            */
  /*                    is_cp_hk_prior_section                                */
  /*           OUTPUT:  eng_data                                              */
  /*           RETURN:  None                                                  */
  /*      ENDFOR                                                              */
  /*      BREAK                                                               */
  /*                                                                          */
  /*    CASE 4:                                                               */
  /*      Set index to 9                                                      */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 4                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*                                                                          */
  /*      Set index to 10                                                     */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 4                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*                                                                          */
  /*      BREAK                                                               */
  /*                                                                          */
  /*    CASE 5:                                                               */
  /*      Set index to 11                                                     */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 5                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*                                                                          */
  /*      Set index to 12                                                     */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 5                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*                                                                          */
  /*      BREAK                                                               */
  /*                                                                          */
  /*    CASE 6:                                                               */
  /*      Set index to 13                                                     */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 6                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*                                                                          */
  /*      BREAK                                                               */
  /*                                                                          */
  /*    CASE 7:                                                               */
  /*      Set index to 14                                                     */
  /*      CALL update_eng_data to update the the appropriate "Telemetry Major */
  /*         Cycle X of 7" Vdata in the eng_data structure with the eng data  */
  /*         values in the eng pkt for major cycle 7                          */
  /*        INPUTS:  index, eng_pkt_2_1, scan_number, eng_data,               */
  /*                 is_cp_hk_prior_section                                   */
  /*        OUTPUT:  eng_data                                                 */
  /*        RETURN:  None                                                     */
  /*                                                                          */
  /*      BREAK                                                               */
  /*                                                                          */
  /*  END_SWITCH                                                              */
  /****************************************************************************/

  cycle_x_of_7 = major_cycle % 8;
  switch(cycle_x_of_7)
    {
    case 0:
      index = 3;
      update_eng_data(index,eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);
      break;

    case 1:
      index = 4;
      update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);
      break;

    case 2:
      index = 5;
      update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);
      break;

    case 3:
      for (index=6; index<=8; index++)
        update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                        is_cp_hk_prior_section);
      break;

    case 4:
      index = 9;
      update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);

      index = 10;
      update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);
      break;

    case 5:
      index = 11;
      update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);

      index = 12;
      update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);
      break;

    case 6:
      index = 13;
      update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);
      break;

    case 7:
      index = 14;
      update_eng_data(index, eng_pkt_2_1,scan_number,eng_data,
                      is_cp_hk_prior_section);
      break;
    }

  return;
}
