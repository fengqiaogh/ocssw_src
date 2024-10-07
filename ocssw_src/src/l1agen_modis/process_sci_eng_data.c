#include "L1A_prototype.h"
#include "hdf.h"
#include "PGS_IO_L0.h"
#include "EN_eng_data.h"


void   process_sci_eng_data (EN_VDATA_TYPE_t   *eng_data,
                             PGSt_IO_L0_Packet *eng_pkt_1_2,
                             uint16             scan_number )

/*
!C******************************************************************************

!Description:  This function processes eng group 1 pkt #2.

!Input Parameters:
               uint16              scan_number   ** The scan number (counting  **
                                                 ** from 1) of the current scan**
                                                 ** within the current granule **

               PGSt_IO_LO_Packet   *eng_pkt_1_2  ** eng group 1's packet #2    **
                                                 ** for the current scan       **

!Output Parameters:
               None

!Input/Output Parameters:
               EN_VDATA_TYPE_t     *eng_data     ** The Vdata array structure  **

Return Values:
               None

Externally Defined:
               EN_VDATA_TYPE_t            (EN_eng_data.h)
               EN_ENG_VDATA_START_INDEX   (EN_eng_data.h)
               EN_NUM_VDATAS              (EN_eng_data.h)
               PGSt_IO_LO_Packet          (PGS_IO.h)

Called By:
               process_eng_packet

Routines Called:
               update_eng_data

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
               (c.f. CDRL Table 30-5D and the Vdata_list file)

!END************************************************************************
*/
{
  /***************************************************************************/
  /*                                                                         */
  /*               Define and Initialize Local Variables                     */
  /*                                                                         */
  /***************************************************************************/

  int          index = EN_ENG_VDATA_START_INDEX;
  int          is_cp_hk_prior = FALSE;


  /***************************************************************************/
  /*                                                                         */
  /*  set 'index' to EN_ENG_VDATA_START_INDEX  (done during declaration)     */
  /*                                                                         */
  /*  set is_cp_hk_prior to False  (done during declaration)                 */
  /*                                                                         */
  /***************************************************************************/

  /***************************************************************************/
  /*                                                                         */
  /*  DO_WHILE ( index < EN_NUM_VDATAS )                                     */
  /*                                                                         */
  /*    CALL update_eng_data to update the eng_data (vdata_array) structure  */
  /*         with the eng data values in the eng pkt                         */
  /*      INPUTS:  index, eng_pkt_1_2, scan_number, eng_data, is_cp_hk_prior */
  /*      OUTPUTS: eng_data                                                  */
  /*      RETURN:  None                                                      */
  /*                                                                         */
  /*    increment index                                                      */
  /*                                                                         */
  /*  END_WHILE                                                              */
  /*                                                                         */
  /***************************************************************************/

  while (index < EN_NUM_VDATAS)
    {
      update_eng_data(index,eng_pkt_1_2,scan_number,eng_data,is_cp_hk_prior);
      index++;
    }

  return;

}
