#ifndef PD_PKT_DATA_H
#define PD_PKT_DATA_H

/*
!C-INC************************************************************************

!Description:  This include file contains information about the kinds of
               data that each type of packet contains, along with the sizes
               of arrays' dimensions.

               All definitions in this file will begin with "PD_".
               Definitions for Day and Night mode packets will begin with 
               "PD_DN_".
               Definitions for Engineering group 1 packet 1 packets will begin
               with "PD_E1P1".
               Definitions for Engineering group 1 packet 2 packets will begin
               with "PD_E1P2".
               Definitions for Engineering group 2 packet 1 packets will begin
               with "PD_E2P1".
               Definitions for Engineering group 2 packet 2 packets will begin
               with "PD_E2P2".


!Input Parameters:
               N/A

!Output Parameters:
               N/A

Return Values: 
               N/A

Externally Defined:  
               None 

Called By:
               N/A

Routines Called:
               N/A

!Revision History:
               Revision 2.0  1997/07/09  15:15
               Tom Johnson/SAIC/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Created include file from Version 1 include files
               modis_pkt_data_pos.h

               Revision 1.0  1996/06/19  16:30 EDT
               Keith Degnan/SAIC/GSC (keith.degnan@gsfc.nasa.gov)
               Created include module

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The ".h" below was specifically written for development in C.
               Any other language choice may require reworking of the ".h"
               before coding can begin.

!END**********************************************************************
*/

#define  PD_PKT_BUF_MAX                              645

#define  PD_NUM_BITS_IN_BYTE                           8 
#define  PD_NUM_BITS_IN_WORD                          16 
#define  PD_NUM_BYTES_IN_WORD                          2 
#define  PD_FIRST_BIT_IN_BYTE                          0 
#define  PD_NUM_BITS_IN_RAD_DATA                      12

#define  PD_PKT_CONTENTS_BYTE_OFFSET                  18 
#define  PD_NUM_ELMTS_IN_DATA_FIELD_NIGHT_PKT        171
#define  PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT          415


/********************************************************/
/*  The following contains information for the Day and  */
/*  Night mode packets.                                 */
/********************************************************/

#define  PD_DN_FIRST_250M_BAND                         1  
#define  PD_DN_LAST_250M_BAND                          2 
#define  PD_DN_BAND_RATIO_250M                         4  
#define  PD_DN_NUM_250M_BANDS                          2
#define  PD_DN_NUM_250M_DETECTORS                     40
#define  PD_DN_NUM_250M_DETECTORS_IN_IFOV              4  

#define  PD_DN_FIRST_500M_BAND                         3 
#define  PD_DN_LAST_500M_BAND                          7  
#define  PD_DN_BAND_RATIO_500M                         2 
#define  PD_DN_NUM_500M_BANDS                          5
#define  PD_DN_NUM_500M_DETECTORS                     20
#define  PD_DN_NUM_500M_DETECTORS_IN_IFOV              2  

#define  PD_DN_FIRST_1KM_DAY_BAND                      8 
#define  PD_DN_LAST_1KM_DAY_BAND                      21  
#define  PD_DN_FIRST_1KM_NIGHT_BAND                   22 
#define  PD_DN_LAST_1KM_NIGHT_BAND                    38  
#define  PD_DN_BAND_RATIO_1KM                          1  
#define  PD_DN_NUM_1KMDAY_BANDS                       14
#define  PD_DN_NUM_1KMDAY_DETECTORS                   10
#define  PD_DN_NUM_1KMNIGHT_BANDS                     17
#define  PD_DN_NUM_1KMNIGHT_DETECTORS                 10

#define  PD_DN_NUM_IFOVS_IN_DAY_PKT                    5  
#define  PD_DN_NUM_IFOVS_IN_NIGHT_PKT                 10  

#define  PD_DN_FIRST_IFOV_DAY_PKT_1                    1
#define  PD_DN_LAST_IFOV_DAY_PKT_1                     5
#define  PD_DN_FIRST_IFOV_DAY_PKT_2                    6
#define  PD_DN_LAST_IFOV_DAY_PKT_2                    10


/********************************************************/
/*  The following contains information for the          */
/*  Engineering group 1 packet 1                        */
/********************************************************/

#define  PD_E1P1_NUM_FPA_DCR_OFFSETS                 550 
#define  PD_E1P1_FPA_DCR_OFFSETS_BYTE_OFFSET          18


/********************************************************/
/*  The following contains information for the          */
/*  Engineering group 1 packet 2                        */
/********************************************************/

#define  PD_E1P2_NUM_EARTH_ENCODER_TIMES              78 
#define  PD_E1P2_EARTH_ENCODER_TIMES_BYTE_OFFSET      18
#define  PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS          40
#define  PD_E1P2_VIEW_SECTOR_DEFINITIONS_BYTE_OFFSET 174
#define  PD_E1P2_NUM_VIEW_SECTOR_ACTUALS              24
#define  PD_E1P2_VIEW_SECTOR_ACTUALS_BYTE_OFFSET     254
#define  PD_E1P2_NUM_SCI_ENG_BYTES                   212 
#define  PD_E1P2_SCI_ENG_BYTE_OFFSET                 302


/********************************************************/
/*  The following contains information for the          */
/*  Engineering group 2 packet 1                        */
/********************************************************/

#define  PD_E2P1_NUM_HK_TELEM_BYTES                  128 
#define  PD_E2P1_CURR_HK_BYTE_OFFSET                  18 
#define  PD_E2P1_PRIOR_HK_BYTE_OFFSET                 82 
#define  PD_E2P1_NUM_SC_ANCIL_WORDS                   64 
#define  PD_E2P1_SC_ANCIL_BYTE_OFFSET                146 
#define  PD_E2P1_NUM_PARAM_BYTES                      40 
#define  PD_E2P1_PARAM_BYTE_OFFSET                   274 


/********************************************************/
/*  The following contains information for the          */
/*  Engineering group 2 packet 2                        */
/********************************************************/

#define  PD_E2P2_NUM_PV_GAINS                        550 
#define  PD_E2P2_PV_GAINS_BYTE_OFFSET                 18


#endif  /* PD_PKT_DATA_H */

