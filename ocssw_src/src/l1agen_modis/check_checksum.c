#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "hdfi.h"
#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "L1A_prototype.h"


PGSt_SMF_status  check_checksum (PH_PACKET_HEADER_t  pkt_header,
                                 uint16              *pkt_contents)


/*
!C************************************************************************

!Description:  This function will evaluate the current packet's check
               sum by calculating the check sum and comparing it to
               the checksum value in the packet.  If there is any 
               discrepancy, the function will return an error message.

!Input Parameters:
               uint16            *pkt_contents ** The current packet, with **
                                               ** the 12-bit words already **
                                               ** unpacked into 16-bit     **
                                               ** integers                 **

               PH_PACKET_HEADER_t  pkt_header  ** Buffer that contains pkt **
                                               ** header info of the pckt  **

!Output Parameters: 
               None

Return Values: 
               MODIS_S_SUCCESS                        (PGS_MODIS_35005.h)
               MODIS_E_CHECKSUM_NOT_VALID             (PGS_MODIS_35005.h)

Externally Defined:  
               PGSt_SMF_status                        (PGS_SMF.h) 
               uint16                                 (hdfi.h)
               PH_PACKET_HEADER_t                     (PH_pkt_hdt.h)
               PH_NUM_12BIT_WORDS_IN_HEADER           (PH_pkt_hdt.h)
               PH_PRI_LONG_PKT_LENGTH                 (PH_pkt_hdt.h)
               PH_PRI_SHORT_PKT_LENGTH                (PH_pkt_hdt.h)
               PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT     (PD_pkt_data.h)

Called by:
               load_eng_data
               process_a_packet
               process_a_scan

Routines Called:
               log_fmt_msg 

!Revision History:
               Revision 2.0  1997/08/14  10:14 EDT 
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Restructered module per version2 development.                

               Revision 1.0  1997/06/18  16:40 EDT 
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits
               None

!Design Notes: The Checksum Algorithm is provided by Santa Barbara
               Remote Sensing (SBRS), and can be found in (TBS).

               The CODE below was developed in C language.    

               The checksum Algorithm was produce by Joe Aucther for Sci
               Test Pkts (6/06/96). It states as follows:
               
               Packet checksums
               ----------------
               FYI, The checksum is computed as follows:
               1: Within each packet, the 12-bit science data (or test data) 
                  is summed into a 16-bit register, with overflows ignored.  
               2: After all the data is summed, the 16-bit result is right 
                  shifted four bits, with zeros inserted into the leading 
                  four bits.  
               3: The 12 bits left in the lower end of the word becomes the 
                  checksum.

!END**************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*               Define and Initialize Local Variables                    */
  /*                                                                        */
  /**************************************************************************/

  char             *routine;        /* Variable to hold routine name            */
  PGSt_SMF_status  returnStatus;    /* SMF-style message returned by function   */
  uint16           calc_checksum;   /* The checksum calculated by this routine  */
  uint16           pkt_checksum;    /* The checksum contained within the packet */
  int              begin_sci_data;  /* Begining of data field in current packet */
  int              end_sci_data;    /* End of data field in current packet      */
  int              i;               /* loop variable                            */
  int              num_16bit_words; /* number of 16-bit words in a packet       */


  routine         = "check_checksum";
  returnStatus    = MODIS_S_SUCCESS;
  calc_checksum   = 0;
  pkt_checksum    = 0;
  begin_sci_data  = PH_NUM_12BIT_WORDS_IN_HEADER;
  end_sci_data    = 0;
  num_16bit_words = 0;


  /**************************************************************************/
  /*                                                                        */
  /* IF pkt_header.pkt_length is equal to PH_PRI_LONG_PKT_LENGTH            */
  /* THEN                                                                   */
  /*   set num_16bit_words equal to PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT      */
  /* ELSE                                                                   */
  /*   set num_16bit_words equal to PD_NUM_ELMTS_IN_DATA_FIELD_NIGHT_PKT    */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* set end_science_data equal to PH_NUM_12BIT_WORDS_IN_HEADER +           */
  /*                                                  num_16bit_words       */
  /*                                                                        */
  /**************************************************************************/

  if (pkt_header.pkt_length == PH_PRI_LONG_PKT_LENGTH) 
    num_16bit_words = PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT;

  else 
    num_16bit_words = PD_NUM_ELMTS_IN_DATA_FIELD_NIGHT_PKT;

  end_sci_data = PH_NUM_12BIT_WORDS_IN_HEADER + num_16bit_words;



  /**************************************************************************/
  /*                                                                        */
  /* Calculate the checksum with the algorithm provided by Joe Aucther      */
  /* for Science Test Packets. Suggested loop: FOR LOOP                     */
  /*                                                                        */
  /* FOR index equal to begin_of_science_data to end_of_science_data        */
  /*   set calc_checksum equal to calc_checksum + pkt_contents[index]       */
  /* END_FOR                                                                */
  /*                                                                        */
  /* shift right calc_checksum 4 times                                      */
  /*                                                                        */
  /**************************************************************************/

  for (i = begin_sci_data; i < end_sci_data; i++)
    calc_checksum = calc_checksum + pkt_contents[i];

  calc_checksum >>= 4;


  /**************************************************************************/ 
  /*                                                                        */
  /* Extract the actual checksum value from the packet and compare it       */
  /* with the calculated value of checksum                                  */
  /*                                                                        */
  /* set pkt_checksum equal to                                              */
  /*        pkt_contents[PH_NUM_12BIT_WORDS_IN_HEADER + num_16bit_words]    */
  /*                                                                        */
  /* IF calc_checksum is not equal to pkt_checksum                          */
  /* THEN                                                                   */
  /*   set returnStatus equal to MODIS_E_CHECKSUM_NOT_VALID                 */
  /*   set routine equal to "check_checksum"                                */
  /*   set msg equal to "checksum value of packet does not match            */
  /*       calculated checksum"                                             */
  /*   CALL log_fmt_msg to report descrepancy in checksum values            */
  /*     INPUTS: returnStatus, routine, msg                                 */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

  pkt_checksum = pkt_contents[PH_NUM_12BIT_WORDS_IN_HEADER + num_16bit_words];

  if (calc_checksum != pkt_checksum) {
    returnStatus = MODIS_E_CHECKSUM_NOT_VALID;
    log_fmt_msg (MODIS_E_CHECKSUM_NOT_VALID, routine, "%d does not match pkt checksum",
                  calc_checksum, pkt_checksum);
  }

  return returnStatus;

 }  /* End of routine check_checksum */
