#include "PGS_IO.h"
#include "hdfi.h"
#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "L1A_prototype.h"

void  unpack_packet_contents (PGSt_IO_L0_Packet   *pkt,
                              PH_PACKET_HEADER_t  *pkt_header,
                              uint16              *pkt_contents)

/*
!C************************************************************************

!Description:  This function extracts the contents contained in the data
               field in the MODIS packet. The contents is extracted from 
               the MODIS packet in 12 bit words and converted it into 16 
               bit words and placed in an array.

!Input Parameters:
               PGSt_IO_L0_Packet   *pkt    ** The current MODIS packet to  **
                                           ** be unpacked                  **

               PH_PACKET_HEADER_t  *pkt_header  ** The current MODIS       **
                                                ** packet's packet header  **

!Output Parameters:     
               uint16              *pkt_contents  ** The unpacked data     **
	                                          ** structure; every      **
                                                  ** 12-bit set is unpacked**
                                                  ** into 16 bits.         **

Return Values: 
               None

Externally Defined:
               PD_NUM_ELMTS_IN_DATA_FIELD_NIGHT_PKT      (PD_pkt_data.h)
               PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT        (PD_pkt_data.h)
               PH_NUM_12BIT_WORDS_IN_HEADER              (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_NIGHT_GROUP               (PH_pkt_hdr.h)
               PH_PACKET_HEADER_t                        (PH_pkt_hdr.h)
               PGSt_IO_L0_Packet                         (PGS_IO.h) 
               uint16                                    (hdfi.h)

Called by:
               load_eng_data
               process_a_packet
               process_a_scan

Routines Called:
               extr_bits                          

!Revision History:
               Revision 2.0  1997/08/14  13:45 EDT
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Recreating module per Version 2 development.

               Revision 1.0  1997/06/18  16:40 EDT 
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The CODE below was developed in C language.    

!END***********************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*              Define and Initialize Local Variables                     */
  /*                                                                        */
  /**************************************************************************/

   int       num_12bit_words;     /* number of unpacked words that will be  */
                                  /*   created from the packet              */
   int       i;                   /* loop variable                          */

   PGSt_IO_L0_Packet  a[3], b;	  /* temporary variables to hold bytes	    */
   uint16    nibble1, nibble2;    /* temporary variables to hold words      */ 
   uint16    f1, f2;		  /* temporary variables to hold words      */
   int 	     cnt;		  /* counter for the 12-bit fields	    */
   int 	     offset;	          /* offset inside the packet		    */ 

   num_12bit_words = 0;

  /**************************************************************************/
  /*                                                                        */
  /* Determine the number of 12 bit words in the data field based on the    */
  /* packet type. All non-night packets have the longer number of elements  */
  /*                                                                        */
  /* IF PH_PACKET_HEADER_t.pkt_type equals PH_SEC_PKT_TYPE_NIGHT_GROUP      */
  /* THEN                                                                   */
  /*   Set num_12bit_words equal to PH_NUM_12BIT_WORDS_IN_HEADER +          */
  /*                        PD_NUM_ELMTS_IN_DATA_FIELD_NIGHT_PKT + 1        */
  /* ELSE                                                                   */
  /*   Set num_12bit_words equal to PH_NUM_12BIT_WORDS_IN_HEADER +          */
  /*                          PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT + 1        */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

  if (pkt_header->pkt_type == PH_SEC_PKT_TYPE_NIGHT_GROUP)
    num_12bit_words = PH_NUM_12BIT_WORDS_IN_HEADER + 
              PD_NUM_ELMTS_IN_DATA_FIELD_NIGHT_PKT + 1;
  else
    num_12bit_words = PH_NUM_12BIT_WORDS_IN_HEADER +
                PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT + 1;


  /**************************************************************************/
  /*                                                                        */
  /* FOR i equal to 0 upto (num_12bit_words/2)                              */
  /*                                                                        */
  /*    read 3 bytes from pkt and place them into the array a[3]            */
  /*                                                                        */
  /**************************************************************************/ 

  cnt = num_12bit_words/2;
  for (i=0; i<cnt; i++) 
    {
     offset = i*3;

     a[0] = *(pkt+offset);
     a[1] = *(pkt+offset+1);
     a[2] = *(pkt+offset+2);
   

  /**************************************************************************/
  /*                                                                        */
  /*    extract low 4 bits from the middle byte                             */
  /*                                                                        */
  /*    extract high 4 bits from the middle byte                            */
  /*                                                                        */
  /**************************************************************************/
   
     b = a[1];
     nibble1 = b & 0x0f;
     nibble2 = (b & 0xf0) >>4;


  /**************************************************************************/
  /*                                                                        */
  /*    put low 4 bits with third byte together                             */
  /*                                                                        */
  /*    put high 4 bits with first byte together                            */
  /*                                                                        */
  /**************************************************************************/

     f2 = a[2] | ((uint16)nibble1)<<8;	
     f1 = nibble2 | ((uint16)a[0])<<4;


  /**************************************************************************/
  /*                                                                        */
  /*    place padded values into pkt_contents[]                             */
  /*                                                                        */
  /* ENDFOR                                                                 */
  /*                                                                        */
  /**************************************************************************/

     pkt_contents[i*2] = f1;
     pkt_contents[i*2+1] = f2;
    }

 }  /* End of routine unpack_packet_contents */
