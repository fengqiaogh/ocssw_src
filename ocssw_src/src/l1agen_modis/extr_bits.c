#include "PGS_MODIS_35005.h"
#include "PGS_SMF.h"
#include "hdfi.h"
#include "L1A_prototype.h"

uint32   extr_bits (uint8  *a, 
                    int    start_bit,
                    int    start_byte,
                    int    num_bits )

/*
!C************************************************************************

!Description:  This function will extract up to 32 bits from a uint8
               (byte) array. This routine assumes that bits are numbered 
               from MSB to LSB in each byte, and that bit streams connect 
               from bit 7 in byte b to bit 0 in byte b + 1.

               The routine will not work for streams larger than 32 bits.
               It won't crash in that case; it just won't return a valid
               value.

!Input Parameters:      
               int8 *a              ** array from which bits are to be **
                                    ** extracted                       **
               int start_byte       ** byte within array that contains **
                                    ** the first bit to be extracted   **
               int start_bit        ** bit within start_byte at witch  **
                                    ** extraction will start           **
               int num_bits         ** number of bits to be extracted  **

!Output Parameters: 
               None

Return Values: 
               Bits extracted from *a, or 0 (zero) if extract cannot
               be performed

Externally Defined:   
               PGSt_SMF_status            (PGS_SMF.h) 
               MODIS_E_TOO_MANY_BITS      (PGS_MODIS_35005.h)
               MODIS_S_SUCCESS            (PGS_MODIS_35005.h)
               uint8                      (hdfi.h)
               uint32                     (hdfi.h)
               pkt_num                    (ext var to hold packet number)

Called By:
               unpack_primary_header
               unpack_secondary_header
               unpack_MODIS_header
               unpack_packet_contents

Routines Called:
               log_fmt_msg

!Revision History:     
               Revision 1.0  1997/09/15  14:35
               Tom Johnson/GSC    (johnson@ltpmail.gsfc.nasa.gov)
               Developed the code for the module based on existing 
               code from beta version and version 1

!Team-unique Header: 
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               None

!END***********************************************************************
*/

{
  /***************************************************************************/
  /*                                                                         */
  /*              Declare and Initialize Local Variables                     */
  /*                                                                         */
  /***************************************************************************/

  char             msg[300];       /* amplifying message                     */
  int              i;              /* loop counter                           */
  uint32           ret;            /* Variable holding extracted bits        */
  int              end_bit;        /* Calculated last bit to be extracted    */
  int              end_byte;       /* Calculated last byte from which bits   */
                                   /*   will be extracted                    */

  char             *routine = "extr_bits";

  /***************************************************************************/


  /***************************************************************************/
  /*                                                                         */
  /*  Set routine to "extr_bits" (done during declaration)                   */
  /*                                                                         */
  /***************************************************************************/

  /***************************************************************************/
  /*                                                                         */
  /*  Set ret to 0                                                           */
  /*  Set end_bit to 0                                                       */
  /*  Set end_byte to 0                                                      */
  /*                                                                         */
  /***************************************************************************/

  ret = 0;
  end_bit = 0;
  end_byte = 0;

  /***************************************************************************/
  /*                                                                         */
  /*  IF num_bits is greater than 32                                         */
  /*  THEN                                                                   */
  /*     Set ret to 0                                                        */
  /*     Set msg to "Number of bits to be extracted: <number>                */
  /*     CALL log_fmt_msg to report that extraction can't be done            */
  /*       INPUT:  MODIS_E_TOO_MANY_BITS, routine, msg                       */
  /*       OUTPUT: None                                                      */
  /*       RETURN: None                                                      */
  /*                                                                         */
  /***************************************************************************/

  if (num_bits > 32) 
    { 
        /** Not able to extract that many bits, so just return error msg**/
     ret = 0;
     sprintf (msg, "Number of bits to be extracted: %d", num_bits);
     log_fmt_msg (MODIS_E_TOO_MANY_BITS, routine, msg);
    }


  /***************************************************************************/
  /*                                                                         */
  /*  ELSE                                                                   */
  /*     Set end_bit to start_bit + num_bits - 1                             */
  /*     Set end_byte to start_byte + end_bit/8                              */
  /*     Set end_bit to end_bit%8                                            */
  /*                                                                         */
  /*     Set ret to head_buff[start_byte] & (Oxff >> start_bit)              */
  /*                                                                         */
  /*     FOR i from start_byte + 1 to end_byte - 1                           */
  /*        Set ret to (ret << 8) | head_buff[i]                             */
  /*     ENDFOR                                                              */
  /*                                                                         */
  /*     IF start_byte is not equal to end_byte                              */
  /*     THEN                                                                */
  /*        Set ret to ((head_buff[end_byte] >> (7 - end_bit)) & Oxff) |     */
  /*           (ret << (end_bit + 1))                                        */
  /*     ELSE                                                                */
  /*       Set ret to ret >> (7 - end_bit)                                   */
  /*     ENDIF                                                               */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  else 
    { /** Calculate the end byte and bit **/
     end_bit = ( start_bit + num_bits - 1 );
     end_byte = ( start_byte + (end_bit / 8) );
     end_bit = ( end_bit % 8 );

              /*******************************************************/
              /**          Build the extracted number               **/
              /*******************************************************/
              /**                 Start byte                        **/
              /**                 ----------                        **/
              /**  shift right the hexadecimal constant "Oxff"      **/
	      /**   (start_bit) times and AND it to a[start_byte],  **/
              /**   then  ret = the result                          **/
              /*******************************************************/
     ret = a[start_byte] & (0xff >> start_bit); 

              /*******************************************************/
              /**                Middle byte                        **/
              /**                -----------                        **/
              /**  shift left ret 8 times and OR it with a[i] and   **/
              /**   ret = the result                                **/
              /*******************************************************/
     for (i = (start_byte + 1); i < end_byte; i++)
       ret = (ret << 8) | a[i];

              /*******************************************************/
              /**                  End byte                         **/
              /**                  --------                         **/
     if (start_byte != end_byte) 
              /**  shift right a[end_byte] (7 - end_bit) times and  **/
              /**   AND it with the hexadecimal constant "0xff",    **/
              /**   then OR this with the result of ret shifted     **/
              /**   left (end_bit + 1) times, then ret = the result **/
        ret = ((a[end_byte] >> (7 - end_bit)) & 0xff) | 
               (ret << (end_bit + 1));
     else
              /**  shift right ret (7 - end_bit) times andi,        **/
              /**   ret = the result                                **/
        ret = ret >> (7 - end_bit);
              /*******************************************************/
    } 


  /***************************************************************************/
  /*                                                                         */
  /*  RETURN ret                                                             */
  /*                                                                         */
  /***************************************************************************/

  return (ret);

} /* End extr_bits */
