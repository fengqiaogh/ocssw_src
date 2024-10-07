#include "L1A_prototype.h"
#include "EN_eng_data.h"
#include "hdfi.h"
#include "PGS_MODIS_35005.h"


void   create_eng_data_vdata_array (char             *eng_data_name,
                                    EN_VDATA_TYPE_t  *eng_data,
                                    uint16           *curr_eng_data_index,
                                    uint16           *curr_field_index)
/*
!C************************************************************************

!Description:  This function is a helper function to parse_eng_data_list().
               An element of an array is created in the eng_data structure. 
               A field called LAST_VALID_SCAN is added as its first field.

!Input Parameters:
               char        *eng_data_name        ** The name of the Vdata 
               					   to create              **

!Output Parameters: 
               None

!Input/Output Parameters:
               EN_VDATA_TYPE_t   *eng_data       ** The eng_data array 
               					   structure (a new Vdata 
               					   is added)              **

               uint16      *curr_field_index     ** The field index (at which
               					   to place the newly created 
               					   field) within the Vdata 
               					   currently being created 
               					   (this is reset and 
               					   incremented)           **

               uint16      *curr_eng_data_index  ** The eng_data index of the 
               					   Vdata currently being 
               					   created (this is 
               					   incremented)           **

Return Values: 
               None

Externally Defined:
               PGS_TRUE                     (PGS_SMF.h)
               PGS_FALSE                    (PGS_SMF.h)
               uint16                       (hdfi.h)
               EN_VDATA_TYPE_t              (EN_eng_data.h)
               PGSt_booean                  (PGS_TYPES.h)
               MODIS_E_NULL_POINTER         (PGS_MODIS_35005.h)

Called By:
               parse_eng_data_list

Routines Called:
               create_eng_data_vdata_array_field
               log_fmt_msg

!Revision History:
               Revision 2.0  1998/10/26  10:03 EST
               John Seaton/SAIC/GSC  (seaton@ltpmail.gsfc.nasangov)
               Added logic not to create LAST_VALID_SCAN field
               for Current/Prior S/C Ancillary Data Vdatas.

               Revision 1.0  1997/07/16  15:58 EDT
               David Catozzi/SAIC/GSC (cato@ltpmail.gsfc.nasa.gov)
               Original design.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               Vdatas are created in the order they appear in the 
               eng_data_list file so that the correct indexes correspond
               to the correct Vdatas (as expected in other routines).

!END************************************************************************
*/

{
  /***************************************************************************/
  /*                                                                         */
  /*              Declare and Initialize Local Variables                     */
  /*                                                                         */
  /***************************************************************************/

  static PGSt_boolean  first_vdata = PGS_TRUE;   /* first vdata (first time  */
                                                 /*  in routine) flag        */

  EN_VDATA_TYPE_t      *vdata;
  char                 *routine = "create_eng_data_vdata_array";

  /***************************************************************************/

  /***************************************************************************/
  /* IF NULL parameters passed in to procedure                               */
  /*  THEN                                                                   */
  /*    CALL log_fmt_msg to print NULL paramaters error message in LogStatus */
  /*      INPUTS: MODIS_E_NULL_POINTER, "create_eng_data_vdata_array", " "   */
  /*      OUTPUTS: None                                                      */
  /*      RETURNS: None                                                      */
  /*                                                                         */
  /*   return                                                                */
  /* ENDIF                                                                   */
  /***************************************************************************/

  if ((eng_data_name == NULL) ||
      (eng_data == NULL) ||
      (curr_eng_data_index == NULL) ||
      (curr_field_index == NULL)) {
    log_fmt_msg(MODIS_E_NULL_POINTER, routine, " ");
    return;
  }

  /***************************************************************************/
  /*                                                                         */
  /*  IF this is NOT the first vdata (i.e. if this isn't                     */
  /*     the first call to this function)                                    */
  /*  THEN                                                                   */
  /*     increment curr_eng_data_index                                       */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  if (first_vdata != PGS_TRUE)
     (*curr_eng_data_index)++;

  first_vdata = PGS_FALSE;


  /***************************************************************************/
  /*                                                                         */
  /*  set curr_field_index to 0                                              */
  /*  set eng_data[curr_eng_data_index].vdata_name to eng_data_name          */
  /*  set eng_data[curr_eng_data_index].num_fields to 0                      */
  /*                                                                         */
  /*  set field_name to "LAST_VALID_SCAN"                                    */
  /*  set num_bits to 0                                                      */
  /*  set start_bit_pos to 0                                                 */
  /*  IF eng_data_name NOT equal to Current S/C Ancillary Data Vdata AND     */
  /*   eng_data_name NOT equal to Prior S/C Ancillary Data Vdata             */
  /*  THEN                                                                   */
  /*    CALL create_eng_data_vdata_array_field to create the first field     */
  /*       'LAST_VALID_SCAN', which is mandatory for all Vdatas.             */
  /*      INPUTS:  field_name, num_bits, start_bit_pos, eng_data,            */
  /*               curr_eng_data_index, curr_field_index                     */
  /*      OUTPUTS: eng_data, curr_field_index                                */
  /*      RETURN:  None                                                      */
  /*                                                                         */
  /***************************************************************************/

  *curr_field_index = 0;

  vdata = &(eng_data[(int)*curr_eng_data_index]);
  strcpy(vdata->vdata_name, eng_data_name);
  vdata->num_fields = 0;
 
  if ((strcmp(eng_data_name, "Current S/C Ancillary Data") != 0) &&
      (strcmp(eng_data_name, "Prior S/C Ancillary Data") != 0))
         create_eng_data_vdata_array_field(EN_LAST_VALID_SCAN, 0, 0, 1, 23, 
                     eng_data,  *curr_eng_data_index, curr_field_index);

}
