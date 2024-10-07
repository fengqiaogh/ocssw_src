#include "L1A_prototype.h"
#include "EN_eng_data.h"
#include "hdfi.h"
#include "PGS_MODIS_35005.h"


void  create_eng_data_vdata_array_field ( char     *field_name,
                                          uint16   num_bits,
                                          uint16   start_bit_pos,
                                          uint16   order,
                                          uint16   type,
                                          EN_VDATA_TYPE_t *eng_data,
                                          uint16   curr_eng_data_index,
                                          uint16   *curr_field_index )

/*
!C************************************************************************

!Description:  This function is a helper function to parse_eng_data_list().
               A field is created in the current structure.  Its value is 
               initialized to EN_INITIAL_FIELD_VALUE (or 
               EN_INITIAL_LAST_VALID_SCAN_VALUE if it is the LAST_VALID_SCAN
               field). 

!Input Parameters:
               char          *field_name         ** The name of the Vdata 
               					    field to create        **

               uint16        num_bits            ** The field's bit length **

               uint16        start_bit_pos       ** The field's start bit 
               					    position within the  
                                        	    eng packet             **

               uint16        curr_eng_data_index ** The eng_data index of 
               					    the Vdata currently being
               					    created                **

!Output Parameters:
               None
      
!Input/Output Parameters:
               EN_VDATA_TYPE_t  *eng_data        ** The eng_data array 
               					    structure (a new field 
                                        	    is added)              **

               uint16        *curr_field_index   ** The field index (at which
               					    to place the newly created
               					    field) within the Vdata 
               					    currently being created 
               					    (this is incremented)  **
Return Values: 
               None

Externally Defined:
               EN_VDATA_TYPE_t                       (EN_eng_data.h)
               EN_FIELD_TYPE_t                       (EN_eng_data.h)
               EN_INITIAL_LAST_VALID_SCAN_VALUE      (EN_eng_data.h)
               EN_INITIAL_FIELD_VALUE                (EN_eng_data.h)

Called By:
               parse_eng_data_list
               create_eng_data_vdata_array

Routines Called:
               log_fmt_msg

!Revision History:
               Revision 2.0 1998/10/26   10:16 EST
               John Seaton/SAIC/GSC  (seaton@ltpmail.gsfc.nasa.gov)   
               Added order and type fields to handle Vdatas with
               different orders and types.

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
               Vdata fields are created in the order they appear in the 
               eng_data_list file so that the correct fields are associated
               with the correct vdatas.

!END************************************************************************
*/


{
  /***************************************************************************/
  /*                                                                         */
  /*              Declare and Initialize Local Variables                     */
  /*                                                                         */
  /***************************************************************************/

  EN_VDATA_TYPE_t  *vdata;
  EN_FIELD_TYPE_t  *field;
  char             *routine = "create_eng_data_vdata_array_field";

  /***************************************************************************/
  /* Check for NULL input paramaters                                         */
  /***************************************************************************/
  if ((field_name == NULL) ||
      (eng_data == NULL) ||
      (curr_field_index == NULL)) {
    log_fmt_msg(MODIS_E_NULL_POINTER, routine, " ");
    return;                                         
  }

  vdata = &(eng_data[curr_eng_data_index]);
  field = &(vdata->field[*curr_field_index]);

  /***************************************************************************/
  /* Set Vdata field members                                                 */
  /***************************************************************************/
  strcpy(field->field_name, field_name);
  field->num_bits = num_bits;
  field->start_bit_pos = start_bit_pos;
  field->order = order;
  field->type = type;


  /***************************************************************************/
  /* Initialize field structure value to 0 ot 65535 if LAST_VALID_SCAN field */
  /***************************************************************************/
  if (strcmp(field_name, EN_LAST_VALID_SCAN) == 0)
    field->value = EN_INITIAL_LAST_VALID_SCAN_VALUE; 
  else
    field->value = EN_INITIAL_FIELD_VALUE;


  /***************************************************************************/
  /*  increment curr_field_index and eng_data[curr_eng_data_index].num_fields*/
  /***************************************************************************/
  (*curr_field_index)++;
  vdata->num_fields++;

}
