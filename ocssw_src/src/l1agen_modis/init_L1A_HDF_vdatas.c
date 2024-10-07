#include "L1A_prototype.h"
#include "VU_vdata_utility.h"
#include "PGS_MODIS_35005.h"
#include "mapi.h"
#include "mapiL1A.h"
#include "hdf.h"
#include "EN_eng_data.h"
#include "SC_scan.h"
#include "PGS_SMF.h"
#include "hntdefs.h"

PGSt_SMF_status   init_L1A_HDF_vdatas (EN_VDATA_TYPE_t  *eng_data,
                                       MODFILE          *L1A_file)

/*
!C*****************************************************************************

!Description:  This function creates in the L1A file all the Vdatas that  
               will be written to by the L1A processing.

!Input Parameters:
               EN_VDATA_TYPE_t    *eng_data       ** The Vdata array 
               					       structure  **
               MODFILE            *L1A_file       ** The address of MODFILE 
               					       structure  ** 

!Output Parameters:
               None

Return Values:
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               FAIL                               (HDF)

Externally Defined:  
               EN_VDATA_TYPE_t                    (EN_eng_data.h)
               EN_NUM_VDATAS                      (EN_eng_data.h)
               EN_MAX_FIELDS_PER_VDATA            (EN_eng_data.h)
               EN_MAX_DATA_TYPE_STRING_LENGTH     (EN_eng_data.h)
               EN_MAX_FIELD_NAME_LENGTH           (EN_eng_data.h)
	       EN_SC_ANCILLARY_VDATA_START	  (EN_eng_data.h)
	       EN_SC_ANCILLARY_VDATA_END	  (EN_eng_data.h)
               PGSt_SMF_status                    (PGS_SMF.h)
               MODIS_E_START_VDATA_ACCESS_TO_FILE (PGS_MODIS_35005.h)
               MODIS_E_CREATE_VDATA               (PGS_MODIS_35005.h)
               global_H_ID                        (level1a)
               DFNT_INT8                          (hntdefs.h)
               DFNT_UINT8                         (hntdefs.h)
               DFNT_INT16                         (hntdefs.h)
               DFNT_UINT16                        (hntdefs.h)
               DFNT_INT32                         (hntdefs.h)
               DFNT_UINT32                        (hntdefs.h)
	       SC_FILL_VALUE			  (SC_scan.h)
               VSFIELDMAX                         
               VU_MAX_NAME_LENGTH                 (VU_vdata_utility.h)
               VU_MAX_DATA_TYPE_STRING_LENGTH     (VU_vdata_utility.h)
               uint16                             (hdfi.h)
               MODIS_E_INVALID_VDATA_TYPE         (PGS_MODIS_35005.h)
               I8                                 (mapi.h)
               UI8                                (mapi.h)
               I16                                (mapi.h)
               UI16                               (mapi.h)
               I32                                (mapi.h)
               UI32                               (mapi.h)

Called By:
               create_L1A_granule

Routines Called:
               create_Vdata
               log_fmt_msg

!Revision History:
               Revision 2.0  1998/10/26  10:09 EST
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Added logic to allow a Vdata to have different size
               data types and different size orders.

               revision 1.0 1997/10/02  17:30:00
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
  /**************************************************************************/
  /*                      Define Global Variables                           */
  /**************************************************************************/

  extern   int32   global_H_ID;                       /* stores HDF file ID */


  /**************************************************************************/
  /*                      Define Local Variables                            */
  /**************************************************************************/

  char              *routine = "init_L1A_HDF_vdatas";
  char              msg[300];
  PGSt_SMF_status   returnStatus = MODIS_S_SUCCESS;
  char              eng_field_names[VSFIELDMAX][VU_MAX_NAME_LENGTH];
  char              eng_data_types[VSFIELDMAX][VU_MAX_DATA_TYPE_STRING_LENGTH];
  uint16            eng_field_order[VSFIELDMAX];
  char              Vdata_name[EN_MAX_VDATA_NAME_LENGTH];
  int               current_index=0;
  int               i;

  /**************************************************************************/
  /* Check for NULL input paramaters                                        */
  /**************************************************************************/

  if ((L1A_file == NULL) || (eng_data == NULL)) {
    log_fmt_msg(MODIS_E_NULL_POINTER, routine, " ");
    return FAIL;
  }

  global_H_ID = (int32)L1A_file->hdf_id;


  /**************************************************************************/
  /* Loop through all vdatas from eng_data structure                        */
  /**************************************************************************/
  while ((current_index < EN_NUM_VDATAS) && 
         (returnStatus == MODIS_S_SUCCESS))
    {
  /**************************************************************************/
  /* Loop through all fields for each Vdata                                 */
  /**************************************************************************/
      for (i=0; i<eng_data[current_index].num_fields; i++)
        {
          strcpy(eng_field_names[i], eng_data[current_index].field[i].field_name);
          switch (eng_data[current_index].field[i].type) {
             case DFNT_INT8   : strcpy(eng_data_types[i], I8); break;
             case DFNT_UINT8  : strcpy(eng_data_types[i], UI8); break;
             case DFNT_INT16  : strcpy(eng_data_types[i], I16); break;
             case DFNT_UINT16 : strcpy(eng_data_types[i], UI16); break;
             case DFNT_INT32  : strcpy(eng_data_types[i], I32); break;
             case DFNT_UINT32 : strcpy(eng_data_types[i], UI32); break;
             default          : sprintf(msg, "\nVdata Name: %s Vdata Field: %s Vdata type: %d", 
                                             eng_data[current_index].vdata_name,
                                             eng_data[current_index].field[i].field_name,
                                             eng_data[current_index].field[i].type);
                                log_fmt_msg(MODIS_E_INVALID_VDATA_TYPE, routine, msg);
                                break;
          }
          eng_field_order[i] = eng_data[current_index].field[i].order;
        }

      if(current_index == EN_SC_ANCILLARY_VDATA_START ||
	 current_index == EN_SC_ANCILLARY_VDATA_END)
	{
	   int j,k;

	   for(j=0; j<eng_data[current_index].num_fields; j++)
	   {
	     switch(eng_data[current_index].field[j].type)
	     {
	     case DFNT_INT8:
	        for(k=0; k<eng_data[current_index].field[j].order; k++)
		   eng_data[current_index].field[j].union_value[k].i8type
		   = (int8)SC_FILL_VALUE;
	        break;

	     case DFNT_UINT8:
	        for(k=0; k<eng_data[current_index].field[j].order; k++)
		   eng_data[current_index].field[j].union_value[k].ui8type 
		   = (uint8)SC_FILL_VALUE;
	        break;

	     case DFNT_INT16:
	        for(k=0; k<eng_data[current_index].field[j].order; k++)
		   eng_data[current_index].field[j].union_value[k].i16type
		   = (int16)SC_FILL_VALUE;
	        break;

	     case DFNT_UINT16:
	        for(k=0; k<eng_data[current_index].field[j].order; k++)
		   eng_data[current_index].field[j].union_value[k].ui16type
		   = (uint16)SC_FILL_VALUE;
	        break;

	     case DFNT_INT32:
	        for(k=0; k<eng_data[current_index].field[j].order; k++)
		   eng_data[current_index].field[j].union_value[k].i32type
		   = (int32)SC_FILL_VALUE;
	        break;

	     case DFNT_UINT32:
	        for(k=0; k<eng_data[current_index].field[j].order; k++)
		   eng_data[current_index].field[j].union_value[k].ui32type
		   = (uint32)SC_FILL_VALUE;
	        break;

	     /* default: Not needed, since valid type is guaranteed */
	     }
	   }
        }
  /**************************************************************************/
  /* Routine to create the Vdata.                                           */
  /**************************************************************************/
      returnStatus = create_Vdata(eng_data[current_index].vdata_name,
                                  eng_field_names,
                                  eng_data_types,
                                  eng_data[current_index].num_fields,
                                  eng_field_order);

  /**************************************************************************/
  /* Error if create_Vdata returns a fail code.                             */
  /**************************************************************************/
      if (returnStatus != MODIS_S_SUCCESS)
        {
          sprintf(msg,"Vdata Name: %s", eng_data[current_index].vdata_name);
          log_fmt_msg(MODIS_E_CREATE_VDATA, routine, msg);
        }

      current_index ++;
    }

  /**************************************************************************/
  /* Set up discarded packets vdata                                         */
  /**************************************************************************/
  if (returnStatus == MODIS_S_SUCCESS)
    {
      memset(eng_field_order, 0, sizeof(eng_field_order));
      eng_field_order[0] = 650;
      strcpy(Vdata_name, M01DISCARDED_PACKETS);
      strcpy(eng_field_names[0], M01DISCARDED_PKTS_FIELD);
      strcpy(eng_data_types[0], I8);
      returnStatus = create_Vdata(Vdata_name,
                                  eng_field_names,
                                  eng_data_types,
                                  1,
                                  eng_field_order);

  /**************************************************************************/
  /* Error if create_Vdata returns a fail code.                             */
  /**************************************************************************/
      if (returnStatus != MODIS_S_SUCCESS)
          log_fmt_msg(MODIS_E_CREATE_VDATA, routine, 
                     "The Discarded Packets Vdata could not be created in the file");

    }

  return (returnStatus);

}
