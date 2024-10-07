#include <ctype.h>
#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_IO.h"
#include "EN_eng_data.h"
#include "PC_pcf_info.h"
#include "MS_misc.h"
#include "PGS_MODIS_35005.h"
#include "hdfi.h"
#include "hdf.h"
#include "PGS_TD.h"


PGSt_SMF_status  parse_eng_data_list (EN_VDATA_TYPE_t  *eng_data)

/*
!C************************************************************************

!Description:  This function parses the eng_data_list_file creating Vdata
               in the eng_data array structure.

!Input Parameters:
               None

!Output Parameters:
               None

!Input/Output Parameters:
               EN_VDATA_TYPE_t    *eng_data    ** The eng_data array  **
	                                       ** structure           **

Return Values:
               MODIS_F_INVALID_ENG_DATA_LIST     (MODIS_35005.h)
               MODIS_S_SUCCESS                   (MODIS_35005.h)
               FAIL                              (hdf.h)

Externally Defined:  
               EN_MAX_FIELD_NAME_LENGTH      (EN_eng_data.h)
               EN_MAX_LINE_LEN               (EN_eng_data.h)
               EN_MAX_VDATA_NAME_LENGTH      (EN_eng_data.h)
               EN_MAX_VDATA_NUM_BITS         (EN_eng_data.h)
               EN_MAX_VDATA_ORDER            (EN_eng_data.h)
               EN_MAX_VDATA_TYPE             (EN_eng_data.h)
               EN_MIN_VDATA_NUM_BITS         (EN_eng_data.h)
               EN_MIN_VDATA_ORDER            (EN_eng_data.h)
               EN_MIN_VDATA_TYPE             (EN_eng_data.h)
               EN_NUM_SCAN_ELEMENTS          (EN_eng_data.h)
               EN_NUM_VDATAS                 (EN_eng_data.h)
               EN_VDATA_TYPE_t               (EN_eng_data.h)
               FAIL                          (hdf.h)
               global_input_pointer          (L1A_prototype.h)
               MODIS_E_GET_CONFIG_FAILED     (PGS_MODIS_35005.h)
               MODIS_E_NULL_POINTER          (PGS_MODIS_35005.h)
               MODIS_E_PGS_IO_GEN_CLOSE      (PGS_MODIS_35005.h)
               MODIS_E_PGS_IO_GEN_OPEN       (PGS_MODIS_35005.h)
               MODIS_F_INVALID_ENG_DATA_LIST (PGS_MODIS_35005.h)
               MODIS_S_SUCCESS               (PGS_MODIS_39501.h)
               MODIS_W_EXCESS_DATA_LIST      (PGS_MODIS_35005.h)
               MS_BLANK                      (MS_misc.h)
               MS_COMMENT                    (MS_misc.h)
               MS_NEW_LINE                   (MS_misc.h)
               PC_L1A_ENG_DATA_LIST_FILE     (PC_pcf_info.h)
               PGS_S_SUCCESS                 (PGS_SMF.h)
               PGS_TRUE                      (PGS_SMF.h)
               PGSd_IO_Gen_Read              (PGS_IO_Gen.h)
               PGSt_IO_Gen_FileHandle        (PGS_IO.h)
               PGSd_PC_VALUE_LENGTH_MAX      (PGS_PC.h)
               PGSt_integer                  (PGS_IO.h)
               PGSt_SMF_status               (PGS_SMF.h)
               PGSt_tag                      (PGS_TYPES.h)

Called By:
               initialize_level1a

Routines Called:
               create_eng_data_vdata_array
               create_eng_data_vdata_array_field
               log_fmt_msg
               PGS_IO_Gen_Close
               PGS_IO_Gen_Open
               PGS_PC_GetUniversalRef
               PGS_SMF_TestSuccessLevel 

!Revision History:
    $Log: parse_eng_data_list.c,v $
    Revision 5.1  2005/08/15 16:27:40  kuyper
    Corrected code to handle properly the possibility that field_name[0]<0.

    Revision 3.2  2002/10/03 19:45:38  vlin
    global_input_pointer removed.
    vlin@saicmodis.com

    Revision 3.1  2002/09/09 20:59:00  vlin
    Updated according to parse_eng_data_list.pdl revision 3.2

               Revision 2.1  2000/07/17
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Added code to work with aqua instrument.

!Team-unique Header:

       This software is developed by the MODIS Science Data Support Team 
       for the National Aeronautics and Space Administration, 
       Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:
       None

Design Notes:
       The input list must have the following format:

The first line of the ENG_DATA_LIST MUST be as follows:
-| 3333 $Revision: 5.1 $ Satellite Instrument (3333 - Auqa) (2222 - Terra)

This line designates the Satellite Instrument that the ENG_DATA_LIST is for.
The line starts as a comment for backward compatibility.

----------------------                          (<== Lines beginning w/ a dash are)
- Comment Line                                  (<== considered comments... )
----------------------                          (<== Comments are only allowed outside)
                                                (<== eng_data definition blocks)
------------------------------------            (<== Start of a eng_data definition block)
eng_data #1 Name                                (<== eng_data name)
------------------------------------            (<== a comment line separates eng_data name and field(s))
field_name_1 num_bits start_bit_pos order type  (<== field 1 )
field_name_2 num_bits start_bit_pos order type  ( etc ) 
   .
   .
   .
field_name_j num_bits start_bit_pos order type
------------------------------------            (<== comment line signals end of fields)

------------------------------------
eng_data #2 Name
------------------------------------
field_name_1 num_bits start_bit_pos order type
field_name_2 num_bits start_bit_pos order type
   .
   .
   .
field_name_k num_bits start_bit_pos order type
------------------------------------

..etc...
------------------------------------

   NOTES:
   1) Field names can't have spaces in them (although eng_data names can). 

   2) An extra field called "LAST_VALID_SCAN" will be added to each
      eng_data as its first field (field 0). Except for the
      Current/Prior S/C Ancillary Data Vdata. They do not need a
      LAST_VALID_SCAN field.

   3) Field.values can be any of the following types:
          int8 uint8 int16 uint16 int32 uint32.

   4) The bit positions listed in the eng_data starting with 
      "Engineering BB data" and below, are bit positions within eng.
      packet 1-2 (so subtract offset 144 from these bit positions to
      reference the corresponding bit positions within the 
      packet->data field)

!END*************************************************************************/

{
  PGSt_SMF_status  returnStatus=MODIS_S_SUCCESS;  
                                /* SMF-style message returned by function */
  PGSt_SMF_status  PGSstatus;   /* SMF-style message returned by function */
  PGSt_SMF_status  fileStatus;  /* SMF-style message returned by function */

  int   scanStatus;               /* stores the status of fscanf routines */
  char  *readStatus;              /* stores lines read from eng_data_list */

  char  *routine = "parse_eng_data_list";  /* pointer to name of function */

  char  msg [300];                       /* array to store error messages */

  PGSt_integer            version_num=1;                /* version number */
  PGSt_IO_Gen_AccessType  access_mode=PGSd_IO_Gen_Read; /* access mode to 
                                                            eng_data_list */
  PGSt_integer            file_version=1;                 /* file version */
  PGSt_IO_Gen_FileHandle  *eng_data_list_handle; /* file handle to 
                                                            eng_data_list */

  unsigned short num_bits;     /* stores the number of bits for the field */
  unsigned short start_bit_pos;/* stores the start bit position for data 
                                              in the packet for the field */
  unsigned short order;                 /* stores the order for the field */
  unsigned short type;   /* stores the type of field (int8, uint8, etc...)*/

  uint16 curr_eng_data_index=0;         /* stores the current vdata index */
  uint16 curr_field_index;  /* stores the current field count for a vdata */
  int    extra_lines;   /* counts the extra lines at end of eng_data_list */
  int    total_lines;   /* counts the total lines read from eng_data_list */
  int    num_vdatas_created;       /* stores the number of vdata created */
  int    number_of_vdata_fields; /* stores the number of fields 
                                                           for each vdata */
  char   line[EN_MAX_LINE_LEN];  /* stores line read in from eng_data_list
                                                                     file */
  char   field_name[EN_MAX_FIELD_NAME_LENGTH]; /* stores vdata field name */
  char   eng_data_name[EN_MAX_VDATA_NAME_LENGTH];    /* stores vdata name */

/**************************************************************************/
/* Check for NULL input parameters                                        */
/**************************************************************************/
  if (eng_data == NULL) {
    log_fmt_msg(MODIS_E_NULL_POINTER, routine, " ");
    return FAIL;
  }

/**************************************************************************/
/* Get Engineering Data List File Name                                    */
/**************************************************************************/
  PGSstatus = PGS_PC_GetUniversalRef(PC_L1A_ENG_DATA_LIST_FILE, 
                                  &version_num,
                                  global_input_pointer[0]);

  if (PGS_SMF_TestSuccessLevel(PGSstatus) != PGS_TRUE)
    log_fmt_msg(MODIS_E_GETCONFIG_FAILED, routine, 
                "unable to retrieve eng data list file name from pcf file");

/**************************************************************************/
/* Open Engineering Data List File                                        */
/**************************************************************************/
  PGSstatus = PGS_IO_Gen_Open(PC_L1A_ENG_DATA_LIST_FILE, 
                              access_mode, 
                              &eng_data_list_handle,
                              file_version);

  if (PGSstatus != PGS_S_SUCCESS)
    {
     log_fmt_msg(MODIS_E_PGS_IO_GEN_OPEN, routine, 
                                      "error opening eng_data_list file");
     returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
    }

/***************************************************************************/
/* Begin Parsing the Engineering Data List File                            */
/***************************************************************************/
  else
    {
     num_vdatas_created = 0;
     total_lines = 0;
     readStatus = fgets(line, sizeof(line), eng_data_list_handle);

     if (readStatus == line) {
         total_lines++;
         scanStatus = sscanf(line, "%s %d %*[$]Revision: %s $ ", field_name, 
                      &eng_data->instrument, (char *)&eng_data->revision);
         /* The %*[$] serves solely to prevent interpretation of this format
	     string as an RCS keyword in this code. It's eng_data->revision
	    which is supposed to match an RCS Revision keyword string. */
         if (scanStatus < 2) {
             log_fmt_msg(MODIS_F_INVALID_ENG_DATA_LIST, routine,
                         "Spacecraft Instrument not valid");
             returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
         }
         if (scanStatus < 3) {
             log_fmt_msg(MODIS_F_INVALID_ENG_DATA_LIST, routine,
                        "No RCS Revision number");
             returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
         }
     }
     else
        returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;

/***************************************************************************/
/* Loop until all vdata have been read, or an error has occurred           */
/***************************************************************************/
     while ((num_vdatas_created < EN_NUM_VDATAS) &&
            (returnStatus == MODIS_S_SUCCESS))
       {
        line[0] = MS_BLANK;

/***************************************************************************/
/* Loop until first non-comment or blank line line has been read           */
/***************************************************************************/
        while (((line[0] == MS_BLANK) ||
                (line[0] == MS_COMMENT) ||
                (line[0] == MS_NEW_LINE)) &&
                (returnStatus == MODIS_S_SUCCESS))
          {
           readStatus = fgets(line, sizeof(line), 
                              eng_data_list_handle);
           if (readStatus == line)
              total_lines++;
           else 
              returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
          }
/****************************************************************************/
/* Create a new vdata with the vdata name = line just read in               */
/****************************************************************************/ 
        if (returnStatus == MODIS_S_SUCCESS)
          {
           strncpy (eng_data_name, line, sizeof(eng_data_name));
           eng_data_name[strlen(line)-1]='\0';
           create_eng_data_vdata_array(eng_data_name, 
                                       eng_data, 
                                       &curr_eng_data_index,
                                       &curr_field_index);
/*****************************************************************************/
/* Read next line, if not a comment, eng_data_list file not correct          */
/*****************************************************************************/
           readStatus = fgets(line, EN_MAX_VDATA_NAME_LENGTH, 
                              eng_data_list_handle);
           if (readStatus == line)
             {
              if (line[0] != MS_COMMENT)
                {
                 returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
                 log_fmt_msg(MODIS_F_INVALID_ENG_DATA_LIST, routine, 
                             "missing separator line in eng_data_list");
                }
             }
           else
             returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
          }
/******************************************************************************/
/* Read the first Vdata field line. Field name, number of bits, starting bit  */
/*  position, order and type.                                                 */
/******************************************************************************/
        if (returnStatus == MODIS_S_SUCCESS)
          {
           total_lines++;
           memset(field_name, '\0', sizeof(field_name));
           number_of_vdata_fields = 0;

           readStatus = fgets(line, sizeof(line),
                              eng_data_list_handle);

           if (readStatus != line) 
              returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;

           scanStatus = sscanf(line, "%s %hu %hu %hu %hu", field_name,
                                  &num_bits, &start_bit_pos, &order, &type);
/******************************************************************************/
/* Read lines from Engineering Data List File which will be Vdata fields      */
/* until a comment line is reached, signaling the end of the field list.     */
/******************************************************************************/
           while ((field_name[0] != MS_COMMENT) &&
                  (returnStatus == MODIS_S_SUCCESS) &&
                  (scanStatus != EOF))
             {
              total_lines++;
              number_of_vdata_fields++;
/******************************************************************************/
/* If first character is not a letter, invalid Vdata field has been read      */
/******************************************************************************/
              if (!isalpha((int)*(unsigned char*)field_name))
                {
                 returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
                 sprintf(msg, "In line %d vdata field name %s in eng_data_list does not begin with a letter", 
                         total_lines, field_name);
                 log_fmt_msg(MODIS_F_INVALID_ENG_DATA_LIST, routine, msg);
                }
/******************************************************************************/
/* Check num_bits, order, and type for valid input ranges.                    */
/******************************************************************************/
              if ((num_bits < EN_MIN_VDATA_NUM_BITS) || (num_bits > EN_MAX_VDATA_NUM_BITS))
                { 
                 returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
                 sprintf(msg, "In line %d the number of bits is invalid for "
                 "vdata field %s in %s. Value: %d Range: %d to %d", total_lines,
                 field_name, eng_data_name, num_bits, EN_MIN_VDATA_NUM_BITS, 
                 EN_MAX_VDATA_NUM_BITS);
                 log_fmt_msg (MODIS_F_INVALID_ENG_DATA_LIST, routine, msg);
                }
              if (scanStatus != EN_NUM_SCAN_ELEMENTS)
               {
                 returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
                 sprintf(msg, "In line %d the number of attributes read for "
                 "vdata field %s in %s was %d.", total_lines, field_name, 
                 eng_data_name, scanStatus);
                 log_fmt_msg (MODIS_F_INVALID_ENG_DATA_LIST, routine, msg);
               }
              if ((order < EN_MIN_VDATA_ORDER) || (order > EN_MAX_VDATA_ORDER))
                {
                 returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
                 sprintf(msg, "In line %d the order number is invalid for "
                 "vdata field %s in %s. Value: %d Range: %d to %d", 
                 total_lines, field_name, eng_data_name, order, 
                 EN_MIN_VDATA_ORDER, EN_MAX_VDATA_ORDER);
                 log_fmt_msg (MODIS_F_INVALID_ENG_DATA_LIST, routine, msg);
                }
              if ((type < EN_MIN_VDATA_TYPE) || (type > EN_MAX_VDATA_TYPE))
                {
                 returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
                 sprintf(msg, "In line %d the type is invalid for vdata field "
                 " %s in %s. Value: %d Range: %d to %d", total_lines, field_name, 
                 eng_data_name, type, EN_MIN_VDATA_TYPE, EN_MAX_VDATA_TYPE);
                 log_fmt_msg (MODIS_F_INVALID_ENG_DATA_LIST, routine, msg);
                }

/*******************************************************************************/
/* Create Vdata Field with data read in from Engineering Data List File, and   */
/* read the next line from Engineering Data List File.                         */
/*******************************************************************************/           
              if (returnStatus == MODIS_S_SUCCESS) 
                 create_eng_data_vdata_array_field(field_name, 
                                                   num_bits, 
                                                   start_bit_pos,
                                                   order,
                                                   type,
                                                   eng_data,
                                                   curr_eng_data_index,
                                                   &curr_field_index);

              readStatus = fgets(line, sizeof(line), eng_data_list_handle);
              if (readStatus != line) 
                 returnStatus = MODIS_F_INVALID_ENG_DATA_LIST;
              else
                 scanStatus = sscanf(line, "%s %hu %hu %hu %hu", field_name,
                                    &num_bits, &start_bit_pos, &order, &type);
          }
        }
        else
           log_fmt_msg (MODIS_F_INVALID_ENG_DATA_LIST, routine, 
                 "Unexpected end of file in eng_data_list file");

        num_vdatas_created++;
       } /* end while < EN_NUM_VDATAS && returnStatus = MODIS_S_SUCCESS */
/********************************************************************************/
/* If all Vdata were created successfully, check to see if any extra lines     */
/*  are in the Engineering Data List File.                                      */
/********************************************************************************/    
     if (returnStatus == MODIS_S_SUCCESS) 
       {
        extra_lines = 0;
        fgets(line, sizeof(line), eng_data_list_handle);
        fileStatus = MODIS_S_SUCCESS;
/********************************************************************************/
/* Continue through Engineering Data List File, counting the number of extra    */
/* lines at the end of the file.                                                */
/********************************************************************************/
        while ((!feof(eng_data_list_handle))&&(!ferror(eng_data_list_handle)))
          {
           total_lines++;
           extra_lines++;
           fileStatus = MODIS_W_EXCESS_DATA_LIST;
           fgets(line, sizeof(line), eng_data_list_handle);
          } 
/*********************************************************************************/
/* If there were extra lines of data in Engineering Data List File, Print warning*/
/*********************************************************************************/
        if (fileStatus == MODIS_W_EXCESS_DATA_LIST)
          {
           sprintf(msg, "There are %d additional lines after the end of the "
           "Engineering Vdata. They were ignored", extra_lines);
           log_fmt_msg(MODIS_W_EXCESS_DATA_LIST, routine, msg);
          }
       }
       else
          log_fmt_msg(MODIS_F_INVALID_ENG_DATA_LIST, routine, 
                      "General Error Reading from ENG_DATA_LIST File");

/*********************************************************************************/
/* Close the Engineering Data List File                                          */
/*********************************************************************************/
     PGSstatus = PGS_IO_Gen_Close(eng_data_list_handle);
     if (PGSstatus != PGS_S_SUCCESS)
        log_fmt_msg(MODIS_E_PGS_IO_GEN_CLOSE, routine, 
                    "error closing eng_data_list_file");

    } /* end else */

  return returnStatus;

}
