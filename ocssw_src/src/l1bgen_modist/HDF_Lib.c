#include "HDF_Lib.h"
#include "PGS_Error_Codes.h"

/*
 * The following global variables are used in all functions within this module
 * for forming short message strings and assigning the error code.
 */

static char errmsgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
static PGSt_SMF_status errcode;

/*
 * The following are commonly used reasons which the low-level function
 * could fail (these are oriented toward operations -- in development, there
 * is always the possibility of a code error).
 */

char *invalidinputfile = "This is most likely due to an invalid file.\n"
                         "(does not meet file specs or is corrupted)";
char *corruptinputfile = "Operationally, this should never occur.\n"
                         "The file may be incomplete or corrupted.";
static char *corruptoutputfile = "Operationally, this should never occur.\n"
                          "The file or file pointer may have become corrupted.";
static char *nowriteoutputfile = "May be out of disk space or file became corrupted.";

/*
 * The following function (prototype is below, function is later in module)
 * is local in scope to this module.
 */
static
PGSt_SMF_status  assign_data_type(char *, int32 *);

PGSt_SMF_status read_attribute (int32   s_id,
                                char    *attr_name,
                                int32   TypeID,
                                void    *buffer)
/*
!C**************************************************************************
!Description:   Reads an attribute into buffer from a HDF file.

!Input Parameters:
     int32   s_id              sd_id for global, or sds_id for local, attrbute
     char    *attr_name        name of the attribute
     int32   TypeID            DFNT Type (used for safety checking only)

!Output Parameters:
     void    *buffer           buffer to hold the value(s) of the attribute.
                               If the attribute is a string, the buffer size
                               in the caller should be at least one more than
                               the string length. 

!Revision History:
 $Log: HDF_Lib.c,v $
 Revision 1.12  2005-01-18 14:34:43-05  ltan
 MOD_PR02_TERRA update to V5.0.4


 Revision 01.11  October 16, 2004  Razor Issue #200
 Casted Int32 variables in sprintf calls to "long" with the 
 format specifier "%ld" for better code portability.  
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 01.10 Feb. 1997
 Changed return value to PGSt_SMF_status.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.00 1996/03/19
 Initial development
 Neal Devine(neal.devine@gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  int32    attr_index  = 0;
  int32    count       = 0;
  intn     res         = 0;
  int32    data_type   = 0;
  char attr_name_buf[MAX_NC_NAME];

  attr_index = SDfindattr(s_id, attr_name);
  if (attr_index == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find attribute \"%s\" in the file.",
            attr_name);
    L1BErrorMsg("read_attribute", errcode, errmsgbuf, "SDfindattr", 0,
                invalidinputfile, False);
    return errcode;
  }

  res = SDattrinfo(s_id, attr_index, attr_name_buf, &data_type, &count);
  if (res == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not get attribute information for \"%s\".",
            attr_name);
    L1BErrorMsg("read_attribute", errcode, errmsgbuf, "SDattrinfo", 0,
                corruptinputfile, False);
    return errcode;
  }

   /*Make sure data_type is as expected*/

  if (data_type != TypeID) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Attribute \"%s\" has data type mismatch.\n"
            "File data type = %ld, expected type = %ld\n",
            attr_name, (long)data_type, (long)TypeID);
    L1BErrorMsg("read_attribute", errcode, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    return errcode;
  }

  res = SDreadattr(s_id, attr_index, buffer);
  if (res == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read attribute \"%s\" from file.",
            attr_name);
    L1BErrorMsg("read_attribute", errcode, errmsgbuf, "SDreadattr", 0,
                corruptinputfile, False);
    return errcode;
  }

  /* 
   * SDreadattr does not add a '\0' at the end of string. It should be
   * added explicitly if the attribute is a string.
   */

   if (data_type == DFNT_CHAR8)
   {
     char *buffer_ptr;
     buffer_ptr = buffer;
     buffer_ptr[count] = '\0';
   } 

  return(MODIS_S_OK);
}


PGSt_SMF_status read_part_sds_rank2 (int32   sd_id,
                                     char    *sds_name,
                                     int32   start0,
                                     int32   start1,
                                     int32   edge0,
                                     int32   edge1,
                                     void    *data)
/*
!C***************************************************************
!Description:        Reads a part of a any data type of 2D sds from 
                     a HDF file. 

!Input Parameters:
     int32   sd_id
     char    *sds_name
     int32   start0
     int32   start1
     int32   edge0
     int32   edge1

!Output Parameters:
     void    *data

!Revision History:
 Revision 01.00 Jan 1999  
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
    This function read part of sds data into a array and remain other part
    of the array. The array size should be the same as the SDS data in the
    HDF file. Otherwise, the location of the read data in the array is 
    meaningless. 
!END********************************************************************
*/
{
  int32   sds_index   = 0;
  int32   sds_id      = 0;
  int32   start[2]    = {0, 0};
  int32   edge[2]     = {0, 0};

  start[0] = start0;
  start[1] = start1;

  edge[0] = edge0;
  edge[1] = edge1;

  sds_index = SDnametoindex (sd_id, sds_name);
  if (sds_index == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("read_part_sds_rank2", errcode, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    return errcode;
  }

  sds_id = SDselect (sd_id, sds_index);
  if (sds_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not open access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_part_sds_rank2", errcode, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    return errcode;
  }

  if (SDreaddata(sds_id, start, NULL, edge, data) == FAIL) {
    SDendaccess(sds_id);
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read data from SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_part_sds_rank2", errcode, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    return errcode;
  }

  if (SDendaccess(sds_id)== FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_part_sds_rank2", errcode, errmsgbuf, "SDendaccess", 0,
                corruptinputfile, False);
    return errcode;
  }

  return (MODIS_S_OK);
}


PGSt_SMF_status read_part_sds_rank3 (int32   sd_id,
                                     char    *sds_name,
                                     int32   start0,
                                     int32   start1,
                                     int32   start2,
                                     int32   edge0,
                                     int32   edge1,
                                     int32   edge2,
                                     void    *data)
/*
!C***************************************************************
!Description:        Reads a part of a any data type of 3D sds from 
                     a HDF file. 

!Input Parameters:
     int32   sd_id
     char    *sds_name
     int32   start0
     int32   start1
     int32   start2
     int32   edge0
     int32   edge1
     int32   edge2

!Output Parameters:
     void    *data

!Revision History:
 Revision 02.10 April 1998
 Changed the interface  
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)

 Revision 01.00 March 1997
 Initial development.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
    This function read part of sds data into a array and remain other part
    of the array. The array size should be the same as the SDS data in the
    HDF file. Otherwise, the location of the read data in the array is 
    meaningless. 
!END********************************************************************
*/
{
  int32   sds_index   = 0;
  int32   sds_id      = 0;
  int32   start[3]    = {0, 0, 0};
  int32   edge[3]     = {0, 0, 0};

  start[0] = start0;
  start[1] = start1;
  start[2] = start2;

  edge[0] = edge0;
  edge[1] = edge1;
  edge[2] = edge2;

  sds_index = SDnametoindex (sd_id, sds_name);
  if (sds_index == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("read_part_sds_rank3", errcode, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    return errcode;
  }

  sds_id = SDselect (sd_id, sds_index);
  if (sds_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not open access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_part_sds_rank3", errcode, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    return errcode;
  }

  if (SDreaddata(sds_id, start, NULL, edge, data) == FAIL) {
    SDendaccess(sds_id);
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read data from SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_part_sds_rank3", errcode, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    return errcode;
  }

  if (SDendaccess(sds_id)== FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_part_sds_rank3", errcode, errmsgbuf, "SDendaccess", 0,
                corruptinputfile, False);
    return errcode;
  }

  return (MODIS_S_OK);
}


PGSt_SMF_status read_sds_rank1 (int32      file_id, 
                                char       *sds_name, 
                                int32      dim, 
                                void       *data)
/*
!C***************************************************************
!Description:        Reads a 1D array data from a HDF file.

!Input Parameters:
      int32     file_id          * HDF file ID *
      char      *sds_name        * sds name in HDF file *
      int32     dim              * dim of sds array *

!Output Parameters:
      void   *data               * 1D array *

!Revision History:
 Revision 01.00 October 1996
 Initial development.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
    Warning: The data type passed in must be the same as the data in hdf file
 Otherwise, the return results will be wrong. 
!END********************************************************************
*/
{
  int32   sds_index = 0;
  int32   sds_id    = 0;
  intn    result    = 0;
  int32   start[1]  = {0};
  int32   edge[1]   = {0};

  edge[0] = dim;

  sds_index = SDnametoindex (file_id, sds_name);
  if (sds_index == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("read_sds_rank1", errcode, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    return errcode;
  }

  sds_id = SDselect (file_id, sds_index);
  if (sds_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not open access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank1", errcode, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    return errcode;
  }

  result = SDreaddata (sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    SDendaccess(sds_id);
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read data from SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank1", errcode, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    return errcode;
  }

  result = SDendaccess (sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank1", errcode, errmsgbuf, "SDendaccess", 0,
                corruptinputfile, False);
    return errcode;
  }

  return (MODIS_S_OK);
}


PGSt_SMF_status read_sds_rank2 (int32       file_id, 
                                char        *sds_name, 
                                int32       dim1, 
                                int32       dim2, 
                                void        *data)
/*
!C***************************************************************
!Description:        Reads a 2D array.

!Input Parameters:
      int32     file_id          * HDF file ID *
      char      *sds_name        * sds name in HDF file *
      int32     dim1             * size of dimension 1 *
      int32     dim2             * size of dimension 2 *

!Output Parameters:
      void      data[][]           

!Revision History:
 Revision 02.10 April   1998
   This function is base on read_sds_rank2_float32()
 Zhenying Gu(zgu@gscmail.gsfc.nasa.gov
  
 Revision 01.00 October 1996
 Initial development.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
    The data type passed in should be the same as the data type in the
 hdf file. Otherwise the results will be wrong.
!END********************************************************************
*/
{
  int32   sds_index   = 0;
  int32   sds_id      = 0;
  intn    result      = 0;
  int32   start[2]    = {0, 0};
  int32   edge[2]     = {0, 0};

  edge[0] = dim1;
  edge[1] = dim2;

  sds_index = SDnametoindex (file_id, sds_name);
  if (sds_index == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("read_sds_rank2", errcode, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    return errcode;
  }

  sds_id = SDselect (file_id, sds_index);
  if (sds_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not open access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank2", errcode, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    return errcode;
  }

  result = SDreaddata (sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    SDendaccess(sds_id);
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read data from SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank2", errcode, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    return errcode;
  }

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank2", errcode, errmsgbuf, "SDendaccess", 0,
                corruptinputfile, False);
    return errcode;
  }

  return (MODIS_S_OK);
}


PGSt_SMF_status read_sds_rank3 (int32     file_id, 
                                char      *sds_name, 
                                int32     dim1, 
                                int32     dim2, 
                                int32     dim3, 
                                void      *data)
/*
!C***************************************************************
!Description:        Reads a 3D sds array.
!Input Parameters:
      int32     file_id          * HDF file ID *
      char      *sds_name        * sds name in HDF file *
      int32     dim1             * size of dimension 1 *
      int32     dim2             * size of dimension 2 *
      int32     dim3             * size of dimension 3 *

!Output Parameters:
      void      data[][][]         

!Revision History:
 Revision 02.10 April   1998
   Change from read float32 *** array to read general type three dimention array. Name changed 
 from read_sds_rank3_float32() to read_sds_rand3_1p. 
 Zhenying Gu(zgu@gscmail.gsfc.nasa.gov)
 
 Revision 01.00 October 1996
 Initial development.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  int32   sds_index   = 0;
  int32   sds_id      = 0;
  intn    result      = 0;
  int32   start[3]    = {0, 0, 0};
  int32   edge[3]     = {0, 0, 0};

  edge[0] = dim1;
  edge[1] = dim2;
  edge[2] = dim3;

  sds_index = SDnametoindex (file_id, sds_name);
  if (sds_index == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("read_sds_rank3", errcode, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    return errcode;
  }

  sds_id = SDselect (file_id, sds_index);
  if (sds_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not open access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank3", errcode, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    return errcode;
  }

  result = SDreaddata (sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    SDendaccess(sds_id);
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read data from SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank3", errcode, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    return errcode;
  }

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank3", errcode, errmsgbuf, "SDendaccess", 0,
                corruptinputfile, False);
    return errcode;
  }

  return (MODIS_S_OK);
}

PGSt_SMF_status read_sds_rank4(int32     file_id, 
                               char      *sds_name, 
                               int32     dim1, 
                               int32     dim2, 
                               int32     dim3, 
                               int32     dim4, 
                               void      *data)
/*
!C*********************************************************************
!Description:        Reads a 4D array data.

!Input Parameters:
      int32     file_id
      char      *sds_name
      int32     dim1
      int32     dim2
      int32     dim3
      int32     dim4
  
!Output Parameters:
      void      data[][][][]

!Revision History:
 Revision 02.10 April 1998
   Change reading to float32 **** to general type 4 dimention array. And change
 the name from read_sds_rank4_float32() to read_sds_rank4_1p.  
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)

 Revision 01.00 October 1996
 Initial development.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:

 !References and Credits
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

 !Design Notes:

!END********************************************************************
*/
{
  int32   sds_index   = 0;
  int32   sds_id      = 0;
  intn    result      = 0;
  int32   start[4]    = {0, 0, 0, 0};
  int32   edge[4]     = {0, 0, 0, 0};
 
  edge[0] = dim1;
  edge[1] = dim2;
  edge[2] = dim3;
  edge[3] = dim4;

  sds_index = SDnametoindex (file_id, sds_name);
  if (sds_index == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("read_sds_rank4", errcode, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    return errcode;
  }

  sds_id = SDselect (file_id, sds_index);
  if (sds_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not open access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank4", errcode, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    return errcode;
  }

  result = SDreaddata (sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    SDendaccess(sds_id);
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read data from SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank4", errcode, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    return errcode;
  }

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rank4", errcode, errmsgbuf, "SDendaccess", 0,
                corruptinputfile, False);
    return errcode;
  }

  return(MODIS_S_OK);
}



PGSt_SMF_status read_vdata (int32   v_id,
                            int32   start_record,
                            int32   records,
                            char    *vname,
                            char    *fname,
                            void    *buffer)
/*
!C**************************************************************************
!Description:   Reads vdata for given vdata_name and field_name from a HDF file. 

!Input Parameters:
     int32   v_id           file id for the vdata interfaces
     int32   start_record   starting record 
     int32   records        number of records to read
     char    *vname         vdata_name
     char    *fname         field_name
     void    *buffer        buffer to hold the data

!Output Parameters:
     void    *buffer        buffer to hold the data

!Revision History:
 Revision 01.00 Aug 1997
 Initial development
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  intn  hdf_return  = FAIL;
  int32 vd_ref      = 0;
  int32 vd_id       = 0;
  int32 n_records   = 0;
  int32 interlace   = 0;

  vd_ref = VSfind(v_id, vname);
  /*VSfind() returns the vdata reference number if successful and 0 otherwise.*/
  if (vd_ref == 0) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find Vdata \"%s\" in the file.",
            vname);
    L1BErrorMsg("read_vdata", errcode, errmsgbuf, "VSfind", 0,
                invalidinputfile, False);
    return errcode;
  }
    
  vd_id  = VSattach(v_id, vd_ref, "r");
  if (vd_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not attach to Vdata \"%s\".",
            vname);
    L1BErrorMsg("read_vdata", errcode, errmsgbuf, "VSattach", 0,
                corruptinputfile, False);
    return errcode;
  }

  hdf_return = VSinquire(vd_id, &n_records, &interlace, NULL, NULL, NULL);
  if (hdf_return == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not get information about Vdata \"%s\".",
            vname);
    L1BErrorMsg("read_vdata", errcode, errmsgbuf, "VSinquire", 0,
                corruptinputfile, False);
    return errcode;
  }

  if (records == -1)
    records = n_records - start_record;

  if (n_records < records + start_record ) {
    errcode = MODIS_F_OUT_OF_RANGE;
    sprintf(errmsgbuf, "Vdata \"%s\" contains too few records based on input arguments.",
            vname);
    L1BErrorMsg("read_vdata", errcode, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    return errcode;
  }

  hdf_return = VSseek (vd_id, start_record);
  if (hdf_return == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "For Vdata \"%s\", could not set current record to desired start_record.",
            vname);
    L1BErrorMsg("read_vdata", errcode, errmsgbuf, "VSseek", 0,
                corruptinputfile, False);
    return errcode;
  }

  hdf_return = VSsetfields (vd_id, fname);
  if (hdf_return == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "For Vdata \"%s\", could not set field to \"%s\".",
            vname, fname);
    L1BErrorMsg("read_vdata", errcode, errmsgbuf, "VSsetfields", 0,
                invalidinputfile, False);
    return errcode;
  }

  hdf_return = VSread (vd_id, buffer, records, interlace);
  if (hdf_return == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read data from Vdata \"%s\".",
            vname);
    L1BErrorMsg("read_vdata", errcode, errmsgbuf, "VSread", 0,
                corruptinputfile, False);
    return errcode;
  }

  hdf_return = VSdetach(vd_id);
  if (hdf_return == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not detach from Vdata \"%s\".",
            vname);
    L1BErrorMsg("read_vdata", errcode, errmsgbuf, "VSdetach", 0,
                corruptinputfile, False);
    return errcode;
  }

  return MODIS_S_OK;
}

static
PGSt_SMF_status  assign_data_type(char *data_type, int32 *number_type)
/*
!C**************************************************************************
!Description:   Reads vdata for given vdata_name and field_name from a HDF file. 

!Input Parameters:
     char *data_type

!Output Parameters:
     int32 *number_type

!Revision History:
 Revision 01.00 April 1998
 This is based on David Catozzi's L1A code.
 Initial development
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  
  if (!data_type || !number_type) {
    errcode = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg("assign_data_type", errcode, "NULL input pointer(s)",
                NULL, 0, "Coding error.", False);
    return errcode;
  }
  else if (strcmp(data_type, "int8") == 0)    *number_type = DFNT_INT8;
  else if (strcmp(data_type, "uint8") == 0)   *number_type = DFNT_UINT8;
  else if (strcmp(data_type, "int16") == 0)   *number_type = DFNT_INT16;
  else if (strcmp(data_type, "uint16") == 0)  *number_type = DFNT_UINT16;
  else if (strcmp(data_type, "int32") == 0)   *number_type = DFNT_INT32;
  else if (strcmp(data_type, "uint32") == 0)  *number_type = DFNT_UINT32;
  else if (strcmp(data_type, "float32") == 0) *number_type = DFNT_FLOAT32;
  else if (strcmp(data_type, "float64") == 0) *number_type = DFNT_FLOAT64;
  else if (strcmp(data_type, "char") == 0)    *number_type = DFNT_CHAR;
  else {
    errcode = MODIS_F_INVALID_ARGUMENT;
    sprintf(errmsgbuf,
            "data_type \"%s\" does not match any allowed type in function.",
            data_type);
    L1BErrorMsg("assign_data_type", errcode, errmsgbuf, NULL, 0,
                "Coding error.", False);
    return errcode;
  }
  return(returnStatus);
}


PGSt_SMF_status write_sds_rank1 (int32     file_id,
                                  char      *sds_name,
                                  char      *dim_name,
                                  int32     dim,
                                  char      *datatype,
                                  void     *data)
/*
!C*********************************************************************
!Description:      Writes a 1D array 

!Input Parameters:
      int32     file_id
      char      *sds_name
      char      *dim_name
      int32     dim
      char      *datatype
      int8      *data

!Output Parameters:

!Revision History:
 Revision 01.00 April. 1998
 Initial development.
 Zhenying Gu (zgu@barebackride.gsfc.nasa.gov)

!Team-unique Header:

 !References and Credits
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

 !Design Notes:
    This function is based on Zhidong hao's write_sds_rank2_float32() and
  David Catozzi's assign_data_type() with modification. 

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  intn     result      = 0;
  int32    sds_id      = 0;
  int32    rank        = 1;
  int32    dim_id      = 0;
  int32    start[1]    = {0};
  int32    edge[1]     = {0};
  int32    dim_size[1] = {0};
  int32    number_type;   

  dim_size[0] = dim;
  edge[0]     = dim;

  returnStatus = assign_data_type(datatype, &number_type);
  if (returnStatus != MODIS_S_OK) {
    L1BErrorMsg("write_sds_rank1", returnStatus, NULL, "assign_data_type", 0,
                NULL, False);
    return returnStatus;
  }

  sds_id = SDcreate(file_id, sds_name, number_type, rank, dim_size);
  if (sds_id == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not create SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("write_sds_rank1", errcode, errmsgbuf, "SDcreate", 0,
                corruptoutputfile, False);
    return errcode;
  }

  dim_id = SDgetdimid(sds_id, 0);
  if (dim_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not retrieve dimension id for SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank1", errcode, errmsgbuf, "SDgetdimid", 0,
                corruptoutputfile, False);
    SDendaccess(sds_id);
    return errcode;
  }

  result = SDsetdimname (dim_id, dim_name);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not set dimension name for SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank1", errcode, errmsgbuf, "SDsetdimname", 0,
                corruptoutputfile, False);
    SDendaccess(sds_id);
    return errcode;
  }

  result = SDwritedata(sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not write data to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank1", errcode, errmsgbuf, "SDwritedata", 0,
                nowriteoutputfile, False);
    SDendaccess(sds_id);
    return errcode;
  }

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank1", errcode, errmsgbuf, "SDendaccess", 0,
                corruptoutputfile, False);
    return errcode;
  }

  return(MODIS_S_OK);
}
 
PGSt_SMF_status write_sds_rank2 (int32     file_id, 
                                 char      *sds_name, 
                                 char      *dim_name1, 
                                 char      *dim_name2, 
                                 int32     dim1, 
                                 int32     dim2, 
                                 char      *datatype, 
                                 void      *data)
/*
!C*********************************************************************
!Description:      Writes a 2D sds

!Input Parameters:
      int32     file_id
      char      *sds_name
      char      *dim_name1
      char      *dim_name2
      int32     dim1
      int32     dim2
      char      *datatype
      void      data[][]

!Output Parameters:

!Revision History:
 Revision 01.00 April 1998
 Initial development.
 Zhenying Gu (zgu@barebackride.gsfc.nasa.gov)

!Team-unique Header:

 !References and Credits
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

 !Design Notes:
    This function is based on Zhidong hao's write_sds_rank2_float32() and
  David Catozzi's assign_data_type() with modification. 
!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int      i;
  int32    sds_id       = 0;
  intn     result       = 0;
  int32    rank         = 2;
  int32    dim_id[2]    = {0, 0};
  int32    start[2]     = {0, 0};
  int32    edge[2];
  int32    dim_size[2];
  char     *dim_name[2];
  int32    number_type;  

  edge[0] = dim1;
  edge[1] = dim2;
  dim_size[0] = dim1;
  dim_size[1] = dim2;
  dim_name[0] = dim_name1;
  dim_name[1] = dim_name2;
 
  returnStatus = assign_data_type(datatype, &number_type);
  if (returnStatus != MODIS_S_OK) {
    L1BErrorMsg("write_sds_rank2", returnStatus, NULL, "assign_data_type", 0,
                NULL, False);
    return returnStatus;
  }

  sds_id = SDcreate(file_id, sds_name, number_type, rank, dim_size);
  if (sds_id == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not create SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("write_sds_rank2", errcode, errmsgbuf, "SDcreate", 0,
                corruptoutputfile, False);
    return errcode;
  }
 
  for (i = 0; i < rank; i++)
  {
    dim_id[i] = SDgetdimid(sds_id, i);
    if (dim_id[i] == FAIL) {
      errcode = MODIS_F_HDF_ERROR;
      sprintf(errmsgbuf, "Could not retrieve dimension #%d id for SDS \"%s\".",
              (i+1), sds_name);
      L1BErrorMsg("write_sds_rank2", errcode, errmsgbuf, "SDgetdimid", 0,
                  corruptoutputfile, False);
      SDendaccess(sds_id);
      return errcode;
    }

    result = SDsetdimname (dim_id[i], dim_name[i]);
    if (result == FAIL) {
      errcode = MODIS_F_HDF_ERROR;
      sprintf(errmsgbuf, "Could not set dimension #%d name for SDS \"%s\".",
              (i+1), sds_name);
      L1BErrorMsg("write_sds_rank2", errcode, errmsgbuf, "SDsetdimname", 0,
                  corruptoutputfile, False);
      SDendaccess(sds_id);
      return errcode;
    }
  }

  result = SDwritedata(sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not write data to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank2", errcode, errmsgbuf, "SDwritedata", 0,
                nowriteoutputfile, False);
    SDendaccess(sds_id);
    return errcode;
  }

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank2", errcode, errmsgbuf, "SDendaccess", 0,
                corruptoutputfile, False);
    return errcode;
  }

  return(MODIS_S_OK);
}

PGSt_SMF_status write_sds_rank3 (int32      file_id, 
                                 char       *sds_name, 
                                 char       *dim_name1, 
                                 char       *dim_name2, 
                                 char       *dim_name3, 
                                 int32      dim1, 
                                 int32      dim2, 
                                 int32      dim3, 
                                 char       *datatype, 
                                 void       *data)
/*
!C*********************************************************************
!Description:     Writes a 3D sds 

!Input Parameters:
      int32     file_id
      char      *sds_name
      char      *dim_name1
      char      *dim_name2
      char      *dim_name3
      int32     dim1
      int32     dim2
      int32     dim3
      char      *datatype
      void      data[][][]

!Output Parameters:

!Revision History:
 Revision 01.00 April 1998
 Initial development.
 Zhenying Gu (zgu@barebackride.gsfc.nasa.gov)

!Team-unique Header:

 !References and Credits
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

 !Design Notes:
    This function is based on Zhidong hao's write_sds_rank3_float32() and
  David Catozzi's assign_data_type() with modification. 
!END********************************************************************
*/
{ 
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int      i;
  int32    sds_id       = 0;
  intn     result       = 0;
  int32    rank         = 3;
  int32    dim_id[3]    = {0, 0, 0};
  int32    start[3]     = {0, 0, 0};
  int32    edge[3];
  int32    dim_size[3];
  char     *dim_name[3];
  int32    number_type;  

  edge[0] = dim1;
  edge[1] = dim2;
  edge[2] = dim3;
  dim_size[0] = dim1;
  dim_size[1] = dim2;
  dim_size[2] = dim3;
  dim_name[0] = dim_name1;
  dim_name[1] = dim_name2;
  dim_name[2] = dim_name3;
  
  returnStatus = assign_data_type(datatype, &number_type);
  if (returnStatus != MODIS_S_OK) {
    L1BErrorMsg("write_sds_rank3", returnStatus, NULL, "assign_data_type", 0,
                NULL, False);
    return returnStatus;
  }

  sds_id = SDcreate(file_id, sds_name, number_type, rank, dim_size);
  if (sds_id == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not create SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("write_sds_rank3", errcode, errmsgbuf, "SDcreate", 0,
                corruptoutputfile, False);
    return errcode;
  }

  for (i = 0; i < rank; i++)
  {
    dim_id[i] = SDgetdimid(sds_id, i);
    if (dim_id[i] == FAIL) {
      errcode = MODIS_F_HDF_ERROR;
      sprintf(errmsgbuf, "Could not retrieve dimension #%d id for SDS \"%s\".",
              (i+1), sds_name);
      L1BErrorMsg("write_sds_rank3", errcode, errmsgbuf, "SDgetdimid", 0,
                  corruptoutputfile, False);
      SDendaccess(sds_id);
      return errcode;
    }

    result = SDsetdimname (dim_id[i], dim_name[i]);
    if (result == FAIL) {
      errcode = MODIS_F_HDF_ERROR;
      sprintf(errmsgbuf, "Could not set dimension #%d name for SDS \"%s\".",
              (i+1), sds_name);
      L1BErrorMsg("write_sds_rank3", errcode, errmsgbuf, "SDsetdimname", 0,
                  corruptoutputfile, False);
      SDendaccess(sds_id);
      return errcode;
    }
  }

  result = SDwritedata(sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not write data to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank3", errcode, errmsgbuf, "SDwritedata", 0,
                nowriteoutputfile, False);
    SDendaccess(sds_id);
    return errcode;
  }

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank3", errcode, errmsgbuf, "SDendaccess", 0,
                corruptoutputfile, False);
    return errcode;
  }

  return(MODIS_S_OK);
}

PGSt_SMF_status write_sds_rank4 (int32     file_id,
                                 char      *sds_name,
                                 char      *dim_name1,
                                 char      *dim_name2,
                                 char      *dim_name3,
                                 char      *dim_name4,
                                 int32     dim1,
                                 int32     dim2,
                                 int32     dim3,
                                 int32     dim4,
                                 char      *datatype,
                                 void      *data)
/*
!C*********************************************************************
!Description:     Writes a 4D sds

!Input Parameters:
      int32     file_id
      char      *sds_name
      char      *dim_name1
      char      *dim_name2
      char      *dim_name3
      char      *dim_name4
      int32     dim1
      int32     dim2
      int32     dim3
      int32     dim4
      char      *datatype
      void      data[][][][]

!Output Parameters:

!Revision History:
 Revision 01.01  April 2002  Razor Issue #183
 Changed definition of number_type to be consistent with other routines.
 Gwyn Fireman, SAIC GSO (fireman@mcst.gsfc.nasa.gov)

 Revision 01.00 April 1998
 Initial development.
 Zhenying Gu (zgu@barebackride.gsfc.nasa.gov)

!Team-unique Header:

 !References and Credits
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

 !Design Notes:
    This function is based on Zhidong hao's write_sds_rank4_float32() and
  David Catozzi's assign_data_type() with modification.
!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int      i;
  int32    sds_id       = 0;
  intn     result       = 0;
  int32    rank         = 4;
  int32    dim_id[4]    = {0, 0, 0, 0};
  int32    start[4]     = {0, 0, 0, 0};
  int32    edge[4];
  int32    dim_size[4];
  char     *dim_name[4];
  int32    number_type;

  edge[0] = dim1;
  edge[1] = dim2;
  edge[2] = dim3;
  edge[3] = dim4;
  dim_size[0] = dim1;
  dim_size[1] = dim2;
  dim_size[2] = dim3;
  dim_size[3] = dim4;
  dim_name[0] = dim_name1;
  dim_name[1] = dim_name2;
  dim_name[2] = dim_name3;
  dim_name[3] = dim_name4;

  returnStatus = assign_data_type(datatype, &number_type);
  if (returnStatus != MODIS_S_OK) {
    L1BErrorMsg("write_sds_rank4", returnStatus, NULL, "assign_data_type", 0,
                NULL, False);
    return returnStatus;
  }

  sds_id = SDcreate(file_id, sds_name, number_type, rank, dim_size);
  if (sds_id == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not create SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("write_sds_rank4", errcode, errmsgbuf, "SDcreate", 0,
                corruptoutputfile, False);
    return errcode;
  }

  for (i = 0; i < rank; i++)
  {
    dim_id[i] = SDgetdimid(sds_id, i);
    if (dim_id[i] == FAIL) {
      errcode = MODIS_F_HDF_ERROR;
      sprintf(errmsgbuf, "Could not retrieve dimension #%d id for SDS \"%s\".",
              (i+1), sds_name);
      L1BErrorMsg("write_sds_rank4", errcode, errmsgbuf, "SDgetdimid", 0,
                  corruptoutputfile, False);
      SDendaccess(sds_id);
      return errcode;
    }

    result = SDsetdimname (dim_id[i], dim_name[i]);
    if (result == FAIL) {
      errcode = MODIS_F_HDF_ERROR;
      sprintf(errmsgbuf, "Could not set dimension #%d name for SDS \"%s\".",
              (i+1), sds_name);
      L1BErrorMsg("write_sds_rank4", errcode, errmsgbuf, "SDsetdimname", 0,
                  corruptoutputfile, False);
      SDendaccess(sds_id);
      return errcode;
    }
  }

  result = SDwritedata(sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not write data to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank4", errcode, errmsgbuf, "SDwritedata", 0,
                nowriteoutputfile, False);
    SDendaccess(sds_id);
    return errcode;
  }

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("write_sds_rank4", errcode, errmsgbuf, "SDendaccess", 0,
                corruptoutputfile, False);
    return errcode;
  }

  return(MODIS_S_OK);
}


PGSt_SMF_status read_sds_rankn (int32     sd_id,
                                char      *sds_name,
                                int32     data_type,
                                int32     rank,
                                int32     *start,
                                int32     *edge,
                                void      *data)
/*
!C**************************************************************************
!Description:
   Read part or all of an n-dimensional SDS into a 1D buffer.  If part of
   the SDS is to be read, then elements within a dimension must be
   consecutive (the "stride" is assumed to be 1 for all dimensions).

!Input Parameters:
   int32  sd_id        File ID for Science Data (SD) access.
   char   *sds_name    SDS name in the HDF file (must be non-NULL)
   int32  data_type    Data type of the SDS (for consistency checking)
   int32  rank         Rank of the SDS (dimensions of start and edge,
                       for consistency checking)
   int32  *start       Starting element index (0-based) for each dimension
   int32  *edge        Number of elements to read in each dimension

!Output Parameters:
   void   *data        One-dimension buffer holding the data.

!Revision History:
   Revision 01.01  October 16, 2004  Razor Issue #200
   Casted Int32 variables in sprintf calls to "long" with the 
   format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.00 November 6, 1999
   Initial development
   Jim Rogers (rogers@mst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   Based on the L1B functions "read_sds_rank1", and 2, and 3, etc.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  int32   sds_index   = 0;
  int32   sds_id      = 0;
  intn    result      = 0;

    /* The following used to retrieve SDS information */

  char    f_sds_name[MAX_NC_NAME];
  int32   f_rank          = 0;
  int32   f_dim_sizes[MAX_VAR_DIMS];
  int32   f_data_type     = 0;
  int32   f_nattrs        = 0;


  /*
   * Check input arguments.
   */

  if (sd_id == FAIL || !sds_name || rank <= 0 || !start || !edge || !data) {
    errcode = MODIS_F_INVALID_ARGUMENT;
    sprintf(errmsgbuf, "One or more input arguments are invalid.");
    L1BErrorMsg("read_sds_rankn", errcode, errmsgbuf, NULL, 0,
                "This is a code defect.", False);
    return errcode;
  }

  /*
   * Convert the name of the SDS to the index in the file.
   */

  sds_index = SDnametoindex (sd_id, sds_name);
  if (sds_index == FAIL) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "Could not find SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg("read_sds_rankn", errcode, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    return errcode;
  }

  /*
   * Open SDS access to the data set.
   */

  sds_id = SDselect (sd_id, sds_index);
  if (sds_id == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not open access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rankn", errcode, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    return errcode;
  }

  /*
   * Get info about the SDS.  Check the type and rank vs. inputs.
   */

  result = SDgetinfo (sds_id, f_sds_name, &f_rank, f_dim_sizes,
                      &f_data_type, &f_nattrs);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not retrieve info for SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rankn", errcode, errmsgbuf, "SDgetinfo", 0,
                corruptinputfile, False);
    return errcode;
  }
  if (data_type != f_data_type) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "SDS \"%s\" has data type mismatch.\n"
            "File data type = %ld, expected type = %ld\n",
            sds_name, (long)f_data_type, (long)data_type);
    L1BErrorMsg("read_sds_rankn", errcode, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    return errcode;
  }
  if (rank != f_rank) {
    errcode = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, "SDS \"%s\" has matrix rank mismatch.\n"
            "File data rank = %ld, expected rank = %ld\n",
            sds_name, (long)f_rank, (long)rank);
    L1BErrorMsg("read_sds_rankn", errcode, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    return errcode;
  }

  /*
   * Read the data.
   */

  result = SDreaddata (sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    SDendaccess(sds_id);
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not read data from SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rankn", errcode, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    return errcode;
  }

  /*
   * End SDS access.
   */

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg("read_sds_rankn", errcode, errmsgbuf, "SDendaccess", 0,
                corruptinputfile, False);
    return errcode;
  }

  return (MODIS_S_OK);
}

PGSt_SMF_status write_sds_rankn (int32     file_id,
                                 char      *sds_name,
                                 int32     data_type,
                                 int32     rank,
                                 int32     *edge,
                                 char      **dim_name,
                                 void      *data)
/*
!C**************************************************************************
!Description:
   Write all of an n-dimensional SDS into an HDF file. (The maximum
   rank is MAXRANK, defined below)

!Input Parameters:
   int32  file_id      File ID for Science Data (SD) access.
   char   *sds_name    SDS name in the HDF file (must be non-NULL)
   int32  data_type    HDF Data type of the SDS
   int32  rank         Rank of the SDS (maximum of MAXRANK)
   int32  *edge        Array holding the size of each dimension.
   char   **dim_name   Array of dimension names.
   void   *data        One-dimension buffer holding the data.

!Output Parameters:
   (none)

!Revision History:
   Revision 01.01  October 16, 2004  Razor Issue #200
   Casted variables of type Int32 in sprintf calls to "long"
   with format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.00 November 22, 1999
   Initial development.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

 !References and Credits
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

 !Design Notes:

!END********************************************************************
*/
{ 
  int      i;
  int32    sds_id       = 0;
  intn     result       = 0;
#define MAXRANK 7
  int32    dim_id[MAXRANK];
  int32    start[MAXRANK];
  char *location = "write_sds_rankn";

  if (rank > MAXRANK)
  {
    errcode = MODIS_F_NOK;
    sprintf(errmsgbuf, "Could not create SDS \"%s\". Rank (%ld) exceeds max (%d)",
            sds_name, (long)rank, MAXRANK);
    L1BErrorMsg(location, errcode, errmsgbuf, NULL, 0,
                "Code bug.", False);
    return errcode;
  }

  for (i = 0; i < rank; i++)
    start[i] = 0;

  sds_id = SDcreate(file_id, sds_name, data_type, rank, edge);
  if (sds_id == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not create SDS \"%s\" in the file.",
            sds_name);
    L1BErrorMsg(location, errcode, errmsgbuf, "SDcreate", 0,
                corruptoutputfile, False);
    return errcode;
  }

  for (i = 0; i < rank; i++)
  {
    dim_id[i] = SDgetdimid(sds_id, i);
    if (dim_id[i] == FAIL) {
      errcode = MODIS_F_HDF_ERROR;
      sprintf(errmsgbuf, "Could not retrieve dimension #%d id for SDS \"%s\".",
              (i+1), sds_name);
      L1BErrorMsg(location, errcode, errmsgbuf, "SDgetdimid", 0,
                  corruptoutputfile, False);
      SDendaccess(sds_id);
      return errcode;
    }

    result = SDsetdimname (dim_id[i], dim_name[i]);
    if (result == FAIL) {
      errcode = MODIS_F_HDF_ERROR;
      sprintf(errmsgbuf, "Could not set dimension #%d name for SDS \"%s\".",
              (i+1), sds_name);
      L1BErrorMsg(location, errcode, errmsgbuf, "SDsetdimname", 0,
                  corruptoutputfile, False);
      SDendaccess(sds_id);
      return errcode;
    }
  }

  result = SDwritedata(sds_id, start, NULL, edge, data);
  if (result == FAIL) {
    errcode = MODIS_F_WRITE_ERROR;
    sprintf(errmsgbuf, "Could not write data to SDS \"%s\".",
            sds_name);
    L1BErrorMsg(location, errcode, errmsgbuf, "SDwritedata", 0,
                nowriteoutputfile, False);
    SDendaccess(sds_id);
    return errcode;
  }

  result = SDendaccess(sds_id);
  if (result == FAIL) {
    errcode = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, "Could not end access to SDS \"%s\".",
            sds_name);
    L1BErrorMsg(location, errcode, errmsgbuf, "SDendaccess", 0,
                corruptoutputfile, False);
    return errcode;
  }

  return(MODIS_S_OK);
}


#include <math.h>        /* for the fabs function */
#include <stdlib.h>      /* for the atoi, atof and atol functions */


PGSt_SMF_status Check_Valid_Range (char  *data_name,
                                   int32 data_type,
                                   char  *a_lb,
                                   char  *a_ub,
                                   char  *a_fillvalue,
                                   int32 count,
                                   void  *buffer)
/*
!C**************************************************************************
!Description:
   Check that all values of an array are within a valid range or are a fill
   value (which could be outside the defined valid range).  If all values
   meet the above described criteria, then return MODIS_S_OK.  If any value
   is found not to meet the above criteria, then write an error message and
   return a value of MODIS_F_OUT_OF_RANGE.  This function only checks
   integer or floating point types, not CHAR8 arrays.

!Input Parameters:
   char    *data_name        Name the data item, for writing to error message
   int32   data_type         DFNT type of data in the array (should be one of
                             the integer or floating types, not CHAR8).
   char    *a_lb             ASCII representation of the lower bound on the
                             values (NULL means no lower bound)
   char    *a_ub             ASCII representation of the upper bound on the
                             valid range of values (NULL means no upper bound)
   char    *a_fillvalue      ASCII representation of the fill value on the
                             valid range of values (NULL means no fill value)
   int32   count             Size of array
   void    *buffer           buffer holding the value(s) of the array

!Output Parameters:
   (none)

!Revision History:

   Revision 1.0.2  March 27, 2003,  Razor Issue #173
   In macros CHECK_FVALS_ON_LB and CHECK_FVALS_ON_UB, remove the length character
   "l" in the type code "%lf" of the double type array element data[i] for ANSI-C
   compliance.
   Liqin Tan, SAIC GSO (ltan@saicmodis.com)

   Revision 1.0.1  January 17, 2002  Razor Issue #172
   Improve portability of code to 64-bit mode.
   Change variable i to type long
   Rework macros so that explicit casts are used based on data type 
   (especially ASSIGN_INT_VALUE)
   Cast atol output to int32
   Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.00 November 6, 1999
   Initial development
   Jim Rogers (rogers@mst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   The lower bound, upper bound and fill values are ASCII representations
   of the actual numbers (e.g., ".011").  atol or atof are used to convert
   to actual values consistent with the data type of the attribute.
   Checking of values is not accomplished on char type attributes.
   The fill value may be any one value that lies outside the valid range.

!END********************************************************************
*/
{
  long i;
  char *fvptr;
  char *unknown_data_name = "(unknown)";
  double del, absd, absdfv;

/***************** MACRO ASSIGN_INT_VALUE(aptr,a,alb,aub)*****************
For one of the integer types holding less than 32 bits precision, check
the input value against the limits of the precision of the variable.
If OK, then assign to the variable (implicit cast occurs).  Note: the
atol operation is assumed done with at least 32 bits of precision.
Global variables used:
  errcode, errmsgbuf
Function variables used:
  i, data_name
Macro variables:
  aptr = pointer to the string holding the value.
  a    = variable to assign
  alb  = lower extremum of data type range
  aub  = upper extremum of data type range
**************************************************************************/

#define ASSIGN_INT_VALUE(aptr, a, alb, aub, type_of_data)                 \
  i = atol(aptr);                                                         \
  if (i < (int32) alb || i > (int32) aub) {                               \
    errcode = MODIS_F_NOK;                                                \
    sprintf(errmsgbuf,                                                    \
            "Checking valid range of \"%s\", the value \"%s\" is not\n"   \
            "a valid number for the precision of the data type.",         \
            data_name, aptr);                                             \
    L1BErrorMsg("Check_Valid_Range", errcode, errmsgbuf, NULL, 0,         \
                "Code defect.", False);                                   \
    return errcode;                                                       \
  }                                                                       \
  else                                                                    \
    {                                                                     \
      switch (type_of_data)                                               \
      { case DFNT_INT8:    a = (int8)    i; break;                        \
        case DFNT_UINT8:   a = (uint8)   i; break;                        \
        case DFNT_INT16:   a = (int16)   i; break;                        \
        case DFNT_UINT16:  a = (uint16)  i; break;                        \
        case DFNT_INT32:   a = (int32)   i; break;                        \
        case DFNT_UINT32:  a = (uint32)  i;                               \
      }                                                                   \
    }
    
/********************* MACRO CHECK_IVALS_ON_LB ************
Loop through one of the integer type arrays and check each value
against the lower bound and the fill value.
Global variables used:
  errcode, errmsgbuf
Function variables used:
  i, data_name, a_lb, data, lb, fillvalue
**********************************************************************/
#define CHECK_IVALS_ON_LB                                            \
for (i = 0; i < count; i++) {                                        \
  if (data[i] < lb && data[i] != fillvalue) {                        \
    errcode = MODIS_F_OUT_OF_RANGE;                                  \
    sprintf(errmsgbuf,                                               \
            "One or more values of \"%s\" is less than the lower bound.\n" \
            "value[%ld] = %ld\nlower bnd = %s\n",                    \
            data_name, i, (long) data[i], a_lb);                     \
    L1BErrorMsg("Check_Valid_Range", errcode, errmsgbuf, NULL, 0,    \
                invalidinputfile, False);                            \
    return errcode;                                                  \
  }                                                                  \
}

/********************* MACRO CHECK_FVALS_ON_LB ************
Loop through one of the floating point type arrays and check each value
against the lower bound and the fill value.
Global variables used:
  errcode, errmsgbuf
Function variables used:
  i, data_name, a_lb, absd, absdfv, del, data, lb, fillvalue
**********************************************************************/
#define CHECK_FVALS_ON_LB                                            \
for (i = 0; i < count; i++) {                                        \
  absd = fabs((double) data[i]); absdfv = fabs((double)fillvalue);   \
  del = fabs((double) (data[i] - fillvalue));                        \
  if (absd > 0) del /= absd; else if(absdfv > 0) del /= absdfv;      \
  if (data[i] < lb && del > 1.e-05) {                                \
    errcode = MODIS_F_OUT_OF_RANGE;                                  \
    sprintf(errmsgbuf,                                               \
            "One or more values of \"%s\" is less than the lower bound.\n" \
            "value[%ld] = %f\nlower bnd = %s\n",                          \
            data_name, i, (double) data[i], a_lb);                         \
    L1BErrorMsg("Check_Valid_Range", errcode, errmsgbuf, NULL, 0,    \
                invalidinputfile, False);                            \
    return errcode;                                                  \
  }                                                                  \
}

/********************* MACRO CHECK_IVALS_ON_UB ************
Loop through one of the integer type arrays and check each value
against the lower bound and the fill value.
Global variables used:
  errcode, errmsgbuf
Function variables used:
  i, data_name, a_ub, data, ub, fillvalue
**********************************************************************/
#define CHECK_IVALS_ON_UB                                            \
for (i = 0; i < count; i++) {                                        \
  if (data[i] > ub && data[i] != fillvalue) {                        \
    errcode = MODIS_F_OUT_OF_RANGE;                                  \
    sprintf(errmsgbuf,                                               \
            "One or more values of \"%s\" is greater than the upper bound.\n" \
            "value[%ld] = %ld\nupper bnd = %s\n",                    \
            data_name, i, (long) data[i], a_ub);                     \
    L1BErrorMsg("Check_Valid_Range", errcode, errmsgbuf, NULL, 0,    \
                invalidinputfile, False);                            \
    return errcode;                                                  \
  }                                                                  \
}

/********************* MACRO CHECK_FVALS_ON_UB ************
Loop through one of the floating point type arrays and check each value
against the lower bound and the fill value.
Global variables used:
  errcode, errmsgbuf
Function variables used:
  i, data_name, a_ub, absd, absdfv, del, data, ub, fillvalue
**********************************************************************/
#define CHECK_FVALS_ON_UB                                            \
for (i = 0; i < count; i++) {                                        \
  absd = fabs((double) data[i]); absdfv = fabs((double)fillvalue);   \
  del = fabs((double) (data[i] - fillvalue));                        \
  if (absd > 0) del /= absd; else if(absdfv > 0) del /= absdfv;      \
  if (data[i] > ub && del > 1.e-05) {                                \
    errcode = MODIS_F_OUT_OF_RANGE;                                  \
    sprintf(errmsgbuf,                                               \
            "One or more values of \"%s\" is greater than the upper bound.\n" \
            "value[%ld] = %f\nupper bnd = %s\n",                              \
            data_name, i, (double) data[i], a_ub);                            \
    L1BErrorMsg("Check_Valid_Range", errcode, errmsgbuf, NULL, 0,    \
                invalidinputfile, False);                            \
    return errcode;                                                  \
  }                                                                  \
}

  /*
   * Check that array exists.
   */

  if (count <= 0 || !buffer) {
    errcode = MODIS_F_INVALID_ARGUMENT;
    sprintf(errmsgbuf, "Input arguments are invalid.  \"count\" is zero "
                       "or \"buffer\" is NULL.");
    L1BErrorMsg("Check_Valid_Range", errcode, errmsgbuf, NULL, 0,
                "Code defect (check that data buffer pointer is assigned).",
                False);
    return errcode;
  }

  /*
   * The data_name should be non-NULL.  If not, set it to a default value
   * (only important for writing error messages).
   */


  if (!data_name) data_name = unknown_data_name;

  /*
   * Check values against lower bound.
   */

  if (a_lb) {

    /*
     * Set the temporary fillvalue pointer (fvptr) to a valid value.
     */

    if (a_fillvalue)
      fvptr = a_fillvalue;
    else
      fvptr = a_lb;

    /*
     * Step through each data type and perform checking.
     */

    if (data_type == DFNT_INT8) {
      int8  *data = (int8 *) buffer;
      int8  lb, fillvalue;
      ASSIGN_INT_VALUE(a_lb,lb,-128,127, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,-128,127, data_type)
      CHECK_IVALS_ON_LB
    }
    else if (data_type == DFNT_UINT8) {
      uint8  *data = (uint8 *) buffer;
      uint8  lb, fillvalue;
      ASSIGN_INT_VALUE(a_lb,lb,0,256, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,0,256, data_type)
      CHECK_IVALS_ON_LB
    }
    else if (data_type == DFNT_INT16) {
      int16  *data = (int16 *) buffer;
      int16  lb, fillvalue;
      ASSIGN_INT_VALUE(a_lb,lb,-32768,32767, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,-32768,32767, data_type)
      CHECK_IVALS_ON_LB
    }
    else if (data_type == DFNT_UINT16) {
      uint16  *data = (uint16 *) buffer;
      uint16  lb, fillvalue;
      ASSIGN_INT_VALUE(a_lb,lb,0,65535, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,0,65535, data_type)
      CHECK_IVALS_ON_LB
    }
    else if (data_type == DFNT_INT32) {
      int32  *data = (int32 *) buffer;
      int32  lb, fillvalue;
      lb = (int32) atol(a_lb);
      fillvalue = (int32) atol(fvptr);
      CHECK_IVALS_ON_LB
    }
    else if (data_type == DFNT_UINT32) {
      uint32  *data = (uint32 *) buffer;
      uint32  lb, fillvalue;
      ASSIGN_INT_VALUE(a_lb,lb,0,2147483647, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,0,2147483647, data_type);
      CHECK_IVALS_ON_LB
    }
    else if (data_type == DFNT_FLOAT32) {
      float32  *data = (float32 *) buffer;
      float32  lb, fillvalue;
      lb = (float32) atof(a_lb);
      fillvalue = (float32) atof(fvptr);
      CHECK_FVALS_ON_LB
    }
    else if (data_type == DFNT_FLOAT64) {
      float64  *data = (float64 *) buffer;
      float64  lb, fillvalue;
      lb = atof(a_lb);
      fillvalue = atof(fvptr);
      CHECK_FVALS_ON_LB
    }
    else {
      return MODIS_S_OK;
    }
  }

  /*
   * Check values against upper bound.
   */

  if (a_ub) {

    /*
     * Set the temporary fillvalue pointer (fvptr) to a valid value.
     */

    if (a_fillvalue)
      fvptr = a_fillvalue;
    else
      fvptr = a_ub;

    /*
     * Step through each data type and perform checking.
     */

    if (data_type == DFNT_INT8) {
      int8  *data = (int8 *) buffer;
      int8  ub, fillvalue;
      ASSIGN_INT_VALUE(a_ub,ub,-128,127, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,-128,127, data_type)
      CHECK_IVALS_ON_UB
    }
    else if (data_type == DFNT_UINT8) {
      uint8  *data = (uint8 *) buffer;
      uint8  ub, fillvalue;
      ASSIGN_INT_VALUE(a_ub,ub,0,256, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,0,256, data_type)
      CHECK_IVALS_ON_UB
    }
    else if (data_type == DFNT_INT16) {
      int16  *data = (int16 *) buffer;
      int16  ub, fillvalue;
      ASSIGN_INT_VALUE(a_ub,ub,-32768,32767, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,-32768,32767, data_type)
      CHECK_IVALS_ON_UB
    }
    else if (data_type == DFNT_UINT16) {
      uint16  *data = (uint16 *) buffer;
      uint16  ub, fillvalue;
      ASSIGN_INT_VALUE(a_ub,ub,0,65535, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,0,65535, data_type)
      CHECK_IVALS_ON_UB
    }
    else if (data_type == DFNT_INT32) {
      int32  *data = (int32 *) buffer;
      int32  ub, fillvalue;
      ub = (int32) atol(a_ub);
      fillvalue = (int32) atol(fvptr);
      CHECK_IVALS_ON_UB
    }
    else if (data_type == DFNT_UINT32) {
      uint32  *data = (uint32 *) buffer;
      uint32  ub, fillvalue;
      ASSIGN_INT_VALUE(a_ub,ub,0,2147483647, data_type);
      ASSIGN_INT_VALUE(fvptr,fillvalue,0,2147483647, data_type);
      CHECK_IVALS_ON_UB
    }
    else if (data_type == DFNT_FLOAT32) {
      float32  *data = (float32 *) buffer;
      float32  ub, fillvalue;
      ub = (float32) atof(a_ub);
      fillvalue = (float32) atof(fvptr);
      CHECK_FVALS_ON_UB
    }
    else if (data_type == DFNT_FLOAT64) {
      float64  *data = (float64 *) buffer;
      float64  ub, fillvalue;
      ub = atof(a_ub);
      fillvalue = atof(fvptr);
      CHECK_FVALS_ON_UB
    }
    else {
      return MODIS_S_OK;
    }
  }
  return MODIS_S_OK;
}

