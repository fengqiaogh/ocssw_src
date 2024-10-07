#ifndef HDF_LIB_H
#define HDF_LIB_H

#include    "Granule.h" /* for SMF_ERROR, include "mfhdf.h", "hdf.h", "SMF_PGS.h"*/

/*
!C-INC**********************************************************************
!Description:  Header file HDF_Lib.h for HDF_Lib.c. 
               Prototypes to read/write science-data-sets (multidimensional 
               arrays), attributes, etc, from/to HDF files.

!Revision History:
$Log: HDF_Lib.h,v $
Revision 1.3  2006-10-30 10:00:11-05  ltan
Changed for ANSI-C compliance. Correction for the generation of code change log.

 Revision 02.11, November 15, 1999
 Added functions:
 read_sds_rankn
 Check_Valid_Range
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 
 Revision 02.10 April 9, 1998
 Removed the prototypes for the following functions:
   read_sds_rank2_float64.c,
   read_sds_rank3_float32.c,
   read_sds_rank4_float32.c,
   read_sds_rank4_int16.c,
   write_sds_rank1_float32.c,
   write_sds_rank1_int8.c,
   write_sds_rank4_float32.c
 These aren't called anywhere in the L1B code, and we want
 to phase these multiple pointer functions out. If any new 
 data structures shouldn't use the multiple pointers.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)
 
 Revision 02.11 July  1998
 Removed read_part_sds_rank3_int16_1p(), added read_part_sds_rank3()
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)

 Revision 02.10 April 1998
 Added read_sds_rank2(), read_sds_rank3(), read_sds_rank4(), write_sds_rank1(), 
       write_sds_rank2(), write_sds_rank3(), write_sds_rank4();
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)
 
 Revision 01.20 Aug 1997
 Added read_vdata().
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.10 March 1997
 Added read_attributes and read_part_sds routines.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.00 Nov. 1996.
 Initial development.
 Zhidong Hao (hao@acrobat.gsfc.nasa.gov)

!References and Credits
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/

/* 
 * The following error messages are currently used in L1B_Tables.c and 
 * possibally used in other files in the future.
 */

extern char *invalidinputfile;
extern char *corruptinputfile;

PGSt_SMF_status read_vdata (int32   v_id,
                            int32   start_record,
                            int32   records,
                            char    *vname,
                            char    *fname,
                            void    *buffer);

PGSt_SMF_status read_attribute (int32   s_id,
                                char    *attr_name,
                                int32   TypeID,
                                void    *buffer);

PGSt_SMF_status read_part_sds_rank2 (int32   sd_id,
                                     char    *sds_name,
                                     int32   start0,
                                     int32   start1,
                                     int32   edge0,
                                     int32   edge1,
                                     void    *data);

PGSt_SMF_status read_part_sds_rank3 (int32   sd_id,
                                     char    *sds_name,
                                     int32   start0,
                                     int32   start1,
                                     int32   start2,
                                     int32   edge0,
                                     int32   edge1,
                                     int32   edge2,
                                     void    *data);

PGSt_SMF_status read_sds_rank1 (int32, char *, int32, void *);

PGSt_SMF_status read_sds_rank2 (int32, char *, int32, int32, void *);

PGSt_SMF_status read_sds_rank3 (int32, char *, int32, int32, int32, void *);

PGSt_SMF_status read_sds_rank4 (int32, char *, int32, int32, int32, int32, void *);

PGSt_SMF_status write_sds_rank1 (int32, char *, char *, int32, char *, void *);

PGSt_SMF_status write_sds_rank2 (int32, char *, char *, char *, int32, int32, char *,  void *);

PGSt_SMF_status write_sds_rank3 (int32, char *, char *, char *, char *, 
                                         int32, int32, int32, char *, void *);
PGSt_SMF_status write_sds_rank4 (int32, char *, char *, char *, char *, char *,
                                 int32, int32, int32, int32, char *, void *);

PGSt_SMF_status read_sds_rankn (int32     sd_id,
                                char      *sds_name,
                                int32     data_type,
                                int32     rank,
                                int32     *start,
                                int32     *edge,
                                void      *data);

PGSt_SMF_status write_sds_rankn (int32     file_id,
                                 char      *sds_name,
                                 int32     data_type,
                                 int32     rank,
                                 int32     *edge,
                                 char      **dim_name,
                                 void      *data);

PGSt_SMF_status Check_Valid_Range (char  *data_name,
                                   int32 data_type,
                                   char  *a_lb,
                                   char  *a_ub,
                                   char  *a_fillvalue,
                                   int32 count,
                                   void  *buffer);

#endif
 





