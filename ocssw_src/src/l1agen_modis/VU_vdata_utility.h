#ifndef   VDATA_UTILITY_H
#define   VDATA_UTILITY_H

/*
!C-INC************************************************************************

!Description:  This header file contains the macros for HDF vdata access.

!Input Parameters:
               N/A

!Output Parameters:
               N/A

Return Values: 
               N/A

Externally Defined:  
               MAX_NC_NAME             (mfhdf.h)
               MAX_VAR_DIMS            (mfhdf.h)

!Revision History:
$Log: VU_vdata_utility.h,v $
Revision 4.2  2003/03/10 16:10:12  vlin
Missing "]" fixed

Revision 4.1  2003/03/05 22:01:39  kuyper
Added Globals that are used in Vdata processing.

James Kuyper Jr. (kuyper@saicmodis.com)

Revision 1.0  1996/09/24  14:45 EDT
David Catozzi/GSC (cato@ltpmail.gsfc.nasa.gov)
Created original header file

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
	       HDF portions developed at the National Center for Supercomputing
	       Applications at the University of Illinois at Urbana-Champaign.

!Design Notes: The ".h" file below was specifically written for development
               in C. Any other language choice may require reworking of the
               ".h" file before coding can begin.

!END**********************************************************************
*/

#include   <string.h>
#include   "PGS_SMF.h"
#include   "mfhdf.h"
 
/***************************************************************/
/**                       MACROS                              **/
/***************************************************************/

#define   VU_MAX_NAME_LENGTH    MAX_NC_NAME    /* (max data set name length) */
#define   VU_MAX_RANK           MAX_VAR_DIMS   /* (max number of dimensions) */

#define   VU_MAX_DIM_NAME_LENGTH           40

#define   VU_MAX_DATA_TYPE_STRING_LENGTH    8

#define   VU_MAX_NUMBER_OF_VDATAS         200  

#define   VU_REPORT                         0   
#define   VU_INCREMENT                      1
#define   VU_DECREMENT                      2
#define   VU_EMPTY                         -1

#define  VU_NEW_VDATA                      -1
#define  VU_MAX_FIELD_NAME_LIST_SIZE     2000   /* true upper bound is 32,000
                                                (c.f. VSsetfields) */


/***************************************************************/
/**                       STRUCTURES                          **/
/***************************************************************/

typedef struct id_table {            /* used to store the Vdata IDs (returned */
    int32 id;                        /* by VSattach) so that functions need   */
    char  name[VU_MAX_NAME_LENGTH];  /* only refer to Vdatas by name          */
  } VU_ID_TABLE; 

/***************************************************************/
/**                         GLOBALS                           **/
/***************************************************************/

extern int32		global_H_ID;
extern VU_ID_TABLE	global_VU_VDATA_ID[VU_MAX_NUMBER_OF_VDATAS];
extern int		global_VU_ID_TABLE_READY;
#endif
