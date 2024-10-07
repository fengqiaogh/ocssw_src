#ifndef PC_PCF_INFO_H
#define PC_PCF_INFO_H
#include "PGS_PC.h"

/*
!C-INC************************************************************************

!Description:  This header file contains the Process Control File Information.

!Input Parameters:  N/A

!Output Parameters: N/A

!Revision History:
 $Log: PC_pcf_info.h,v $
 Revision 5.1  2004/09/30 18:54:46  seaton
 Updated to run Collection 5 L1A code.
 seaton@saicmodis.com

 Revision 4.2  2002/12/16 20:24:37  vlin
 global_L0_logical added and macro "PC_L1A_PROCESSINGENVIRONMENT" defined.
 vlin@saicmodis.com

 Revision 4.1  2002/08/28 18:51:45  vlin
 PCF_CONFIG_t declaration added.

               Revision 2.0  2000/07/21
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Added #defines for instrument data

               Revision 1.0  1997/09/24  14:45 EDT
               Timi Adelekan/GSC/SAIC (adelekan@ltpmail.gsfc.nasa.gov)
               Created include file from Version 1 include files:
               modis_init.h,

!Team-unique Header:

        This software is developed by the MODIS Science Data Support Team 
        for the National Aeronautics and Space Administration, 
        Goddard Space Flight Center, under contract NAS5-32373.

 Design Notes: 

        The ".h" file below was specifically written for development
        in C. Any other language choice may require reworking of the
        ".h" file before coding can begin.

!END**********************************************************************
*/

#define	PC_PRIOR_L0_PCF_ID		599001
#define	PC_CURRENT_L0_PCF_ID		599002
#define PC_L1A_ENG_DATA_LIST_FILE	599003

#define PC_L1A_FIRST_GRAN_START_TIME	10258
#define PC_L1A_LAST_GRAN_STOP_TIME	10259
#define PC_L1A_LUT_REVISION             599004  /* RCS Revision number */

#define PC_L1A_GRAN_TIME_LENGTH	        503000
#define PC_L1A_PGE_VERSION              800500
#define PC_L1A_PROCESSING_ENVIRONMENT   800550
#define PC_L1A_REPROCESS_ACTUAL         800600
#define PC_L1A_REPROCESS_PLANNED        800605
#define PC_L1A_SCAN_RATE		504000
#define PC_L1A_VERSION  		505000

#define PC_L1A_PROD_PCF_ID              500100
#define PC_MCF_PCF_ID                   500500
#define PC_INSTRUMENT                   800510
#define PC_L1A_PROCESSINGENVIRONMENT    "ProcessingEnvironment"
#define INSTRUMENT_TERRA                "AM1M"
#define INSTRUMENT_AQUA                 "PM1M"
#define TERRA_SCID                      42
#define AQUA_SCID                       154


extern  PGSt_PC_Logical  global_L0_logical;
                         /* L0 logical file unit number */

typedef struct {
        PGSt_double first_gran_start_time;
        PGSt_double last_gran_stop_time;
        PGSt_double gran_time_length;
        PGSt_double scan_rate;
        char localversionid[PGSd_PC_VALUE_LENGTH_MAX];
        char lutrevision[PGSd_PC_VALUE_LENGTH_MAX];
        char pgeversion[PGSd_PC_VALUE_LENGTH_MAX];
        char processingenvironment[PGSd_PC_VALUE_LENGTH_MAX];
        PGSt_tag instrument;
        char reprocessingactual[PGSd_PC_VALUE_LENGTH_MAX];
        char reprocessingplanned[PGSd_PC_VALUE_LENGTH_MAX];
} PCF_CONFIG_t;

#endif /* PC_PCF_INFO_H */
