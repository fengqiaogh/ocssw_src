#ifndef FNAMES_H
#define FNAMES_H
/***********************************************************************
!C-INC

!Description:
   Header file containing macros for PCF logical unit numbers.

!Revision History:
   $Log: FNames.h,v $
   Revision 1.10  2006-10-30 10:01:29-05  ltan
   Defined PGE-2_VERSION_LUN with the value of 800500. Changed for ANSI_C complianec. Correction for the generation of code change log.
 
   Revision 01.06, October 31, 2003  Razor Issue #195
   Added LUN for Processing Center
   Liqin Tan, SAIC GSO (ltan@saicmodis.com)

   Revision 01.05, March 21, 2002  Razor Issue #181
   Added LUN for Processing Environment
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)  
  
   Revision 01.04, November 12, 2001 (Razor Issue #169)
   Added LUN for Writing Night mode 250m and 500m data, WRITE_NIGHT_HIRES_LUN
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.03, November 6, 2001  (Razor Issue #167)
   Added LUN for MCST LUT Version Number, MCST_LUT_VERSION_LUN
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.02, November 21, 2000
   Added LUN for Satellite Instrument
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.01 1996/04/05
   Update to match Version 1 Design Document
   Joan Baden (baden@highwire.gsfc.nasa.gov)
   John Hannon(hannon@highwire.gsfc.nasa.gov)

   Revision 01.00 1993
   Initial development
   Geir E. Kvaran(geir@highwire.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   
!END********************************************************************
*/

/*---------------------------------------
   PGS reference numbers 
     L1A:    500,000 - 599,999
     Geo:    600,000 - 699,999
     L1B:    700,000 - 799,999
----------------------------------------*/

#define MCF_FILE_QKM        700250
#define MCF_FILE_HKM        700251
#define MCF_FILE_1KM        700252
#define MCF_FILE_OBC        700253

#define L1B_EV_250M_FILE                 700000
#define L1B_EV_500M_FILE                 700001
#define L1B_EV_1000M_FILE                700002
#define L1B_OBC_FILE                     700010

#define REFLECTIVE_TABLES_FILE           700050
#define EMISSIVE_TABLES_FILE             700060
#define QA_TABLES_FILE                   700070

#define LEADING_L1A_GRANULE              500000
#define FIRST_L1A_GRANULE                500001
#define TRAILING_L1A_GRANULE             500002

#define GEOLOCATION_FILE                 600000

#define PGE02_VERSION_LUN                800500
#define SATELLITE_INSTRUMENT_LUN         800510
#define PROCESSING_ENVIRONMENT_LUN       800550
#define REPROCESSING_PLANNED_LUN         800600
#define REPROCESSING_ACTUAL_LUN          800605
#define MCST_LUT_VERSION_LUN             800610
#define WRITE_NIGHT_HIRES_LUN            800615
#define PROCESSING_CENTER_LUN            800620


#endif
