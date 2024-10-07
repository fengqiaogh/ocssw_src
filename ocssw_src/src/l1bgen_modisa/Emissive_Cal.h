
#ifndef EMISSIVE_CAL_H
#define EMISSIVE_CAL_H


#include "L1B_Tables.h" /*use emiss_tables_t and it also includes Granule.h, in 
                          which the L1B_granule_t and L1B_Scan_t ... are used
                          in Emissive_Cal.c */
#include "L1B_Setup.h"  /*use L1B_Scan_Metadata_t data type */
/*
!C-INC***********************************************************************
!Description:   Public header file for externals in Emissive_Cal.c

!Revision History:
 $Log: Emissive_Cal.h,v $
 Revision 1.3  2006-10-27 11:01:33-04  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.

 Revision 01.00 Nov 1998  
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST) for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/

PGSt_SMF_status  Emissive_Cal (int16                S,
                               int16                MS,
                               L1B_granule_t        *L1B_Gran,
                               L1A_Scan_t           *L1A_Scan,
                               L1B_Scan_t           *L1B_Scan,
                               Preprocess_Emiss_t   *PP_emiss,
                               emiss_tables_t       *emiss_tables,
                               QA_tables_t          *QA_tables,
                               QA_Data_t            *QA,
                               L1B_Scan_Metadata_t  *L1B_Scan_Meta);

#endif









