#ifndef GRANULEP_H
#define GRANULEP_H

#include "Granule.h"

/*
!C-INC***********************************************************************
!Description:  Private header file, to be used in Granule.c only. Contains 
               function prototypes.

!Revision History:
 $Log: GranuleP.h,v $
 Revision 1.3  2006-10-30 10:00:11-05  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.

   Revision 01.02, October 29, 2000
   Moved Determine_Split_Scans and Get_Split_Scan_Indexes from L1B.c and
   removed adjusting the QA counters. (issue 142)
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.01 Feb 2, 1999
   Separated out parts belonging to L1B_Setup as part of the
   process of separating L1B_Setup out into its own module.
   Deleted obsolete (unused) declaration.
   Jim Rogers (rogers@msct.gsfc.nasa.gov)

   Revision 01.00 Nov. 1998
   This is part of original Granule.h, L1B.h
   Initial development
   Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

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

/*---------------------------------------------------------------
  Prototypes
---------------------------------------------------------------*/

PGSt_SMF_status  Compute_Aggregates (int16     scale,
                                     int16     line_dim,
                                     int16     pixel_dim,
                                     uint16    *SI_in,
                                     uint8     *UI_in,
                                     uint16    *SI_out,
                                     uint8     *UI_out,
                                     int8      *SU_out);

PGSt_SMF_status Write_L1B_SI_UI (int16 S, L1B_Scan_t *L1B_Scan, int16 R);

#endif





