/*  BEGIN THIRD_PARTY_COPYRIGHT  */
/*==============================================================
The SZIP Science Data Lossless Compression Program is Copyright (C) 2001
Science & Technology Corporation @ UNM.  All rights released.
Copyright (C) 2003 Lowell H. Miles and Jack A. Venbrux.  Licensed to ICs, LLC,
for distribution by the University of Illinois' National Center for
Supercomputing Applications as a part of the HDF data storage and retrieval
file format and software library products package.  All rights reserved.
Do not modify or use for other purposes.

SZIP implements an extended Rice adaptive lossless compression algorithm for
sample data.  The primary algorithm was developed by R. F. Rice at Jet
Propulsion Laboratory.

SZIP embodies certain inventions patented by the National Aeronautics &
Space Administration.  United States Patent Nos. 5,448,642, 5,687,255,
and 5,822,457 have been licensed to ICs, LLC, for distribution with the
HDF data storage and retrieval file format and software library products.
All rights reserved.

Revocable, royalty-free, nonexclusive sublicense to use SZIP decompression
software routines and underlying patents is hereby granted by ICs, LLC, to
all users of and in conjunction with HDF data storage and retrieval file
format and software library products.

Revocable, royalty-free, nonexclusive sublicense to use SZIP compression
software routines and underlying patents for non-commercial, scientific use
only is hereby granted by ICs, LLC, to users of and in conjunction with HDF
data storage and retrieval file format and software library products.

For commercial use license to SZIP compression software routines and
underlying patents please contact ICs, LLC, at:
        ICs, LLC
        2600-A E. Seltice Way #234,
        Post Falls, ID 83854

        Phone: (208) 755-8990
        Fax: (208) 773-5747.

*************Begin NPOESS Specific Permission *******************
NPOESS Permissions: Email correspondence, dated February 6, 2007, from
Dr. Joseph J. Feeley, Chief Executive Officer, of ICs, LLC
2040 Warren Wagon Road, PO Box 2236, McCall ID 83638, USA, grants:
1. Raytheon Company is granted Royalty free permission to use this code
file (as modified) for the compression and decompression of all NPOESS data
packets.
2. Same permission is granted to the NPOESS Integrated Program Office,
all NPOESS contractors, and all U.S. Government agencies using NPOESS data.
  **
3. NPOESS data users may use SZIP decompression software routines and
underlying patents only in conjunction with decompressing NPOESS data in
accordance with and as limited by the above Third Party License.

4. All other rights reserved.  Do not modify or use for other purposes.
Further information available at http://www.hdfgroup.org/doc.resource/SZIP/
or by contacting ICs, LLC, at support@ICs4chips.com
*************End NPOESS Specific Permission *******************

==============================================================================*/
/*  END THIRD_PARTY_COPYRIGHT  */

/* BEGIN GOVERNMENT RIGHTS */
/**************************************************************************
* This code was partially developed under NASA, and other
* U.S. Government contracts.
*
* Patent/copyrights are retained as described in THIRD_PARTY_COPYRIGHT above.
*
* Please also see "LIMITATIONS" below.
*
****************************************************************************/
/* END GOVERNMENT RIGHTS */

/* BEGIN ITAR */
/* END ITAR CONTROL */

/* BEGIN CATEGORY */
/**************************************************************************
* CATEGORY: 1 (see SEN for definitions)
*
* rice.cpp
*
* Refer to Software Standards and Practices Manual (SSPM) for
* coding standards.
*
* NOTE: This code is exempt from NPOESS SSPM because it is Third Party,
* copyrighted code.  (Please see THIRD_PARTY_COPYRIGHT above)
*
**************************************************************************/
/* END CATEGORY */

/**************************************************************************
*
* NAME: rice.h
*
* REVISION/EVENT HISTORY:
* DATE        PR#      AUTHOR            Build    DESCRIPTION
* ---------   ---      ------            -----    -----------
* 29AUG2006            Miles              1.5     received source from L. Miles
*                                                at ICs Corp. Source is same as
*                                                U. of Illinois HDF Library.
* 25SEP2006            Arbeiter           1.5     inserted modifications
*                                               developed by Jennings, SBRS
* 27SEP2006            Arbeiter           1.5     added NPOESS headings
* 16FEB2007            Arbeiter           1.5     added NPOESS specific
*                                               permissions to 3rd Party
*                                               copyright Notices.
* 16APR2007            Arbeiter           1.5     Updated header comments
* 21DEC2007            Arbeiter           1.5     Updated header comments
*                                               and LIMITATIONS notes
* 31JAN2008            Arbeiter          1.5x1    Limatations updated
* 21NOV2008            Arbeiter          1.5x1    comment udpates
* ---------   ---      ------            -----    -----------
*
* REFERENCES: None.
*
*
* LIMITATIONS:
* 1. This code is not generalized Rice Compress/Uncompress.  It will only
* work with the VIIRS sensor data.
*
* 2. Use of this software for data from satellites other than NPP or NPOESS
* is prohibited.
*
* 3. This code has only been tested with the numbers of pixels in the VIIRS
* aggregation or compression zones, or the number of calibration pixels:
* 16, 48, 64, 96, 368, 488, 592, 640, 736, 760, 784, 1184, 1280, and 1776.
* Results with other numbers of pixels may be unreliable.
*
* 4. Raytheon has only tested this code with settings for 15 bits per pixel,
* 8 pixels per block, and option mask for NN, MSB, RAW, and CHIP.  Results with
* other settings may be unreliable. (These are the settings used by
* the VIIRS sensor for all compressions.)
*
* NOTES (MISCELLANEOUS) SECTION:
* Raytheon Company, Bellevue, Nebraska
*  2007 Raytheon Company.  All Rights Reserved.
*  2008 Raytheon Company.  All Rights Reserved.
*
* Comment marking "MJJ" is code written by Michael J. Jennings,
* of Raytheon SBRS, Goleta, California, 04/28/2004.
*
* Comment marking "RGA" is code written by Randolph G. Arbeiter,
* of Raytheon IIS, Omaha, Nebraska, September 2006.
*
*******************************************************************************/

/* ## start code added for VIIRS - RGA 14SEP2006  */
#ifndef _RICE_H
#define _RICE_H

#ifdef __cplusplus
extern "C"     /* prevents name dithering by C++ compiler RGA */
{
#endif
/* ## end code added for VIIRS - RGA 14SEP2006 */

#if !defined(__MWERKS__) 
typedef int boolean;
#endif

#define VIIRS_BLOCKS_PER_REFERENCE  128                            /* RGA */

#define FALSE   0
#define TRUE   1

#define eq(a, b) (!strcmp((a), (b)))
#define eqn(a, b, n) (!strncmp((a), (b), (n)))
#define MIN(x,y) ((x)<(y)? (x): (y)) 

#define EC_MODE   0
#define NN_MODE   1

#define DEFAULT_BITS_PER_PIXEL         8
#define DEFAULT_BLOCKS_PER_SCANLINE    32
#define DEFAULT_PIXELS_PER_BLOCK    16
#define DEFAULT_PIXELS_PER_SCANLINE  (DEFAULT_BLOCKS_PER_SCANLINE)*(DEFAULT_PIXELS_PER_BLOCK)

#define MAX_EXT2                   7 
#define MAX_EXT2_SUM              (MAX_EXT2*(MAX_EXT2+1)/2 + MAX_EXT2)

#define MAX_COMMAND_LINE_FILES      256
#define MAX_FILENAME_SIZE          256

#define MAX_ZERO_BLOCKS                64   /*** Must be a power of two ***/

  //#define MAX_BLOCKS_PER_SCANLINE      128
#define MAX_BLOCKS_PER_SCANLINE      256
#define MAX_PIXELS_PER_BLOCK        32
#define MAX_PIXELS_PER_SCANLINE     (MAX_BLOCKS_PER_SCANLINE)*(MAX_PIXELS_PER_BLOCK)

#define ID_ZERO         -1
#define ID_LOW           0
#define ID_FS            1
#define ID_K1          2
#define ID_K2          3
#define ID_K3          4
#define ID_K4          5
#define ID_K5          6
#define ID_K6          7
#define ID_K7          8
#define ID_K8          9
#define ID_K9         10
#define ID_K10         11
#define ID_K11         12
#define ID_K12         13
#define ID_K13         14
#define ID_K14         15
#define ID_K15         16
#define ID_K16         17
#define ID_K17         18
#define ID_K18         19
#define ID_K19         20
#define ID_K20         21
#define ID_K21         22
#define ID_K22         23
#define ID_K23         24
#define ID_DEFAULT      31

#define ID_DEFAULT1       7
#define ID_DEFAULT2      15
#define ID_DEFAULT3      31

#define K_FACTOR      1

#define FILE_DATA   1
#define MEMORY_DATA   2
/* This doesn't work when scaline is set to max value 4K EIP 2004/8/05
#define INPUT_BUFFER_SIZE    16*1024
#define OUTPUT_BUFFER_SIZE  16*1024
*/
//#define INPUT_BUFFER_SIZE    2*MAX_PIXELS_PER_SCANLINE*4
//#define OUTPUT_BUFFER_SIZE      2*MAX_PIXELS_PER_SCANLINE*4 

//JMG  10/17/14
#define INPUT_BUFFER_SIZE   65536*2
#define OUTPUT_BUFFER_SIZE  65536*2

static struct
   {
   int bits_per_pixel[8];
   int pixels_per_block[8];
   int pixels_per_block_mult[16];
   int scanlines_per_file[128];
   } short_header =
{
   {7, 8, 9, 10, 12, 14, 15, 16},
   {8, 10, 12, 16, 18, 20, 24, 32},
   {1, 2, 4, 6, 8, 10, 12, 16, 18, 20, 24, 32, 34, 36, 40, 48},
   {   1,    2,    4,    8,   16,   32,   64,  128,  256,  512, 1024, 2048, 4096, 
        3,    6,   12,   24,   48,   96,  192,  384,  768, 1536, 3072, 
        5,   10,   20,   40,   80,  160,  320,  640, 1280, 2560, 
        9,   18,   36,   72,  144,  288,  576, 1152, 2304, 
       17,   34,   68,  136,  272,  544, 1088, 2176,
       25,   50,   75,  100,  125,  150,  175,  200,  225,  250,
     275,  300,  325,  350,  375,  400,  425,  450,  475,  500,
     550,  600,  650,  700,  750,  800,  850,  900,  950, 1000,
    1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500,
    1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000,
    2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
     3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000}
};

#define ALLOW_K13_OPTION_MASK        1
#define CHIP_OPTION_MASK             2 
#define EC_OPTION_MASK              4
#define LSB_OPTION_MASK                8
#define MSB_OPTION_MASK             16
#define NN_OPTION_MASK             32
#define OVERWRITE_OPTION_MASK       64
#define RAW_OPTION_MASK            128
#define KEEP_IMAGE_OPTION_MASK      256
#define KEEP_COMPRESSED_OPTION_MASK   512

#define SZIP_PROGRAM_NAME      "szip"
#define SUNZIP_PROGRAM_NAME      "sunzip"

#define MEMORY_ERROR      (-2)
#define PARAM_ERROR         (-4)
#define NO_ENCODER_ERROR   (-5)

/* ## start code added for VIIRS - RGA 14SEP2006  */
#ifdef __cplusplus
}
#endif

#endif
/* ## end code added for VIIRS - RGA 14SEP2006  */
