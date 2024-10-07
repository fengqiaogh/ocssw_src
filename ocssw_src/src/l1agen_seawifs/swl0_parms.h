#ifndef _SWL0_PARMS_H
#define _SWL0_PARMS_H

#define L01VERSION "4.5"
#define FATAL_ERROR 1

#define NBANDS     8        /* Number of spectral bands           */
#define L0LEN      21504    /* L0 frame length in bytes           */
#define MAXFRAMES  16000    /* Maximum number of frames per file  */
#define MAXSCENES  100      /* Maximum number of scenes per file  */
#define MAXTILTS   20       /* Maximum number of tilts per scene  */
#define MINORBVEC  0        /* Minimum number of GPS vectors fit  */
#define MINFRAMES  50       /* Minimum number of frames per scene */
#define MAXORBTAB  1000     /* Maximum orbit table entries        */

#define  GAC       0        /* L0 file type code for GAC          */
#define  LAC       0        /* L0 file type code for LAC          */
#define  HRPT      1        /* L0 file type code for HRPT         */
#define  GACTYPE   15       /* Minor frame type code for GAC      */
#define  LACTYPE   0        /* Minor frame type code for LAC      */
#define  LUNTYPE   1        /* Minor frame type code for LUN Cal  */
#define  SOLTYPE   2        /* Minor frame type code for SOL Cal  */
#define  IGCTYPE   3        /* Minor frame type code for IGN Cal  */
#define  TDITYPE   4        /* Minor frame type code for TDI Cal  */
#define  NPIXLAC   1285     /* Number of LAC pixels per scan      */
#define  NPIXGAC   248      /* Number of GAC pixels per scan      */
#define  SPIXLAC   1        /* First LAC pixel in  scan           */
#define  SPIXGAC   147      /* First GAC pixel in  scan           */
#define  IPIXLAC   1        /* LAC pixel increment in  scan       */
#define  IPIXGAC   4        /* GAC pixel increment in  scan       */

#define  O_SCID    3        /* Byte-offset to SCID field          */
#define  O_TIME    7        /* Byte-offset to frame timetag       */
#define  O_SOH     15       /* Byte-offset to SOH Block           */
#define  O_SGA     15       /* Byte-offset to SOH/SGA             */
#define  O_SAC     185      /* Byte-offset to SOH/SAC             */
#define  O_SAA     481      /* Byte-offset to SOH/SAA             */
#define  O_DATAGAC 790      /* Byte-offset to GAC image block     */
#define  O_DATALAC 878      /* Byte-offset to LAC image block     */
#define  O_INSTGAC 20950    /* Byte-offset to GAC instrument tlm  */
#define  O_INSTLAC 790      /* Byte-offset to LAC instrument tlm  */

#define  S_DATAGAC 4032     /* Byte-size of GAC image block       */
#define  S_DATALAC 20624    /* Byte-size of LAC image block       */

#define DTGAC      (10./3.) /* Nominal GAC time difference (sec)  */
#define DTLAC      (1./6.)  /* Nominal LAC time difference (sec)  */
#define CLOCKVAR   0.010    /* Nominal S/C clock variance  (sec)  */

#define DELSCENELAC 5       /* Sec thresh to break LAC scenes     */
#define DELSCENEGAC 3000    /* Sec thresh to break GAC scenes     */

#define SEADAS 0
#define SDPS   1


/*
 * NOTES: 
 * 
 * The image block includes GTDI-PIX, START-PIX, DARK-PIX, IMAGE-PIX,
 * and STOP-PIX.
 *
 */
#endif








