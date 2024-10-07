#ifndef _SWL0_STRUC_H
#define _SWL0_STRUC_H

#include <stdio.h>
#include "swl0_parms.h"
#include "swl0_types.h"

typedef struct swl0ctl_struct {
    INT16 fileType;
    INT16 timerangeFactor;
    INT16 maxBitErrors;
    INT16 gainSetting;
    INT16 stopTimeDelta;
    INT16 env;
    char *stationInfoFile;
    char *progname;
} swl0ctl;

typedef struct swl0hdr_struct {
    BYTE id[5]; /* Should be ascii CWIF\0                    */
    BYTE type; /* 0 = LAC/GAC, 1 = HRPT                     */
    BYTE fill1[310]; /* Not used                                  */
    INT32 numbits; /* Number of sync bits checked by formatter  */
    INT32 errbits; /* Number of error bits counted by formatter */
    BYTE fill2[8]; /* Not used                                  */
    INT32 startTime; /* Time of AOS in seconds since 1/1/70 00:00 */
    INT32 stopTime; /* Time of LOS in seconds since 1/1/70 00:00 */
    INT32 numlines; /* Number of frames captured by formatter    */
    BYTE fill3[168]; /* Not used                                  */
} swl0hdr;

typedef struct frameqc_struct {
    FLOAT64 time; /* Minor frame time (secs since 1/1/70 00:00)*/
    FLOAT64 timeShift; /* Minor frame time shift (secs)             */
    BYTE mnftype; /* Minor frame type (0-4,15)                 */
    BYTE mnfnum; /* Minor frame number (1-3)                  */
    BYTE tRanError; /* 1 => time out of range                    */
    BYTE tSeqError; /* 1 => time not increasing with frame type  */
    BYTE tDifError; /* 1 => time difference not expected multiple*/
    BYTE tShfError; /* 1 => time shifted relative to multiple    */
    BYTE scidError; /* 1 => mnfnum or mnftype out of range       */
    BYTE sgaError; /* 1 => bit errors in sga header             */
    BYTE sacError; /* 1 => bit errors in sac header             */
    BYTE saaError; /* 1 => bit errors in saa header             */
    BYTE bitError; /* 1 => bit errors in sync or image bits     */
    INT16 numBits; /* number of bits checked                    */
    INT16 errBits; /* number of bit errors found                */
    BYTE maxErrBits; /* 1 => max bit errors exceeded              */
    INT32 pixVariance; /* total of pix-to-pix differences, absolute */
} frameqc;

typedef struct swl0indx_struct {
    char l0file[FILENAME_MAX]; /* Full path name of L0 file       */
    BYTE type; /* File type (HRPT=1, LAC/GAC=0)             */
    INT16 nrecs; /* Total number of records in L0 file        */
    INT16 srec; /* First valid frame                         */
    INT16 erec; /* Last valid frame                          */
    INT16 ngac; /* Number of GAC records                     */
    INT16 nlac; /* Number of GAC records                     */
    INT16 nlun; /* Number of lunar cal records               */
    INT16 nsol; /* Number of solar cal records               */
    INT16 nigc; /* Number of intergain cal records           */
    INT16 ntdi; /* Number of TDI cal records                 */
    INT16 nund; /* Number frames with undefined typecode     */
    INT16 timeError; /* Number of Timing Errors                   */
    INT16 scidError; /* Number of SCID errors                     */
    INT16 sohError; /* Number of SOH header errors               */
    FLOAT64 bitRateRep; /* Bit rate reported by frame formatter      */
    FLOAT64 bitRateComp; /* Bit rate computed from image start pixel  */
    INT32 timeAOS; /* Time of AOS in seconds since 1/1/70 00:00 */
    INT32 timeLOS; /* Time of LOS in seconds since 1/1/70 00:00 */
    FLOAT64 timeFirst; /* Time of first valid scanline              */
    FLOAT64 timeLast; /* Time of last  valid scanline              */
    frameqc *rec; /* Frame-by-frame quality info               */
} swl0indx;

typedef struct swl0scene_struct {
    char l0file[FILENAME_MAX]; /* Full path name of L0 file       */
    BYTE type; /* File type (HRPT=1, LAC/GAC=0)             */
    BYTE mnftype; /* Minor frame type (LAC=0, SOL=2, GAC=15)   */
    INT16 nrec; /* Number of frames in scene                 */
    INT16 nscan; /* Number of scans  in scene                 */
    INT16 srec; /* First L0 frame number of scene            */
    INT16 crec; /* Central L0 frame number of scene          */
    INT16 erec; /* Last L0 frame number of scene             */
    FLOAT64 stime; /* Start time in secs since 1/1/70 00:00     */
    FLOAT64 ctime; /* Central time in secs since 1/1/70 00:00   */
    FLOAT64 etime; /* End time in secs since 1/1/70 00:00       */
    INT32 orbnum; /* Orbit number at center of scene           */
    INT16 indx[MAXFRAMES]; /* L0 file frame numbers for this scene*/
    INT16 qual[MAXFRAMES]; /* Quality index per scene frame       */
    /*                                                               */
    /* The rest of this stuff can only be filled after navigation    */
    /*                                                               */
    INT32 center_scan_line; /* Scene center-most scan w valid nav  */
    FLOAT32 node_lon; /* Longitude of dsc node crossing (deg)      */
    FLOAT64 node_time; /* Time of dsc node xing in secs since 1/1/70*/
    INT32 ntilts; /* Number of tilt states in scene            */
    INT16 tilt_flags[MAXTILTS]; /* 2=aft, 3=changing, 1=fwd   */
    INT16 tilt_ranges[MAXTILTS][2]; /* Scan num range per state   */
    FLOAT32 tilt_lons[MAXTILTS][2][2]; /* Lon of tilt period corners */
    FLOAT32 tilt_lats[MAXTILTS][2][2]; /* Lat of tilt period corners */
    FLOAT32 center_lat;
    FLOAT32 center_lon;
    FLOAT32 center_solz;
    FLOAT32 upper_left_lat;
    FLOAT32 upper_left_lon;
    FLOAT32 upper_right_lat;
    FLOAT32 upper_right_lon;
    FLOAT32 lower_left_lat;
    FLOAT32 lower_left_lon;
    FLOAT32 lower_right_lat;
    FLOAT32 lower_right_lon;
    FLOAT32 northern_lat;
    FLOAT32 southern_lat;
    FLOAT32 western_lon;
    FLOAT32 eastern_lon;
    FLOAT32 start_center_lat;
    FLOAT32 start_center_lon;
    FLOAT32 end_center_lat;
    FLOAT32 end_center_lon;
    char start_node[11];
    char end_node[11];

} swl0scene;

typedef struct time_table_struct {
    char time1[17];
    char time2[17];
    INT32 type;
    INT32 action;
    INT32 delmsec;
} timeTab;

typedef struct swl0mindx_struct {
    char filename[FILENAME_MAX];
    INT32 filenum;
    INT32 framenum;
    frameqc *qual;
} swl0mindx;


#endif
