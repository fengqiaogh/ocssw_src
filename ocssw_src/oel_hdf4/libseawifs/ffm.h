/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.5/L012_Util/util/data/ffm.h,v 4.13 1995/12/07 15:51:28 seawifsd Exp seawifsd $
$Log: ffm.h,v $
Revision 4.13  1995/12/07 15:51:28  seawifsd
added macro definitions about on-board recorder sizes used with
new MOPS command scheduling. The values are not completely settled
because different sets of values were shown up in MOPS schedule files.

Revision 4.12  1995/08/03 20:46:04  seawifsd
added definition of FFM_RAW_LEN.

Revision 4.11  1995/05/04 15:10:17  seawifsd
changed the way to set macros PCT_MSEC_ERR, GAC_MSEC_ERR, and
LAC_MSEC_ERR so that it can be set in the compiling time without
edit the file.

Revision 4.10  1995/01/17 19:57:55  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:14:13  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:46:20  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:18  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.3  1994/11/03 19:42:07  frank
changed FFM_GAP_MAX from 150 to 30(minorframes). Original 150 is using the
unit of 'GAC scan line'.
changed typo on 'LAC_FFM_GAP_MAX'(from LAC_FFM_GAC_MAX)

Revision 1.1.1.2  1994/11/03 19:39:37  frank
defined GAC_FFM_MSEC_INC, LAC_FFM_MSEC_INC, PCT_MSEC_ERR, GAC_MSEC_ERR
LAC_MSEC_ERR, FFM_GAC_MAX, GAC_FFM_GAP_MAX, LAC_FFM_GAP_MAX, G
GAC_MSEC_GAP_MAX and LAC_MSEC_GAP_MAX.

Revision 1.1.1.1  1994/10/04 15:21:18  frank
added definition of GAC_PIX_START,GAC_PIX_SUB, LAC_PIX_START, and LAC_PIX_SUB
added DARK_OFF and used the macro in the definition of LAC_DARK_OFF and GAC_DARK_OFF.

Revision 1.2  1994/05/10 18:47:25  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:29:55  seawifst
Initial revision


 */


#ifndef FFM_
#define FFM_

#define byte unsigned char

/* short integer, number of synch bits used for error rate */
#define QF1_OFF 0
#define QF1_LEN 2
#define QF1_SCALE_FACTOR 1/5
/* byte, number of bit errors in synch bits */
#define QF2_OFF 2
#define QF2_LEN 1

#define ID_OFF 3
#define ID_LEN 4
#define TAG_OFF 7
#define TAG_LEN 8
#define SOH_OFF 15
#define SOH_LEN 775
#define TLM_LEN 88
#define TDI_LEN 16
#define START_LEN 16
#define DARK_LEN 16
#define STOP_LEN 16

#define DARK_OFF 32
/* LAC */
#define LAC_TLM_OFF 790 
#define LAC_TLM_LEN TLM_LEN 
#define LAC_OFF 878
#define LAC_LEN 20624
#define  LAC_TDI_OFF 0
#define  LAC_TDI_LEN TDI_LEN
#define  LAC_START_OFF 16
#define  LAC_START_LEN START_LEN
#define  LAC_DARK_OFF DARK_OFF
#define  LAC_DARK_LEN DARK_LEN
#define  LAC_IMAGE_OFF 48
#define  LAC_IMAGE_LEN 20560
#define  LAC_STOP_OFF 20608
#define  LAC_STOP_LEN STOP_LEN
/* END LAC */

/* GAC */
#define GAC_OFF 790
#define GAC_LEN 4032
#define  GAC_TDI_OFF 0
#define  GAC_TDI_LEN TDI_LEN
#define  GAC_START_OFF 16
#define  GAC_START_LEN START_LEN
#define  GAC_DARK_OFF DARK_OFF
#define  GAC_DARK_LEN DARK_LEN
#define  GAC_IMAGE_OFF 48
#define  GAC_IMAGE_LEN 3968
#define  GAC_STOP_OFF 4016
#define  GAC_STOP_LEN STOP_LEN

#define GAC_TLM_OFF 20950
#define GAC_TLM_LEN TLM_LEN
/*
        For each GAC minorframe, following equations show how to access
        each GAC image and the corresponding Inst/Anc Telemetry.
        for(i=0; i < 5 ; i++) {
          memcpy(image,&rec[GAC_OFF+i*GAC_LEN+GAC_IMAGE_OFF],GAC_IMAGE_LEN);
          memcpy(tlm,&rec[GAC_TLM_OFF+i*GAC_TLM_LEN],GAC_TLM_LEN);
        }
 */

#define GAC_SPR_OFF 21390
#define GAC_SPR_LEN 112

/* END GAC */
#define SPR_OFF 21502
#define SPR_LEN 2
/* total length */
#define FFM_LEN 21504
#define FFMRECLEN FFM_LEN
#define FFM_RAW_LEN 13860

/* RECORDER_SIZE was changed from 62500 to 65536 because the 62500 is	*/
/* using MB unit instead of million bytes used in other places		*/
/* 65536 = 62500 * 1.024 * 1.024					*/
#define RECORDER_SIZE 65536
#define RECORDER2_LIMIT 59990
#define RECORDER_TLM_SIZE (RECORDER_SIZE - RECORDER2_LIMIT)
#define ALLOCATE_GAC_SIZE(x) ((x * FFM_RAW_LEN + 500) / 1000 - RECORDER_SIZE)
#define ALLOCATE_LAC_SIZE(x) ((x * FFM_RAW_LEN + 500) / 1000)

/* FFM header block size */
#define FFMHDRLEN 512
/* Frame formatter missing frames */
#define FF_MISSING_FRAMES 1

#define GAC_PER_FFM 5

#define BANDS 8

#define PIXEL_BLEN BANDS*2
#define PIXEL_WLEN BANDS
/*
   define byte offset from the beginning offset of each gac/lac segment
   which is equal to the size for START_SYNC and DARK_RESTORE and should
   be equal to 32 bytes.
 */
#define SCI_DATA_BOFF  2*PIXEL_BLEN
#define SCI_DATA_WOFF  2*PIXEL_WLEN
#define START_SYNC_OFF  0
#define DARK_RESTORE_OFF 1
#define GAC_PIXEL_NUM  248
#define LAC_PIXEL_NUM  1285
#define MAX_PIXEL_VALUE  1023
#define MIN_PIXEL_VALUE  0
#define SATURATED_VALUE  1023

#define LAC_PIX_START 1
#define GAC_PIX_START 147
#define LAC_PIX_SUB 1
#define GAC_PIX_SUB 4

typedef struct ffm_hdr_Struc {
    byte hdr_byte[FFMHDRLEN];
} ffm_hdr_Type;

typedef union ffm_rec_Struc {
    byte rec_byte[FFMRECLEN];
} ffm_rec_Type;

#define MSEC_PER_DAY 86400000
#define GAC_MSEC_INC (1000.0 * 4.0/6.0)
#define LAC_MSEC_INC (1000.0/6.0)
#define GAC_FFM_MSEC_INC ((GAC_MSEC_INC) * (GAC_PER_FFM))
#define LAC_FFM_MSEC_INC (LAC_MSEC_INC)
/*
   Might not using the percentage tolerance. 1% is 1.667 msec for LAC
   and 6.667 msec for GAC.
 */
#ifndef PCT_MSEC_ERR
#define PCT_MSEC_ERR 1
#endif /* !PCT_MSEC_ERR */
#ifndef GAC_MSEC_ERR
#define GAC_MSEC_ERR 1
#endif /* !GAC_MSEC_ERR */
#ifndef LAC_MSEC_ERR
#define LAC_MSEC_ERR 1
#endif /* !LAC_MSEC_ERR */
/*
   Maximum number of minorframe gap are the same(for now) for both GAC
   and LAC. But they can be set to different number.
 */
#define FFM_GAP_MAX 30
#define GAC_FFM_GAP_MAX FFM_GAP_MAX
#define LAC_FFM_GAP_MAX FFM_GAP_MAX
#define GAC_MSEC_GAP_MAX ((GAC_FFM_GAP_MAX) * (GAC_FFM_MSEC_INC))
#define LAC_MSEC_GAP_MAX ((LAC_FFM_GAP_MAX) * (LAC_FFM_MSEC_INC))

#define REVERSE_MSEC_GAP 1.0
#ifndef ERR_UTIME
#define ERR_UTIME  REVERSE_MSEC_GAP
#endif /* !ERR_UTIME */

/*
   this value is used to be filled into the 2 spare byte at the end of
   each minorframe. OFFSET = SPR_OFF
   Value 89(59h) has bit pattern of 01011001
 */
#define FILL_FRAME_PAT 89
#define REAL_FRAME_PAT 0

#define GOOD_FFM 0
#define FILL_FFM 1

/*
   define sub-item within each FFM block
 */

typedef struct ffmStruct {
    char buf1[GAC_TLM_OFF];
    short int TLM1[44];
    short int TLM2[44];
    short int TLM3[44];
    short int TLM4[44];
    short int TLM5[44];
    char buf2[114];
} gacffmType;
#endif /* FFM_ */
