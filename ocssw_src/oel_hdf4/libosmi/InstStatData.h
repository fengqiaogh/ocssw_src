/*
 * Instrument Status Datafile header definition
 *
 * Copyright @1998 Datron/Transco, Inc.
 *
 * Originator: P.G. Heffernan
 *
 * "%Z% %M% %I% - %E%"
 *
 * Description:  This include file provides a set
 *	of structures which identify the data
 *	contained in an Instrument Status Data (ISD) file
 *	produced by the KOMPSAT Mission Control Station.
 *	The structure identifies the data contained
 *	in one record contained in an ISD per the
 *	KOMPSAT MCS-IRPE Interface Control Document.
 *
 * Version	Date	Programmer	Description
 * Orig		10/8/98	PG Heffernan	Original
 * 1.1          12/1/99 WD Benton	ISD format change, fix interpretation
 *					of OSMI Band limits
 *
 */
#ifndef ISD_H_DEFINED
#define ISD_H_DEFINED

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Library Defines
 */

/* Sizes */
#define ONE_BIT_SIZE 2
#define TWO_BIT_SIZE 4
#define THREE_BIT_SIZE 8
#define FOUR_BIT_SIZE 16
#define FIVE_BIT_SIZE 32
#define SIX_BIT_SIZE 64
#define SEVEN_BIT_SIZE 128
#define EIGHT_BIT_SIZE 256

/* Time */
#define NUM_SEC_PER_MIN 60.0
#define NUM_MIN_PER_HOUR 60.0
#define NUM_HOUR_PER_DAY 24.0
#define FRAC_TIME_DIV 16777216.0

/* ISD Specifics */
#define INST_STAT_DATA_FILE_WC "INSTR_%04s%02s%02s*_S.DAT"
#define ISD_RECORD_SIZE 204
#define ISD_SPACECRAFT_TIME_OFFSET 0
#define ISD_GPS_TIME_SIZE 26
#define ISD_GPS_TIME_OFFSET 11

/*
 * Function Prototypes
 */

/*
 * Structure Declarations
 */
struct isd_spacecraft_time {
    int quarter_sec_week; /* Quarter Seconds of Week, 0-2,419,199 */
    int weeks; /* Weeks 0-4,095 */
};

struct isd_glob_pos_syst_time {
    int year; /* Self explanatory */
    int month;
    int day;
    int hour;
    int minute;
    int second;
    int frac_second;
};

struct inst_stat_data_misc {
    char cover_position[14]; /* Cover position, No units */
    int eoc_prim_pwr_stat; /* EOC Primary Power Status, 0:Off, 1:On */
    int eoc_red_pwr_stat; /* EOC Redundant Power Status, 0:Off, 1:On */
    int eoc_pri_htr_pwr_stat; /* EOC Primary Heater Power Status, 0:Off, 1:On */
    int eoc_red_htr_pwr_stat; /* EOC Redundant Heater Power Status, 0:Off, 1:On */
    int lrc_htr_pwr_stat; /* LRC Heater Power Status, 0:Off, 1:On */
    int lrc_pwr_stat; /* LRC Power Status, 0:Off, 1:On */
    char mirror_position[14]; /* Mirror Position */
    int sps_pri_pwr_stat; /* SPS Primary Power Status, 0:Off, 1:On */
    int sps_red_pwr_stat; /* SPS Redundant Power Status, 0:Off, 1:On */
};

struct inst_stat_data_osmi {
    int band_1_lo; /* OSMI (LRC) Band 1 Low Edge Frequency Value */
    int band_1_hi; /* OSMI (LRC) Band 1 High Edge Frequency Value */
    int band_2_lo; /* OSMI (LRC) Band 2 Low Edge Frequency Value */
    int band_2_hi; /* OSMI (LRC) Band 2 High Edge Frequency Value */
    int band_3_lo; /* OSMI (LRC) Band 3 Low Edge Frequency Value */
    int band_3_hi; /* OSMI (LRC) Band 3 High Edge Frequency Value */
    int band_4_lo; /* OSMI (LRC) Band 4 Low Edge Frequency Value */
    int band_4_hi; /* OSMI (LRC) Band 4 High Edge Frequency Value */
    int band_5_lo; /* OSMI (LRC) Band 5 Low Edge Frequency Value */
    int band_5_hi; /* OSMI (LRC) Band 5 High Edge Frequency Value */
    int band_6_lo; /* OSMI (LRC) Band 6 Low Edge Frequency Value */
    int band_6_hi; /* OSMI (LRC) Band 6 High Edge Frequency Value */
    int gain; /* OSMI (LRC) Sensor Gain */
    int offset1; /* OSMI (LRC) Sensor Offset 1 */
    int offset2; /* OSMI (LRC) Sensor Offset 2 */
    int offset3; /* OSMI (LRC) Sensor Offset 3 */
    int offset4; /* OSMI (LRC) Sensor Offset 4 */
    int temp; /* OSMI (LRC) Sensor Temperature */
};

struct inst_stat_data_earth {
    float magfield_x; /* Earth Magnetic Field X Component */
    float magfield_y; /* Earth Magnetic Field Y Component */
    float magfield_z; /* Earth Magnetic Field Z Component */
};

struct inst_stat_data {
    struct isd_spacecraft_time sc_time;
    struct isd_glob_pos_syst_time gps_time;
    struct inst_stat_data_misc isd_misc;
    struct inst_stat_data_osmi isd_osmi;
    struct inst_stat_data_earth isd_earth;
};
#ifdef __cplusplus
}
#endif

#endif /* ISD_H_DEFINED */
