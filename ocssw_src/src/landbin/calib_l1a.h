/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/L1A.V4.4/l1a.h,v 4.15 1995/12/07 19:17:36 seawifsd Exp seawifsd $
$Log: l1a.h,v $
Revision 4.15  1995/12/07 19:17:36  seawifsd
added snode,enode,startclat,startclon,endclat,endclon and
made it IO_SPEC_V4.4 and PROD_SPEC_V2.8 compliant.

Revision 4.14  1995/05/12 13:53:27  seawifsd
1. code conformed to specifications IO_SPEC_V43, PROD_SPEC_V27,
and NONPROD_SPEC_V12(no funtionality changes, just eliminate
code related to versions before IO_SPEC_V42)

Revision 4.13  1995/04/03 16:05:59  seawifsd
fixed typo(GENERA -> GENERAL).

Revision 4.12  1995/02/15 18:45:01  seawifsd
modified definition of the structures meta_l1a and rec_l1a to
add scsol_z, entry_year, entry_day, and csol_z.
rename macro LAKSHMI to EXTRA_FOR_BROWSE.

Revision 4.11  1995/01/27 18:20:47  seawifsd
changed the order of mirror from mirror[bands][sides] to mirror[sides][bands]

Revision 4.10  1995/01/17 20:03:19  seawifsd
Jan. 17, 1994 V4.1

Revision 4.1  1995/01/17 14:15:38  seawifsd
Jan. 9, 1994, 4.0

Revision 3.4  1994/12/15 16:02:14  seawifsd
cleanup code that was IO_SPEC_V10 specific.

Revision 3.3  1994/11/08 18:45:56  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:03:54  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.2  1994/11/01 15:58:49  frank
added 'infiles' in the default meta_l1a structure because the requirement
of L1A browse routines.

Revision 1.1.1.1  1994/10/05 20:08:10  frank
modified to fit the I/O spec v 3.3 and product spec. 1.0.

Revision 1.2  1994/05/10 18:44:42  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:23:43  seawifst
Initial revision

Removed all code that was only compiled when GENERAL was defined.
Norman Kuring		6-Nov-1996

Added Calibration Vgroup members to struct meta_l2Struct.
Norman Kuring		6-Nov-1996

Removed all code that was only compiled when EVERYTHING was defined.
Norman Kuring		7-Nov-1996

Put back all the GENERAL and EVERYTHING code.  Heavy sigh.
Norman Kuring		25-Nov-1996

W. Robinson, GSC, 22 Apr 97  add the cal_mod_struc definition
 */


#ifndef L1A_H
#define L1A_H

#ifdef GENERAL
#define get_l1a_open get_l1a_open_all
#define get_l1a_record get_l1a_record_all
#define get_l1a_close get_l1a_close_all
#define put_l1a_open put_l1a_open_all
#define put_l1a_record put_l1a_record_all
#define put_l1a_close put_l1a_close_all
#define set_l1aget_buffer_blksize set_l1aget_buffer_blksize_all
#define set_l1aput_buffer_blksize set_l1aput_buffer_blksize_all
#endif /* GENERAL */

#include "navigation.h"
#include "level_1a_index.h"

#ifndef MAX_HDF
#define MAX_HDF 16
#endif

#ifndef MAX_HDF_L1AGET
#define MAX_HDF_L1AGET 4
#endif

typedef struct meta_l1aStruct {
    char *product_name; /* ATTR Product name(file name)		*/
#ifdef GENERAL
    char *title; /* ATTR title				*/
#endif /* GENERAL */
    char *data_center; /* ATTR data_center, processing center	*/
    char *station; /* ATTR station				*/
    float station_lat; /* ATTR station latitude		*/
    float station_lon; /* ATTR station longitude		*/
    char *mission; /* ATTR mission				*/
    char *mission_char; /* ATTR Mission Characteristics		*/
    char *sensor; /* ATTR sensor				*/
    char *sensor_char; /* ATTR instrumentInformation		*/
#ifdef EVERYTHING
    char *dtype; /* ATTR Data type			*/
#endif /* EVERYTHING */
#ifdef GENERAL
    char *replace; /* ATTR Replacement flag (""/fname)	*/
#endif /* GENERAL */
#define ISSUE "this might need to be changed to GENERAL or EVERYTHING"
    char *sw_id; /* ATTR Software ID			*/
#ifdef GENERAL
    char *ptime; /* ATTR Processing time			*/
#endif /* GENERAL */
    char *infiles; /* ATTR Input files			*/
#ifdef GENERAL
    char *proc_con; /* ATTR Processing control		*/
    char *proc_log; /* ATTR Processing log			*/
#endif /* GENERAL */
    char *stime; /* ATTR Start time			*/
    char *etime; /* ATTR End time			*/
    char *ctime; /* ATTR scene center time		*/
    char *ntime; /* ATTR Node crossing time		*/
#ifdef EVERYTHING
    short syear; /* ATTR Start year			*/
    short sday; /* ATTR Start day			*/
    int smsec; /* ATTR	Start millisec			*/
    short eyear; /* ATTR End year			*/
    short eday; /* ATTR	End day				*/
    int emsec; /* ATTR End millisec			*/
#endif /* EVERYTHING */
    char *snode; /* ATTR Start Node			*/
    char *enode; /* ATTR End Node			*/
    int orbnum; /* ATTR orbit number			*/
    char *norad1; /* ATTR NORAD elements, first line	*/
    char *norad2; /* ATTR NORAD elements, second line	*/
#ifdef EVERYTHING
    int nsamp; /* ATTR pixel per scan line		*/
    int nrec; /* ATTR Number of scan line		*/
#endif /* EVERYTHING */
    int pix_start; /* ATTR LAC Pixel Start Number		*/
    int pix_sub; /* ATTR LAC Pixel Subsampling		*/
    int ncrec; /* ATTR scene center scan line		*/
    int nfrec; /* ATTR number of filled scan line	*/
    int ff_mis; /* ATTR FF missing frames		*/
    int sd_mis; /* ATTR SDPS missing frames		*/
#ifdef GENERAL
    /* ATTR Gain 1 saturated pixels					*/
    int satg1[BANDS_DIMS_1A];
    /* ATTR Gain 2 saturated pixels					*/
    int satg2[BANDS_DIMS_1A];
    /* ATTR Gain 1 non-saturated pixels				*/
    int nsatg1[BANDS_DIMS_1A];
    /* ATTR Gain 2 non-saturated pixels				*/
    int nsatg2[BANDS_DIMS_1A];
    /* ATTR Zero pixels						*/
    int zeroes[BANDS_DIMS_1A];
    /* ATTR Mean gain 1 radiance					*/
    float meanr1[BANDS_DIMS_1A];
    /* ATTR Mean gain 2 radiance					*/
    float meanr2[BANDS_DIMS_1A];
#endif /* GENERAL */
    char *lat_units; /* ATTR Latitude units			*/
    char *lon_units; /* ATTR Longitude units			*/
    float northlat; /* ATTR Northernmost latitude		*/
    float southlat; /* ATTR	Southernmost latitude		*/
    float westlon; /* ATTR	Westernmost longitude		*/
    float eastlon; /* ATTR	Easternmost longitude		*/
    float startclat; /* ATTR Start Center Latitude		*/
    float startclon; /* ATTR Start Center Longtitude		*/
    float endclat; /* ATTR End Center Latitude		*/
    float endclon; /* ATTR End Center Longtitude		*/
    float nodel; /* ATTR Orbit node longitude		*/
    int ntilts; /* MFSD Sensor Tilt			*/
    /* MFSD tilt indicators						*/
    short tilt_flags[MTILT_DIMS_1A];
    /* MFSD scan-line range of scene tilt				*/
    short tilt_ranges[MTILT_DIMS_1A][LTILT_DIMS_1A];
    /* MFSD lat of tilt-range scan-line end pt			*/
    float tilt_lats[MTILT_DIMS_1A][LTILT_DIMS_1A][PTILT_DIMS_1A];
    /* MFSD lon of tilt-range scan-line end pt			*/
    float tilt_lons[MTILT_DIMS_1A][LTILT_DIMS_1A][PTILT_DIMS_1A];
    /* Calibration Vgroup */
    short entry_year;
    short entry_day;
    short ref_year;
    short ref_day;
    short ref_minute;
    float mirror[SIDES_DIMS_1A][BANDS_DIMS_1A];
    double t_const[BANDS_DIMS_1A];
    double t_linear[BANDS_DIMS_1A];
    double t_quadratic[BANDS_DIMS_1A];
    float cal_offs[BANDS_DIMS_1A];
    float counts[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
    float rads[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
} meta_l1aType;

typedef struct rec_l1aStruct {
#ifdef EVERYTHING
    int *msec; /* MFSD scan-line time, msec	*/
    unsigned char *eng_qual; /* MFSD engr data quality flags	*/
    unsigned char *s_flags; /* MFSD scan-line quality flags	*/
#endif /* EVERYTHING */
#ifdef GENERAL
    short int *s_satp; /* MFSD saturated pixels per band*/
    short int *s_zerop; /* MFSD zero pixels per band	*/
#endif /* GENERAL */
    float *slat; /* MFSD Scan start pixel latitude*/
    float *slon; /* MFSD Scan start pixel longitude*/
    float *clat; /* MFSD Scan center pixel lat	*/
    float *clon; /* MFSD Scan center pixel lon	*/
    float *csol_z; /* MFSD Scan center solar zenith*/
    float *elat; /* MFSD Scan end pixel latitude	*/
    float *elon; /* MFSD Scan end pixel longitude*/
    float *tilt; /* MFSD tilt angle for scan line*/
#ifdef GENERAL
    short int *sc_id; /* MFSD Spacecraft ID		*/
    short int *sc_ttag; /* MFSD	Spacecraft time tag	*/
    unsigned char *sc_soh; /* MFSD spacecraft SOH data	*/
    short int *inst_tlm; /* MFSD SeaWiFS instrument Telemetry*/
#endif /* GENERAL */
#ifdef EVERYTHING
    short int *l1a_data; /* MFSD scan-line data		*/
#endif /* EVERYTHING */
#ifdef GENERAL
    short int *start_syn; /* MFSD Start synch pixel	*/
    short int *stop_syn; /* MFSD Stop sync pixel		*/
#endif /* GENERAL */
    short int *dark_rest; /* MFSD Dark restore pixel	*/
    short int *gain; /* MFSD Band gain		*/
    short int *tdi; /* MFSD Band TDI		*/
    float *inst_ana; /* MFSD Inst analog telemetry	*/
    unsigned char *inst_dis; /* MFSD Inst discrete telemetry	*/
#ifdef GENERAL
    float *sc_ana; /* MFSD Spacecraft analog tlm	*/
    unsigned char *sc_dis; /* MFSD Spacecraft discrete tlm	*/
#endif /* GENERAL */
    short *scan_temp; /* MFSD digitized scan temp.	*/
    short *side; /* MFSD mirror side		*/
    navblockType nav;
} rec_l1aType;

/* out of band correction */
#define OOB_OFF   0
#define OOB_DEFAULT_METHOD 1
#define OXYGEN_CORR_FACTOR 1.12f
extern int out_band_corr(float *radiances, float oxygen_factor, int nsample);

/* calibration modification structure definition */
typedef struct cal_mod_def {
    int flag; /* use flag: 0 - use cal file values
                                   1 - use input gain, cal file offset
                                   2 - use cal file gain, input offset
                                   3 - use input gain and cal  */
    double gain[8]; /* calibration gain */
    double offset[8]; /* calibration offset */
} cal_mod_struc;

#endif /* L1A_H */
