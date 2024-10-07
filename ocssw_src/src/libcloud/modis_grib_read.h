
#ifndef __MODIS_GRIB_READ_H
#define __MODIS_GRIB_READ_H
/*
!C-INC**********************************************************************

!Description:
    Include file for the read ancillary data routines.

!Input Parameters:  N/A
!Output Parameters: N/A

!Revision History:

!Team-unique Header:

    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and credits:
    N. Devine (devine@gsfc.nasa.gov)
    W. Wolf   (walter.wolf@sse.wisc.edu)

!Design Notes:

!END*************************************************************
*/

/* Structure for date & hours */
typedef struct {
  int   year;
  int   month;
  int   day;
  int   hour;
} grib_date_t;

/* Structure for specifying grib records to read.  See modis_grib_read(). */
typedef struct {
  int   rec_num;
  char  *name;
  float *data;
  int   nx;
  int   ny;
  grib_date_t date;
} grib_record_t;


/*  Definitions */

/*     GRIB 2 business   */
/*     end GRIB 2 business  */

#ifndef INT2
#define INT2(a,b)   ((1-(int) ((unsigned) (a & 0x80) >> 6)) \
                   * (int) (((a & 0x7f) << 8) + b))
#endif

#define BDS_LEN(bds)            ((int) ((bds[0]<<16)+(bds[1]<<8)+bds[2]))

#define BDS_MoreFlags(bds)      ((bds[3] & 16) != 0)
#define BDS_UnusedBits(bds)     ((int) (bds[3] & 15))

#define BDS_BinScale(bds)       INT2(bds[4],bds[5])

#define BDS_RefValue(bds)       (ibm2flt(bds+6))
#define BDS_NumBits(bds)        ((int) bds[10])

#define BDS_DataStart(bds)      ((int) (11 + BDS_MoreFlags(bds)*3))

/* breaks if BDS_NumBits(bds) == 0 */
#define BDS_NValues(bds)        (((BDS_LEN(bds) - BDS_DataStart(bds))*8 - \
                                BDS_UnusedBits(bds)) / BDS_NumBits(bds))

/* undefined value -- if bitmap */
#define UNDEFINED               9.999e20

/* version 1.2 of grib headers  w. ebisuzaki */

#define BMS_LEN(bms)            ((bms) == NULL ? 0 : \
                                 (bms[0]<<16)+(bms[1]<<8)+bms[2])
#define BMS_UnusedBits(bms)     ((bms) == NULL ? 0 : bms[3])
#define BMS_StdMap(bms)         ((bms) == NULL ? 0 : ((bms[4]<<8) + bms[5]))
#define BMS_bitmap(bms)         ((bms) == NULL ? NULL : (bms)+6)
#define BMS_nxny(bms)           ((((bms) == NULL) || BMS_StdMap(bms)) \
        ? 0 : (BMS_LEN(bms)*8 - 48 - BMS_UnusedBits(bms)))

/* version 1.3 of grib headers  w. ebisuzaki */
/* this version is incomplete */
#ifndef INT3
#define INT3(a,b,c) ((1-(int) ((unsigned) (a & 0x80) >> 6)) \
                   * (int) (((a & 127) << 16)+(b<<8)+c))
#endif
#ifndef INT2
#define INT2(a,b)   ((1-(int) ((unsigned) (a & 0x80) >> 6)) \
                   * (int) (((a & 127) << 8) + b))
#endif

#define GDS_LEN(gds)            ((int) ((gds[0]<<16)+(gds[1]<<8)+gds[2]))

#define GDS_NV(gds)             (gds[3])
#define GDS_DataType(gds)       (gds[5])

#define GDS_LatLon(gds)         (gds[5] == 0)
#define GDS_Gaussian(gds)       (gds[5] == 4)

#define GDS_LatLon_nx(gds)      ((int) ((gds[6] << 8) + gds[7]))
#define GDS_LatLon_ny(gds)      ((int) ((gds[8] << 8) + gds[9]))

/* index of NV and PV */
#define GDS_PL(gds)             ((gds[4] == 255) ? -1 : \
                                 (int) gds[3] * 4 + (int) gds[4] - 1)



/*
 * ************* local prototypes *************
 */

int modis_grib_setup (char* inlun, char* outlun, grib_record_t *grib_records,
                      int numofrecords, char *errmsg);
int modis_grib_read (char* lun, int version,
	             grib_record_t *grib_records, int N_grib_records,
                     char *errmsg );

unsigned char *seek_grib(FILE *file, long *pos, long *len_grib,
        unsigned char *buffer, unsigned int buf_len);

int read_grib(FILE *file, long pos, long len_grib, unsigned char *buffer);

int Name_Matches_PDS( char *name, char *pds );

int parse_grib_message( unsigned char *buffer,
                        unsigned char **pds,
                        unsigned char **gds,
                        unsigned char **bms,
                        unsigned char **bds,
                        int *nx, int *ny, long *nxny, char *errmsg);

void BDS_unpack(float *flt, unsigned char *bits, unsigned char *bitmap,
        int n_bits, int n, double ref, double scale);

double ibm2flt(unsigned char *ibm);

double int_power(double x, int y);

int GDS_grid(unsigned char *gds, int *nx, int *ny, long int *nxny);

#define TRUE (0==0)
#define FALSE (1==0)

/* version 3.4 of grib headers  w. ebisuzaki */
/* this version is incomplete */

#define PDS_Len1(pds)           (pds[0])
#define PDS_Len2(pds)           (pds[1])
#define PDS_Len3(pds)           (pds[2])
#define PDS_LEN(pds)            ((int) ((pds[0]<<16)+(pds[1]<<8)+pds[2]))
#define PDS_Vsn(pds)            (pds[3])
#define PDS_Center(pds)         (pds[4])
#define PDS_Model(pds)          (pds[5])
#define PDS_Grid(pds)           (pds[6])
#define PDS_HAS_GDS(pds)        ((pds[7] & 128) != 0)
#define PDS_HAS_BMS(pds)        ((pds[7] & 64) != 0)
#define PDS_PARAM(pds)          (pds[8])
#define PDS_L_TYPE(pds)         (pds[9])
#define PDS_LEVEL1(pds)         (pds[10])
#define PDS_LEVEL2(pds)         (pds[11])

#define PDS_KPDS5(pds)          (pds[8])
#define PDS_KPDS6(pds)          (pds[9])
#define PDS_KPDS7(pds)          ((int) ((pds[10]<<8) + pds[11]))

/* this requires a 32-bit default integer machine */
#define PDS_Field(pds)          ((pds[8]<<24)+(pds[9]<<16)+(pds[10]<<8)+pds[11])

#define PDS_Year(pds)           (pds[12])
#define PDS_Month(pds)          (pds[13])
#define PDS_Day(pds)            (pds[14])
#define PDS_Hour(pds)           (pds[15])
#define PDS_Minute(pds)         (pds[16])
#define PDS_ForecastTimeUnit(pds)       (pds[17])
#define PDS_P1(pds)             (pds[18])
#define PDS_P2(pds)             (pds[19])
#define PDS_TimeRange(pds)      (pds[20])
#define PDS_NumAve(pds)         ((int) ((pds[21]<<8)+pds[22]))
#define PDS_NumMissing(pds)     (pds[23])
#define PDS_Century(pds)        (pds[24])
#define PDS_Subcenter(pds)      (pds[25])
#define PDS_DecimalScale(pds)   INT2(pds[26],pds[27])
#define PDS_Year4(pds)          (pds[12] + 100*(pds[24] - (pds[12] != 0)))


/* Error messages added by N. Devine */
#define GRIB_SUCCESS                     0
#define GRIB_MISSING_END_SECTION        -1
#define GRIB_MISSING_GDS_NUM_DATAPOINTS  1
#define GRIB_BAD_SEEK                    8
#define GRIB_MALLOC_ERROR               11
#define GRIB_INPUT_ERROR                12
#define GRIB_ERROR_IN_CALLED            13
#define GRIB_BAD_NUM_VS_NAME            14

#define VERSION "wgrib v1.5.0b10 (5-07-96) Wesley Ebisuzaki"

#define CHECK_GRIB

/*
 * MSEEK = I/O buffer size for seek_grib
 */

#define MSEEK 1024
#define BUFF_ALLOC0     40000

#ifndef min
#define min(a,b)  ((a) < (b) ? (a) : (b))
#define max(a,b)  ((a) < (b) ? (b) : (a))
#endif


#endif
