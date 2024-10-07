#include <stdint.h>
/* o3_toms.h
 * specific definitions needed for o3_toms program
 */
#define TOMS_TYP_UNKNOWN 0
#define TOMS_TYP_N7TOMS 1
#define TOMS_TYP_EPTOMS 2
#define TOMS_TYP_OMITOMS 3

struct toms_txt_info_struc_def {
    int32_t toms_typ; /* toms file type */
    int32_t year; /* year and doy of data */
    int32_t doy;
    int32_t nlon; /* # longitudes in grid and start, end, delta */
    float slon, elon, del_lon;
    int32_t nlat; /* # latitudes etc. */
    float slat, elat, del_lat;
    char gen_str[11]; /* date raw product was generated */
    char node_time[17]; /* ascending node time */
    int16_t *datarr; /* the array of ozone  */
};

typedef struct toms_txt_info_struc_def toms_txt_info_struc;
/*
 *  prototypes
 */
int rd_toms_ascii(char *, toms_txt_info_struc *);
int fill_smooth(int16_t *, char *, int32_t, int32_t);
int field_extend(float *, char *, int, int, int, float*);
int mk_ker(float, float **, int *, int *, int *, int *);
