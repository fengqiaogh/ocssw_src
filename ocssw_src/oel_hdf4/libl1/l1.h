#ifndef _L1_H
#define _L1_H

#include "l1_input.h"
#include "filehandle.h"
#include "l2_flags.h"

#include <clo.h>
#include <genutils.h>
#include <timeutils.h>
#include <libnav.h>

#include <stdint.h>
#include <stdbool.h>
#include "uncertainty.h"

// Metadata standard strings
#define   INSTITUTION "NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group"
#define   LICENSE "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
#define   NAMING_AUTHORITY "gov.nasa.gsfc.sci.oceandata"
#define   KEYWORDS_VOCABULARY "NASA Global Change Master Directory (GCMD) Science Keywords"
#define   STDNAME_VOCABULARY "CF Standard Name Table v36"
#define   CREATOR_NAME "NASA/GSFC/OBPG"
#define   CREATOR_EMAIL "data@oceancolor.gsfc.nasa.gov"
#define   CREATOR_URL "https://oceandata.sci.gsfc.nasa.gov"
#define   PROJECT "Ocean Biology Processing Group (NASA/GSFC/OBPG)"
#define   PUBLISHER_NAME "NASA/GSFC/OBPG"
#define   PUBLISHER_EMAIL "data@oceancolor.gsfc.nasa.gov"
#define   PUBLISHER_URL "https://oceandata.sci.gsfc.nasa.gov"

#ifndef MIN
#define MIN(a,b)    (((a)<(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b)    (((a)>(b)) ? (a) : (b))
#endif

#ifndef PI
#define PI  3.141592654
#endif
#define RADEG 57.29577951

#define OFF            0
#define ON             1
#define NO             0
#define YES            1

#define SUCCESS        0
#define FATAL_ERROR    1
#define LONLAT_ERROR 110
#define MAXPIX     10000

#define BANDW 10

#define BT_LO       -1000
#define BT_HI        1000

#define SOLZNIGHT      90.0
#define SOLZNIGHTA     80.0
#define GLINT_MIN  0.0001


//VIIRS aggregation zone pixel indexes (M-band)
#define AGZONE0    0
#define AGZONE1  640
#define AGZONE2 1008
#define AGZONE3 1600
#define AGZONE4 2192
#define AGZONE5 2560

// evalmask
#define STDPROC          0   /* evalmask bit definitions                 */
#define OLDAERMOD        1   /* init to old aerosol models               */
#define MODCLOUD         2   /* enables MODIS/MERIS cloud mask algorithm */
#define MODCIRRUS       16   /* enables MODIS cirrus mask                */
#define NEWSENSINFO     32   /* use test sensor info file                */
#define NEWRAYTAB       64   /* use test rayleigh tables                 */
#define NEWAERTAB      128   /* use test aerosol tables                  */
#define NEWPOLTAB      256   /* use test polarization tables             */
#define MSKMODMIR1    1024   /* mask modis mirror-side 1 (navfail)       */
#define MSKMODMIR2    2048   /* mask modis mirror-side 2 (navfail)       */
#define SSTMODS       4096   /* reserved for testing SST changes         */
#define ALTSENSORINFO 8192   /* use .alt sensor infor file in eval       */
#define TRANSSPHER   32768   /* enables spherical path geom for dtran    */

#define XCALRVS          1
#define XCALPOL          2
#define XCALOLI          4 //Sudipta added for OLI SCA based XCAL

#ifdef __cplusplus
extern "C" {
#endif


/* Notice: any changes to this structure may require modifications to the */
/* following routines: alloc_l1.c, cpl1rec.c, l1subpix.c.                 */

typedef struct geom_struc_def {
    float *senz;
    float *sena;
    float *csenz;
    float *solz;
    float *sola;
    float *csolz;
    float *delphi;
    float *scattang;
} geom_struc;

typedef struct anc_add_struc_def {
    int32_t nlvl;   /* number of profile levels */
    float *prof_temp;
    float *prof_rh;
    float *prof_height;
    float *prof_q;
    float *prof_o3;
} anc_struc;

// ancillary aerosol from GMAO MERRA model
typedef struct anc_aer_struc_def {
    float *black_carbon_ext;
    float *black_carbon_scat;
    float *dust_ext;
    float *dust_scat;
    float *sea_salt_ext;
    float *sea_salt_scat;
    float *sulphur_ext;
    float *sulphur_scat;
    float *organic_carbon_ext;
    float *organic_carbon_scat;
    float *total_aerosol_ext;
    float *total_aerosol_scat;
    float *total_aerosol_angstrom;
} anc_aer_struc;

// for the cloud parameters, now just the sfc albedos
typedef struct cld_dat_struc_def {
    float *sfc_albedo_659;
    float *sfc_albedo_858;
    float *sfc_albedo_1240;
    float *sfc_albedo_1640;
    float *sfc_albedo_2130;
    //  added for the cloud height need
    float *cth_alb_init;
    float *cth_alb_unc_init;
} cld_struc;

typedef struct cld_rad_dat_struct_def
{
    // par needed data
    float **taucld;
    float **cfcld;
    float *timecldrange;
    size_t ntimes;
}
cld_rad_struct;
typedef struct l1_struct {
    int32_t length; /* number of bytes allocated to data block */
    int32_t npix; // number of pixels

    int32_t iscan; // number of lines
    int32_t detnum; // detector index for multidetector sensors (ie. MODIS)
    int32_t mside; // mirror side 0 or 1

    double scantime; // time of scan in unix time
    double fsol;

    bool is_l2; /**< Lt values are actually (above water?) reflectance; skip atmocor */

    /* scan attributes */

    float tilt;
    float alt; //altitude of sensor

    /* All parameters below are scan-length dependent */

    /* sensor band-pass-specific data */


    char *data; /* points to start of variable-length data block */

    int32_t *nobs;
    float *lon;
    float *lat;
    float *solz; // solz[pix] is solar zenith angle in degrees
    float *sola; // sola[pix] is solar azimuth angle in degrees
    float *senz; // senz[pix] is sensor zenith angle in degrees
    float *sena; // sena[pix] is sensor azimuth angle in degrees
    float *Lt; // Lt[pix][band]

    float *Ltir; // Ltir[pix][IRband]
    float *Bt;

    float *delphi;
    float *csolz; // csolz[pix] is cos(solz) in radians
    float *csenz; // csenz[pix] is cos(senz) in radians
    int32_t *pixnum; // pixnum[pix] is pixel index from non-extracted L1 file
    char *slot; /**< slot number                                */
    float *alpha;
    float *scattang;

    float *ws;
    float *wd;
    float *mw;
    float *zw;
    float *pr;
    float *oz;
    float *wv;
    float *rh;
    float *no2_tropo;
    float *no2_strat;
    float *no2_frac;
    float *sfcp;
    float *sfcrh;
    float *sfct;
    float *icefr;
    float *height;
    float *dem;
    // TODO: can get rid of this.  Only used in setanc.c
    short *ancqc;


    short *ssttype; /* per pixel - reference type or climatology */

    int32_t *flags;
    char *mask; // this group of params is the flags expanded into a byte
    char *hilt;
    char *cloud;
    char *glint;
    char *land;
    char *swater;
    char *ice;
    char *solzmax;
    char *senzmax;
    char *stlight;
    char *absaer;
    char *navfail;
    char *navwarn;

    char *filter;

    float *t_h2o;
    float *t_o2;
    float *t_sol;
    float *t_sen;
    float *rhof;
    float *tLf;
    float *Lr;
    float *L_q;
    float *L_u;
    float *polcor;
    float *dpol;
    float *TLg;
    float *rhos;
    float *glint_coef;
    float *cloud_albedo;
    float *aerindex;
    float *sstref;
    float *sssref;
    float *sw_n;
    float *sw_a;
    float *sw_bb;
    float *sw_a_avg;
    float *sw_bb_avg;
    float *rho_cirrus;

    double *tg_sol;
    double *tg_sen;
    double *tg;

    // TODO: move MERIS L1 to private_data pointer in l1rec
    /* for MERIS L1 */
    int32_t *pixdet; /* detector index of pixel */
    float *radcor; /* smile correction */


    float *Fo;

    // TODO: this needs to go into private_data pointer in filehandle
    /* for VIIRS unaggregated and superscan */
    int16_t scn_fmt; /* scan format of data, 0 std, else unaggregated */
    float margin_s; /* extra scan margin beyond actual samples */

    filehandle *l1file;

    // pointer to data needed by specific readers so far just meris, hawkeye
    void *private_data;
     
    // geometry per band
    geom_struc *geom_per_band;

    // added ancillary data, for CHIMAERA profiles, etc
    anc_struc *anc_add;

    // ancillary aerosol information from MERRA-2
    anc_aer_struc *anc_aerosol;

    // cloud processing data
    cld_struc *cld_dat;
    //cloud RAD data for PAR
    cld_rad_struct *cld_rad;
    //uncertainty record
    uncertainty_t *uncertainty;

} l1str;

void l1_input_init();
void l1_input_delete(l1_input_t *input);
void l1_add_options(clo_optionList_t* list);
void l1_read_default_files(clo_optionList_t *list, filehandle *l1file, const char *ifile);
void l1_load_options(clo_optionList_t *list, filehandle *l1file);
void l1_get_input_params(filehandle *l1file, char *input_parms);
void l1_get_input_files(filehandle *l1file, char *input_files);

void bindex_set(int32_t wave[], int nwave, int dwave_vswir);
int bindex_get(int32_t wave);
int bindex_get_555(void);
int windex(float wave, float twave[], int ntwave);
int invbindx(int band, int32_t *bindx, int nbands);

void radiance2bt(l1str *l1rec, int resolution);
void flag_bowtie_deleted(l1str *l1rec, size_t ipix, int extract_offset);

int get_f0_neckel(int32_t wl, int32_t width, float *f0);
int get_f0_thuillier(int32_t wl, int32_t width, float *f0);
void get_f0_thuillier_ext(int32_t wl, int32_t width, float *f0);

int init_geom_per_band(l1str *);
int geom_per_band_deriv(l1str *);
int destroy_geom_per_band(geom_struc *);

void init_l1(l1str *l1rec);
int32_t alloc_l1(filehandle *l1file, l1str *l1rec);
void free_l1(l1str *l1rec);

int openl1(filehandle *l1file);
int readl1(filehandle *l1file, int32_t recnum, l1str *l1rec);
int readl1_lonlat(filehandle *l1file, int32_t recnum, l1str *l1rec);
void closel1(filehandle *l1file);
void closel1_generic(filehandle *l1file);

int openl1_write(filehandle *l1file);
int writel1(filehandle *l1file, int32_t recnum, l1str *l1rec);
int openl1_read_hdf(filehandle *l1file);
int openl1_write_hdf(filehandle *l1file);
int writel1_hdf(filehandle *l1file, int32_t recnum, l1str *l1rec);
void closel1_hdf(filehandle *l1file);

int l1subpix(filehandle *l1file, l1str *l1rec);

void l1_mask_set(l1str *l1rec, int32_t ip);
int setflags(l1str *l1rec);
void setflagbits_l1(int level, l1str *l1rec, int32_t ipix);

char get_cldmask(l1str *l1rec, int32_t ip);
int modis_cloud_mask(l1str *l1rec, int32_t ip);
int get_sdps_cld_mask( l1str *, int32_t, char *);

#ifdef __cplusplus
}
#endif

#endif 
