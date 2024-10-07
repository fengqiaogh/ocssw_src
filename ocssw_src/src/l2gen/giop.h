#ifndef _GIOP_H
#define _GIOP_H


/* optimization routines */
#define AMOEBA    0
#define LEVMARQ   1
#define SVDFIT    3
#define SVDSIOP   4

/* adg function types */
#define ADGTAB    0
#define ADGS      1
#define ADGSQAA   2
#define ADGSOBPG  3
#define ADGSIOPTAB 4

/* bbp function types */
#define BBPTAB     0
#define BBPS       1
#define BBPSHAL    2
#define BBPSQAA    3
#define BBPSPML    4
#define BBPSCIOTTI 5
#define BBPSMM01   6
#define BBPSLAS    7
#define BBPLAS     8
#define BBPLASFIX  9
#define BBPQAAFIX 10
#define BBPLH     11
#define BBPCHL    12

/* aph function types */
#define APHTAB     0
#define APHGAUSS   1
#define APHBRICAUD 2
#define APHCIOTTI  3


/* acdom function type*/
#define ACDOMTAB  0
#define ACDOMNONE 1

/* anap function type*/
#define ANAPTAB  0
#define ANAPNONE 1

/* bbph function type*/
#define BBPHTAB  0
#define BBPHNONE 1

/* bbnap function type*/
#define BBNAPTAB  0
#define BBNAPNONE 1


/* rrs function types */
#define RRSGRD    0
#define RRSFOQ    1

/*uncertainty input option*/
#define URRSNONE 0
#define URRSCALC 1
#define URRSABS  2
#define URRSREL  3

typedef struct giop_ctl_str {
    int npar; /* # optimized params   */
    int maxiter; /* max iterations       */

    int fit_opt; /* optimization routine */
    int aph_opt; /* aph function type    */
    int adg_opt; /* adg function type    */
    int bbp_opt; /* bbp function type    */
    int acdom_opt; /* acdom function type    */
    int anap_opt; /* anap function type    */
    int bbph_opt; /* bbph function type    */
    int bbnap_opt; /* bbnap function type    */
    int rrs_opt; /* rrs function type    */
    int wt_opt; /* do we have input wts */
    int urrs_opt; /*relative or absolute uncertainties option*/

    int nwave; /* # fit wavelengths    */
    float *wave; /* fit wavelengths  [NBANDS]    */
    int *bindx; /* index to sensor wave [NBANDS]*/
    float *aw; /* aw per fit wave  [NBANDS]    */
    float *bbw; /* bbw per fit wave  [NBANDS]   */
    float *wts; /* input Rrs wts  [NBANDS]       */

    float *Rrs_a; /* input Rrs uncertainties  [ALL BANDS]       */
    float *uRrs_a; /* input Rrs uncertainties  [ALL BANDS]       */

    int siopIdx; /*Optimal set of IOPs for svd_siop implementation*/

    float grd[2]; /* RRSGRD param         */
    float *foq; /* RRSFOQ variable  [NBANDS]    */
    float aph_s; /* APHS spectral param  */
    float adg_s; /* ADGS spectral param  */
    float uadg_s; /* ADGS spectral param uncertainty */
    float bbp_s; /* BBPS spectral param  */
    float ubbp_s; /* BBPS spectral param  uncertainty*/
    float aph_w; /* APHS referance wave  */
    float adg_w; /* ADGS reference wave  */
    float bbp_w; /* BBPS reference wave  */

    float chl; /* Input chlorophyll    */
    float uchl; /* Input chlorophyll uncertainty */

    char aph_tab_file[FILENAME_MAX];
    char uaph_tab_file[FILENAME_MAX];
    int aph_tab_nw; /* elements in aph tab */
    float *aph_tab_w; /* aph tab wavelengths */
    float **aph_tab_s; /* aph tab values      */
    float **uaph_tab_s; /* aph tab uncertainty values      */
    int aph_nvec; /* number of eigenvect */

    char adg_tab_file[FILENAME_MAX];
    char uadg_tab_file[FILENAME_MAX];
    int adg_tab_nw; /* elements in adg tab */
    float *adg_tab_w; /* adg tab wavelengths */
    float **adg_tab_s; /* adg tab values      */
    float **uadg_tab_s; /* adg tab uncertainty values      */
    int adg_nvec; /* number of eigenvect */

    char anap_tab_file[FILENAME_MAX];
    char uanap_tab_file[FILENAME_MAX];
    int anap_tab_nw; /* elements in anap tab */
    float *anap_tab_w; /* anap tab wavelengths */
    float **anap_tab_s; /* anap tab values      */
    float **uanap_tab_s; /* anap tab uncertainty values      */
    int anap_nvec; /* number of eigenvect */

    char acdom_tab_file[FILENAME_MAX];
    char uacdom_tab_file[FILENAME_MAX];
    int acdom_tab_nw; /* elements in acdom tab */
    float *acdom_tab_w; /* acdom tab wavelengths */
    float **acdom_tab_s; /* acdom tab values      */
    float **uacdom_tab_s; /* acdom tab uncertainty values      */
    int acdom_nvec; /* number of eigenvect */

    char bbp_tab_file[FILENAME_MAX];
    char ubbp_tab_file[FILENAME_MAX];
    int bbp_tab_nw; /* elements in bbp tab */
    float *bbp_tab_w; /* bbp tab wavelengths */
    float **bbp_tab_s; /* bbp tab values      */
    float **ubbp_tab_s; /* bbp tab values      */
    int bbp_nvec; /* number of eigenvect */

    char bbph_tab_file[FILENAME_MAX];
    char ubbph_tab_file[FILENAME_MAX];
    int bbph_tab_nw; /* elements in bbph tab */
    float *bbph_tab_w; /* bbph tab wavelengths */
    float **bbph_tab_s; /* bbph tab values      */
    float **ubbph_tab_s; /* bbph tab uncertainty values      */    
    int bbph_nvec; /* number of eigenvect */

    char bbnap_tab_file[FILENAME_MAX];
    char ubbnap_tab_file[FILENAME_MAX];
    int bbnap_tab_nw; /* elements in bbnap tab */
    float *bbnap_tab_w; /* bbnap tab wavelengths */
    float **bbnap_tab_s; /* bbnap tab values      */
    float **ubbnap_tab_s; /* bbnap tab uncertainty values      */
    int bbnap_nvec; /* number of eigenvect */

    double *par; /* starting params [NBANDS]     */
    double *len; /* characteristic length [NBANDS] */

} giopstr;

void set_giop_aph_s(giopstr *g, float aph_s);
void set_giop_adg_s(giopstr *g, float adg_s);
void set_giop_bbp_s(giopstr *g, float bbp_s);
void set_giop_aph_t(giopstr *g, float aphw[], float aphs[], int aphn);

float* giop_get_chl_pointer();
float* giop_get_adg_pointer();
float* giop_get_aph_pointer();
float* giop_get_bbp_pointer();
float** giop_get_fitpar_pointer();
float* giop_get_mrrs_pointer();

void run_giop(l2str *l2rec);

#endif
