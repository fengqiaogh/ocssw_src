/* ========================================================================= */
/* module anc_acq.c - functions to read alternate ancillary data             */
/*     Intended to initially do the MET and OZ with ECMWF data               */
/*                                                                           */
/* Written By: W. Robinson, SAIC, Aug, 2013                                  */
/* W. Robinson, SAIC, 3 Oct 24  to separate the FRSNOICE to be used only     */
/*     land pixels and FRSEAICE only over water                              */
/*                                                                           */
/* ========================================================================= */
#include "l12_proto.h"
#include "met_cvt.h"
#include "anc_acq.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
/* status of the 3 ancillary files: MET1, 2, 3
 or OZONE1, 2, 3 - also set it up so the 1,2,3 mean which files are usable */
#define ANC_STAT_1T 1 /* one time used */
#define ANC_STAT_2T_END                       \
    4 /* 2 different anc times in end of list \
         = files[0] = [1] != [2] */
#define ANC_STAT_2T_START                                        \
    2                   /*2 different anc times in start of list \
                         = files[0] ! [1] = [2] */
#define ANC_STAT_3T 3   /* all 3 files with different times */
#define ANC_STAT_CLIM 0 /* a climatology in  use */
#define NPRM 7
#define ANCBAD -999.
#define OZ_KG_M2_TO_DU 1. / 2.1415e-5
#define USE_PMSL                                                                    \
    1                         /* choice to use surface pressure (0)  or use the MSL \
                                 (mean sea level) pressure (1) with appropriate     \
                                 adjustment for land height */
#define ANC_SRC_TYP_ECMWF 0   /* types of anc file data */
#define ANC_SRC_TYP_STD_HDF 1 /* old HDF from NCEP (met) TOMS (oz) */
#define ANC_SRC_TYP_OLCI 2    /* OLCI tie point at data points */
#define ANC_SRC_TYP_GMAO 3    /* GMAP FP-IT or FP forecast */
#define ANC_SRC_TYP_BAD -1    /* a bad or undefined type */
#define FILE_CHECK(file_name)             if (access(file_name, F_OK) || access(file_name, R_OK)) { \
printf("--Error-- : Input file '%s' does not exist or cannot open.\n.See line %d in %s\n", \
file_name, __LINE__, __FILE__); \
exit(FATAL_ERROR); \
}

// WDR carry ice over land, water separately
enum out_nam { ZW, MW, PR, WV, RH, SFCP, SFCRH, SFCT, ICEFR_WTR, ICEFR_LND };

enum out_prof { TPROF, RHPROF, HPROF, QPROF, O3PROF };

enum out_cldrad { TAUCLD, CFCLD };

enum out_aerosol {
    BCEXTTAU,
    BCSCATAU,
    DUEXTTAU,
    DUSCATAU,
    SSEXTTAU,
    SSSCATAU,
    SUEXTTAU,
    SUSCATAU,
    OCEXTTAU,
    OCSCATAU,
    TOTEXTTAU,
    TOTSCATAU,
    TOTANGSTR
};

struct met_sto_str_d {
    float s_lon;         /* start longitude for a grid */
    float lon_step;      /* longitude incriment */
    int nlon;            /* # longitude points  */
    int nlat;            /* # latitude points  */
    float e_lon;         /* end longitude for a grid */
    float s_lat;         /* start latitude for a grid */
    float lat_step;      /* latitude incriment */
    float e_lat;         /* end latitude for a grid */
    double data_time[3]; /* time (in Julian days and fraction)
     for the data 1, 2, 3 */
    int anc_f_stat;      /* status of the met data */
    float *data[3];      /*  storage for MET1, 2, 3 */
};

typedef struct met_sto_str_d met_sto_str;
static met_sto_str met_sto[NPRM];

static int proc_land;

int anc_acq_init(instr *input, l1str *l1rec, int32_t *anc_id)
/*******************************************************************

   anc_acq_init

   purpose: Identify the incoming MET files and if netcdf, set up the
   ancillary data so it can be accessed by anc_acq_line
   For now, only ECMWF netcdf files can be processed which contain
   only 1 time and (at least) the parameters listed in prm_nm

   Returns type: int - 0 - good -1 any trouble checking input anc files

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     instr *           input            I      program inputs
     l1str *           l1rec            I      level 1 record structure
     int32_t *         anc_id           O      size 2 ID of anc data being
                                               handled: for [met, ozone]
                                               0 - ECMWF files, 1 - NCEP/TOMS
                                               files, -1 - bad

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       12-Aug-2013     Original development
     W. Robinson, SAIC 22 Feb 2016     enhance to handle OLCI tie point
                                       ancillary data

 *******************************************************************/
{
    char *files[3]; /* the 3 MET file names */
    /* ECMWF parameter names in the 2 arrays below */
    /*  name   descrip   ends up as value in met_sto
     sp     sfc pressure           pressure at surface (not currently used)
     tcwv   precip water           precip water
     msl    MSL pressure           seal lvl pressure
     u10    10 m u wind            uwind meridional S->N (+)
     v10    10 m v wind            vwind zonal W->E (+)
     t2m    2m temp           -    Relative humidity
     d2m    2m dewpoint temp  /
     tco3   total Ozone            ozone
     */
    char *prm_nm_met[] = {"sp", "tcwv", "msl", "u10", "v10", "t2m", "d2m"};
    char *prm_nm_oz[] = {"tco3"};
    int32_t n_prm_met = 7, n_prm_oz = 1, sto_ix;
    static char file_olci_str[FILENAME_MAX];
    char *file_olci = (char *)0; /* for olci file name, 0 if not olci */
    /*
     *  determine if rad data ia OLCI and get tie point met file name if so
     */
    if (l1rec->l1file->format == FT_OLCI) {
        strcpy(file_olci_str, l1rec->l1file->name);
        char *ptr = strrchr(file_olci_str, '/');
        if (ptr) {
            *ptr = 0;
            strcat(file_olci_str, "/tie_meteo.nc");
        } else {
            strcpy(file_olci_str, "tie_meteo.nc");
        }
        file_olci = file_olci_str;
    }
    proc_land = input->proc_land; /* carry to the anc_acq_lin routine */
    /*
     *  get a general classification of the ancillary data from the 1st file
     *  name for met and ozone
     */
    anc_id[0] = anc_acq_ck(input->met1, file_olci);
    anc_id[1] = anc_acq_ck(input->ozone1, file_olci);
    if ((anc_id[0] == -1) || (anc_id[1] == -1))
        return -1;
    else {
        if (anc_id[0] == ANC_SRC_TYP_ECMWF) {
            /*  do ecmwf met data full check and init */
            files[0] = input->met1;
            files[1] = input->met2;
            files[2] = input->met3;
            sto_ix = 0; /* location of met in storage struct */

            anc_id[0] = anc_acq_ecmwf_init(files, prm_nm_met, n_prm_met, sto_ix);
        }
        /* do olci set-up in the read: anc_acq_lin_olci */
        if (anc_id[1] == ANC_SRC_TYP_ECMWF) {
            /*  do ecmwf ozone data full check and init */
            files[0] = input->ozone1;
            files[1] = input->ozone2;
            files[2] = input->ozone3;
            sto_ix = 6; /* location of oz in storage struct */

            anc_id[1] = anc_acq_ecmwf_init(files, prm_nm_oz, n_prm_oz, sto_ix);
        }
    }
    if ((anc_id[0] == -1) || (anc_id[1] == -1))
        return -1;
    else
        return 0;
}

int32_t anc_acq_ck(char *file, char *file_olci)
/*******************************************************************

   anc_acq_ck

   purpose: Identify the incoming MET or OZONE file #1 as either netcdf OLCI,
   netcdf ECMWF, or not netcdf, which means it is a std older file format

   Returns type: int ancillary data file source type, see file top for defs:
     form:  ANC_SRC_TYP_<type>

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     char *            file             I      1st anc file name
     char *            file_olci        I      if null, data type is not OLCI
                                               else, the name of the olci met
                                               tie point data file

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       22-Feb-2016     Original development

 *******************************************************************/
{
    size_t attlen = 0;
    int status, ncid, fstat, tie_pt;
    char *str_attr, lcl_fil[FILENAME_MAX];
    /*
     *  make sure 1st file has something
     */
    fstat = ANC_SRC_TYP_BAD;
    if ((file == NULL) || (file[0] == 0)) {
        fprintf(stderr, "-E- %s %d: First anc file name is null\n", __FILE__, __LINE__);
        return -1;
    }
    /*
     *  First, any keyword substitutions
     *  If OLCI data and file name is special OLCI designation: OLCI_TIE_METEO,
     *  set the file to the olci met tie point file
     */
    strcpy(lcl_fil, file);
    if ((file_olci != (char *)0) && (strcmp(upcase(lcl_fil), "OLCI_TIE_METEO") == 0)) {
        fprintf(stderr,
                "-I- %s %d: Ancillary file will be taken from OLCI metadata tie point meteo file: %s\n",
                __FILE__, __LINE__, file_olci);
        strcpy(file, file_olci);
    }
    /*
     *  next, see if file opens with netcdf
     */
    // first check if it is an HDF4 file
    if (Hishdf(file)) {
        return ANC_SRC_TYP_STD_HDF;
    }

    if (nc_open(file, 0, &ncid) != NC_NOERR) {
        // File is not netcdf, thus assume it's a "standard" (HDF4) file
        return ANC_SRC_TYP_STD_HDF;
    } else {
        /*
         *  for the OLCI data only, see if the file is tie point data
         */
        tie_pt = 0;
        if (file_olci != (char *)0) /* data type is OLCI */ {
            status = nc_inq_attlen(ncid, NC_GLOBAL, "title", &attlen);
            if (status == NC_NOERR) {
                str_attr = (char *)calloc(attlen + 1, sizeof(char));
                status = nc_get_att_text(ncid, NC_GLOBAL, "title", str_attr);
                if (status != NC_NOERR) {
                    fprintf(stderr, "-I- %s %d: nc_get_att_string error, file: %s\n", __FILE__, __LINE__,
                            file);
                    return ANC_SRC_TYP_BAD;
                }

                if (strcmp(str_attr, "OLCI Level 1b Product, Tie-Point Meteo Data Set") == 0) {
                    free(str_attr);
                    /* should free the space with  nc_free_string() as said in nc_get_att_string call */
                    fstat = ANC_SRC_TYP_OLCI;
                    tie_pt = 1;
                } else
                    fstat = ANC_SRC_TYP_BAD;
            }
        }
        /*
         *  look to see if the file is a GMAO file
         */
        if (tie_pt == 0) {
            status = nc_inq_attlen(ncid, NC_GLOBAL, "title", &attlen);
            if (status == NC_NOERR) {
                /* Check for any of the GMAO ids set up */
                str_attr = (char *)calloc(attlen + 1, sizeof(char));
                status = nc_get_att_text(ncid, NC_GLOBAL, "title", str_attr);
                if (status != NC_NOERR) {
                    fprintf(stderr, "-I- %s %d: nc_get_att_string error, file: %s\n", __FILE__, __LINE__,
                            file);
                    return ANC_SRC_TYP_BAD;
                }
                if (strstr(str_attr, "GMAO") != 0) {
                    fstat = ANC_SRC_TYP_GMAO;
                }
                free(str_attr);
            }
        }
        status = nc_close(ncid);
        return fstat;
    }
    /*
     *  data type id not olci or GMAO, but file is netcdf, so only
     *  thing left is ecmwf
     */
    return ANC_SRC_TYP_ECMWF;
}

int32_t anc_acq_f_stat(char **files, char prioritize_files, int32_t n_anc)
/*******************************************************************
   anc_acq_f_stat

   purpose: find the status of the input set of ancillary files.  Enlarged
            from the treatment for the MET, OZONE 3-file set and added treatment
            for 1 and 2 files

   Returns type: int32_t - the status of the files sent in: ANC_STAT_...
                         (see declarations at anc_acq start)
                          -1 for bad

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     char **           files            I      list of files of the ancillary
                                               data
     char              prioritize_files I      0 - do nothing, 1 - for 3 files
                                               with only 2 that are different,
                                               make sure it is always the 1st
                                               and 2nd files
     int32_t           n_anc            I      Number of files to consider

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       2 May 2018      original development

 *******************************************************************/
{
    int32_t f12_mtch, f23_mtch, anc_f_stat;

    anc_f_stat = -1;
    if (n_anc == 3) {
        if ((files[1] == NULL) || (files[2] == NULL) || (files[1][0] == 0) || (files[2][0] == 0))
            anc_f_stat = ANC_STAT_CLIM;
        else {
            f12_mtch = 0;
            if (strcmp(files[0], files[1]) == 0)
                f12_mtch = 1;

            f23_mtch = 0;
            if (strcmp(files[1], files[2]) == 0)
                f23_mtch = 1;

            if ((strcmp(files[0], files[2]) == 0) && (f12_mtch == 0)) {
                printf("%s, %d E: ANC1 and 3 match while ANC2 different\n", __FILE__, __LINE__);
                return -1;
            }
            if ((f12_mtch == 1) && (f23_mtch == 1))
                anc_f_stat = ANC_STAT_1T;
            else if ((f12_mtch == 1) && (f23_mtch == 0))
                anc_f_stat = ANC_STAT_2T_END;
            else if ((f12_mtch == 0) && (f23_mtch == 1))
                anc_f_stat = ANC_STAT_2T_START;
            else
                anc_f_stat = ANC_STAT_3T;
        }
        /*  get the 2 different file names in 1st, 2nd slots */
        if (prioritize_files && (anc_f_stat == ANC_STAT_2T_END)) {
            anc_f_stat = ANC_STAT_2T_START;
            files[1] = files[2];
        }
    } else if (n_anc == 2) {
        if ((files[0] == files[1]) || files[1][0] == 0)
            anc_f_stat = ANC_STAT_1T;
        else
            anc_f_stat = ANC_STAT_2T_START;
    } else {
        anc_f_stat = ANC_STAT_1T;
    }
    /*  END of file set identification  */
    if (prioritize_files) {
        switch (anc_f_stat) {
            case ANC_STAT_1T:
                return 1;
                break;
            case ANC_STAT_2T_START:
            case ANC_STAT_2T_END:
                return 2;
                break;
            case ANC_STAT_3T:
                return 3;
                break;
        }
    }
    return anc_f_stat;
}

int32_t anc_acq_lin_met(l1str *l1rec)
/*******************************************************************

   anc_acq_lin_met

   purpose: process the met data for the l1rec, for now, just do the GMAO

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1str *           l1rec           I/O     all L1 line info

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       7 May 2018      Original development
      W. Robinson       3 Oct 2024      adapt to use ice over water only for 
                                        water pixels and same for land

 *******************************************************************/
{
    static int32_t firstcall = 0;
    static int32_t ntim_int; /* # needed times for interp arrays */
    static float r2d = OEL_RADEG;
    static int32_t anc_id = -1;
    static gen_int_str *met_int, *met_tim;

    char *files[3]; /* the 3 MET file names */
    static char file_olci_str[FILENAME_MAX];
    char *file_olci = (char *)0; /* for olci file name, 0 if not olci */
    int32_t n_met = 10;           /* # met parms finally needed, with ice 
                                     over land and water separately */
    int32_t nlvl_int = 1;        /* # levels in interpolation 1 in the 2-d case */
    int32_t ilvl = 0, nlvl = 1;  /* for met, only 1 level */
    int32_t itim, ilon, npix, iprm, t_interp, data_ix[2];
    int32_t field_skip;   /* To skip water ice over land and land ice over water */
    float val, wt_t1, uwnd, vwnd, unc, u_unc, v_unc, ws_2;
    double l_time, anc_times[3], lat, lon, last_time;

    /* initialize  for the met processing */
    if (firstcall == 0) {
        firstcall = 1;
        /*
         *  determine if rad data is OLCI and get tie point met file name if so
         */
        if (l1rec->l1file->format == FT_OLCI) {
            strcpy(file_olci_str, l1rec->l1file->name);
            char *ptr = strrchr(file_olci_str, '/');
            if (ptr) {
                *ptr = 0;
                strcat(file_olci_str, "/tie_meteo.nc");
            } else {
                strcpy(file_olci_str, "tie_meteo.nc");
            }
            file_olci = file_olci_str;
        }
        /*  get the files from l1rec->met1,2,3  */
        files[0] = input->met1;
        files[1] = input->met2;
        files[2] = input->met3;

        /*  although done previously (for now), do again to order the
            important files */
        if ((ntim_int = anc_acq_f_stat(files, 1, 3)) == -1)
            return 1;

        /*  set up the interpolation structure array */
        /*  note met_int( ntim_int, nlvl_int, n_met ) with met running fastest */
        if ((met_int = (gen_int_str *)malloc(n_met * nlvl_int * ntim_int * sizeof(gen_int_str))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate met interpolation array\n", __FILE__, __LINE__);
            return 1;
        }
        /*  Check the type of the data, call anc_acq_ck
         *   (not needed as this is only GMAO
         */
        if ((anc_id = anc_acq_ck(files[0], file_olci)) == -1)
            return 1;

        /* if GMAO type, set up the GMAO met data */
        if (anc_id == ANC_SRC_TYP_GMAO) {
            /* loop thru times / files from 1 - ntim_int */
            last_time = 0;
            for (itim = 0; itim < ntim_int; itim++) {
                met_tim = met_int + itim * n_met;
                /* call anc_acq_gmao_met_prep() *to get all parms for that file */
                if (anc_acq_gmao_met_prep(files[itim], met_tim) != 0)
                    return 1;
                if (last_time > met_tim[0].anc_time) {
                    fprintf(stderr, "-E- %s %d: met file times are out of sequence\n", __FILE__, __LINE__);
                    return 1;
                }
                last_time = met_tim[0].anc_time;
            }
        } /* END GMAO */
        else {
            /* report an error - something is wrong - no other types now */
            fprintf(stderr, "-E- %s %d: Error reading ancillary file: %s\n", __FILE__, __LINE__, files[0]);
            return -1;
        }
    } /* END INITIALIZE */

    /*
     * Start getting the data for the line in l1rec
     *  get the time of the current line
     */
    l_time = l1rec->scantime;
    npix = l1rec->npix;

    /* go through the parameters and interpolate in space and time */
    for (iprm = 0; iprm < n_met; iprm++) {
        /* find the times and weighting to use for this line */
        for (itim = 0; itim < ntim_int; itim++)
            anc_times[itim] = met_int[iprm + n_met * itim].anc_time;

        /* get the proper times to use and weights */
        if (anc_acq_fnd_t_interp(l_time, anc_times, ntim_int, &t_interp, data_ix, &wt_t1) != 0)
            return 1;

        /* interpolate for 1 or 2 times */
        for (ilon = 0; ilon < npix; ilon++) {
            lon = l1rec->lon[ilon];
            lat = l1rec->lat[ilon];
            // gsl gets cranky with bad inputs, so...
            if (lon < -180.0 || lon > 180.0 || lat < -90.0 || lat > 90.0 || isnan(lat) || isnan(lon)) {
                continue;
            }
            // only evaluate ice over land/water at land/water pixels
            field_skip = 0;
            if( ( (iprm == ICEFR_LND ) && (!l1rec->land[ilon]) ) ||
                ( (iprm == ICEFR_WTR ) && (l1rec->land[ilon]) ) )
               field_skip = 1;
            if(field_skip == 0 ) {
                if (anc_acq_eval_pt(met_int, iprm, ilvl, lat, lon, t_interp, data_ix, wt_t1, ntim_int, nlvl,
                  n_met, &val, &unc) != 0) {
                     fprintf(stderr, "-E- %s %d: Error interpolating to file: %s\n", __FILE__, __LINE__, files[0]);
                    return -1;
                }
            }
            /*  fill in the proper ancillary data slot */
            switch (iprm) {
                case ZW:
                    l1rec->zw[ilon] = val;
                    if (l1rec->uncertainty)
                        l1rec->uncertainty->dzw[ilon] = unc;
                    break;
                case MW: /* the mwind comes second, so then combinations cn be done */
                    l1rec->mw[ilon] = val;
                    if (l1rec->uncertainty)
                        l1rec->uncertainty->dmw[ilon] = unc;
                    /* compute the wind speed, angle with the 2 components */
                    uwnd = l1rec->zw[ilon];
                    vwnd = l1rec->mw[ilon];
                    if (l1rec->uncertainty) {
                        u_unc = l1rec->uncertainty->dzw[ilon];
                        v_unc = l1rec->uncertainty->dmw[ilon];
                    }
                    if (input->windspeed != -2000)
                        l1rec->ws[ilon] = sqrt(pow(uwnd, 2.) + pow(vwnd, 2.));
                    if (input->windangle != -2000)
                        l1rec->wd[ilon] = atan2f(-uwnd, -vwnd) * r2d;
                    /* deal with the uncertainties in speed, direction */
                    uwnd = fabsf(uwnd);
                    vwnd = fabsf(vwnd);
                    ws_2 = uwnd * uwnd + vwnd * vwnd;
                    if (l1rec->uncertainty) {
                        if ((uwnd + vwnd) > 0.05 * (u_unc + v_unc)) {
                            l1rec->uncertainty->dws[ilon] =
                                sqrt((uwnd * uwnd * u_unc * u_unc + vwnd * vwnd * v_unc * v_unc) / ws_2);
                            l1rec->uncertainty->dwd[ilon] =
                                sqrt(vwnd * vwnd * u_unc * u_unc + uwnd * uwnd * v_unc * v_unc) / ws_2;
                            if (l1rec->uncertainty->dwd[ilon] > OEL_PI)
                                l1rec->uncertainty->dwd[ilon] = OEL_PI;
                        } else {
                            l1rec->uncertainty->dws[ilon] = sqrt(0.5 * (u_unc * u_unc + v_unc * v_unc));
                            l1rec->uncertainty->dwd[ilon] = OEL_PI;
                        }
                        l1rec->uncertainty->dwd[ilon] *= r2d;
                    }
                    break;
                case PR:
                    if (input->pressure != -2000) {
                        if (proc_land && (l1rec->height[ilon] != 0))
                            val *= exp(-l1rec->height[ilon] / 8434.);
                        l1rec->pr[ilon] = val;
                        if (l1rec->uncertainty)
                            l1rec->uncertainty->dpr[ilon] = unc;
                    }
                    break;
                case WV: /* to make the kg m^-2 to g cm^-2 */
                    if (input->watervapor != -2000)
                        l1rec->wv[ilon] = val / 10.;
                    if (l1rec->uncertainty)
                        l1rec->uncertainty->dwv[ilon] = unc / 10.;
                    break;
                case RH:
                    if (input->relhumid != -2000)
                        l1rec->rh[ilon] = val;
                    if (l1rec->uncertainty)
                        l1rec->uncertainty->drh[ilon] = unc;
                    break;
                case SFCP:
                    l1rec->sfcp[ilon] = val;
                    break;
                case SFCRH:
                    l1rec->sfcrh[ilon] = val;
                    break;
                case SFCT:
                    l1rec->sfct[ilon] = val;
                    break;
                // for these last 2 fill icefr based on underlying sfc type
                case ICEFR_WTR:
                    if( !l1rec->land[ilon] ) {
                        l1rec->icefr[ilon] = val;
                        if (val > input->ice_threshold)
                            l1rec->ice[ilon] = ON;
                    }
                    break;
                case ICEFR_LND:
                    if( l1rec->land[ilon] ) {      
                        l1rec->icefr[ilon] = val;
                        if (val > input->ice_threshold)
                            l1rec->ice[ilon] = ON;
                    }
                    break;
            }
        }
    }
    /*  and end */
    return 0;
}
int32_t anc_acq_lin_aerosol(l1str *l1rec) {
    static int32_t firstcall = 0;
    static int32_t ntim_int; /* # needed times for interp arrays */
    static gen_int_str *aer_int, *aer_tim;

    char *files[3];     /* the 3 MET aerosol file names */
    int32_t n_met = 13; /* # met parms finally needed */
    int32_t itim, ilon, npix, iprm, t_interp, data_ix[2];
    int32_t ilvl = 0, nlvl = 1; /* only 1 level */

    float val, wt_t1, unc;
    double l_time, anc_times[3], lat, lon, last_time;
    anc_aer_struc *anc_aerosol;

    npix = l1rec->npix;
    /* initialize  for the met processing */
    if (firstcall == 0) {
        firstcall = 1;
        /*  get the files from l1rec->anc_profile1,2,3  */
        files[0] = input->anc_aerosol1;
        files[1] = input->anc_aerosol2;
        files[2] = input->anc_aerosol3;

        printf("\nOpening ancillary aerosol files.\n");
        printf("  anc_aerosol1 = %s\n", input->anc_aerosol1);
        printf("  anc_aerosol2 = %s\n", input->anc_aerosol2);
        printf("  anc_aerosol3 = %s\n", input->anc_aerosol3);
        printf("\n");

        /*  although done previously (for now), do again to order the
            important files */
        if ((ntim_int = anc_acq_f_stat(files, 1, 3)) == -1)
            return 1;

        /*  for a climatology (ntim_int = 0) just say that climatology
            profies not implemented yet */
        if (ntim_int == ANC_STAT_CLIM) {
            fprintf(stderr, "-E- %s %d: Ancillary met profile climatology not implemented yet\n", __FILE__,
                    __LINE__);
            return 1;
        }

        /*  set up the interpolation structure array */
        if ((aer_int = (gen_int_str *)malloc(n_met * nlvl * ntim_int * sizeof(gen_int_str))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate profile interpolation array\n", __FILE__,
                    __LINE__);
            return 1;
        }

        /* loop thru times / files from 1 - ntim_int */
        last_time = 0;
        for (itim = 0; itim < ntim_int; itim++) {
            aer_tim = aer_int + itim * n_met;
            /* call anc_acq_gmao_met_prep() *to get all parms for that file */
            if (anc_acq_gmao_aer_prep(files[itim], aer_tim) != 0)
                return 1;
            if (last_time > aer_tim[0].anc_time) {
                fprintf(stderr, "-E- %s %d: met file times are out of sequence\n", __FILE__, __LINE__);
                return 1;
            }
            last_time = aer_tim[0].anc_time;
        }
    } /* END INITIALIZE */

    /*  set up the storage for the profile data if needed */
    anc_aerosol = l1rec->anc_aerosol;
    if (anc_aerosol == NULL) {
        if (init_anc_aerosol(l1rec) != 0)
            return 1;
    }
    /*
     * Start getting the data for the line in l1rec
     *  get the time of the current line
     */
    l_time = l1rec->scantime;
    npix = l1rec->npix;
    /* go through the parameters and interpolate in space and time */
    for (iprm = 0; iprm < n_met; iprm++) {
        /* find the times and weighting to use for this line */
        for (itim = 0; itim < ntim_int; itim++)
            anc_times[itim] = aer_int[iprm + n_met * itim].anc_time;

        /* get the proper times to use and weights */
        if (anc_acq_fnd_t_interp(l_time, anc_times, ntim_int, &t_interp, data_ix, &wt_t1) != 0)
            return 1;

        /* interpolate for 1 or 2 times */
        for (ilon = 0; ilon < npix; ilon++) {
            lon = l1rec->lon[ilon];
            lat = l1rec->lat[ilon];
            // gsl gets cranky with bad inputs, so...
            if (lon < -180.0 || lon > 180.0 || lat < -90.0 || lat > 90.0 || isnan(lat) || isnan(lon)) {
                continue;
            }
            if (anc_acq_eval_pt(aer_int, iprm, ilvl, lat, lon, t_interp, data_ix, wt_t1, ntim_int, nlvl,
                                n_met, &val, &unc) != 0) {
                fprintf(stderr, "-E- %s %d: Error interpolating to file: %s\n", __FILE__, __LINE__, files[0]);
                return -1;
            }
            /*  fill in the proper ancillary data slot */
            switch (iprm) {
                case BCEXTTAU:
                    l1rec->anc_aerosol->black_carbon_ext[ilon] = val;
                    break;
                case BCSCATAU:
                    l1rec->anc_aerosol->black_carbon_scat[ilon] = val;
                    break;
                case DUEXTTAU:
                    l1rec->anc_aerosol->dust_ext[ilon] = val;
                    break;
                case DUSCATAU:
                    l1rec->anc_aerosol->dust_scat[ilon] = val;
                    break;
                case SSEXTTAU:
                    l1rec->anc_aerosol->sea_salt_ext[ilon] = val;
                    break;
                case SSSCATAU:
                    l1rec->anc_aerosol->sea_salt_scat[ilon] = val;
                    break;
                case SUEXTTAU:
                    l1rec->anc_aerosol->sulphur_ext[ilon] = val;
                    break;
                case SUSCATAU:
                    l1rec->anc_aerosol->sulphur_scat[ilon] = val;
                    break;
                case OCEXTTAU:
                    l1rec->anc_aerosol->organic_carbon_ext[ilon] = val;
                    break;
                case OCSCATAU:
                    l1rec->anc_aerosol->organic_carbon_scat[ilon] = val;
                    break;
                case TOTEXTTAU:
                    l1rec->anc_aerosol->total_aerosol_ext[ilon] = val;
                    break;
                case TOTSCATAU:
                    l1rec->anc_aerosol->total_aerosol_scat[ilon] = val;
                    break;
                case TOTANGSTR:
                    l1rec->anc_aerosol->total_aerosol_angstrom[ilon] = val;
                    break;
            }
        }
    }
    /*  and end */
    return 0;
}

int32_t anc_acq_lin_prof(l1str *l1rec)
/*******************************************************************

   anc_acq_lin_prof

   purpose: process the met profile data for the l1rec, from the GMAO
     FP-IT data

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1str *           l1rec           I/O     all L1 line info

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       7 May 2018      Original development

 *******************************************************************/
{
    static int32_t firstcall = 0;
    static int32_t ntim_int; /* # needed times for interp arrays */
    static gen_int_str *prof_int, *prof_tim;

    char *files[3];                     /* the 3 MET profile file names */
    int32_t n_met = 5;                  /* # met parms finally needed */
    int32_t nlvl = l1rec->l1file->nlvl; /* # levels in interpolation 42 for
                             GMAO, set in filehandle_init.c */
    int32_t itim, ilon, ilvl, npix, iprm, loc, t_interp, data_ix[2];
    float val, wt_t1, unc;
    double l_time, anc_times[3], lat, lon, last_time;
    anc_struc *anc_add;

    npix = l1rec->npix;
    /* initialize  for the met processing */
    if (firstcall == 0) {
        firstcall = 1;
        /*  get the files from l1rec->anc_profile1,2,3  */
        files[0] = input->anc_profile1;
        files[1] = input->anc_profile2;
        files[2] = input->anc_profile3;

        printf("\nOpening meteorological profile files.\n");
        printf("  anc_profile1 = %s\n", input->anc_profile1);
        printf("  anc_profile2 = %s\n", input->anc_profile2);
        printf("  anc_profile3 = %s\n", input->anc_profile3);
        printf("\n");

        /*  although done previously (for now), do again to order the
            important files */
        if ((ntim_int = anc_acq_f_stat(files, 1, 3)) == -1)
            return 1;

        /*  for a climatology (ntim_int = 0) just say that climatology
            profies not implemented yet */
        if (ntim_int == ANC_STAT_CLIM) {
            fprintf(stderr, "-E- %s %d: Ancillary met profile climatology not implemented yet\n", __FILE__,
                    __LINE__);
            return 1;
        }

        /*  set up the interpolation structure array */
        /*  note prof_int( ntim_int, nlvl, n_met ) with met running fastest */
        if ((prof_int = (gen_int_str *)malloc(n_met * nlvl * ntim_int * sizeof(gen_int_str))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate profile interpolation array\n", __FILE__,
                    __LINE__);
            return 1;
        }

        /* get GMAO profile set-up  */
        /* loop thru times / files from 1 - ntim_int */
        last_time = 0;
        for (itim = 0; itim < ntim_int; itim++) {
            prof_tim = prof_int + n_met * nlvl * itim;
            /* call anc_acq_gmao_prof_prep() *to get all parms for that file */
            if (anc_acq_gmao_prof_prep(files[itim], prof_tim, nlvl) != 0)
                return 1;
            if (last_time > prof_tim[0].anc_time) {
                fprintf(stderr, "-E- %s %d: met file times are out of sequence\n", __FILE__, __LINE__);
                return 1;
            }
            last_time = prof_tim[0].anc_time;
        }
    } /* END INITIALIZE */

    /*  set up the storage for the profile data if needed */
    anc_add = l1rec->anc_add;
    if (anc_add == NULL) {
        if (init_anc_add(l1rec) != 0)
            return 1;
    }

    /*
     * Start getting the data for the line in l1rec
     *  get the time of the current line
     */
    l_time = l1rec->scantime;

    /* go through the parameters, levels and interpolate in space and time */
    for (iprm = 0; iprm < n_met; iprm++) {
        for (ilvl = 0; ilvl < nlvl; ilvl++) {
            /* find the times and weighting to use for this line */
            for (itim = 0; itim < ntim_int; itim++) {
                loc = iprm + n_met * (ilvl + nlvl * itim);
                anc_times[itim] = prof_int[loc].anc_time;
            }

            /* get the proper times to use and weights */
            if (anc_acq_fnd_t_interp(l_time, anc_times, ntim_int, &t_interp, data_ix, &wt_t1) != 0)
                return 1;

            /* interpolate for 1 or 2 times */
            for (ilon = 0; ilon < npix; ilon++) {
                lon = l1rec->lon[ilon];
                lat = l1rec->lat[ilon];
                // gsl gets cranky with bad inputs, so...
                if (lon < -180.0 || lon > 180.0 || lat < -90.0 || lat > 90.0 || isnan(lat) || isnan(lon)) {
                    continue;
                }
                if (anc_acq_eval_pt(prof_int, iprm, ilvl, lat, lon, t_interp, data_ix, wt_t1, ntim_int, nlvl,
                                    n_met, &val, &unc) != 0) {
                    fprintf(stderr, "-E- %s %d: Error interpolating file: %s\n", __FILE__, __LINE__,
                            files[0]);
                    return -1;
                }
                /*  fill in the proper ancillary data slot */
                loc = ilvl + nlvl * ilon;
                switch (iprm) {
                    case TPROF:
                        l1rec->anc_add->prof_temp[loc] = val;
                        break;
                    case RHPROF:
                        l1rec->anc_add->prof_rh[loc] = val;
                        // set l1rec->rh to level 4 (index 3, 925mb), if valid and no user defined relhumid was passed
                        if ((input->relhumid != -2000) && (ilvl == 3) && (val != BAD_FLT))
                            l1rec->rh[ilon] = val;
                        break;
                    case HPROF:
                        l1rec->anc_add->prof_height[loc] = val;
                        break;
                    case QPROF:
                        l1rec->anc_add->prof_q[loc] = val;
                        break;
                    case O3PROF:
                        l1rec->anc_add->prof_o3[loc] = val;
                        break;
                }
            }
        }
    }
    /*  and end */
    return 0;
}

/**
 * @brief 
 * Interpolates RAD layers
 * @param l1rec 
 * @return int32_t 
 */
int32_t anc_acq_lin_rad(l1str *l1rec) {
    static int32_t firstcall = 0;
    static size_t ntim_int, ntime_step; /* # needed times for interp arrays */
    static gen_int_str *rad_int, *rad_tim;
    static float time_hours[72];
    char *files[3];   /* the 3 RAD  file names */
    size_t n_rad = 2; /* # rad parms finally needed */
    size_t npix = l1rec->npix;
    /* initialize  for the rad processing */
    if (firstcall == 0) {
        firstcall = 1;
        size_t file_count = 1;
        /*  get the files from l1rec->anc_profile1,2,3  */
        files[0] = input->cld_rad1;
        files[1] = input->cld_rad2;
        files[2] = input->cld_rad3;
        printf("\nOpening ancillary RAD files.\n");
        printf("  rad1 = %s\n", input->cld_rad1);
        printf("  rad2 = %s\n", input->cld_rad2);
        printf("  rad3 = %s\n", input->cld_rad3);
        printf("\n");
        size_t file_offset = 0;
        {
            char isodate1[32], isodate2[32], isodate3[32];
            memset(isodate1, '\0', 32);
            memset(isodate2, '\0', 32);
            memset(isodate3, '\0', 32);
            FILE_CHECK(input->cld_rad1)
            FILE_CHECK(input->cld_rad2)
            FILE_CHECK(input->cld_rad3)
            // opening the first file
            idDS rad_file = openDS(input->cld_rad1);
            {
                int dim_id;
                int status = nc_inq_dimid((rad_file).fid, "time", &dim_id);
                if (status == 0)
                    status = nc_inq_dimlen((rad_file).fid, dim_id, &ntime_step);
                if (status != 0) {
                    printf("--Error--: Dimension %s could not be read. See line %d in file %s.\nExiting...\n ", "time",__LINE__, __FILE__);
                    exit(EXIT_FAILURE);
                }
            }
            nc_get_att_text(rad_file.fid, NC_GLOBAL, "time_coverage_start", isodate1);
            endDS(rad_file);  // closing the first file
            // opening the second file
            rad_file = openDS(input->cld_rad2);
            nc_get_att_text(rad_file.fid, NC_GLOBAL, "time_coverage_start", isodate2);
            endDS(rad_file);  // closing the second file
            // opening the third file
            rad_file = openDS(input->cld_rad3);
            nc_get_att_text(rad_file.fid, NC_GLOBAL, "time_coverage_start", isodate3);
            endDS(rad_file);  // closing the third file

            double st_time1 = isodate2unix((const char *)isodate1);
            double st_time2 = isodate2unix((const char *)isodate2);
            double st_time3 = isodate2unix((const char *)isodate3);
            if (strcmp(isodate1, isodate2) != 0) {
                if (st_time1 > st_time2) {
                    printf("-Error-: Wrong time ordering for GMAO RAD inputs 1 and 2\n");
                    exit(EXIT_FAILURE);
                }
                file_count++;
            } else {
                file_offset++;  // start from the second file  because the first file is the same
            }
            if (strcmp(isodate3, isodate2) != 0) {
                if (st_time2 > st_time3) {
                    printf("-Error-: Wrong time ordering for GMAO RAD ancillary inputs 2 and 3\n");
                    exit(EXIT_FAILURE);
                }
                file_count++;
            }
            ntim_int = ntime_step * file_count;  // or 2 or 1
            // need to get ntim_int and ntime_step
        }
        /*  set up the interpolation structure array */
        if ((rad_int = (gen_int_str *)malloc(n_rad * ntim_int * sizeof(gen_int_str))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate profile interpolation array\n", __FILE__,
                    __LINE__);
            return 1;
        }
        for (int32_t ifile = 0; ifile < file_count; ifile++) {
            if (anc_acq_gmao_rad_prep(files[ifile + file_offset], rad_int, ifile, n_rad, ntime_step) != 0)
                return 1;
        }
        {
            float time_range[ntim_int];
            float zero_time_offset;
            switch (file_count) {
                case 1:
                    zero_time_offset = rad_int[0].anc_time;
                    break;
                case 2:
                    if (file_offset == 1)
                        zero_time_offset = rad_int[0].anc_time;
                    else
                        zero_time_offset = (rad_int + ntime_step * n_rad)[0].anc_time;
                    break;
                case 3:
                    zero_time_offset = (rad_int + ntime_step * n_rad)[0].anc_time;
                    break;
                default:
                    zero_time_offset = rad_tim[0].anc_time;
                    break;
            }
            for (int32_t ifile = 0; ifile < file_count; ifile++) {
                for (size_t itim = 0; itim < ntime_step; itim++) {
                    size_t time_slice_start = ntime_step * ifile;
                    gen_int_str *rad_tim = rad_int + (itim + time_slice_start) * n_rad;

                    time_range[time_slice_start + itim] = (rad_tim[0].anc_time - zero_time_offset) / 3600.0;
                }
            }

            for (size_t itim = 0; itim < ntim_int; itim++) {
                time_hours[itim] = time_range[itim];
            }
        }
    }
    if (l1rec->cld_rad == NULL)
        init_anc_cld_rad(l1rec, ntim_int, time_hours);
    /* go through the parameters and interpolate in space */
    for (size_t iprm = 0; iprm < n_rad; iprm++) {             // size_t iprm = 0; iprm < n_rad; iprm++
        for (size_t ip = 0; ip < npix; ip++) {                // size_t ip = 0; ip < npix; ip++
            for (size_t itim = 0; itim < ntim_int; itim++) {  // size_t itim = 0; itim < ntim_int; itim++
                float lon = l1rec->lon[ip];
                float lat = l1rec->lat[ip];
                float val = BAD_FLT;
                // bounds check
                if (lon < -180.0 || lon > 180.0 || lat < -90.0 || lat > 90.0 || isnan(lat) || isnan(lon))
                    continue;
                if (anc_rad_eval_pt(rad_int, iprm, itim, n_rad, lat, lon, &val) != 0) {
                    fprintf(stderr, "-E- %s %d: Error interpolating from file: %s\n", __FILE__, __LINE__,
                            files[itim % ntime_step]);
                    return -1;
                }
                // int32_t iscan = l1rec->iscan;
                switch (iprm) {
                    case TAUCLD:
                        l1rec->cld_rad->taucld[ip][itim] = val;
                        break;
                    case CFCLD:
                        l1rec->cld_rad->cfcld[ip][itim] = val;
                        break;
                    default:
                        break;
                }
            }
        }
    }

    return 0;
}
int32_t init_anc_aerosol(l1str *l1rec) {
    if ((l1rec->anc_aerosol = (anc_aer_struc *)malloc(sizeof(anc_aer_struc))) == NULL) {
        fprintf(stderr, "-E- %s %d: Unable to allocate additional ancillary data storage 1\n", __FILE__,
                __LINE__);
        return 1;
    }
    if (((l1rec->anc_aerosol->black_carbon_ext = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->black_carbon_scat = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->dust_ext = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->dust_scat = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->sea_salt_ext = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->sea_salt_scat = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->sulphur_ext = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->sulphur_scat = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->organic_carbon_ext = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->organic_carbon_scat = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->total_aerosol_ext = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->total_aerosol_scat = (float *)malloc(l1rec->npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_aerosol->total_aerosol_angstrom = (float *)malloc(l1rec->npix * sizeof(float))) ==
         NULL)) {
        fprintf(stderr, "-E- %s %d: Unable to allocate additional ancillary data storage 2\n", __FILE__,
                __LINE__);
        return 1;
    }
    return 0;
}

int32_t init_anc_cld_rad(l1str *l1rec, size_t times_dim, const float *time_range) {
    if ((l1rec->cld_rad = (cld_rad_struct *)malloc(sizeof(cld_rad_struct))) == NULL) {
        fprintf(stderr, "-E- %s %d: Unable to allocate additional ancillary data storage for RAD data 1\n",
                __FILE__, __LINE__);
        return 1;
    }
    const size_t npix = l1rec->npix;
    l1rec->cld_rad->cfcld = (float **)malloc(sizeof(float *) * npix);
    l1rec->cld_rad->taucld = (float **)malloc(sizeof(float *) * npix);
    l1rec->cld_rad->ntimes = times_dim;
    l1rec->cld_rad->timecldrange = (float *)malloc(sizeof(float) * times_dim);
    for (size_t itim = 0; itim < times_dim; itim++) {
        l1rec->cld_rad->timecldrange[itim] = time_range[itim];
    }
    for (size_t ip = 0; ip < npix; ip++) {
        if (((l1rec->cld_rad->cfcld[ip] = (float *)malloc(sizeof(float) * times_dim)) == NULL) ||
            ((l1rec->cld_rad->taucld[ip] = (float *)malloc(sizeof(float) * times_dim)) == NULL)) {
            fprintf(stderr,
                    "-E- %s %d: Unable to allocate additional ancillary data storage for RAD data 2\n",
                    __FILE__, __LINE__);
            return 1;
        }
    }
    l1rec->cld_rad->taucld[0][0] = -100;
    return 0;
}

int init_anc_add(l1str *l1rec)
/*-----------------------------------------------------------------------
 init_anc_add

 purpose: general storage allocation for added ancillary parameters

 Returns 0 if all checks are OK
 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 l1str *           l1rec            I      record containing the band-
                                        dependent geometry fields and
                                        the sizes in pixels, bands

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 Wayne Robinson    22 June 2018    Original development

 -----------------------------------------------------------------------*/
{
    int32_t npix, nlvl = 42;

    npix = l1rec->npix;

    if ((l1rec->anc_add = (anc_struc *)malloc(sizeof(anc_struc))) == NULL) {
        fprintf(stderr, "-E- %s %d: Unable to allocate additional ancillary data storage 1\n", __FILE__,
                __LINE__);
        return 1;
    }
    if (((l1rec->anc_add->prof_temp = (float *)malloc(nlvl * npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_add->prof_rh = (float *)malloc(nlvl * npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_add->prof_height = (float *)malloc(nlvl * npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_add->prof_q = (float *)malloc(nlvl * npix * sizeof(float))) == NULL) ||
        ((l1rec->anc_add->prof_o3 = (float *)malloc(nlvl * npix * sizeof(float))) == NULL)) {
        fprintf(stderr, "-E- %s %d: Unable to allocate additional ancillary data storage 2\n", __FILE__,
                __LINE__);
        return 1;
    }
    l1rec->anc_add->nlvl = nlvl;
    return 0;
}

int32_t anc_acq_lin_oz(l1str *l1rec)
/*******************************************************************

   anc_acq_lin_met

   purpose: process the ozone data for the l1rec, for now, just do the GMAO

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1str *           l1rec           I/O     all L1 line info

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       12 June 2018     Original development

 *******************************************************************/
{
    static int32_t firstcall = 0;
    static int32_t ntim_int; /* # needed times for interp arrays */
    static int32_t anc_id = -1;
    static gen_int_str *oz_int, *oz_tim;

    char *files[3]; /* the 3 OZONE file names */
    static char file_olci_str[FILENAME_MAX];
    char *file_olci = (char *)0; /* for olci file name, 0 if not olci */
    int32_t n_oz = 1;            /* # met parms finally needed */
    int32_t nlvl_int = 1;        /* # levels in interpolation 1 in the 2-d case */
    int32_t ilvl = 0, nlvl = 1;  /* for oz, only 1 level */
    int32_t itim, ilon, npix, iprm, t_interp, data_ix[2];
    float val, wt_t1, unc;
    double l_time, anc_times[3], lat, lon, last_time;

    /* initialize  for the ozone processing */
    if (firstcall == 0) {
        firstcall = 1;
        /*
         *  determine if rad data is OLCI and get tie point met file name if so
         */
        if (l1rec->l1file->format == FT_OLCI) {
            strcpy(file_olci_str, l1rec->l1file->name);
            char *ptr = strrchr(file_olci_str, '/');
            if (ptr) {
                *ptr = 0;
                strcat(file_olci_str, "/tie_meteo.nc");
            } else {
                strcpy(file_olci_str, "tie_meteo.nc");
            }
            file_olci = file_olci_str;
        }
        /*  get the files from l1rec->met1,2,3  */
        files[0] = input->ozone1;
        files[1] = input->ozone2;
        files[2] = input->ozone3;

        /*  although done previously (for now), do again to order the
            important files */
        if ((ntim_int = anc_acq_f_stat(files, 1, 3)) == -1)
            return 1;

        /*  set up the interpolation structure array */
        /*  note oz_int( ntim_int, nlvl_int, n_oz ) with oz running fastest */
        if ((oz_int = (gen_int_str *)malloc(n_oz * nlvl_int * ntim_int * sizeof(gen_int_str))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate met interpolation array\n", __FILE__, __LINE__);
            return 1;
        }
        /*  Check the type of the data, call anc_acq_ck
         *   (not needed as this is only GMAO
         */
        if ((anc_id = anc_acq_ck(files[0], file_olci)) == -1)
            return 1;

        /* if GMAO type, set up the GMAO met data */
        if (anc_id == ANC_SRC_TYP_GMAO) {
            /* loop thru times / files from 1 - ntim_int */
            last_time = 0;
            for (itim = 0; itim < ntim_int; itim++) {
                oz_tim = oz_int + itim * n_oz;
                /* call anc_acq_gmao_oz_prep() *to get all parms for that file */
                if (anc_acq_gmao_oz_prep(files[itim], oz_tim) != 0)
                    return 1;
                if (last_time > oz_tim[0].anc_time) {
                    fprintf(stderr, "-E- %s %d: ozone file times are out of sequence\n", __FILE__, __LINE__);
                    return 1;
                }
                last_time = oz_tim[0].anc_time;
            }
        } /* END GMAO */
        else {
            /* report an error - something is wrong - no other types now */
            fprintf(stderr, "-E- %s %d: Error reading file: %s\n", __FILE__, __LINE__, files[0]);
            return -1;
        }
    } /* END INITIALIZE */

    /*
     * Start getting the data for the line in l1rec
     *  get the time of the current line
     */
    l_time = l1rec->scantime;
    npix = l1rec->npix;

    /* go through the parameters and interpolate in space and time */
    for (iprm = 0; iprm < n_oz; iprm++) {
        /* find the times and weighting to use for this line */
        for (itim = 0; itim < ntim_int; itim++)
            anc_times[itim] = oz_int[iprm + n_oz * itim].anc_time;

        /* get the proper times to use and weights */
        if (anc_acq_fnd_t_interp(l_time, anc_times, ntim_int, &t_interp, data_ix, &wt_t1) != 0)
            return 1;

        /* interpolate for 1 or 2 times */
        if (input->ozone != -2000) {
            for (ilon = 0; ilon < npix; ilon++) {
                lon = l1rec->lon[ilon];
                lat = l1rec->lat[ilon];
                // gsl gets cranky with bad inputs, so...
                if (lon < -180.0 || lon > 180.0 || lat < -90.0 || lat > 90.0 || isnan(lat) || isnan(lon)) {
                    continue;
                }
                if (anc_acq_eval_pt(oz_int, iprm, ilvl, lat, lon, t_interp, data_ix, wt_t1, ntim_int, nlvl,
                                    n_oz, &val, &unc) != 0) {
                    fprintf(stderr, "-E- %s %d: Error interpolating file: %s\n", __FILE__, __LINE__,
                            files[0]);
                    return -1;
                }
                /*  fill in the proper ozone data slot converted from DU to
                    cm thickness*/
                l1rec->oz[ilon] = val / 1000.;
                if (l1rec->uncertainty)
                    l1rec->uncertainty->doz[ilon] = unc / 1000.;
            }
        }
    }
    /*  and end */
    return 0;
}

/**
 * @brief 
 * 
 * @param file - nc RAD filename
 * @param rad_int - interpolation data structure for RAD
 * @param ifile - file index
 * @param nrad - number of RAD variables (2 for now)
 * @param ntime_step - number of hours in a RAD file
 * @return int32_t 
 */


int32_t anc_acq_gmao_rad_prep(char *file, gen_int_str *rad_int, int32_t ifile, int32_t nrad,
                              int32_t ntime_step) {
    int32_t nlat, nlon, ntime, iv, nv;
    unsigned char *qual;
    double *lat_coord, *lon_coord, *ddata;
    int *time;
    double start_time;
    static float *data = 0;
    const char *rad_vars[] = {"TAUTOT", "CLDTOT"};
    const size_t n_rad_vars = nrad;
    for (size_t iprm = 0; iprm < n_rad_vars; iprm++) {
        if (anc_acq_read_gmao_rad(file, rad_vars[iprm], &data, &qual, &start_time, &ntime, &nlon, &nlat,
                                  &time, &lon_coord, &lat_coord) != 0)
            return 1;

        nv = nlon * nlat * ntime;
        if ((ddata = (double *)malloc(nv * sizeof(double))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate array ddata\n", __FILE__, __LINE__);
            return 1;
        }
        /* make the grid into an interpolation object
        have grid: data  lat_coord, lon_coord as scales nlon, nlat as lengths
       -> met_int[iprm].int_id  */
        size_t time_slice_start = ntime_step * ifile;
        if (ntime_step != ntime) {
            printf("Error reading time dimension\n");
            exit(EXIT_FAILURE);
        }
        for (iv = 0; iv < nv; iv++)
            ddata[iv] = data[iv];
        for (size_t itim = 0; itim < ntime_step; itim++)

        {
            gen_int_str *rad_tim = rad_int + (itim + time_slice_start) * n_rad_vars;
            rad_tim[iprm].accel_lat = gsl_interp_accel_alloc();
            rad_tim[iprm].accel_lon = gsl_interp_accel_alloc();
            rad_tim[iprm].int_id = gsl_spline2d_alloc(gsl_interp2d_bilinear, nlon, nlat);
            rad_tim[iprm].lat_coord = lat_coord;
            rad_tim[iprm].lon_coord = lon_coord;

            if (gsl_spline2d_init(rad_tim[iprm].int_id, rad_tim[iprm].lon_coord, rad_tim[iprm].lat_coord,
                                  ddata + itim * nlon * nlat, nlon, nlat) != 0) {
                fprintf(stderr, "-E- %s %d: GSL 2-D initialization failed, file: %s\n", __FILE__, __LINE__,
                        file);
                return 1;
            }

            /*  delete the parameter data array here
                and assign the qc, lon and lat coords to the int struct element
             */
            rad_tim[iprm].anc_time = start_time + time[itim] * 60.0;  // needs a proper conversion
            rad_tim[iprm].qual = qual;
            rad_tim[iprm].nlat = nlat;
            rad_tim[iprm].nlon = nlon;
        }
        free(data);
        free(ddata);
    }
    return 0;
}
int32_t anc_acq_gmao_met_prep(char *file, gen_int_str *met_int)
/*******************************************************************

   anc_acq_gmao_met_prep

   purpose: set up the interpolation objects and qc arrays for all the
     GMAO met parms for a file

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      GMAO single-level file
      gen_int_str *     met_int         I/O     array of interpolation
                                                structures for the met data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       7 May 2018      Original development
      W. Robinson       16 Jan 2019     adapt to read the T10M product for
                                        sfc T
      W. Robinson       3 Oct 2024      changes to put land/water ice over land/water pixels

 *******************************************************************/
{
    /* list of GMAO groups, parm names and special processing value */
    /* note that for the second RH, use the the 1st RH */
    char *ob_gmao_prm_nm[] = {"U10M", "V10M", "SLP", "TQV", "PS", "PS", "PS", "T10M", "FRSEAICE", "FRSNO"};
    /* also will need for rh: met, T10M and met, QV10M */
    /* also will need for icefr: lnd_ice, FRSNO (future) */
    int32_t n_raw_gmao = 10;
    int32_t iprm, nlat, nlon, nlvl, iv, nv;
    float *data2, *data3, *data_rh, p_lcl;
    unsigned char *qual, *qual2;
    double time, *lat_coord, *lon_coord, *ddata;
    double *lat_coord2, *lon_coord2;
    static float *data = 0;

    /* out prm            GMAO grp(s)  GMAO name(s)  do
       zwind (10 m) zw         met     U10M          as-is
       mwind (10 m) mw         met     V10M          as-is
       pressure (at MSL) pr    met     SLP           cvt Pa ->HPa
       water_vapor (PW) wv     met     TQV           as-is
       humidity (1000 mb) rh   met     PS            convert to RH with:
                               met     T10M
                               met     QV10M
       sfc pressure sfcp       met     PS            cvt Pa ->HPa
       RH sfc       sfcrh                            take above rh for now
       T surface    sfct       met     T10M           as-is
       ice fraction icefr      ocn_ice FRSEAICE      as-is
                               lnd_ice FRSNO         as-is
     */

    /* loop for output parms and get the data array(s) and do any special
       processing to make final param in final units */
    for (iprm = 0; iprm < n_raw_gmao; iprm++) {
        /* get the GMAO array */
        if (anc_acq_read_gmao(file, ob_gmao_prm_nm[iprm], &data, &qual, &time, &nlon,
                              &nlat, &nlvl, &lon_coord, &lat_coord) != 0)
            return 1;

        nv = nlon * nlat;
        if ((ddata = (double *)malloc(nv * sizeof(double))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate array ddata\n", __FILE__, __LINE__);
            return 1;
        }
        /* convert parm to needed type and combine parms to get field needed if
          necessary */
        switch (iprm) {
            case ZW:
            case MW:
            case WV:
                break;
            case PR:
                for (iv = 0; iv < nv; iv++) {
                    if (*(qual + iv) == 0) {
                        p_lcl = *(data + iv);
                        p_lcl = met_cvt_p_cvt(p_lcl, MET_UNITS__P_PA, MET_UNITS__P_HPA);
                        *(data + iv) = p_lcl;
                    }
                }
                break;
            case RH:
                if (anc_acq_read_gmao(file, "T10M", &data2, &qual2, &time, &nlon, &nlat, &nlvl,
                                      &lon_coord2, &lat_coord2) != 0)
                    return 1;
                free(lat_coord2);
                free(lon_coord2);
                for (iv = 0; iv < nv; iv++)
                    *(qual + iv) = *(qual + iv) | *(qual2 + iv);
                free(qual2);

                if (anc_acq_read_gmao(file, "QV10M", &data3, &qual2, &time, &nlon, &nlat, &nlvl,
                                      &lon_coord2, &lat_coord2) != 0)
                    return 1;
                free(lat_coord2);
                free(lon_coord2);
                for (iv = 0; iv < nv; iv++)
                    *(qual + iv) = *(qual + iv) | *(qual2 + iv);
                free(qual2);

                if ((data_rh = (float *)malloc(nlon * nlat * sizeof(float))) == NULL) {
                    fprintf(stderr, "-E- %s, %d: malloc failure\n", __FILE__, __LINE__);
                    return 1;
                }
                met_cvt_q_to_rh(nv, data, MET_UNITS__P_PA, data2, MET_UNITS__T_K, data3, MET_UNITS__Q_KG_KG,
                                data_rh);

                free(data2);
                free(data3);
                for (iv = 0; iv < nv; iv++)
                    *(data + iv) = *(data_rh + iv);
                free(data_rh);
                break;
            case SFCP:
                for (iv = 0; iv < nv; iv++) {
                    if (*(qual + iv) == 0) {
                        p_lcl = *(data + iv);
                        p_lcl = met_cvt_p_cvt(p_lcl, MET_UNITS__P_PA, MET_UNITS__P_HPA);
                        *(data + iv) = p_lcl;
                    }
                }
                break;
            case SFCRH:
                if (anc_acq_read_gmao(file, "T10M", &data2, &qual2, &time, &nlon, &nlat, &nlvl,
                                      &lon_coord2, &lat_coord2) != 0)
                    return 1;
                free(lat_coord2);
                free(lon_coord2);
                for (iv = 0; iv < nv; iv++)
                    *(qual + iv) = *(qual + iv) | *(qual2 + iv);
                free(qual2);

                if (anc_acq_read_gmao(file, "QV10M", &data3, &qual2, &time, &nlon, &nlat, &nlvl,
                                      &lon_coord2, &lat_coord2) != 0)
                    return 1;
                free(lat_coord2);
                free(lon_coord2);
                for (iv = 0; iv < nv; iv++)
                    *(qual + iv) = *(qual + iv) | *(qual2 + iv);
                free(qual2);

                if ((data_rh = (float *)malloc(nlon * nlat * sizeof(float))) == NULL) {
                    fprintf(stderr, "-E- %s, %d: malloc failure\n", __FILE__, __LINE__);
                    return 1;
                }
                met_cvt_q_to_rh(nv, data, MET_UNITS__P_PA, data2, MET_UNITS__T_K, data3, MET_UNITS__Q_KG_KG,
                                data_rh);

                free(data2);
                free(data3);
                for (iv = 0; iv < nv; iv++)
                    *(data + iv) = *(data_rh + iv);
                free(data_rh);
                break;
            case SFCT:
                for (iv = 0; iv < nv; iv++) {
                    if (*(qual + iv) == 0) {
                        p_lcl = *(data + iv);
                        p_lcl = met_cvt_t_cvt(p_lcl, MET_UNITS__T_K, MET_UNITS__T_C);
                        *(data + iv) = p_lcl;
                    }
                }
                break;
            case ICEFR_WTR:
                /* nothing to do with the ice over water */
                break;
            case ICEFR_LND:
                /* the over-water values are missing - we'll set then to 0 ice instead */
                for( iv=0; iv < nv; iv++ ) {
                    if( *(qual + iv) != 0 ) {
                        *(data + iv) = 0.;
                        *(qual + iv) = 0;
                    }
                }
                break;
            default:
                fprintf(stderr, "-E- %s %d: Unknown output identifier: %d\n", __FILE__, __LINE__, iprm);
        }
        /* make the grid into an interpolation object
        have grid: data  lat_coord, lon_coord as scales nlon, nlat as lengths
       -> met_int[iprm].int_id  */

        met_int[iprm].accel_lat = gsl_interp_accel_alloc();
        met_int[iprm].accel_lon = gsl_interp_accel_alloc();
        met_int[iprm].int_id = gsl_spline2d_alloc(gsl_interp2d_bilinear, nlon, nlat);
        met_int[iprm].lat_coord = lat_coord;
        met_int[iprm].lon_coord = lon_coord;

        for (iv = 0; iv < nv; iv++)
            ddata[iv] = data[iv];
        if (gsl_spline2d_init(met_int[iprm].int_id, met_int[iprm].lon_coord, met_int[iprm].lat_coord, ddata,
                              nlon, nlat) != 0) {
            fprintf(stderr, "-E- %s %d: GSL 2-D initialization failed, file: %s\n", __FILE__, __LINE__, file);
            return 1;
        }

        /*  delete the parameter data array here
            and assign the qc, lon and lat coords to the int struct element
         */
        met_int[iprm].anc_time = time;
        met_int[iprm].qual = qual;
        met_int[iprm].nlat = nlat;
        met_int[iprm].nlon = nlon;
        free(data);
        free(ddata);
    }
    return 0;
}

int32_t anc_acq_gmao_prof_prep(char *file, gen_int_str *prof_int, int32_t nlvl_expect)
/*******************************************************************

   anc_acq_gmao_prof_prep

   purpose: set up the interpolation objects and qc arrays for all the
     GMAO profile parms for a file (at 1 time)

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      GMAO single-level file
      gen_int_str *     prof_int        I/O     array of interpolation
                                                structures for the met
                                                profile data
      int32_t           nlvl_expect      I      expected # profile levels

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       20 June 2018    Original development

 *******************************************************************/
{
    /* list of GMAO groups, parm names and special processing value */
    /* note that for the second RH, use the the 1st RH */
    char *ob_gmao_prm_nm[] = {"T", "RH", "H", "QV", "O3"};
    int32_t n_raw_gmao = 5;
    int32_t iprm, nlat, nlon, nlvl, ilvl, loc, iv, nv, ntot;
    float p_lcl;
    unsigned char *qual;
    double time, *lat_coord, *lon_coord, *ddata;
    static float *data = 0;

    /* out prm            GMAO grp(s)  GMAO name(s)  do
       temp profile  prof_temp NONE   T             convert K -> C
       RH profile    prof_rh   NONE   RH            convert fraction to %
       height profile  prof_h  NONE   H             as-is
       specific humidity profile  prof_q  NONE   QV   as-is
       ozone profile  prof_o3  NONE   O3            as-is (unless asked)
     */

    /* loop for output parms and get the data array(s) and do any special
       processing to make final param in final units */
    for (iprm = 0; iprm < n_raw_gmao; iprm++) {
        /* get the GMAO array */
        if (anc_acq_read_gmao(file, ob_gmao_prm_nm[iprm], &data, &qual, &time, &nlon, &nlat, &nlvl,
                              &lon_coord, &lat_coord) != 0)
            return 1;

        /* verify standard # levels  */
        if (nlvl != nlvl_expect) {
            fprintf(stderr, "-E- %s %d: unexpected # profile levels: %d were read from file: %s\n", __FILE__,
                    __LINE__, nlvl, file);
            return 1;
        }
        /* convert temperature parm to degrees C  */
        if (iprm == TPROF) {
            ntot = nlon * nlat * nlvl;
            for (iv = 0; iv < ntot; iv++) {
                if (*(qual + iv) == 0) {
                    p_lcl = *(data + iv);
                    p_lcl = met_cvt_t_cvt(p_lcl, MET_UNITS__T_K, MET_UNITS__T_C);
                    *(data + iv) = p_lcl;
                }
            }
        } /* make the RH in % */
        else if (iprm == RHPROF) {
            ntot = nlon * nlat * nlvl;
            for (iv = 0; iv < ntot; iv++) {
                if (*(qual + iv) == 0) {
                    *(data + iv) *= 100.;
                }
            }
        }
        /*  loop through each level to make separate interpolation objects */
        nv = nlon * nlat;
        if ((ddata = (double *)malloc(nv * sizeof(double))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate array ddata\n", __FILE__, __LINE__);
            return 1;
        }
        for (ilvl = 0; ilvl < nlvl; ilvl++) {
            loc = iprm + n_raw_gmao * ilvl;
            /* make the grid into an interpolation object
              have grid: data  lat_coord, lon_coord as scales nlon, nlat as lengths
             -> prof_int[loc].int_id  */

            prof_int[loc].accel_lat = gsl_interp_accel_alloc();
            prof_int[loc].accel_lon = gsl_interp_accel_alloc();
            prof_int[loc].int_id = gsl_spline2d_alloc(gsl_interp2d_bilinear, nlon, nlat);
            prof_int[loc].lat_coord = lat_coord;
            prof_int[loc].lon_coord = lon_coord;

            for (iv = 0; iv < nv; iv++)
                ddata[iv] = data[iv + nv * ilvl];
            if (gsl_spline2d_init(prof_int[loc].int_id, prof_int[loc].lon_coord, prof_int[loc].lat_coord,
                                  ddata, nlon, nlat) != 0) {
                fprintf(stderr, "-E- %s %d: GSL 2-D initialization failed, file: %s\n", __FILE__, __LINE__,
                        file);
                return 1;
            }
            prof_int[loc].anc_time = time;
            prof_int[loc].qual = qual + nv * ilvl;
            prof_int[loc].nlat = nlat;
            prof_int[loc].nlon = nlon;
        }

        /*  delete the parameter data array here
            and assign the qc, lon and lat coords to the int struct element
         */
        free(data);
        free(ddata);
    }
    return 0;
}

int32_t anc_acq_gmao_oz_prep(char *file, gen_int_str *oz_int)
/*******************************************************************

   anc_acq_gmao_oz_prep

   purpose: set up the interpolation objects and qc arrays for all the
     GMAO ozone parms for a file

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      GMAO single-level file
      gen_int_str *     oz_int          I/O     array of interpolation
                                                structures for the ozone data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       12 Jun 2018     Original development

 *******************************************************************/
{
    /* Only the GMAO TO3 is needed from the 'met' file */
    int32_t nlat, nlon, nlvl, iv, nv;
    int32_t iprm = 0;
    unsigned char *qual;
    double time, *lat_coord, *lon_coord, *ddata;
    static float *data = 0;

    /* read the ozone */
    if (anc_acq_read_gmao(file, "TO3", &data, &qual, &time, &nlon, &nlat, &nlvl, &lon_coord,
                          &lat_coord) != 0)
        return 1;

    nv = nlon * nlat;
    if ((ddata = (double *)malloc(nv * sizeof(double))) == NULL) {
        fprintf(stderr, "-E- %s %d: Unable to allocate array ddata\n", __FILE__, __LINE__);
        return 1;
    }
    /* make the grid into an interpolation object
       have grid: data  lat_coord, lon_coord as scales nlon, nlat as lengths
       -> oz_int[iprm].int_id  */

    oz_int[iprm].accel_lat = gsl_interp_accel_alloc();
    oz_int[iprm].accel_lon = gsl_interp_accel_alloc();
    oz_int[iprm].int_id = gsl_spline2d_alloc(gsl_interp2d_bilinear, nlon, nlat);
    oz_int[iprm].lat_coord = lat_coord;
    oz_int[iprm].lon_coord = lon_coord;

    for (iv = 0; iv < nv; iv++)
        ddata[iv] = data[iv];
    if (gsl_spline2d_init(oz_int[iprm].int_id, oz_int[iprm].lon_coord, oz_int[iprm].lat_coord, ddata, nlon,
                          nlat) != 0) {
        fprintf(stderr, "-E- %s %d: GSL 2-D initialization failed, file: %s\n", __FILE__, __LINE__, file);
        return 1;
    }
    /*  delete the parameter data array here
     and assign the qc, lon and lat coords to the int struct element
     */
    oz_int[iprm].anc_time = time;
    oz_int[iprm].qual = qual;
    oz_int[iprm].nlat = nlat;
    oz_int[iprm].nlon = nlon;
    free(data);
    free(ddata);
    /*  and end */
    return 0;
}

int32_t anc_acq_gmao_aer_prep(char *file, gen_int_str *aer_int) {
    /* Only the GMAO TO3 is needed from the 'met' file */
    int32_t nlat, nlon, nlvl, iv, nv;
    int32_t iprm = 0;
    unsigned char *qual;
    double time, *lat_coord, *lon_coord, *ddata;
    static float *data = 0;

    /* list of GMAO groups, parm names  */
    char *ob_gmao_prm_nm[] = {"BCEXTTAU",  "BCSCATAU",  "DUEXTTAU", "DUSCATAU", "SSEXTTAU",
                              "SSSCATAU",  "SUEXTTAU",  "SUSCATAU", "OCEXTTAU", "OCSCATAU",
                              "TOTEXTTAU", "TOTSCATAU", "TOTANGSTR"};

    int32_t n_raw_gmao = 13;
    /* read the aerosol parameters */
    for (iprm = 0; iprm < n_raw_gmao; iprm++) {
        if (anc_acq_read_gmao(file, ob_gmao_prm_nm[iprm], &data, &qual, &time, &nlon, &nlat, &nlvl,
                              &lon_coord, &lat_coord) != 0)
            return 1;

        nv = nlon * nlat;
        if ((ddata = (double *)malloc(nv * sizeof(double))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate array ddata\n", __FILE__, __LINE__);
            return 1;
        }
        /* make the grid into an interpolation object
           have grid: data  lat_coord, lon_coord as scales nlon, nlat as lengths
           -> aer_int[iprm].int_id  */

        aer_int[iprm].accel_lat = gsl_interp_accel_alloc();
        aer_int[iprm].accel_lon = gsl_interp_accel_alloc();
        aer_int[iprm].int_id = gsl_spline2d_alloc(gsl_interp2d_bilinear, nlon, nlat);
        aer_int[iprm].lat_coord = lat_coord;
        aer_int[iprm].lon_coord = lon_coord;

        for (iv = 0; iv < nv; iv++)
            ddata[iv] = data[iv];
        if (gsl_spline2d_init(aer_int[iprm].int_id, aer_int[iprm].lon_coord, aer_int[iprm].lat_coord, ddata,
                              nlon, nlat) != 0) {
            fprintf(stderr, "-E- %s %d: GSL 2-D initialization failed, file: %s\n", __FILE__, __LINE__, file);
            return 1;
        }
        /*  delete the parameter data array here
         and assign the qc, lon and lat coords to the int struct element
         */
        aer_int[iprm].anc_time = time;
        aer_int[iprm].qual = qual;
        aer_int[iprm].nlat = nlat;
        aer_int[iprm].nlon = nlon;
        free(data);
        free(ddata);
        /*  and end */
    }
    return 0;
}

int32_t anc_rad_eval_pt(gen_int_str *rad_int, int32_t iprm, int32_t itim, int32_t nrad, float lat, float lon,
                        float *val) {
    gen_int_str *rad_tim = rad_int + (itim)*nrad;
    double val_get;
    gsl_interp_accel *accel_lon = rad_tim[iprm].accel_lon;
    gsl_interp_accel *accel_lat = rad_tim[iprm].accel_lat;
    if (gsl_spline2d_eval_e(rad_tim[iprm].int_id, lon, lat, accel_lon, accel_lat, &val_get) != 0) {
        fprintf(stderr, "-E- %s %d: gsl_spline2d_eval_e error: %d, %d \n", __FILE__, __LINE__, itim, iprm);
        return -1;
    }
    size_t nlat = rad_tim[iprm].nlat;
    size_t nlon = rad_tim[iprm].nlon;
    /*  make sure the value is OK and make bad if not */
    size_t ilat = gsl_interp_accel_find(accel_lat, rad_tim[iprm].lat_coord, nlat, lat);
    size_t ilon = gsl_interp_accel_find(accel_lon, rad_tim[iprm].lon_coord, nlon, lon);
    if ((rad_tim[iprm].qual[ilon + nlon * ilat] == 1) || (rad_tim[iprm].qual[ilon + 1 + nlon * ilat] == 1) ||
        (rad_tim[iprm].qual[ilon + nlon * (ilat + 1)] == 1) ||
        (rad_tim[iprm].qual[ilon + 1 + nlon * (ilat + 1)] == 1))
        val_get = BAD_FLT;
    *val = (float)val_get;
    return 0;
}

int32_t anc_acq_eval_pt(gen_int_str *met_int, int32_t iprm, int32_t ilvl, float lat, float lon,
                        int32_t t_interp, int32_t *data_ix, float wt_t1, int32_t ntim_int, int32_t nlvl,
                        int32_t nprm, float *final_val, float *unc)
/*******************************************************************

   anc_acq_eval_pt

   purpose: find a parameter for a certain lat, lon by interpolating in
     space and time

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      gen_int_str *     met_int          I      array of all interpolation
                                                info for parameters, levels,
                                                times
      int32_t           iprm             I      parameter index
      int32_t           ilvl             I      level index
      float             lat              I      latitude
      float             lon              I      longitude
      int32_t           t_interp         I      time interpolation indicator:
                                                1 - interpolate, 0 - none
      int32_t           data_ix          I      indicies of first and possibly
                                                2nd times to interpolate
      float             wt_t1            I      weight to place on 1st time
      int32_t           iprm             I      parameter to work on
      int32_t           ilvl             I      level to work on
      int32_t           ntim_int         I      total # times available
      int32_t           nlvl             I      # levels available
      int32_t           nprm             I      # parameters available
      float *           final_val        O      final interpolated value
      float *           unc              O      Uncertainty in value
                                                (absolute diff in values at
                                                2 times used in interpolation)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       7 May 2018      Original development

 *******************************************************************/
{
    float wt_t2;
    gsl_interp_accel *accel_lat, *accel_lon;
    int32_t met_ptr, tim_ix, itim, ntim;
    double val;
    int32_t ilat, ilon, nlon, nlat;

    wt_t2 = 1. - wt_t1;
    /*  the # times is just the interpolation indicator + 1 */
    ntim = t_interp + 1;
    *unc = 0.;
    for (itim = 0; itim < ntim; itim++) {
        tim_ix = data_ix[itim];
        met_ptr = iprm + nprm * (ilvl + nlvl * tim_ix);
        accel_lon = met_int[met_ptr].accel_lon;
        accel_lat = met_int[met_ptr].accel_lat;
        nlat = met_int[met_ptr].nlat;
        nlon = met_int[met_ptr].nlon;

        if (gsl_spline2d_eval_e(met_int[met_ptr].int_id, lon, lat, accel_lon, accel_lat, &val) != 0) {
            fprintf(stderr, "-E- %s %d: gsl_spline2d_eval_e error: %d, %d, %d\n", __FILE__, __LINE__, itim,
                    iprm, ilvl);
            return -1;
        }
        /*  make sure the value is OK and make bad if not */
        ilat = gsl_interp_accel_find(accel_lat, met_int[met_ptr].lat_coord, nlat, lat);
        ilon = gsl_interp_accel_find(accel_lon, met_int[met_ptr].lon_coord, nlon, lon);

        if ((met_int[met_ptr].qual[ilon + nlon * ilat] == 1) ||
            (met_int[met_ptr].qual[ilon + 1 + nlon * ilat] == 1) ||
            (met_int[met_ptr].qual[ilon + nlon * (ilat + 1)] == 1) ||
            (met_int[met_ptr].qual[ilon + 1 + nlon * (ilat + 1)] == 1))
            val = BAD_FLT;
        /*  for the  case where time interpolation is done, get the other
            time's interpolated value */
        if (itim == 1) {
            if (*final_val == BAD_FLT) {
                *final_val = (val == BAD_FLT) ? BAD_FLT : val;
            } else if (val == BAD_FLT) {
                *final_val = (*final_val == BAD_FLT) ? BAD_FLT : *final_val;
            } else {
                *unc = fabsf(*final_val - (float)val);
                *final_val = *final_val * wt_t1 + val * wt_t2;
            }
        } else {
            *final_val = val;
        }
    }
    if (*final_val == BAD_FLT)
        *unc = BAD_FLT;
    return 0;
}

int32_t anc_acq_fnd_t_interp(double s_time, double *anc_time, int32_t anc_f_stat, int32_t *t_interp,
                             int32_t *data_ix, float *wt)
/*******************************************************************

   anc_acq_fnd_t_interp

   purpose: determine the anc times to use and their weight from the anc
     times and the scan time

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      double            s_time           I      scan time
      double *          anc_time         I      times for ancillary data
      int32_t           anc_f_stat       I      ancillary time status
                                                (ANC_STAT...) with 1-3 being
                                                the number of active anc times
                                                and 0, >3 other possible
                                                configurations
      int32_t *         t_interp         O      interpolation flag: 0 - none
                                                1 - use the 2 files pointed to
      int32_t *         data_ix          O      indecies of 1st, possibly 2nd
                                                anc time to use
      float *           wt               O      for t_interp = 1, the weights
                                                of time 1 to use for
                                                interpolation. weight for 2nd
                                                time = 1. - wt
  NOTE that there could be a anc_t_warn set (1) if we still want to use the
  early definitions for bad anc.  I'll comment those cases but not impliment

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       7 May 2018      Original development

 *******************************************************************/
{
    int32_t data1_ix = 0, data2_ix = 0;

    if (anc_f_stat == ANC_STAT_1T) {
        /* one time only available */
        *t_interp = 0;
        data1_ix = 0;
    } else if (anc_f_stat == ANC_STAT_2T_START) {
        /* 2 different times in the 1, 2 positions */
        if (s_time < anc_time[0]) {
            /* scan is earlier than 1st time */
            *t_interp = 0;
            data1_ix = 0;
            /* anc_t_warn = 1 */
        } else if (s_time > anc_time[1]) {
            /* use anc 2 only */
            *t_interp = 0;
            data1_ix = 1;
        } else {
            /* in-between anc 1, 2 use data 0, 1 and time interpolate */
            *t_interp = 1;
            data1_ix = 0;
            data2_ix = 1;
            *wt = (anc_time[1] - s_time) / (anc_time[1] - anc_time[0]);
        }
    } else if (anc_f_stat == ANC_STAT_2T_END) {
        if (s_time < anc_time[0]) {
            /* outside on the low end, use 1st t only */
            *t_interp = 0;
            data1_ix = 0;
        } else if (s_time > anc_time[2]) {
            /* beyond the high end, use 2nd time only */
            *t_interp = 0;
            data1_ix = 2;
            /* anc_t_warn = 1 */
        } else {
            /* between the MET 1 and 3 */
            *t_interp = 1;
            data1_ix = 0;
            data2_ix = 2;
            *wt = (anc_time[2] - s_time) / (anc_time[2] - anc_time[0]);
        }
    } else if (anc_f_stat == ANC_STAT_3T) {
        if (s_time < anc_time[0]) {
            *t_interp = 0;
            data1_ix = 0;
            /* anc_t_warn = 1 */
        } else if (s_time > anc_time[2]) {
            *t_interp = 0;
            data1_ix = 2;
            /* anc_t_warn = 1 */
        } else if (s_time < anc_time[1]) {
            /* between times 0 and 1 */
            *t_interp = 1;
            data1_ix = 0;
            data2_ix = 1;
            *wt = (anc_time[1] - s_time) / (anc_time[1] - anc_time[0]);
        } else {
            /* what's left: between data 1 and 2 */
            *t_interp = 1;
            data1_ix = 1;
            data2_ix = 2;
            *wt = (anc_time[2] - s_time) / (anc_time[2] - anc_time[1]);
        }
    } else {
        /* this should not happen at this time - a status that is either
         climatology or not defined */
        printf("%s, %d: Undefined anc_f_stat - should not happen\n", __FILE__, __LINE__);
        return -1;
    }
    *data_ix = data1_ix;
    *(data_ix + 1) = data2_ix;

    return 0;
}
/**
 * @brief 
 *  READ data from a RAD file
 * @param file - file name
 * @param var_name - variable name
 * @param data - data to save
 * @param qa - quality flags
 * @param start_time -start time
 * @param ntime - number of hours in RAD file
 * @param nlon - size of lon 
 * @param nlat - size of lat
 * @param time - time array (size of ntime)
 * @param lon_coord - lon coordinates
 * @param lat_coord - lat coordinates
 * @return int32_t 
 */
int32_t anc_acq_read_gmao_rad(char *file, const char *var_name, float **data, unsigned char **qa,
                              double *start_time, int32_t *ntime, int32_t *nlon, int32_t *nlat, int **time,
                              double **lon_coord, double **lat_coord) {
    static char lcl_fil[FILENAME_MAX] = "";
    static int ncid;
    int var_id, var_id_lat, var_id_lon, var_id_time;
    int32_t loc;
    int32_t nv, nlon_pre, il, ip;
    float *data_tmp, fillv, missv;
    /* if input file is different, close old file (if open) and open new file */
    if (strcmp(file, lcl_fil) != 0) {
        if (lcl_fil[0] != 0) {
            if (nc_close(ncid) != NC_NOERR) {
                fprintf(stderr, "-E- %s %d: nc_close error, file: %s\n", __FILE__, __LINE__, file);
                return 1;
            }
        }
        if (nc_open(file, 0, &ncid) != NC_NOERR) {
            fprintf(stderr, "-E- %s %d: file: %s is not netcdf, not acceptable GMAO file\n", __FILE__,
                    __LINE__, file);
            return -1;
        }
        strcpy(lcl_fil, file);
    }
    /*  get the main dimensions */
    if (((*nlat = ncio_dim_siz(ncid, "lat")) == -1) || ((nlon_pre = ncio_dim_siz(ncid, "lon")) == -1) ||
        ((*ntime = ncio_dim_siz(ncid, "time")) == -1)) {
        fprintf(stderr, "-E- %s %d: file: %s error reading lat, lon and time dimension size\n", __FILE__,
                __LINE__, file);
        return 1;
    }
    *nlon = nlon_pre + 1;
    nv = *nlat * *nlon * *ntime;
    char isodate[32];
    memset(isodate, '\0', 32);
    nc_get_att_text(ncid, NC_GLOBAL, "time_coverage_start", isodate);
    *start_time = isodate2unix((const char *)isodate);
    /* set up data and lon, lat arrays using the sizes */
    if (((*data = (float *)malloc(nv * sizeof(float))) == NULL) ||
        ((*lon_coord = (double *)malloc(*nlon * sizeof(double))) == NULL) ||
        ((*lat_coord = (double *)malloc(*nlat * sizeof(double))) == NULL) ||
        ((*time = (int *)malloc(*ntime * sizeof(int))) == NULL) ||
        ((*qa = (unsigned char *)calloc(nv, sizeof(char))) == NULL) ||
        ((data_tmp = (float *)malloc(nlon_pre * *nlat * *ntime * sizeof(float))) == NULL)) {
        fprintf(stderr, "-E- %s %d: file: %s error allocating gmao parameter allocation\n", __FILE__,
                __LINE__, file);
        return 1;
    }
    if ((nc_inq_varid(ncid, var_name, &var_id) != NC_NOERR) ||
        (nc_inq_varid(ncid, "lat", &var_id_lat) != NC_NOERR) ||
        (nc_inq_varid(ncid, "lon", &var_id_lon) != NC_NOERR) ||
        (nc_inq_varid(ncid, "time", &var_id_time) != NC_NOERR)) {
        fprintf(stderr, "-E- %s %d: file: %s error setting an id for product: %s\n", __FILE__, __LINE__, file,
                var_name);
        return 1;
    }
    if (((nc_get_var_double(ncid, var_id_lat, *lat_coord)) != NC_NOERR) ||
        ((nc_get_var_double(ncid, var_id_lon, *lon_coord)) != NC_NOERR) ||
        ((nc_get_var_int(ncid, var_id_time, *time)) != NC_NOERR)) {
        fprintf(stderr, "-E- %s %d: file: %s error reading the scales and parameter\n", __FILE__, __LINE__,
                file);
        return 1;
    }
    /*
     *  Get the data and its bad values
     */
    if ((nc_get_att_float(ncid, var_id, "_FillValue", &fillv) != NC_NOERR) ||
        (nc_get_att_float(ncid, var_id, "missing_value", &missv) != NC_NOERR)) {
        printf("%s, %d: nc_get_att_float error for parm %s\n", __FILE__, __LINE__, var_name);
        return -1;
    }

    if (nc_get_var_float(ncid, var_id, data_tmp) != NC_NOERR) {
        printf("%s, %d: nc_get_var_float error for parm %s\n", __FILE__, __LINE__, var_name);
        return -1;
    }
    /*
     *  read the dataset and replace either fill or missing with bad value
     */

    for (int32_t it = 0; it < *ntime; it++) {
        for (il = 0; il < *nlat; il++) {
            for (ip = 0; ip < *nlon; ip++) {
                loc = ip + *nlon * (il + *nlat * it);
                if (ip < nlon_pre) {
                    *(*data + loc) = *(data_tmp + ip + nlon_pre * (il + *nlat * it));
                } else {
                    *(*data + loc) = *(data_tmp + nlon_pre * (il + *nlat * it));
                }
                if ((*(*data + loc) == missv) || (*(*data + loc) == fillv)) {
                    *(*data + loc) = BAD_FLT;
                    *(*qa + loc) = 1;
                }
            }
        }
    }
    free(data_tmp);
    (*lon_coord)[nlon_pre] = 180.;
    return 0;
}
int32_t anc_acq_read_gmao(char *file, char *ds_name,
                          float **data, unsigned char **qa, double *time, int32_t *nlon, int32_t *nlat,
                          int32_t *nlvl, double **lon_coord, double **lat_coord)
/*******************************************************************
   anc_acq_read_gmao

   purpose: read in the (OBPG modified or not) GMAO parameter from the
     specified file.  Also, add a column so longitude runs -180 -> 180

   Returns type: int32_t - status 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      GMAO single-level file
      char *            ds_name          I      name of dataset to read
      float **          data             O      array of parameter data,
                                                BAD_FLT for bad values
      unsigned char **  qa               O      quality of data 0 good, 1 bad
      double *          time             O      time for this data
      int32_t *         nlon             O      # longitudes
      int32_t *         nlat             O      # latitudes
      int32_t *         nlvl             O      # levels
      double **         lon_coord        O      array of longitude values
      double **         lat_coord        O      array of latitude values

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       7 May 2018      Original development

  Note that the level coordinates are not returned (though they could easily be)
  - the GMAO data has a fairly fixed definition of what the pressure levels is

 *******************************************************************/
{
    static char lcl_fil[FILENAME_MAX] = "";
    static int ncid;
    int var_id, var2_id, var3_id, dim_id;
    int32_t loc;
    int32_t nv, ilvl, nlon_pre, il, ip;
    size_t tlvl;
    float *data_tmp, fillv, missv;

    /* if input file is different, close old file (if open) and open new file */
    if (strcmp(file, lcl_fil) != 0) {
        if (lcl_fil[0] != 0) {
            if (nc_close(ncid) != NC_NOERR) {
                fprintf(stderr, "-E- %s %d: nc_close error, file: %s\n", __FILE__, __LINE__, file);
                return 1;
            }
        }
        if (nc_open(file, 0, &ncid) != NC_NOERR) {
            fprintf(stderr, "-E- %s %d: file: %s is not netcdf, not acceptable GMAO file\n", __FILE__,
                    __LINE__, file);
            return -1;
        }
        strcpy(lcl_fil, file);
    }
    /*  get the main dimensions */
    if (((*nlat = ncio_dim_siz(ncid, "lat")) == -1) || ((nlon_pre = ncio_dim_siz(ncid, "lon")) == -1)) {
        fprintf(stderr, "-E- %s %d: file: %s error reading lat, lon dimension size\n", __FILE__, __LINE__,
                file);
        return 1;
    }
    /* separate level treatment */
    if (nc_inq_dimid(ncid, "lev", &dim_id) != NC_NOERR) {
        *nlvl = 1;
    } else {
        if (nc_inq_dimlen(ncid, dim_id, &tlvl) != NC_NOERR) {
            fprintf(stderr, "-E- %s %d: file: %s error reading level dimension size\n", __FILE__, __LINE__,
                    file);
            return 1;
        }
        *nlvl = tlvl;
    }
    *nlon = nlon_pre + 1;
    nv = *nlat * *nlon * *nlvl;

    /* get the time for this data */
    char isodate[32];
    memset(isodate, '\0', 32);
    nc_get_att_text(ncid, NC_GLOBAL, "time_coverage_start", isodate);
    *time = isodate2unix((const char *)isodate);

    /* set up data and lon, lat arrays using the sizes */
    if (((*data = (float *)malloc(nv * sizeof(float))) == NULL) ||
        ((*lon_coord = (double *)malloc(*nlon * sizeof(double))) == NULL) ||
        ((*lat_coord = (double *)malloc(*nlat * sizeof(double))) == NULL) ||
        ((*qa = (unsigned char *)calloc(nv, sizeof(char))) == NULL) ||
        ((data_tmp = (float *)malloc(nlon_pre * *nlat * *nlvl * sizeof(float))) == NULL)) {
        fprintf(stderr, "-E- %s %d: file: %s error allocating gmao parameter allocation\n", __FILE__,
                __LINE__, file);
        return 1;
    }

    /* read the data and the lon, lat values */
    if ((nc_inq_varid(ncid, ds_name, &var_id) != NC_NOERR) ||
        (nc_inq_varid(ncid, "lat", &var2_id) != NC_NOERR) ||
        (nc_inq_varid(ncid, "lon", &var3_id) != NC_NOERR)) {
        fprintf(stderr, "-E- %s %d: file: %s error setting an id for product: %s\n", __FILE__, __LINE__, file,
                ds_name);
        return 1;
    }

    if (((nc_get_var_double(ncid, var2_id, *lat_coord)) != NC_NOERR) ||
        ((nc_get_var_double(ncid, var3_id, *lon_coord)) != NC_NOERR)) {
        fprintf(stderr, "-E- %s %d: file: %s error reading the scales and parameter\n", __FILE__, __LINE__,
                file);
        return 1;
    }
    /*
     *  Get the data and its bad values
     */
    if ((nc_get_att_float(ncid, var_id, "_FillValue", &fillv) != NC_NOERR) ||
        (nc_get_att_float(ncid, var_id, "missing_value", &missv) != NC_NOERR)) {
        printf("%s, %d: nc_get_att_float error for parm %s\n", __FILE__, __LINE__, ds_name);
        return -1;
    }
    /*
     *  read the dataset and replace either fill or missing with bad value
     */
    if (nc_get_var_float(ncid, var_id, data_tmp) != NC_NOERR) {
        printf("%s, %d: nc_get_var_float error for parm %s\n", __FILE__, __LINE__, ds_name);
        return -1;
    }
    for (ilvl = 0; ilvl < *nlvl; ilvl++) {
        for (il = 0; il < *nlat; il++) {
            for (ip = 0; ip < *nlon; ip++) {
                loc = ip + *nlon * (il + *nlat * ilvl);
                if (ip < nlon_pre) {
                    *(*data + loc) = *(data_tmp + ip + nlon_pre * (il + *nlat * ilvl));
                } else {
                    *(*data + loc) = *(data_tmp + nlon_pre * (il + *nlat * ilvl));
                }
                if ((*(*data + loc) == missv) || (*(*data + loc) == fillv)) {
                    *(*data + loc) = BAD_FLT;
                    *(*qa + loc) = 1;
                }
            }
        }
    }

    free(data_tmp);
    (*lon_coord)[nlon_pre] = 180.;
    /*
     *  All is set up, return the data
     */
    return 0;
}

/*
int32_t anc_acq_fin()
 *******************************************************************
   anc_acq_fin

   purpose: finish the ancillary processing

   Returns type: int32_t - return status 0 is OK

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     char **           files            I      list of files of the ancillary

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       17 May 2018     original development

 *******************************************************************/

/*
Mainly, free the int objects, lat, lon arrays in the interpolation struct and
also the interpolation structure
 */

int32_t anc_acq_ecmwf_init(char **files, char **prm_nm, int n_prm, int32_t sto_ix)
/*******************************************************************

   anc_acq_ecmwf_init

   purpose: Identify the incoming MET or OZONE file set and if netcdf,
   set up the ancillary data so it can be accessed by anc_acq_line.
   Otherwise, do nothing and have getanc.c process the NCEP/TOMS data.
   For now, only ECMWF netcdf files can be processed which contain
   only 1 time and (at least) the parameters listed in prm_nm

   Returns type: int ancillary data identification: 0 - ECMWF data,
   1 non-ECMWF data, -1 any trouble checking input anc files

     Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     char **           files            I      3 file name set, either MET
     or OZONE source
     char **           prm_nm           I      array of parameter names to read
     from the ECMWF
     int               n_prm            I      # parameters in prm_nm.  Note
                                               that the met parms (n_prm >1)
                                               have 7 fields to read but the
                                               t(2 m) and td(2 m) combine to
                                               make RH so the storage for met
                                               is 1 less (nprm_sto)
     int32_t           sto_ix           I      position in storage array to
                                               place result
     int32_t           anc_typ          O      type of anc data 0 - ECMWF,
                                               1 - non ECMWF, -1 - problem

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       24-Sep-2013     Original development

 *******************************************************************/
{
    int ids, iprm, npix, nlin, ilin, ipix;
    int t_days, t_hrs, lon_gt_180, ird;
    int dstpix, npix0, nlin0, status;
    int npix_ext, nlin_ext, ntim, f12_mtch, f23_mtch, anc_f_stat;
    int ncid, iprm_sto, nprm_sto;
    float s_lon, lon_step, e_lon, s_lat, lat_step, e_lat, time;
    float prm_bad[] = {ANCBAD, ANCBAD, ANCBAD, ANCBAD, ANCBAD, ANCBAD, ANCBAD, ANCBAD};
    float *base_data, *lat, *lon, *comp1, *comp2;
    double data_time;
    int64_t jd1900;

    nprm_sto = (n_prm > 1) ? n_prm - 1 : n_prm;
    /*
     *  identify the type of files set from their existance and name matching
     */
    anc_f_stat = anc_acq_f_stat(files, 0, 3);

    f12_mtch = 0;
    if ((anc_f_stat == ANC_STAT_2T_END) || (anc_f_stat == ANC_STAT_1T))
        f12_mtch = 1;
    f23_mtch = 0;
    if ((anc_f_stat == ANC_STAT_2T_START) || (anc_f_stat == ANC_STAT_1T))
        f23_mtch = 1;

    if (anc_f_stat == ANC_STAT_CLIM) {
        /*
         printf( "%s, %d I: Assuming standard (not ECMWF) climatology\n",
         __FILE__, __LINE__ );
         */
        return 1;
    }
    jd1900 = jd4713bc_get_jd(1900, 1, 1);
    for (ids = 0; ids < 3; ids++) {
        /*
         *  ECMWF is only anc data in netcdf format for now, so if it opens
         *  it should be ECMWF format
         */
        if ((ids > 0) && (anc_f_stat == ANC_STAT_1T))
            break;
        if ((ids == 1) && (f12_mtch == 1))
            continue;
        if ((ids == 2) && (f23_mtch == 1))
            break;

        if (Hishdf(files[ids]))
            status = NC2_ERR;
        else
            status = nc_open(files[ids], 0, &ncid);

        if (status != NC_NOERR) {
            if (ids == 0) {
                /*
                 printf(
                 "%s, %d E: file: %s ECMWF is not openable\n",
                 __FILE__, __LINE__, files[ids] );
                 */
                return -1;
            } else {
                printf("%s, %d: nc_open failed on file: %s\n", __FILE__, __LINE__, files[ids]);
                printf("       Mismatch or bad ECMWF file\n");
                return -1;
            }
        }

        /*
         *  get the basic information for MET1
         */
        if (((npix0 = ncio_dim_siz(ncid, "longitude")) == -1) ||
            ((nlin0 = ncio_dim_siz(ncid, "latitude")) == -1) || ((ntim = ncio_dim_siz(ncid, "time")) == -1)) {
            printf("%s, %d: ncio_dim_siz error reading longitude, latitude or time datasets\n", __FILE__,
                   __LINE__);
            return -1;
        }
        if (ids == 0) {
            npix = npix0;
            nlin = nlin0;
        } else {
            if ((npix != npix0) || (nlin != nlin0)) {
                printf("%s, %d: mismatch in size of MET array[%d]\n", __FILE__, __LINE__, ids);
                return -1;
            }
        }
        /*
         * for now, if more than 1 time, we can't proces it
         */
        if (ntim > 1) {
            printf("%s, %d: Number of times > 1, can't deal with at this time\n", __FILE__, __LINE__);
            return -1;
        }
        npix_ext = npix + 2;
        nlin_ext = nlin + 2;
        /*
         *  allocate storage in the structures for the data
         */
        for (iprm = 0; iprm < nprm_sto; iprm++) {
            if ((met_sto[iprm + sto_ix].data[ids] = (float *)malloc(npix_ext * nlin_ext * sizeof(float))) ==
                NULL) {
                printf("%s, %d: malloc failed for data[%d] in met_sto %d\n", __FILE__, __LINE__, ids, iprm);
                return -1;
            }
        }
        /*
         *  for 1st dataset, make array to read data into initially and
         *  arrays for lat, lon
         */
        if (ids == 0) {
            if (((base_data = (float *)malloc(npix * nlin * sizeof(float))) == NULL) ||
                ((comp1 = (float *)malloc(npix * nlin * sizeof(float))) == NULL) ||
                ((comp2 = (float *)malloc(npix * nlin * sizeof(float))) == NULL)) {
                printf("%s, %d: malloc failed for base_data or comp arrays\n", __FILE__, __LINE__);
                return -1;
            }

            if ((lat = (float *)malloc(nlin * sizeof(float))) == NULL) {
                printf("%s, %d: malloc failed for latitude\n", __FILE__, __LINE__);
                return -1;
            }
            if ((lon = (float *)malloc(npix * sizeof(float))) == NULL) {
                printf("%s, %d: malloc failed for longitude\n", __FILE__, __LINE__);
                return -1;
            }
            if (ncio_grab_f_ds(ncid, "latitude", lat) != 0) {
                printf("%s, %d: ncio_grab_f_ds failed on latitude\n", __FILE__, __LINE__);
                return -1;
            }
            if (ncio_grab_f_ds(ncid, "longitude", lon) != 0) {
                printf("%s, %d: ncio_grab_f_ds failed on longitude\n", __FILE__, __LINE__);
                return -1;
            }
            /*
             *  from the latitude, longitude arrays, determine the nav properties
             *  ECMWF longitudes go 0 -> 360 and we do -180 -> 180
             */
            s_lat = lat[0];
            e_lat = lat[nlin - 1];
            lat_step = lat[1] - lat[0];

            lon_gt_180 = -1;
            for (ipix = 0; ipix < npix; ipix++) {
                if (lon[ipix] > 180.) {
                    s_lon = lon[ipix] - 360.;
                    lon_gt_180 = ipix;
                    e_lon = lon[ipix - 1]; /* if lon_gt_180 = 0, need to upgrade this */
                    break;
                }
            }
            lon_step = lon[1] - lon[0];
        }
        /*
         *  get the time from the dataset and convert to julian days
         */
        if (ncio_grab_f_ds(ncid, "time", &time) != 0) {
            printf("%s, %d: error reading the time in\n", __FILE__, __LINE__);
            return -1;
        }

        t_days = time / 24;
        t_hrs = (int)time % 24;
        data_time = (double)(jd1900 + t_days) + (double)t_hrs / 24.;
        /*
         *  for all params, read the data, put lon in range -180 - 180
         *  and add extra layer to make interpolation easier
         *  The RH is made from T and Td in else below
         */
        ird = 0;
        for (iprm = 0; iprm < nprm_sto; iprm++) {
            iprm_sto = iprm + sto_ix;
            if (ird != 5) {
                if (ncio_grab_stdsclf_ds(ncid, prm_nm[ird], prm_bad[ird], base_data) != 0) {
                    printf("%s, %d: ncio_grab_stdsclf_ds failed on %s\n", __FILE__, __LINE__, prm_nm[ird]);
                    return -1;
                }
                ird++;
            } else {
                if ((ncio_grab_stdsclf_ds(ncid, prm_nm[ird], prm_bad[ird + sto_ix], comp1) != 0) ||
                    (ncio_grab_stdsclf_ds(ncid, prm_nm[ird + 1], prm_bad[ird + sto_ix + 1], comp2) != 0)) {
                    printf("%s, %d: ncio_grab_stdsclf_ds failed on %s or %s\n", __FILE__, __LINE__,
                           prm_nm[ird], prm_nm[ird + 1]);
                    return -1;
                }
                ird = ird + 2;
                /*
                 *  for the td and t at 2 m, we need to make a RH
                 */
                if (met_cvt_ttd_to_rh(npix * nlin, comp1, MET_UNITS__T_K, comp2, MET_UNITS__T_K, base_data) !=
                    0) {
                    printf("met_cvt_ttd_to_rh had an error\n");
                    printf("%s, %d: met_cvt_ttd_to_rh failure\n", __FILE__, __LINE__);
                    return -1;
                }
            }
            /*  rotate to -180 -> 180 */
            for (ilin = 0; ilin < nlin; ilin++) {
                for (ipix = 0; ipix < npix; ipix++) {
                    dstpix = ipix - lon_gt_180; /* put in with lon -180 -> 180 */
                    if (dstpix < 0)
                        dstpix += npix;
                    *(met_sto[iprm_sto].data[ids] + dstpix + 1 + (ilin + 1) * npix_ext) =
                        *(base_data + ipix + npix * ilin);
                }
            }
            /*  now, the extra boarder: lat, then lon */
            /* for lat, repeat the nearest value */
            for (ipix = 0; ipix < npix; ipix++) {
                *(met_sto[iprm_sto].data[ids] + ipix + 1) =
                    *(met_sto[iprm_sto].data[ids] + ipix + 1 + npix_ext);
                *(met_sto[iprm_sto].data[ids] + ipix + 1 + (nlin + 1) * npix_ext) =
                    *(met_sto[iprm_sto].data[ids] + ipix + 1 + nlin * npix_ext);
            }
            /* for lon, use the opposite side value */
            for (ilin = 0; ilin < nlin_ext; ilin++) {
                *(met_sto[iprm_sto].data[ids] + ilin * npix_ext) =
                    *(met_sto[iprm_sto].data[ids] + npix + ilin * npix_ext);
                *(met_sto[iprm_sto].data[ids] + npix + 1 + ilin * npix_ext) =
                    *(met_sto[iprm_sto].data[ids] + 1 + ilin * npix_ext);
            }
            /*  put in the controls found above  */
            met_sto[iprm_sto].s_lon = s_lon;
            met_sto[iprm_sto].lon_step = lon_step;
            met_sto[iprm_sto].nlon = npix;
            met_sto[iprm_sto].e_lon = e_lon;
            met_sto[iprm_sto].s_lat = s_lat;
            met_sto[iprm_sto].lat_step = lat_step;
            met_sto[iprm_sto].nlat = nlin;
            met_sto[iprm_sto].e_lat = e_lat;
            met_sto[iprm_sto].data_time[ids] = data_time;
            met_sto[iprm_sto].anc_f_stat = anc_f_stat;
        }
        /*
         *  close the dataset
         */
        nc_close(ncid);
    }
    return 0;
}

int anc_acq_lin(int32_t anc_class, l1str *l1rec)
/*******************************************************************

   anc_acq_lin

   purpose: get proper ancillary parameters for a particular line of
   points

   Returns type: int - 0 if good, else -1

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     int32_t           anc_class        I      anc data class to access the
                                               correct stored grids: 0 for
                                               MET grids and 1 for ozone grid
     l1str *           l1rec           I/O     structure with information
                                               for the line
   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       16 Aug 2013     original development

 *******************************************************************/
{
    double l_time;
    float data_val, uwnd, vwnd, data_val1, data_val2, dx, dy, lon_frac, lat_frac;
    float trg_lon, trg_lat, wt_t1, wt_t2, s_lon, s_lat;
    float *data;
    int iprm, xbox_st, ybox_st, nx, ny, t_interp, data1_ix, data2_ix, anc_f_stat;
    int npix, ipix, sto_st, sto_en;
    static float r2d = OEL_RADEG;
    /*
     *  find places in the met_sto structure to look in
     */
    if (anc_class == 0) {
        sto_st = 0;
        sto_en = 5;
    } else {
        sto_st = 6;
        sto_en = 6;
    }
    /*
     *  get the time of the current line
     */
    int16_t year, month, day;
    double sec;
    unix2ymds(l1rec->scantime, &year, &month, &day, &sec);
    l_time = (double)jd4713bc_get_jd((int32_t)year, (int32_t)month, (int32_t)day);
    l_time += sec / 86400.0;

    npix = l1rec->npix;
    /*
     *  for this line, decide which of the 3 anc files will be needed based on
     *  the line's time and the ancillary times
     */
    /*   ***** In this set up, all grids are the same, so only one
     determination is needed
     */
    anc_f_stat = met_sto[sto_st].anc_f_stat;
    if (anc_f_stat == ANC_STAT_1T) {
        /* use data[0] only */
        t_interp = 0;
        data1_ix = 0;
        /*  further along, when interpolating, use met_sto[0].data[data1_ix]
         to access the data */
    } else if (anc_f_stat == ANC_STAT_2T_START) {
        /* 2 different times in the 1, 2 positions */
        if (l_time < met_sto[sto_st].data_time[0]) {
            printf("%s, %d: data time is before the ancillary data start time\n", __FILE__, __LINE__);
            return -1;
        } else if (l_time > met_sto[sto_st].data_time[1]) {
            /* use MET2 only */
            t_interp = 0;
            data1_ix = 1;
        } else {
            /* in-between MET 1, 2 use data 0, 1 and time interpolate */
            t_interp = 1;
            data1_ix = 0;
            data2_ix = 1;
            wt_t1 = (met_sto[sto_st].data_time[1] - l_time) /
                    (met_sto[sto_st].data_time[1] - met_sto[sto_st].data_time[0]);
            wt_t2 = (l_time - met_sto[sto_st].data_time[0]) /
                    (met_sto[sto_st].data_time[1] - met_sto[sto_st].data_time[0]);
        }
    } else if (anc_f_stat == ANC_STAT_2T_END) {
        if (l_time < met_sto[sto_st].data_time[0]) {
            /* outside on the low end, use data[0] */
            t_interp = 0;
            data1_ix = 0;
        } else if (l_time > met_sto[sto_st].data_time[2]) {
            /* beyond the high end, Can't use end time alone */
            printf("%s, %d: data time is after the ancillary data end time\n", __FILE__, __LINE__);
            return -1;
        } else {
            /* between the MET 1 and 3 */
            t_interp = 1;
            data1_ix = 0;
            data2_ix = 2;
            wt_t1 = (met_sto[sto_st].data_time[2] - l_time) /
                    (met_sto[sto_st].data_time[2] - met_sto[sto_st].data_time[0]);
            wt_t2 = (l_time - met_sto[sto_st].data_time[0]) /
                    (met_sto[sto_st].data_time[2] - met_sto[sto_st].data_time[0]);
        }
    } else if (anc_f_stat == ANC_STAT_3T) {
        if (l_time < met_sto[sto_st].data_time[0]) {
            printf("%s, %d: data time is before the ancillary data start time\n", __FILE__, __LINE__);
            return -1;
        } else if (l_time > met_sto[sto_st].data_time[2]) {
            printf("%s, %d: data time is after the ancillary data end time\n", __FILE__, __LINE__);
            return -1;
        } else if (l_time < met_sto[sto_st].data_time[1]) {
            /* between data 0 and 1 */
            t_interp = 1;
            data1_ix = 0;
            data2_ix = 1;
            wt_t1 = (met_sto[sto_st].data_time[1] - l_time) /
                    (met_sto[sto_st].data_time[1] - met_sto[sto_st].data_time[0]);
            wt_t2 = (l_time - met_sto[sto_st].data_time[0]) /
                    (met_sto[sto_st].data_time[1] - met_sto[sto_st].data_time[0]);
        } else {
            /* what's left: between data 1 and 2 */
            t_interp = 1;
            data1_ix = 1;
            data2_ix = 2;
            wt_t1 = (met_sto[sto_st].data_time[2] - l_time) /
                    (met_sto[sto_st].data_time[2] - met_sto[sto_st].data_time[1]);
            wt_t2 = (l_time - met_sto[sto_st].data_time[1]) /
                    (met_sto[sto_st].data_time[2] - met_sto[sto_st].data_time[1]);
        }
    } else {
        /* this should not happen at this time - a status that is either
         climatology or not defined */
        printf("%s, %d: Undefined anc_f_stat - should not happen\n", __FILE__, __LINE__);
        return -1;
    }
    /*
     *  this found if time interpolation is needed and the and grids to use.
     *  next, for each pixel, find the bounding grid box and weights
     *  for bi-linear interpolation to be applied to each parameter
     *
     *  AGAIN note that since all parameters are from the same size grid,
     *  we can compute this information once.
     */
    dx = met_sto[sto_st].lon_step;
    dy = met_sto[sto_st].lat_step;
    s_lon = met_sto[sto_st].s_lon;
    s_lat = met_sto[sto_st].s_lat;
    nx = met_sto[sto_st].nlon;
    ny = met_sto[sto_st].nlat;
    for (ipix = 0; ipix < npix; ipix++) {
        trg_lat = l1rec->lat[ipix];
        trg_lon = l1rec->lon[ipix];

        /*
         xbox_st =
         MAX( MIN( (INT)( ( trg_lon - s_lon + dx / 2. ) / dx ), nx + 1 ), 0 );
         ybox_st =
         MAX( MIN( (INT)( ( trg_lat - s_lat + dy / 2. ) / dy ), ny + 1 ), 0 );
         x_dist = xbox_st * dx + s_lon - dx / 2;
         y_dist = ybox_st * dy + s_lat - dy / 2;

         I think below is correct for data at the grid points
         */
        xbox_st = MAX(MIN((int)((trg_lon - s_lon + dx) / dx), nx + 1), 0);
        ybox_st = MAX(MIN((int)((trg_lat - s_lat + dy) / dy), ny + 1), 0);

        lon_frac = (trg_lon - s_lon) / dx - (float)(xbox_st - 1);
        lat_frac = (trg_lat - s_lat) / dy - (float)(ybox_st - 1);

        for (iprm = sto_st; iprm < (sto_en + 1); iprm++) {
            data = met_sto[iprm].data[data1_ix];
            data_val1 = bilin_interp(data, xbox_st, (nx + 2), ybox_st, lon_frac, lat_frac);

            if (t_interp == 1) {
                data = met_sto[iprm].data[data2_ix];
                data_val2 = bilin_interp(data, xbox_st, (nx + 2), ybox_st, lon_frac, lat_frac);

                /*
                 *  do time interpolation
                 */
                if (data_val1 < ANCBAD + 1) {
                    if (data_val2 < ANCBAD + 1)
                        data_val = ANCBAD;
                    else
                        data_val = data_val2;
                } else
                    data_val = wt_t1 * data_val1 + wt_t2 * data_val2;
            } else
                data_val = data_val1;
            /*
             *  place this interpolated value in proper l1rec slot
             */
            switch (iprm) {
                case 0: /*  sfc press */
                    /*  Currently, no use for this, but it may be better than what
                     is done with MSL pressure to take it to a height above sea
                     level.  USE_PMSL of 0 will use this.  Note that pressure on Mt
                     Everest is nominally 337 mb, so enlarge range accordingly */
                    if (USE_PMSL == 1)
                        break;
                    else {
                        if (input->pressure != -2000) {
                            data_val = (data_val < ANCBAD + 1) ? ANCBAD : data_val / 100.;
                            if (data_val < 0)
                                l1rec->pr[ipix] = 1013.25;
                            else if (data_val < 250.)
                                l1rec->pr[ipix] = 250.;
                            else if (data_val > 1100.)
                                l1rec->pr[ipix] = 1100.;
                            else
                                l1rec->pr[ipix] = data_val;
                        }
                    }
                    break;
                case 1: /*  precip water  */
                    /* need to make from kg m^-2 into g cm^-2 */
                    if (input->watervapor != -2000)
                        l1rec->wv[ipix] = (data_val < ANCBAD + 1) ? 0. : data_val / 10.;
                    break;
                case 2: /*  sea level pressure  */
                    /* need to make from pascals to hectopascals (millibars) */
                    if (USE_PMSL == 0)
                        break;
                    else {
                        if (input->pressure != -2000) {
                            data_val = (data_val < ANCBAD + 1) ? ANCBAD : data_val / 100.;
                            if (data_val < 0)
                                l1rec->pr[ipix] = 1013.25;
                            else if (data_val < 900.)
                                l1rec->pr[ipix] = 900.;
                            else if (data_val > 1100.)
                                l1rec->pr[ipix] = 1100.;
                            else
                                l1rec->pr[ipix] = data_val;
                        }

                        /* if processing land, adjust pressure for terrain height */
                        if (proc_land && l1rec->height[ipix] != 0.0)
                            l1rec->pr[ipix] *= exp(-l1rec->height[ipix] / 8434);
                    }
                    break;
                case 3: /*  u wind, zonal W-E */
                    uwnd = (data_val < ANCBAD + 1) ? 0. : data_val;
                    l1rec->zw[ipix] = uwnd;
                    break;
                case 4: /*  v wind, meridional S-N */
                    vwnd = (data_val < ANCBAD + 1) ? 0. : data_val;
                    l1rec->mw[ipix] = vwnd;
                    if (input->windspeed != -2000)
                        l1rec->ws[ipix] = sqrt(pow(uwnd, 2.) + pow(vwnd, 2.));
                    if (input->windangle != -2000)
                        l1rec->wd[ipix] = atan2f(-uwnd, vwnd) * r2d;
                    break;
                case 5: /*  rel humidity % */
                    if (input->relhumid != -2000)
                        l1rec->rh[ipix] = (data_val < ANCBAD + 1) ? 0. : data_val;
                    break;
                case 6: /*  ozone  */
                    if (input->ozone != -2000) {
                        /*  convert from kg m^-2 to DU (units of 10 um STP OZ thickness)
                         to cm of thickness */
                        l1rec->oz[ipix] = (data_val < ANCBAD + 1) ? 0. : data_val * OZ_KG_M2_TO_DU / 1000.;
                    }
                    break;
            }
        }
    }
    return 0;
}

int anc_acq_lin_olci(int anc_class, char *file, l1str *l1rec)
/*******************************************************************

   anc_acq_lin_olci

   purpose: Read a line of ancillary data from the OLCI tie point meteo
   file and fill the met or ozone fields of teh l1 structure

   Returns type: status - 0 if all is good

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     int               anc_class        I      class of ancillary: 0 - met,
                                               1 - ozone
     char *            file             I      ancillary file name
     l1str *           l1rec           I/O     structure to fill with data

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson, SAIC 24 Feb 2016     original development

 *******************************************************************/
{
    static int ncid[2], varid[5]; /* there are 2 file ids for the
   met and ozone tie point files, and 4
   met dataset or variable ids:
   0 - wind - horizontal_wind
   1 - rh - humidity
   2 - msl pressure - sea_level_pressure
   3 - precip water - total_column_water_vapor
   4 - ozone - total_ozone
   */
    static float fill[5], valid_min[5], valid_max[5];
    static int32_t firstcalls[2] = {1, 1}, pix_smp[2];
    static int32_t npix, *qual, *qual_met, *qual_oz;
    static int32_t tie_epix, spix, epix, dpix;
    static size_t tie_nlin[2], tie_npix[2];
    static float *tie_data, *tie_met, *tie_oz;
    static float r2d = OEL_RADEG;
    int32_t anc_field_per_parm[5] = {1, 2, 1, 1, 1};
    int32_t anc_class_n_ds[2] = {4, 1}, class_off, n_ds_prm, ptr_prm;
    int32_t n_field_per_parm, ifld, ipx, px_tie1, px_tie2;
    int nid, lin_smp, ids, stat, ipx_dat;
    size_t start[3], count[3];
    char *anc_prm_nm[] = {"humidity", "horizontal_wind", "sea_level_pressure", "total_columnar_water_vapour",
                          "total_ozone"};
    float *ar, fr_dist, w1, w2;

    int32_t iscan = l1rec->iscan;
    /*
     *  do the initialization on the first calls for met and ozone files
     */
    if (anc_class == 0)
        class_off = 0;
    else
        class_off = 4;

    n_ds_prm = anc_class_n_ds[anc_class];

    if (firstcalls[anc_class]) {
        npix = l1rec->npix;
        spix = l1_input->spixl - 1;
        epix = l1_input->epixl - 1;
        dpix = l1_input->dpixl;
        /*
         * open the file and check that the sampling is 64 in pixls, 1 in lines
         */
        if (nc_open(file, NC_NOWRITE, (ncid + anc_class)) != NC_NOERR) {
            fprintf(stderr, "-E- %s %d: Unable to open OLCI tie point anc file: %s\n", __FILE__, __LINE__,
                    file);
            return -1;
        }
        if (nc_get_att_int(ncid[anc_class], NC_GLOBAL, "ac_subsampling_factor", &pix_smp[anc_class]) !=
            NC_NOERR) {
            fprintf(stderr,
                    "-E- %s %d: Unable to read column sampling attrib from OLCI tie point anc file: %s\n",
                    __FILE__, __LINE__, file);
            return -1;
        }
        if (nc_get_att_int(ncid[anc_class], NC_GLOBAL, "al_subsampling_factor", &lin_smp) != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s %d: Unable to read row sampling attrib from OLCI tie point anc file: %s\n",
                    __FILE__, __LINE__, file);
            return -1;
        }
        /*
         *  get the tie dataset dim sizes
         */
        if (nc_inq_dimid(ncid[anc_class], "tie_columns", &nid) != NC_NOERR) {
            fprintf(stderr, "-E- %s %d: Unable to get OLCI tie point meteo # pixels\n", __FILE__, __LINE__);
            return -1;
        }
        if (nc_inq_dimlen(ncid[anc_class], nid, &tie_npix[anc_class]) != NC_NOERR) {
            fprintf(stderr, "-E- %s %d: Unable to get OLCI tie point meteo # pixels\n", __FILE__, __LINE__);
            return -1;
        }

        if (nc_inq_dimid(ncid[anc_class], "tie_rows", &nid) != NC_NOERR) {
            fprintf(stderr, "-E- %s %d: Unable to get OLCI tie point meteo # lines\n", __FILE__, __LINE__);
            return -1;
        }
        if (nc_inq_dimlen(ncid[anc_class], nid, &tie_nlin[anc_class]) != NC_NOERR) {
            fprintf(stderr, "-E- %s %d: Unable to get OLCI tie point meteo # lines\n", __FILE__, __LINE__);
            return -1;
        }
        /*
         *  If not the expected values, leave
         */
        if (lin_smp != 1) {
            fprintf(stderr, "-E- %s %d: OLCI and tie point line sampling: %d, not = 1\n", __FILE__, __LINE__,
                    lin_smp);
            return -1;
        }
        tie_epix = pix_smp[anc_class] * (tie_npix[anc_class] - 1) + 1;
        if (tie_epix < epix) {
            fprintf(stderr, "-E- %s %d: tie point range out to pixel %d is < data range of %d\n", __FILE__,
                    __LINE__, tie_epix, epix);
            fprintf(stderr, "tie point sampling: %d and # pixels: %d\n", pix_smp[anc_class],
                    (int)tie_npix[anc_class]);
            return -1;
        }
        /*
         *  warn if the # lines in tie and dataset don't match
         *  UNFORTUNATELY, the # scans is not set in that data area at this time
         olci_dat = (olci_t *) ( l1rec->l1file->private_data );
         nlin = olci_dat->nscan;
         if( tie_nlin[anc_class] != nlin )
         fprintf(stderr,
         "-W- %s %d: OLCI tie point and radiance data have unequal # lines\n",
         __FILE__, __LINE__ );
         */
        /*
         * set up storage for the tie point data and quality line
         */
        if ((tie_data = (float *)malloc(tie_npix[anc_class] * sizeof(float))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate tie point data array\n", __FILE__, __LINE__);
            return -1;
        }
        /*  assign address for tie data storage */
        if (anc_class == 0)
            tie_met = tie_data;
        else
            tie_oz = tie_data;

        if ((qual = (int32_t *)malloc(tie_npix[anc_class] * sizeof(int32_t))) == NULL) {
            fprintf(stderr, "-E- %s %d: Unable to allocate tie point quality array\n", __FILE__, __LINE__);
            return -1;
        }
        if (anc_class == 0)
            qual_met = qual;
        else
            qual_oz = qual;
        /*
         *  loop through the parameters and set the access id, get fill value,
         *  valid min and max
         */
        for (ids = 0; ids < n_ds_prm; ids++) {
            ptr_prm = ids + class_off;
            if (nc_inq_varid(ncid[anc_class], anc_prm_nm[ptr_prm], &varid[ptr_prm]) != NC_NOERR) {
                fprintf(stderr, "-E- %s %d: Can't find id for OLCI anc tie point dataset: %s\n", __FILE__,
                        __LINE__, anc_prm_nm[ptr_prm]);
                return -1;
            }
            if (nc_get_att_float(ncid[anc_class], varid[ptr_prm], "_FillValue", &fill[ptr_prm]) != NC_NOERR) {
                fprintf(stderr, "-E- %s %d: Can't get _FillValue for OLCI anc tie point dataset: %s\n",
                        __FILE__, __LINE__, anc_prm_nm[ptr_prm]);
                return -1;
            }

            if (nc_get_att_float(ncid[anc_class], varid[ptr_prm], "valid_min", &valid_min[ptr_prm]) !=
                NC_NOERR) {
                fprintf(stderr, "-E- %s %d: Can't get valid_min for OLCI anc tie point dataset: %s\n",
                        __FILE__, __LINE__, anc_prm_nm[ptr_prm]);
                return -1;
            }
            if (nc_get_att_float(ncid[anc_class], varid[ptr_prm], "valid_max", &valid_max[ptr_prm]) !=
                NC_NOERR) {
                fprintf(stderr, "-E- %s %d: Can't get valid_max for OLCI anc tie point dataset: %s\n",
                        __FILE__, __LINE__, anc_prm_nm[ptr_prm]);
                return -1;
            }
        }
        firstcalls[anc_class] = 0;
    } /* end of init portion */
    /*
     *  read the line of each parameter
     */
    start[1] = 0;
    count[0] = 1;

    for (ids = 0; ids < n_ds_prm; ids++) {
        ptr_prm = ids + class_off;
        /* note that the description of the horizontal wind in the tie point
         file does not indicate the u, v components.  from all descriptions
         of wind components that ecmwf has (the source), they mention u,
         then v.  We can only assume that the 1st index is u (zonal)
         followed by v (meridional)   */
        switch (ptr_prm) {
            case 0:
                ar = l1rec->rh;
                break;
            case 1:
                ar = l1rec->zw;
                break;
            case 2:
                ar = l1rec->pr;
                break;
            case 3:
                ar = l1rec->wv;
                break;
            case 4:
                ar = l1rec->oz;
                break;
        }
        /*
         *  set the start and count to read the current scan
         */
        start[0] = iscan;
        count[1] = tie_npix[anc_class];
        count[2] = 0;
        start[2] = 0;

        tie_data = (anc_class == 0) ? tie_met : tie_oz;
        qual = (anc_class == 0) ? qual_met : qual_oz;

        n_field_per_parm = anc_field_per_parm[ptr_prm];
        for (ifld = 0; ifld < n_field_per_parm; ifld++) {
            if (ids == 1) {
                count[2] = 1;
                if (ifld == 1) {
                    start[2] = 1; /* for 2nd wind field */
                    ar = l1rec->mw;
                }
            }
            /*  make sure we have not exceeded the tie point # lines */
            if (iscan < tie_nlin[anc_class]) {
                /*  read the line  */
                if ((stat = nc_get_vara_float(ncid[anc_class], varid[ptr_prm], start, count, tie_data)) !=
                    NC_NOERR) {
                    fprintf(stderr, "-E- %s %d: Can't read OLCI anc tie point line, dataset: %s\n", __FILE__,
                            __LINE__, anc_prm_nm[ptr_prm]);
                    return -1;
                }
                /*  set a quality for all the tie points and convert to expected units */
                for (ipx = 0; ipx < tie_npix[anc_class]; ipx++) {
                    if ((*(tie_data + ipx) == fill[ptr_prm]) || (*(tie_data + ipx) < valid_min[ptr_prm]) ||
                        (*(tie_data + ipx) > valid_max[ptr_prm]))
                        *(qual + ipx) = 1;
                    else {
                        *(qual + ipx) = 0;
                        if (ptr_prm == 3)
                            *(tie_data + ipx) *= 0.1; /* for pw conversion
kg m-2 -> g cm-2 */
                        if (ptr_prm == 4)
                            *(tie_data + ipx) *= 46.698; /* for oz
conversion kg m-2 -> cm std atmos */
                    }
                }
                /* interpolate  */
                for (ipx = 0; ipx < npix; ipx++) {
                    /* need to find pixel in un-subsetted line */
                    ipx_dat = ipx * dpix + spix;
                    px_tie1 = ipx_dat / pix_smp[anc_class];
                    px_tie2 = px_tie1 + 1;
                    /*  for very last pixel, we need to make this adjustment*/
                    if (px_tie2 > (tie_npix[anc_class] - 1)) {
                        px_tie1 -= 1;
                        px_tie2 -= 1;
                    }
                    if ((*(qual + px_tie1) == 1) || (*(qual + px_tie2) == 1)) {
                        /*  fill with filler */
                        *(ar + ipx) = anc_miss_fill(ptr_prm);
                        l1rec->flags[ipx] |= ATMWARN;
                    } else {
                        /*  interpolate */
                        fr_dist = (float)(ipx_dat - px_tie1 * pix_smp[anc_class]) / (float)pix_smp[anc_class];
                        *(ar + ipx) = tie_data[px_tie1] * (1. - fr_dist) + tie_data[px_tie2] * fr_dist;
                    }
                }
            } else {
                /* place a fill value and set flag to atmwarn for whole line */
                for (ipx = 0; ipx < npix; ipx++) {
                    *(ar + ipx) = anc_miss_fill(ptr_prm);
                    l1rec->flags[ipx] |= ATMWARN;
                }
            }
        }
        /*
         *  set up wind speed, direction
         */
        if (ids == 1) {
            for (ipx = 0; ipx < npix; ipx++) {
                w1 = l1rec->zw[ipx];
                w2 = l1rec->mw[ipx];
                if (input->windspeed != -2000)
                    l1rec->ws[ipx] = sqrt(w1 * w1 + w2 * w2);
                if (input->windangle != -2000)
                    l1rec->wd[ipx] = atan2f(-w1, w2) * r2d;
            }
        }
        /*
         *  adjust the pressure over land
         */
        if (ids == 2) {
            for (ipx = 0; ipx < npix; ipx++) {
                if (l1rec->height[ipx] != 0.0)
                    l1rec->pr[ipx] *= exp(-l1rec->height[ipx] / 8434);
            }
        }
    }
    return 0;
}

float anc_miss_fill(int32_t prod_ix)
/*******************************************************************

   anc_miss_fill

   purpose: return a ancillary fill value for a product type
   using the same fillers as getanc.c/interpolate()

   Returns type: float of the fill value

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     int23_t           prod_ix          I      product number being done

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       16 Aug 2013     original development

 *******************************************************************/
{
    float fill;

    switch (prod_ix) {
        case 0:
            fill = 80.;
            break; /* humidity % */
        case 1:
            fill = 6.;
            break; /* wind m s-1 */
        case 2:
            fill = 1013.;
            break; /* msl pressure hPa */
        case 3:
            fill = 50.;
            break; /* pw in kg m-2  */
        case 4:
            fill = 0.36;
            break; /* ozone in cm at std atmos */
    }
    return fill;
}

float bilin_interp(float *data, int xbox_st, int nx, int ybox_st, float xfrac, float yfrac)
/*******************************************************************

   bilin_interp

   purpose: quick bi-linear interpolation.

   Returns type: float of interpolated result

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     float *           data             I      2-d grid of data
     int               xbox_st          I      x (longitude) index of
                                               grid box to interpolate
     int               nx               I      # pixels in x
     int               ybox_st          I      x (latitude) index of
                                               grid box to interpolate
     float             xfrac            I      fractional grid box distance
                                               from xbox_st to the point
     float             yfrac            I      fractional grid box distance
                                               from ybox_st to the point

   Modification history:
   Programmer        Date            Description of change
   ----------        ----            ---------------------
   W. Robinson       16 Aug 2013     original development

 *******************************************************************/
{
    float data_val;

    if ((*(data + xbox_st + nx * ybox_st) < (ANCBAD + 1)) ||
        (*(data + xbox_st + nx * (ybox_st + 1)) < (ANCBAD + 1)) ||
        (*(data + (xbox_st + 1) + nx * ybox_st) < (ANCBAD + 1)) ||
        (*(data + (xbox_st + 1) + nx * (ybox_st + 1)) < (ANCBAD + 1)))
        data_val = ANCBAD;
    else
        data_val = (1 - xfrac) * (1 - yfrac) * *(data + xbox_st + nx * ybox_st) +
                   (1 - xfrac) * yfrac * *(data + xbox_st + nx * (ybox_st + 1)) +
                   xfrac * (1 - yfrac) * *(data + (xbox_st + 1) + nx * ybox_st) +
                   xfrac * yfrac * *(data + (xbox_st + 1) + nx * (ybox_st + 1));
    return data_val;
}

int64_t jd4713bc_get_jd(int32_t year, int32_t month, int32_t day)
/*******************************************************************

   jd4713bc_get_jd

   purpose: get the julian day (from 4713 BC) for the year,
   month and day of month
   taken from the idl jd routine

   Returns type: int64_t - the julian date

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     int32_t           year             I      standard 4 digit year
     int32_t           month            I      month of the year
     int32_t           day              I      day of month

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       5 Aug 2013     original development

 *******************************************************************/
{
    int64_t lyear, lmonth, lday, jday;

    lyear = (int64_t)year;
    lmonth = (int64_t)month;
    lday = (int64_t)day;
    jday = (367 * lyear - 7 * (lyear + (lmonth + 9) / 12) / 4 + 275 * lmonth / 9 + lday + 1721014);
    /*
     *  this additional step is only needed if you expect to work on dates
     *  outside March 1, 1900 to February 28, 2100
     */
    jday = jday + 15 - 3 * ((lyear + (lmonth - 9) / 7) / 100 + 1) / 4;
    return jday;
}

int jd4713bc_get_date(int64_t jd, int32_t *year, int32_t *month, int32_t *day)
/*******************************************************************

   jd4713bc_get_date

   purpose: get the year, month, day from julian date
   (from 4713 BC)
   taken from the idl jddate routine

   Returns type: int - no set value now

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     int64_t           jd               I      julian date
     int32_t *         year             O      standard 4 digit year
     int32_t *         month            O      month
     int32_t *         day              O      day of month

   Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     W. Robinson       5 Aug 2013     original development

 *******************************************************************/
{
    int64_t v1, v2, v3, v4;
    v1 = jd + 68569;
    v2 = 4 * v1 / 146097;
    v1 = v1 - (146097 * v2 + 3) / 4;
    v3 = 4000 * (v1 + 1) / 1461001;
    v1 = v1 - 1461 * v3 / 4 + 31;
    v4 = 80 * v1 / 2447;
    *day = v1 - 2447 * v4 / 80;
    v1 = v4 / 11;
    *month = v4 + 2 - 12 * v1;
    *year = 100 * (v2 - 49) + v3 + v1;
    return 0;
}
