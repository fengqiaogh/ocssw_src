/***************************************************************** 
 * File:        ancnrt_2p5.c
 *
 * Purpose:     create HDF ancillary real time datafile from 
 *              NCEP (formerly NMC) (unpacked GRIB format) W, P, PW, H input
 *		or a TOVS ozone unpacked GRIB file from NMC.
 *              For temporary use, a 2.5 degree grid version for MET only
 *
 * Description: this program will create a HDF with the components:
 *
 *          Met data:
 *              HDF DAAC compliant Annotation object  
 *              1-D Latitude SDS
 *              1-D Longitude SDS
 *              2-D data parameter SDS (z_wind)
 *              2-D Q/C flag           (z_wind_QC)
 *              2-D data parameter SDS (m_wind)
 *              2-D Q/C flag           (m_wind_QC)
 *              2-D data parameter SDS (press)
 *              2-D Q/C flag           (press_QC)
 *              2-D data parameter SDS (rel_hum)
 *              2-D Q/C flag           (rel_hum_QC)
 *              2-D data parameter SDS (p_water)
 *              2-D Q/C flag           (p_water_QC)
 *
 *          Ozone data:
 *              2-D data parameter SDS (ozone)
 *              2-D Q/C flag           (ozone_QC)
 *
 * Input parms: 
 *    char *annotfile - input of DAAC compliant annotations ("fillenv.met")
 *    char *inz_wind  - input to process ("z_wind.dat")
 *    char *inm_wind  - input to process ("m_wind.dat")
 *    char *inpress   - input to process ("press.dat")
 *    char *inp_water - input to process ("p_water.dat")
 *    char *inrel_hum   - input to process ("rel_hum.dat")
 *
 *    z_wind, m_wind, press, p_water, rel_hum are replaced by the ozone 
 *        parameter only for the ozone dataset
 *
 *    char *outdir    - output directory to write file
 *
 * Output parms:none            
 *
 *       output name is derived from GRIB converted file header. 
 *       Example name: "N199312300_MET_NCEPR_6h.hdf"
 *                     "S199312300_TOVS.OZONE" for ozone file
 *
 * Returns:     run status    
 *
 * Local vars:  numerous variables for storing HDF annotations and SDSs.
 *
 * Subs called:         
 *              wrtattr    - write HDF annotation to file
 *              fill_annot - read HDF annotations from file
 *              wrtsds     - write SDS from array to file
 *              pexit      - print fatal errors
 *              check_usage- check command line args
 *              readgrib  - FORTRAN routine to read binary 
 *                              (GRIB extracted) files
 *
 * History:     none            
 *
 * Note:        This program is meant to generate SeaWiFS style HDF 
 *              ancillary real-time datafiles.
 *
 * Author:      Brian D. Schieber, GSC, 4/93
 *
 * Modification history:        
 *       BDS, 4/19/93 - HDF missing flag removed after discussions
 *                      with Jim F. and M.Darzi.
 *                    - NCEP data in Pascals, divided by 100 => millibars.
 *       BDS, 4/30/93 - Modified GRIB read WPHLONSZ to remove "world wrapped"
 *                      extra Lon. Will skip last longitude in HDF write. 
 *       BDS, 10/4/93 - Modified for new HDF design and to produce and 
 *                      output file which holds all three met fields (WPH).
 *       BDS, 11/8/93 - Modified to create z/m wind SDS rather than windspeed.
 *       BDS, 11/28/93- Modified MAX_EAST to 177.5 (from 180.0).
 *       BDS, 12/6/93 - Assuming no missing values, 0.0's passed through.
 *       BDS, 5/19/95 - Support new SeaWiFS v 2.7 product specs.
 *       BDS, 8/21/96 - Renamed 'perror' to 'pexit' to avoid conflict with
 *                      HDF4.0.
 *       BDS, 9/18/96 - Changed annot structure.  Modify to support new 1x1
 *                      GDAS1 files.
 *       BDS, 9/19/96 - Precipitable water has units kg/m^2
 *       W. Robinson, GSC, 4 Feb 97 put rh back in
 *       W. Robinson, GSC, 18 Jun 98  rename to ancnrt to include added 
 *                    role of processing the TOVS ozone data
 *       W. Robinson, GSC, 18 Dec 98 for the ozone output, write the
 *                    start and end times to be compatible with o3nrt
 *       W. Robinson, SAIC, 23 Jan 2008  to add ability to read the binary 
 *             files derived from the GRIB 2 format files (just a few 
 *             differences, see readgrib2 in readgrib.c)
 *       W. Robinson, SAIC, 18 Sep 2009  modify readgrib2 call args- rgmode
 *       W. Robinson, SAIC, 16 Feb 2010  modify readgrib2 call args npix, nlin
 *       W. Robinson, SAIC, 24 Jun 2014  updates for new met file name
 *       W. Robinson, SAIC, 24 Jun 2014  Formal 2.5 degree grid version
 *       W. Robinson, SAIC, 4 Aug 2014  use wind and height at standard 
 *              heights to make the 10 m wind for the R2 data
 *       W. Robinson, SAIC, 6 Apr 2015  as this is a Reanalysis 2 only version
 *             and never does GRIB2, remove calls for readgrib2
 *****************************************************************/
#include <mfhdf.h>
#include "ancil.h"
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include "ancnrt_proto.h"

/*
 * Met data specific settings
 */

#define VGROUPNAME   "Geophysical Data"

#define BIN_METH      2
#define REGISTRATION  CENTER
#define VSIZE         2.5
#define HSIZE         2.5
#define MAX_NORTH    90.0
#define MAX_SOUTH   -90.0
#define MAX_WEST   -180.0
#define MAX_EAST    177.50
#define WPHLATSZ    73                /* latitude  of GRIB files */
#define WPHLONSZ    144                /* long of MODIFIED GRIB files */
#define OZLATSZ   180     /* TOVS output # lines  */
#define OZLONSZ   288     /* TOVS output made to this size  */
#define ANC_MET  0   /* for ancillary data type (anctyp) say this is met data */
#define ANC_OZONE 1  /*  say this is ozone data */
#define GRIB_MODE_1 1  /*  the GRIB modes, GRIB 1 old, GRIB 2 post jan 08 */
#define GRIB_MODE_2 2
#define NLEV 6

int main(int argc, char *argv[]) {
    int i, j, k, l, anctyp;
    int rank;
    int result = 0;
    int array_size = 0;
    int numannarr = 0;
    int year = 0, month = 0, day = 0, hour = 0, msec = 0,
            julval = 0, pyear = 0;
    int32 shape[2];
    int npix = WPHLONSZ, nlin = WPHLATSZ;

    float32 datarr1[WPHLATSZ][WPHLONSZ];
    float32 datarr2[WPHLATSZ][WPHLONSZ];
    float32 datarr3[WPHLATSZ][WPHLONSZ];
    float32 datarr4[WPHLATSZ][WPHLONSZ];
    float32 datarr5[WPHLATSZ][WPHLONSZ];

    int32 ilev, nsamp, isamp, ihgt, ipix, ilin;
    float32 u_cube[ npix * nlin * NLEV ], v_cube[ npix * nlin * NLEV ];
    float32 hgt_cube[ npix * nlin * (NLEV + 1) ], u_lo, u_hi, v_lo, v_hi;
    float32 hgt_agl, hgt_agl2, corr, frac;
    char nms_uwind[NLEV][MAXNAMELNG], nms_vwind[NLEV][MAXNAMELNG];
    char nms_hgt[ NLEV + 1 ][MAXNAMELNG];
    char message[MAXNAMELNG], inz_wind[MAXNAMELNG], inm_wind[MAXNAMELNG];
    char inpress[MAXNAMELNG], inp_water[MAXNAMELNG], inrel_hum[MAXNAMELNG];
    char in_oz[MAXNAMELNG], longname[MAXNAMELNG], *sdsname;
    char *units;
    char *datafmt, outfile[MAXNAMELNG];
    char outfilename[MAXNAMELNG], outdir[MAXNAMELNG];
    char annotfile[MAXNAMELNG], source_name[100];
    int n_opt_arg;
    float latstep, lonstep; /* lat/lon incr. */
    int latsz, lonsz; /* lat/lon size */
    float nmostlat, smostlat, wmostlon, emostlon; /* lat/lon corners */
    time_t t; /* processing time */
    struct tm *localtm; /* processing time */
    struct annotation *xannot;
    int grib_mode;
    char grib2_t_str[300];

    /*
     * HDF datafile variables
     */

    int32 sdfid, fid, gridid;
    int32 datatype;

    /*
     * data type array pointers
     */

    float32 *float_SDSz_wind, *float_SDSm_wind, *float_SDSpress;
    float32 *float_SDSp_water, *float_SDSrel_hum;
    int16 *float_SDS_oz, *oz_newsiz;
    int8 *int8_SDSdataQC;

    /* 
     * external functions used 
     */

    int readgrib(), julian_(), wrtattr(), startHDF();
    int addAttr(), setSDSref();
    int32 wrtsds();
    void deattachHDFGrid(), pexit();
    int8 check_usage();
    struct annotation * fillenv_annot();
    float wind_hgt_corr(float);

    /*
     * ------- check command line arguments and set args  ------------------
     */

    if ((check_usage(argc, argv, &anctyp, &n_opt_arg, source_name,
            &grib_mode, grib2_t_str)) != 0)
        exit(-1);
    strcpy(annotfile, argv[ n_opt_arg ]);
    if (anctyp == ANC_MET) {
        for (ilev = 0; ilev < NLEV; ilev++) {
            strcpy(nms_uwind[ilev], argv[ 1 + ilev + n_opt_arg ]);
            strcpy(nms_vwind[ilev], argv[ 1 + ilev + NLEV + n_opt_arg ]);
            strcpy(nms_hgt[ilev], argv[ 1 + ilev + 2 * NLEV + n_opt_arg ]);
        }
        strcpy(nms_hgt[NLEV], argv[ 1 + 3 * NLEV + n_opt_arg ]);
        /*  end new wind input */
        strcpy(inpress, argv[ 2 + 3 * NLEV + n_opt_arg ]);
        strcpy(inp_water, argv[ 3 + 3 * NLEV + n_opt_arg ]);
        strcpy(inrel_hum, argv[ 4 + 3 * NLEV + n_opt_arg ]);
        if (argc - n_opt_arg < 8) strcpy(outdir, "./");
        else strcpy(outdir, argv[ 5 + 3 * NLEV + n_opt_arg ]);
    } else {
        strcpy(in_oz, argv[ 1 + n_opt_arg ]);
        if (argc - n_opt_arg < 3) strcpy(outdir, "./");
        else strcpy(outdir, argv[ 2 + n_opt_arg ]);
    }

    /*
     * ------- Read each binary file, extract data array and time  ----------
     */

    if (anctyp == ANC_MET) {
        if (grib_mode == GRIB_MODE_2) {
            strcpy(message, "GRIB 2 files not handled");
            pexit(message);
        } else {
            /*
             *  get the u, v, height and find the 10 m wind
             */
            for (ilev = 0; ilev < NLEV; ilev++) {
                result = readgrib(nms_uwind[ilev], npix, nlin,
                        (u_cube + ilev * npix * nlin), &year, &month, &day, &hour);
                if (result == 1) {
                    strcpy(message, "readgrib fail, file: ");
                    strcat(message, nms_uwind[ilev]);
                    pexit(message);
                }
                /*  */
                result = readgrib(nms_vwind[ilev], npix, nlin,
                        (v_cube + ilev * npix * nlin), &year, &month, &day, &hour);
                if (result == 1) {
                    strcpy(message, "readgrib fail, file: ");
                    strcat(message, nms_vwind[ilev]);
                    pexit(message);
                }
                /*  */
                result = readgrib(nms_hgt[ilev], npix, nlin,
                        (hgt_cube + ilev * npix * nlin), &year, &month, &day, &hour);
                if (result == 1) {
                    strcpy(message, "readgrib fail, file: ");
                    strcat(message, nms_hgt[ilev]);
                    pexit(message);
                }
            }
            result = readgrib(nms_hgt[NLEV], npix, nlin,
                    (hgt_cube + NLEV * npix * nlin), &year, &month, &day, &hour);
            if (result == 1) {
                strcpy(message, "readgrib fail, file: ");
                strcat(message, nms_hgt[NLEV]);
                pexit(message);
            }
        }
        /*
         *  Now, derive the wind at 10 m using the wind at the 1st height that is
         *  more than 10 m above the surface.  For a surface -10 -> 10 m above 
         *  surface, use that wind mixed in to the wind at the next higher level
         *  corrected to 10 m
         */
        nsamp = npix * nlin;
        for (isamp = 0; isamp < nsamp; isamp++) {
            for (ihgt = 0; ihgt < NLEV; ihgt++) {
                ilin = isamp / npix;
                ipix = isamp % npix;
                hgt_agl = *(hgt_cube + isamp + nsamp * ihgt) -
                        *(hgt_cube + isamp + nsamp * NLEV);
                if ((hgt_agl < -10.) && (ihgt == NLEV - 1)) {
                    /* if we are at the top and still < -10 m, just use the vals at top */
                    datarr1[ilin][ipix] = *(u_cube + isamp + nsamp * ihgt);
                    datarr2[ilin][ipix] = *(v_cube + isamp + nsamp * ihgt);
                    break;
                }
                if ((hgt_agl >= -10.) && (hgt_agl <= 10.)) {
                    if (ihgt == NLEV - 1) {
                        /* if at the top level, just use the top wind */
                        datarr1[ilin][ipix] = *(u_cube + isamp + nsamp * ihgt);
                        datarr2[ilin][ipix] = *(v_cube + isamp + nsamp * ihgt);
                        break;
                    } else {
                        /*
                         * if around the level, mix wind at that level with corrected 
                         * wind from next level up
                         */
                        u_lo = *(u_cube + isamp + nsamp * ihgt);
                        v_lo = *(v_cube + isamp + nsamp * ihgt);
                        frac = (10. - hgt_agl) / 20.;
                        hgt_agl2 = *(hgt_cube + isamp + nsamp * (ihgt + 1)) -
                                *(hgt_cube + isamp + nsamp * NLEV);
                        corr = wind_hgt_corr(hgt_agl2);
                        u_hi = corr * *(u_cube + isamp + nsamp * ihgt);
                        v_hi = corr * *(v_cube + isamp + nsamp * ihgt);
                        datarr1[ilin][ipix] = u_hi * frac + u_lo * (1. - frac);
                        datarr2[ilin][ipix] = v_hi * frac + v_lo * (1. - frac);
                        break;
                    }
                } else if (hgt_agl > 10.) {
                    /* just correct the wind above down to 10 m */
                    corr = wind_hgt_corr(hgt_agl);
                    datarr1[ilin][ipix] = *(u_cube + isamp + nsamp * ihgt) * corr;
                    datarr2[ilin][ipix] = *(v_cube + isamp + nsamp * ihgt) * corr;
                    break;
                }
            }
        }
        /*
         *  on with the pressure, PW, and RH
         */
        result = readgrib(inpress, npix, nlin, datarr3, &year, &month, &day,
                &hour);

        if (result == 1) {
            strcpy(message, "readgrib  file: ");
            strcat(message, inpress);
            pexit(message);
        }

        result = readgrib(inp_water, npix, nlin, datarr4, &year, &month,
                &day, &hour);

        if (result == 1) {
            strcpy(message, "readgrib  file: ");
            strcat(message, inp_water);
            pexit(message);
        }

        result = readgrib(inrel_hum, npix, nlin, datarr5, &year, &month,
                &day, &hour);

        if (result == 1) {
            strcpy(message, "readgrib  file: ");
            strcat(message, inrel_hum);
            pexit(message);
        }
    } else {
        result = readgrib(in_oz, npix, nlin, datarr1, &year, &month,
                &day, &hour);

        if (result == 1) {
            strcpy(message, "readgrib  file: ");
            strcat(message, in_oz);
            pexit(message);
        }
    }

    /* 
     *** This used to dump the GRIB array to stdout for a data check ***
      dlat = 2.5;
      dlon = 2.5;
      for (i = 0; i < WPHLATSZ; i++) {
         for (j = 0; j < WPHLONSZ; j++) {
            printf("\nLat: %f", 90.0 - (i * dlat));
            printf(" Lon: %f", (-180.0 + (j * dlon)));
            printf(" %f ", datarr1[i][j]/100.0);
         }
      }
     */

    /*
     * Create outfile name
     * Note that for grib 2 files, whole year is available already
     */
    if (grib_mode == GRIB_MODE_1) {
        if (year < 70) year = year + 2000; /* 2 digit 19xx or 20xx yrs to 4 */
        else year = year + 1900;
    }

    if (anctyp == ANC_MET) {
        julval = julian_(&day, &month, &year);
        sprintf(outfilename, "N%04d%03d%02d_MET_%s_6h.hdf",
                year, julval, hour, source_name);
        sprintf(outfile, "%s/%s",
                outdir, outfilename);
        /* For SDPS' benefit! */
        printf("%s\n", outfile);
    } else {
        julval = julian_(&day, &month, &year);
        sprintf(outfilename, "S%04d%03d00%03d23_%s.OZONE",
                year, julval, julval, source_name);
        sprintf(outfile, "%s/%s",
                outdir, outfilename);
        /* For SDPS' benefit! */
        printf("%s+%04d%03d000000000+%04d%03d235959999\n", outfile, year, julval,
                year, julval);
    }


    /*
     * Write HDF annotations, lat and lon SDSs
     */

    /**** check for existing HDF  ***/
    /*  WDR remove and let overwrite work
     if (!access(outfile, F_OK)) pexit ("....output file exists");
     */

    if ((numannarr = count_annot(annotfile)) == 0)
        pexit("no annotations found");

    if ((xannot = (struct annotation *)
            malloc(sizeof (struct annotation) * numannarr))
            == NULL) pexit("malloc Annotation");

    xannot = fillenv_annot(annotfile);

    /*
     * Determine processing time
     */

    (void) time(&t);
    localtm = localtime(&t);
    pyear = localtm->tm_year;
    if (pyear < 70) pyear = pyear + 2000; /* 2 digit 19xx or 20xx yrs to 4 */
    else pyear = pyear + 1900;

    if (hour > 0) msec = hour * 60 * 60 * 1000;
    else msec = 0;

    /*
     * -------  assign metadata values to local descriptive variables ----
     * -------  insert dates and other values to metadata array  ---------
     * -------  follow order of fillenv data file  -----------------------
     */

    for (i = 0; i < numannarr; i++) {

        /* Insert values to descr array element */

        if (!strcmp(xannot[i].label, "Product Name"))
            sprintf(xannot[i].descr, "%s", outfilename);

        else if (!strcmp(xannot[i].label, "Processing Time"))
            sprintf(xannot[i].descr, "%04d%03d%02d%02d%02d%03d",
                pyear, (localtm->tm_yday + 1), localtm->tm_hour, localtm->tm_min,
                localtm->tm_sec, 0);

        else if (!strcmp(xannot[i].label, "Input Files"))
            if (anctyp == ANC_MET)
                sprintf(xannot[i].descr, "%s %s %s %s %s",
                    inz_wind, inm_wind, inpress, inrel_hum, inp_water);
            else
                sprintf(xannot[i].descr, "%s", in_oz);

        else if (!strcmp(xannot[i].label, "Processing Control"))
            if (anctyp == ANC_MET)
                sprintf(xannot[i].descr, "%s %s %s %s %s %s %s %s",
                    argv[0], annotfile, inz_wind, inm_wind, inpress,
                    inrel_hum, inp_water, outdir);
            else
                sprintf(xannot[i].descr, "%s %s %s %s",
                    argv[0], annotfile, in_oz, outdir);

        else if (!strcmp(xannot[i].label, "Start Time"))
            if (anctyp == ANC_MET)
                sprintf(xannot[i].descr, "%04d%03d%02d0000000", year, julval, hour);
            else
                sprintf(xannot[i].descr, "%04d%03d000000000", year, julval);

        else if (!strcmp(xannot[i].label, "End Time"))
            if (anctyp == ANC_MET)
                sprintf(xannot[i].descr, "%04d%03d%02d0000000", year, julval, hour);
            else
                sprintf(xannot[i].descr, "%04d%03d230000000", year, julval);

        else if (!strcmp(xannot[i].label, "Start Year"))
            sprintf(xannot[i].descr, "%04d", year);

        else if (!strcmp(xannot[i].label, "Start Day"))
            sprintf(xannot[i].descr, "%03d", julval);

        else if (!strcmp(xannot[i].label, "Start Millisec"))
            if (anctyp == ANC_MET)
                sprintf(xannot[i].descr, "%08d", msec);
            else
                sprintf(xannot[i].descr, "00000000");

        else if (!strcmp(xannot[i].label, "End Year"))
            sprintf(xannot[i].descr, "%04d", year);

        else if (!strcmp(xannot[i].label, "End Day"))
            sprintf(xannot[i].descr, "%03d", julval);

        else if (!strcmp(xannot[i].label, "End Millisec"))
            if (anctyp == ANC_MET)
                sprintf(xannot[i].descr, "%08d", msec);
            else
                sprintf(xannot[i].descr, "82800000");

            /* WDR for ozone, modify the data type and data source*/
            /* This should come from the metadata file
                  else if ( ( anctyp == ANC_OZONE ) && 
                       ( !strcmp(xannot[i].label, "Data Type")) )
                     sprintf(xannot[i].descr, "Ozone" );

                  else if ( ( anctyp == ANC_OZONE ) && 
                       ( !strcmp(xannot[i].label, "Data Source")) )
                     sprintf(xannot[i].descr, "TOVS" );

                  else if ( ( anctyp == ANC_OZONE ) &&
                       ( !strcmp(xannot[i].label, "Data Source Desc")) )
                     sprintf(xannot[i].descr, "NOAA TIROS Operational Vertical Sounder" );
             */

        else if ((anctyp == ANC_OZONE) &&
                (!strcmp(xannot[i].label, "Temporal Resolution")))
            sprintf(xannot[i].descr, "12");

            /* Extract values from descr array element for program use */

        else if (!strcmp(xannot[i].label, "Northernmost Latitude"))
            sscanf(xannot[i].descr, "%f", &nmostlat);

        else if (!strcmp(xannot[i].label, "Southernmost Latitude"))
            sscanf(xannot[i].descr, "%f", &smostlat);

        else if (!strcmp(xannot[i].label, "Westernmost Longitude"))
            sscanf(xannot[i].descr, "%f", &wmostlon);

        else if (!strcmp(xannot[i].label, "Easternmost Longitude"))
            sscanf(xannot[i].descr, "%f", &emostlon);

        else if (!strcmp(xannot[i].label, "Latitude Step"))
            sscanf(xannot[i].descr, "%f", &latstep);

            /*
             *  WDR now a variable longitude step, # rows and columns
             */
        else if (!strcmp(xannot[i].label, "Longitude Step")) {
            if (anctyp == ANC_OZONE)
                lonstep = 1.25;
            else
                lonstep = 2.5;
            sprintf(xannot[i].descr, "%f", lonstep);
        }
        else if (!strcmp(xannot[i].label, "Number of Rows")) {
            if (anctyp == ANC_OZONE)
                latsz = 180;
            else
                latsz = nlin;
            sprintf(xannot[i].descr, "%d", latsz);
        }
        else if (!strcmp(xannot[i].label, "Number of Columns")) {
            if (anctyp == ANC_OZONE)
                lonsz = 288;
            else
                lonsz = npix;
            sprintf(xannot[i].descr, "%d", lonsz);
        }            /*
      *  SW point is different for the ozone
      */
        else if (!strcmp(xannot[i].label, "SW Point Latitude") &&
                anctyp == ANC_OZONE)
            sprintf(xannot[i].descr, "%f", -89.5);

        else if (!strcmp(xannot[i].label, "SW Point Longitude") &&
                anctyp == ANC_OZONE)
            sprintf(xannot[i].descr, "%f", -179.375);

    }

    /*
     * Create HDF file
     */

    if ((result = startHDF(outfile, &sdfid, &fid, DFACC_CREATE)) != 0)
        pexit("Fatal error starting HDF file");

    /*
     * Write attribute array to HDF file
     */

    if ((result = wrtattr(sdfid, xannot, numannarr)) != 0) pexit("wrtattr");
    free(xannot);

    /*
     * -------- Allocate space for 2D data arrays  -----------------------
     */

    rank = 2;

    if (anctyp == ANC_MET) {
        shape[0] = WPHLATSZ; /* lat */
        shape[1] = WPHLONSZ; /* lon */
        array_size = shape[0] * shape[1];

        if ((int8_SDSdataQC =
                (int8 *) calloc(array_size, sizeof (int8))) == NULL)
            pexit("calloc int8_SDSdataQC");

        if ((float_SDSz_wind =
                (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
            pexit("malloc float_SDSz_wind");

        if ((float_SDSm_wind =
                (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
            pexit("malloc float_SDSm_wind");

        if ((float_SDSpress =
                (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
            pexit("malloc float_SDSpress");

        if ((float_SDSp_water =
                (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
            pexit("malloc float_SDSp_water");

        if ((float_SDSrel_hum =
                (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
            pexit("malloc float_SDSrel_hum");

        /*
         * Process z/m_winds, press, rel_hum, p_water
         *  Shift all arrays so that 0 lon is centered, divide pressure by 100.
         */

        l = 0;
        for (j = 0; j < shape[0]; j++) { /* lats */
            for (k = (shape[1] / 2); k < shape[1]; k++) { /* lons */
                float_SDSz_wind[l] = datarr1[j][k];
                float_SDSm_wind[l] = datarr2[j][k];
                float_SDSpress[l] = datarr3[j][k] / 100.0;
                float_SDSp_water[l] = datarr4[j][k];
                float_SDSrel_hum[l] = datarr5[j][k];
                l++; /* index counter */
            } /* for k */

            for (k = 0; k < (shape[1] / 2); k++) { /* lons */
                float_SDSz_wind[l] = datarr1[j][k];
                float_SDSm_wind[l] = datarr2[j][k];
                float_SDSpress[l] = datarr3[j][k] / 100.0;
                float_SDSp_water[l] = datarr4[j][k];
                float_SDSrel_hum[l] = datarr5[j][k];
                l++; /* index counter */
            } /* for k */
        } /* for j */
    } else {
        shape[0] = OZLATSZ; /* lat */
        shape[1] = OZLONSZ; /* lon */
        array_size = shape[0] * shape[1];

        if ((int8_SDSdataQC =
                (int8 *) calloc(array_size, sizeof (int8))) == NULL)
            pexit("calloc int8_SDSdataQC");

        if ((float_SDS_oz =
                (int16 *) malloc(sizeof (int16) * WPHLATSZ * WPHLONSZ)) == NULL)
            pexit("malloc float_SDS_oz");

        if ((oz_newsiz =
                (int16 *) malloc(sizeof (int16) * array_size)) == NULL)
            pexit("malloc oz_newsiz");

        /*
         *  Process ozone the same way here 
         */
        l = 0;
        for (j = 0; j < WPHLATSZ; j++) { /* lats */
            for (k = (WPHLONSZ / 2); k < WPHLONSZ; k++) { /* lons */
                if (datarr1[j][k] < 0 || datarr1[j][k] > SHRT_MAX)
                    float_SDS_oz[l] = 0;
                else
                    float_SDS_oz[l] = (int16) datarr1[j][k];
                l++; /* index counter */
            } /* for k */

            for (k = 0; k < (WPHLONSZ / 2); k++) { /* lons */
                if (datarr1[j][k] < 0 || datarr1[j][k] > SHRT_MAX)
                    float_SDS_oz[l] = 0;
                else
                    float_SDS_oz[l] = (int16) datarr1[j][k];
                l++; /* index counter */
            } /* for k */
        } /* for j */

        /*
         *  interpolate the array to the traditional TOVS ozone size
         */
        resize_oz((short *) float_SDS_oz, WPHLONSZ, WPHLATSZ,
                OZLONSZ, OZLATSZ, (short *) oz_newsiz);
    }

    /* 
     * ----------------- Create Vgroup  ----------------------------
     */

    gridid = setupGrid(fid, VGROUPNAME);

    if (anctyp == ANC_MET) {
        /* 
         * ----------------- Write z_wind SDS ----------------------------
         */


        sdsname = "z_wind";
        strcpy(longname, "Zonal wind at 10 m");
        units = "m sec^-1";
        datafmt = "float32";
        datatype = DFNT_FLOAT32;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                float_SDSz_wind, gridid)) != 0)
            pexit("SDSinFile z_wind");

        free(float_SDSz_wind);

        /*
         * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
         */

        sdsname = "z_wind_QC";
        strcpy(longname, "Zonal wind at 10 m Q/C flag");
        units = "";
        datafmt = "int8";
        datatype = DFNT_INT8;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                int8_SDSdataQC, gridid)) != 0)
            pexit("SDSinFile z_wind_QC");


        /* 
         * ----------------- Write m_wind SDS ----------------------------
         */

        sdsname = "m_wind";
        strcpy(longname, "Meridional wind at 10 m");
        units = "m sec^-1";
        datafmt = "float32";
        datatype = DFNT_FLOAT32;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                float_SDSm_wind, gridid)) != 0)
            pexit("SDSinFile m_wind");

        free(float_SDSm_wind);

        /*
         * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
         */

        sdsname = "m_wind_QC";
        strcpy(longname, "Meridional wind at 10 m Q/C flag");
        units = "";
        datafmt = "int8";
        datatype = DFNT_INT8;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                int8_SDSdataQC, gridid)) != 0)
            pexit("SDSinFile m_wind_QC");

        /* 
         * ----------------- Write pressure SDS ----------------------------
         */


        sdsname = "press";
        strcpy(longname, "Atmospheric pressure at mean sea level");
        units = "millibars";
        datafmt = "float32";
        datatype = DFNT_FLOAT32;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                float_SDSpress, gridid)) != 0)
            pexit("SDSinFile press");

        free(float_SDSpress);

        /*
         * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
         */

        sdsname = "press_QC";
        strcpy(longname, "Atmospheric pressure at mean sea level Q/C flag");
        units = "";
        datafmt = "int8";
        datatype = DFNT_INT8;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                int8_SDSdataQC, gridid)) != 0)
            pexit("SDSinFile press_QC");

        /*
         * ----------------- Write humidity SDS ----------------------------
         */

        sdsname = "rel_hum";
        strcpy(longname, "Relative humidity at 1000 mb");
        units = "percent";
        datafmt = "float32";
        datatype = DFNT_FLOAT32;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                float_SDSrel_hum, gridid)) != 0)
            pexit("SDSinFile rel_hum");

        free(float_SDSrel_hum);

        /*
         * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
         */

        sdsname = "rel_hum_QC";
        strcpy(longname, "Relative humidity at 1000 mb Q/C flag");
        units = "";
        datafmt = "int8";
        datatype = DFNT_INT8;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                int8_SDSdataQC, gridid)) != 0)
            pexit("SDSinFile rel_hum_QC");

        /*
         * ----------------- Write precipitable water SDS ----------------------------
         */

        sdsname = "p_water";
        strcpy(longname, "Precipitable water");
        units = "kg m^-2";
        datafmt = "float32";
        datatype = DFNT_FLOAT32;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                float_SDSp_water, gridid)) != 0)
            pexit("SDSinFile p_water");

        free(float_SDSp_water);

        /*
         * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
         */

        sdsname = "p_water_QC";
        strcpy(longname, "Precipitable water Q/C flag");
        units = "";
        datafmt = "int8";
        datatype = DFNT_INT8;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                int8_SDSdataQC, gridid)) != 0)
            pexit("SDSinFile p_water_QC");
    } else {
        /*
         * ----------------- Write ozone SDS ----------------------------   */

        sdsname = "ozone";
        strcpy(longname, "Total ozone");
        units = "Dobson units";
        datafmt = "int16";
        datatype = DFNT_INT16;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                oz_newsiz, gridid)) != 0)
            pexit("SDSinFile ozone");

        free(float_SDS_oz);

        /*
         * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
         */

        sdsname = "ozone_QC";
        strcpy(longname, "Total ozone Q/C flag");
        units = "";
        datafmt = "int8";
        datatype = DFNT_INT8;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                int8_SDSdataQC, gridid)) != 0)
            pexit("SDSinFile ozone_QC");
    }

    /*
     * free storage for QC array (used for all products)
     */

    free(int8_SDSdataQC);

    /*
     * deattach HDF grid (used for all products)
     */

    deattachHDFgrid(gridid);

    /*
     * close HDF structures
     */

    if ((result = closeHDFstructs(sdfid, fid)) != 0) pexit("closeHDFstructs");
    return 0;
} /* main */

int8 check_usage(int argc, char *argv[], int *anctyp, int *n_opt_arg,
        char *source_name, int *grib_mode, char *grib2_t_str)
/*****************************************************************
   File:        check_usage
  
   Purpose:     check command line arguments
  
   Description: check command line arguements and report proper usage
                on argument count error or no argumnts.
  
   Returns:     success (0) or failure(-1)

      Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               argc             I      arg count
      char *[]          argv             I      input arg list
      int *             anctyp           O      type of anc data, either
                                                ANC_MET or ANC_OZONE
      int *             n_opt_arg        O      # option args
      char *            source_name      O      An alternate source (as in gbl 
                                              attrib source) name
      int *             grib_mode        O      Type of grib file the binary
                                              file came from: GRIB_MODE_1 for
                                              a GRIB 1 file, GRIB_MODE_2 for a
                                              GRIB 2 file (NCEP std as of
                                              Jan 2008)
      char *            grib2_t_str       O      string with time of the GRIB 2
                                                data - used for setting file
                                                name and attributes
  
   Note:        This routine is custom for each program that uses it since
                the arguments are different of course.
  
   Author:      Brian D. Schieber, GSC, 10/93
  
   Mods:           10/4/93 BDS, new HDF inputs handle all 3 NCEP parms
            W. Robinson, GSC, 4 Feb 97  changes for 5th param:rel_hum
            W. Robinson, GSC, 18 Jun 98  add flexability to either accept
                7 or 8 inputs for met data or 3 or 4 for ozone data
            W. Robinson, SAIC, 12 May 2005  have optional argument -s
                for an optional source name
            W. Robinson, SAIC, 24 Jan 2008  add a grib mode -2 for grib 2 
               format files
 *****************************************************************/ {
    extern char *optarg;
    extern int optind, optopt;
    int c, opt_source_found, errflg;

    /*
     *  get the options first, one for the source name
     */
    *grib_mode = GRIB_MODE_1;
    opt_source_found = 0;
    errflg = 0;
    *n_opt_arg = 0;
    while ((c = getopt(argc, argv, ":s:2:")) != -1) {
        switch (c) {
        case 's':
            strcpy(source_name, optarg);
            opt_source_found = 1;
            break;
        case '2':
            strcpy(grib2_t_str, optarg);
            *grib_mode = GRIB_MODE_2;
            break;
        case ':': /* -s without operand */
            fprintf(stderr,
                    "Option -%c requires an operand\n", optopt);
            errflg++;
            break;
        case '?':
            fprintf(stderr, "Unrecognized option: -%c\n", optopt);
            errflg++;
        }
    }
    if (!errflg) {
        *n_opt_arg = optind;
        /*
         *  for regular arguments
         *  for the R2, met will have 6 files for u,v wind and height, a height
         *  at sfc and then p, pw, rh and the making 22 files and the command, 
         *  metadata, making 24 
         *  
         */
        if (argc - optind == 24 || argc - optind == 25) {
            *anctyp = ANC_MET;
            if (!opt_source_found) strcpy(source_name, "NCEPR");
        } else if (argc - optind == 3 || argc - optind == 4) {
            *anctyp = ANC_OZONE;
            if (!opt_source_found) strcpy(source_name, "TOVS");
        } else {
            errflg++;
        }
    }
    if (errflg) {
        printf("\n\nUsage:\n");
        printf("\t%s <metadata> <fu 1-6> <fv 1-6> <fh 1-7> <fp>\n", argv[0]);
        printf("\t         <fpw> <frh> [outdir] [options]\n");
        printf("  OR (for ozone input):\n");
        printf("\t%s [-s<sname>] <metadata> <oz_file> [outdir]\n", argv[0]);
        printf("\nWhere:\n");
        printf("\t<metadata>:   DAAC-style HDF metadata (ASCII file)\n");
        printf("\t            usually fillenv_met.ancnrt\n");
        printf("\t<fu 1-6> is 6 u wind grid binary files\n");
        printf("\t<fv 1-6> is 6 v wind grid binary files\n");
        printf("\t<fh 1-7> is 7 height files, telling the height in meters\n");
        printf("\t         for the winds, going from 1000 mb up (to 500)\n");
        printf("\t         The 7th height is the surface height\n");
        printf("\t<fp>, <fpw> <frh> are for the pressure, precip water,\n");
        printf("\t         and relative humidity\n");
        printf("\tFor ozone data, only 1 input file (oz_file):\n");
        printf("\t (tovs_oz.dat)\n");
        printf("\toutdir:     Optional output directory ('.' default)\n");
        printf("\t   Options:\n");
        printf("\t-s <sname>   Optional source name to output file with:\n");
        printf("\t   for MET data: SYYYYDDD00_<sname>.MET, default sname = NCEPR\n");
        printf("\t   for OZONE data: SYYYYDDD00DDD23_<sname>.OZONE, default sname = TOVS\n");
        printf(
                "\t-2 <GRIB str>   Specify the -2 with the <GRIB str> to process\n");
        printf(
                "\t        binary files derived from a GRIB 2 format file.  The\n");
        printf(
                "\t        binary files are derived using wgrib2 with the -bin \n");
        printf(
                "\t        option.  The <GRIB str> is a required date / time\n");
        printf(
                "\t        identification string gotten using wgrib -t\n");
        printf(
                "\t         (see $SWFAPP/scripts/met_ingest2.scr for full example)\n");
        printf("\n");
        exit(-1);
    }
    return (0);
}

float wind_hgt_corr(float hgt)
/*****************************************************************

  wind_hgt_corr will return a correction factor to multiply by the wind at 
    height hgt to get the wind at 10 m height

  It uses the power law relation from 
   Hsu, S. A., Eric A. Meindl, and David B. Gilhousen, 1994: Determining
   the Power-Law Wind-Profile Exponent under Near-Neutral Stability Conditions
   at Sea, Applied Meteorology, Vol. 33, No. 6, June 1994.
  with the power reduced to 0.09

     Returns type: float - the multiplier to apply to the wind

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             hgt              I      height in meters above the 
                                                surface at which the wind is 
                                                measured.

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       5 Aug 2014      Original development

 *******************************************************************/ {
    if (hgt < 10.) return 1.;
    return (float) pow(10.0 / (double) hgt, 0.09);
}
