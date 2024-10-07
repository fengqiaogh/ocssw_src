/***************************************************************** 
 * File:        ancnrt.c
 *
 * Purpose:     create HDF ancillary real time datafile from 
 *              NCEP (formerly NMC) (unpacked GRIB format) W, P, PW, H input
 *              or the NCEP above PLUS O3
 *		or a TOVS ozone unpacked GRIB file from NMC.
 *
 * Description: this program will create a HDF with the components:
 *
 *          Met data (along with a '_QC' dataset):
 *              2-D data parameter SDS (z_wind)
 *              2-D data parameter SDS (m_wind)
 *              2-D data parameter SDS (press)
 *              2-D data parameter SDS (rel_hum)
 *              2-D data parameter SDS (p_water)
 *
 *          Ozone data (along with a '_QC' dataset):
 *              2-D data parameter SDS (ozone)
 *          Note that the NCEP O3 can also be put into a met file, if supplied
 *
 * Input parms: 
 *    char *annotfile - input of DAAC compliant annotations ("fillenv.met")
 *    a set of binary file names, 3 forms:
 *      ANC_MET - 5 files (products)
 *        wind, 10 m agl (z_wind.dat, m_wind.dat), 1000 mb pressure (press.dat),
 *        Column water vapor (p_water.dat), syc relative humidity (rel_hum.dat)
 *      ANC_OZONE - 1 files (products)
 *        Column ozone (ozone.dat)
 *      ANC_METOZ - 6 files (products)
 *        all products in ANC_MET and ANC_OZONE, from NCEP
 *    char *outdir    - output directory to write file
 *
 * Output file names
 *       Example name: "N199312300_MET_NCEP_6h.hdf"
 *                     "S199312300_TOVS.OZONE" for ozone file
 *                     "N201412303_MET_NCEP_f12_360x181.hdf"
 *                         (proposed for the forecast with different grid size)
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
 *       W. Robinson, SAIC, 18 Mar 2015  more grid size flexability, and 
 *           use of forecast data
 *****************************************************************/

#include <mfhdf.h>
#include "ancil.h"
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include "ancnrt_proto.h"
#include <timeutils.h>

#define VGROUPNAME   "Geophysical Data"
/*  default grid sizes for older datasets */
#define WPH_DFT_LATSZ    181                /* latitude  of GRIB files */
#define WPH_DFT_LONSZ    360                /* long of MODIFIED GRIB files */
#define OZLATSZ   180     /* fixed TOVS output # lines to cvt to */
#define OZLONSZ   288     /* fixed TOVS output made to this size  */
#define ANC_MET  0   /* for ancillary data type (anctyp) say this is met data */
#define ANC_OZONE 1  /*  say this is ozone data */
#define ANC_METOZ 2  /*  combined met and ozone from NCEP  */
#define ANC_METOZ2 3  /*  for test 2 extra parms of RH at 1000, RH at 10 m */
#define GRIB_MODE_1 1  /*  the GRIB modes, GRIB 1 old, GRIB 2 post jan 08 */
#define GRIB_MODE_2 2

int main(int argc, char *argv[]) {
    int i, j, k, l, h_fcst, anctyp, iprm, rank;
    int result = 0, numannarr = 0;
    int rgmode = 1, pyear = 0;
    int year = 0, month, day, doy = 0, hour = 0, msec = 0;
    int32 shape[2];
    int npix = WPH_DFT_LONSZ, nlin = WPH_DFT_LATSZ, npix_n, nlin_n, n_parm;
    double usec, vsec = 0.;
    float *float_out;
    int8 *qc_out;

    float32 *datarr; /* dynamically allocated storage for all met fields */
    char in_fnames[8][MAXNAMELNG];

    char message[MAXNAMELNG];
    char *longname, *sdsname;
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
    int ipoff, npix_out, nlin_out;
    /*
     *  info for SDSes of data and QC
     */
    /*
      char *metoz_nm[6] = { "z_wind", "m_wind", "press", "p_water", "rel_hum",
        "ozone" };
      char *metoz_ln_nm[6] = { "Zonal wind at 10 m", "Meridional wind at 10 m",
        "Atmospheric pressure at mean sea level", "Precipitable water", 
        "Relative humidity at 1000 mb", "Total ozone" };
      char *metoz_units[6] = { "m sec^-1", "m sec^-1", "millibars", "kg m^-2", 
        "percent", "Dobson units" };
      char *qc_nm[6] = { "z_wind_QC", "m_wind_QC", "press_QC", "p_water_QC",
        "rel_hum_QC", "ozone_QC" };
      char *qc_ln_nm[6] = { "Zonal wind at 10 m Q/C flag", 
        "Meridional wind at 10 m Q/C flag", 
        "Atmospheric pressure at mean sea level Q/C flag",
        "Precipitable water Q/C flag",
        "Relative humidity at 1000 mb Q/C flag",
        "Total ozone Q/C flag" };
      float fac_lst[6] = { 1., 1., 0.01, 1., 1., 1. };
     */
    char *metoz_nm[8] = {"z_wind", "m_wind", "press", "p_water", "rel_hum",
        "ozone", "rel_hum_2m_agl", "rel_hum_d30mb"};
    char *metoz_ln_nm[8] = {"Zonal wind at 10 m", "Meridional wind at 10 m",
        "Atmospheric pressure at mean sea level", "Precipitable water",
        "Relative humidity at 1000 mb", "Total ozone",
        "Relative humidity at 2 m AGL", "Relative humidity at 0-30 mb AGL"};
    char *metoz_units[8] = {"m sec^-1", "m sec^-1", "millibars", "kg m^-2",
        "percent", "Dobson units", "percent", "percent"};
    char *qc_nm[8] = {"z_wind_QC", "m_wind_QC", "press_QC", "p_water_QC",
        "rel_hum_QC", "ozone_QC", "rel_hum_2m_agl_QC", "rel_hum_d30mb_QC"};
    char *qc_ln_nm[8] = {"Zonal wind at 10 m Q/C flag",
        "Meridional wind at 10 m Q/C flag",
        "Atmospheric pressure at mean sea level Q/C flag",
        "Precipitable water Q/C flag",
        "Relative humidity at 1000 mb Q/C flag",
        "Total ozone Q/C flag",
        "Relative humidity at 2 m AGL Q/C flag",
        "Relative humidity at 0-30 mb AGL Q/C flag"};
    float fac_lst[8] = {1., 1., 0.01, 1., 1., 1., 1., 1.};
    /*
     * HDF datafile variables
     */
    int32 sdfid, fid, gridid;
    int32 datatype;
    /*
     * data type array pointers
     */
    int16 *float_SDS_oz, *oz_newsiz;
    /* 
     * external functions used 
     */
    int readgrib(), julian_(), wrtattr(), startHDF();
    int addAttr(), setSDSref(), grib2_t();
    int32 wrtsds();
    void deattachHDFGrid(), pexit();
    int8 check_usage();
    struct annotation * fillenv_annot();

    /*
     * ------- check command line arguments and set args  ------------------
     */

    if ((check_usage(argc, argv, &anctyp, &n_opt_arg, source_name,
            &grib_mode, grib2_t_str)) != 0)
        exit(-1);
    strcpy(annotfile, argv[ n_opt_arg ]);
    if (anctyp == ANC_METOZ) {
        strcpy(in_fnames[0], argv[ 1 + n_opt_arg ]);
        strcpy(in_fnames[1], argv[ 2 + n_opt_arg ]);
        strcpy(in_fnames[2], argv[ 3 + n_opt_arg ]);
        strcpy(in_fnames[3], argv[ 4 + n_opt_arg ]);
        strcpy(in_fnames[4], argv[ 5 + n_opt_arg ]);
        strcpy(in_fnames[5], argv[ 6 + n_opt_arg ]);
        strcpy(outdir, argv[ 7 + n_opt_arg ]);
    } else if (anctyp == ANC_METOZ2) {
        strcpy(in_fnames[0], argv[ 1 + n_opt_arg ]);
        strcpy(in_fnames[1], argv[ 2 + n_opt_arg ]);
        strcpy(in_fnames[2], argv[ 3 + n_opt_arg ]);
        strcpy(in_fnames[3], argv[ 4 + n_opt_arg ]);
        strcpy(in_fnames[4], argv[ 5 + n_opt_arg ]);
        strcpy(in_fnames[5], argv[ 6 + n_opt_arg ]);
        strcpy(in_fnames[6], argv[ 7 + n_opt_arg ]);
        strcpy(in_fnames[7], argv[ 8 + n_opt_arg ]);
        strcpy(outdir, argv[ 9 + n_opt_arg ]);
    } else if (anctyp == ANC_MET) {
        strcpy(in_fnames[0], argv[ 1 + n_opt_arg ]);
        strcpy(in_fnames[1], argv[ 2 + n_opt_arg ]);
        strcpy(in_fnames[2], argv[ 3 + n_opt_arg ]);
        strcpy(in_fnames[3], argv[ 4 + n_opt_arg ]);
        strcpy(in_fnames[4], argv[ 5 + n_opt_arg ]);
        strcpy(outdir, argv[ 6 + n_opt_arg ]);
    } else {
        strcpy(in_fnames[0], argv[ 1 + n_opt_arg ]);
        strcpy(outdir, argv[ 2 + n_opt_arg ]);
    }
    /*
     *  get grib2-specific information
     */
    if (grib_mode == GRIB_MODE_2) {
        if (grib2_t(grib2_t_str, &year, &doy, &hour, &npix_n, &nlin_n, &h_fcst)
                != 0) exit(-1);
        if ((npix_n > 0) && (nlin_n > 0)) {
            npix = npix_n;
            nlin = nlin_n;
        }
    }
    /*
     *  the ozone (TOAST) is at a fixed 1.25x1 degree grid so -
     *  and also reserve storage for all the data in the dararr array
     */
    n_parm = 1;
    if (anctyp == ANC_MET) n_parm = 5;
    if (anctyp == ANC_METOZ) n_parm = 6;
    if (anctyp == ANC_METOZ2) n_parm = 8;

    datarr = (float *) malloc(n_parm * npix * nlin * sizeof ( float));
    /*
     * ------- Read each binary file, extract data array and time  ----------
     */
    for (iprm = 0; iprm < n_parm; iprm++) {
        ipoff = iprm * npix * nlin;
        if (grib_mode == GRIB_MODE_2)
            result = readgrib2(in_fnames[iprm], npix, nlin, rgmode,
                (float *) (datarr + ipoff));
        else {
            result = readgrib(in_fnames[iprm], npix, nlin, (datarr + ipoff),
                    &year, &month, &day, &hour);
            /*  need to expand year out for the old grib year  */
            if (year < 70) year = year + 2000; /* 2 digit 19xx or 20xx yrs to 4 */
            else year = year + 1900;
            usec = ymds2unix(year, month, day, vsec);
            unix2yds(usec, (int16_t *) & year, (int16_t *) & doy, &vsec);
        }
        if (result == 1) {
            strcpy(message, "readgrib  file: ");
            strcat(message, in_fnames[iprm]);
            pexit(message);
        }
    }

    /* 
     *** This used to dump the GRIB array to stdout for a data check ***
      float dlat = 2.5;
      float dlon = 2.5;
      for (i = 0; i < WPH_DFT_LATSZ; i++) {
         for (j = 0; j < WPH_DFT_LONSZ; j++) {
            printf("\nLat: %f", 90.0 - (i * dlat));
            printf(" Lon: %f", (-180.0 + (j * dlon)));
            printf(" %f ", datarr1[i][j]/100.0);
         }
      }
     */

    /*
     * Create outfile name
     */
    /*  if( anctyp == ANC_METOZ )  */
    if ((anctyp == ANC_METOZ) || (anctyp == ANC_METOZ2)) {
        sprintf(outfilename, "N%04d%03d%02d_MET_%s_%4.4dx%4.4d_f%3.3d.hdf",
                year, doy, hour, source_name, npix, nlin, h_fcst);
        sprintf(outfile, "%s/%s",
                outdir, outfilename);
        /* For SDPS' benefit! */
        printf("%s\n", outfile);
    } else if (anctyp == ANC_MET) {
        sprintf(outfilename, "N%04d%03d%02d_MET_%s_6h.hdf",
                year, doy, hour, source_name);
        sprintf(outfile, "%s/%s",
                outdir, outfilename);
        /* For SDPS' benefit! */
        printf("%s\n", outfile);
    } else {
        sprintf(outfilename, "S%04d%03d00%03d23_%s.OZONE",
                year, doy, doy, source_name);
        sprintf(outfile, "%s/%s",
                outdir, outfilename);
        /* For SDPS' benefit! */
        printf("%s+%04d%03d000000000+%04d%03d235959999\n", outfile, year, doy,
                year, doy);
    }


    /*
     * Write HDF annotations
     */

    /**** check for existing HDF  ***/
    /*  WDR disable unless some complaints about behavior change 
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
    if (pyear < 90) pyear = pyear + 2000; /* 2 digit 19xx or 20xx yrs to 4 */
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

        else if (!strcmp(xannot[i].label, "Title")) {
            if ((anctyp == ANC_METOZ) || (anctyp == ANC_METOZ2))
                strcpy(xannot[i].descr, "NCEP Forecast Ancillary Data");
        } else if (!strcmp(xannot[i].label, "Data Type")) {
            if ((anctyp == ANC_METOZ) || (anctyp == ANC_METOZ2))
                strcpy(xannot[i].descr, "Meteorological and Ozone");
        }
        else if (!strcmp(xannot[i].label, "Processing Time"))
            sprintf(xannot[i].descr, "%04d%03d%02d%02d%02d%03d",
                pyear, (localtm->tm_yday + 1), localtm->tm_hour, localtm->tm_min,
                localtm->tm_sec, 0);

        else if (!strcmp(xannot[i].label, "Input Files"))
            if (anctyp == ANC_MET)
                sprintf(xannot[i].descr, "%s %s %s %s %s",
                    in_fnames[0], in_fnames[1], in_fnames[2], in_fnames[3],
                    in_fnames[4]);
            else if (anctyp == ANC_METOZ)
                sprintf(xannot[i].descr, "%s %s %s %s %s %s",
                    in_fnames[0], in_fnames[1], in_fnames[2], in_fnames[3],
                    in_fnames[4], in_fnames[5]);
            else if (anctyp == ANC_METOZ2)
                sprintf(xannot[i].descr, "%s %s %s %s %s %s %s %s",
                    in_fnames[0], in_fnames[1], in_fnames[2], in_fnames[3],
                    in_fnames[4], in_fnames[5], in_fnames[6], in_fnames[7]);
            else
                sprintf(xannot[i].descr, "%s", in_fnames[0]);

        else if (!strcmp(xannot[i].label, "Processing Control"))
            if (anctyp == ANC_MET)
                sprintf(xannot[i].descr, "%s %s %s %s %s %s %s %s",
                    argv[0], annotfile, in_fnames[0], in_fnames[1], in_fnames[2],
                    in_fnames[3], in_fnames[4], outdir);
            else if (anctyp == ANC_METOZ)
                sprintf(xannot[i].descr, "%s %s %s %s %s %s %s %s %s",
                    argv[0], annotfile, in_fnames[0], in_fnames[1], in_fnames[2],
                    in_fnames[3], in_fnames[4], in_fnames[5], outdir);
            else if (anctyp == ANC_METOZ)
                sprintf(xannot[i].descr, "%s %s %s %s %s %s %s %s %s %s %s",
                    argv[0], annotfile, in_fnames[0], in_fnames[1], in_fnames[2],
                    in_fnames[3], in_fnames[4], in_fnames[5], in_fnames[6],
                    in_fnames[7], outdir);
            else
                sprintf(xannot[i].descr, "%s %s %s %s",
                    argv[0], annotfile, in_fnames[0], outdir);

        else if (!strcmp(xannot[i].label, "Start Time"))
            if (anctyp != ANC_OZONE)
                sprintf(xannot[i].descr, "%04d%03d%02d0000000", year, doy, hour);
            else
                sprintf(xannot[i].descr, "%04d%03d000000000", year, doy);

        else if (!strcmp(xannot[i].label, "End Time"))
            if (anctyp != ANC_OZONE)
                sprintf(xannot[i].descr, "%04d%03d%02d0000000", year, doy, hour);
            else
                sprintf(xannot[i].descr, "%04d%03d230000000", year, doy);

        else if (!strcmp(xannot[i].label, "Start Year"))
            sprintf(xannot[i].descr, "%04d", year);

        else if (!strcmp(xannot[i].label, "Start Day"))
            sprintf(xannot[i].descr, "%03d", doy);

        else if (!strcmp(xannot[i].label, "Start Millisec"))
            if (anctyp != ANC_OZONE)
                sprintf(xannot[i].descr, "%08d", msec);
            else
                sprintf(xannot[i].descr, "00000000");

        else if (!strcmp(xannot[i].label, "End Year"))
            sprintf(xannot[i].descr, "%04d", year);

        else if (!strcmp(xannot[i].label, "End Day"))
            sprintf(xannot[i].descr, "%03d", doy);

        else if (!strcmp(xannot[i].label, "End Millisec"))
            if (anctyp != ANC_OZONE)
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
        else if ((anctyp == ANC_METOZ) &&
                (!strcmp(xannot[i].label, "Temporal Resolution")))
            sprintf(xannot[i].descr, "0");

            /* Extract values from descr array element for program use */

        else if (!strcmp(xannot[i].label, "Northernmost Latitude"))
            sscanf(xannot[i].descr, "%f", &nmostlat);

        else if (!strcmp(xannot[i].label, "Southernmost Latitude"))
            sscanf(xannot[i].descr, "%f", &smostlat);

        else if (!strcmp(xannot[i].label, "Westernmost Longitude"))
            sscanf(xannot[i].descr, "%f", &wmostlon);

        else if (!strcmp(xannot[i].label, "Easternmost Longitude"))
            sscanf(xannot[i].descr, "%f", &emostlon);

        else if ((!strcmp(xannot[i].label, "Latitude Step")) &&
                ((anctyp == ANC_METOZ) || (anctyp == ANC_METOZ2))) {
            latstep = 180. / (float) (nlin - 1);
            sprintf(xannot[i].descr, "%f", latstep);
        }            /*
      *  WDR now a variable longitude step, # rows and columns
      */
        else if (!strcmp(xannot[i].label, "Longitude Step")) {
            if (anctyp == ANC_OZONE)
                lonstep = 1.25;
            else if ((anctyp == ANC_METOZ) || (anctyp == ANC_METOZ2))
                lonstep = 360. / (float) npix;
            else
                lonstep = 1.;
            sprintf(xannot[i].descr, "%f", lonstep);
        }
        else if (!strcmp(xannot[i].label, "Number of Rows")) {
            if (anctyp == ANC_OZONE)
                latsz = 180;
            else if (anctyp == ANC_METOZ)
                latsz = nlin;
            else if (anctyp == ANC_METOZ2)
                latsz = nlin;
            else
                latsz = 181;
            sprintf(xannot[i].descr, "%d", latsz);
        }
        else if (!strcmp(xannot[i].label, "Number of Columns")) {
            if (anctyp == ANC_OZONE)
                lonsz = 288;
            else if (anctyp == ANC_METOZ)
                lonsz = npix;
            else if (anctyp == ANC_METOZ2)
                lonsz = npix;
            else
                lonsz = 360;
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

    /* 
     * ----------------- Create Vgroup  ----------------------------
     */
    gridid = setupGrid(fid, VGROUPNAME);
    npix_out = (anctyp == ANC_OZONE) ? OZLONSZ : npix;
    nlin_out = (anctyp == ANC_OZONE) ? OZLATSZ : nlin;

    if ((oz_newsiz =
            (int16 *) malloc(sizeof (int16) * npix_out * nlin_out)) == NULL)
        pexit("malloc oz_newsiz");

    if ((float_out = (float *) malloc(npix_out * nlin_out * sizeof (float)))
            == NULL) pexit("malloc float_out");
    if ((qc_out = (int8 *) calloc(npix * nlin, sizeof (int8)))
            == NULL) pexit("malloc qc_out");
    shape[0] = nlin_out;
    shape[1] = npix_out;
    if (anctyp != ANC_OZONE) {
        for (iprm = 0; iprm < n_parm; iprm++) {
            /*
             *  Shift data from 0-360 lon to -180 -> 180 in arrays
             */
            ipoff = iprm * npix_out * nlin_out;
            shift_180((datarr + ipoff), npix_out, nlin_out, *(fac_lst + iprm),
                    float_out);
            /*
             *  Write each data and qc SDS
             */
            sdsname = metoz_nm[iprm];
            longname = metoz_ln_nm[iprm];
            units = metoz_units[iprm];
            datafmt = "float32";
            datatype = DFNT_FLOAT32;

            if (SDSinFile(sdsname, longname, units, datafmt, datatype, sdfid,
                    rank, shape, float_out, gridid) != 0)
                pexit(strcat("SDSinFile ", metoz_nm[iprm]));

            sdsname = qc_nm[iprm];
            longname = qc_ln_nm[iprm];
            units = "";
            datafmt = "int8";
            datatype = DFNT_INT8;

            if (SDSinFile(sdsname, longname, units, datafmt, datatype, sdfid,
                    rank, shape, qc_out, gridid) != 0)
                pexit(strcat("SDSinFile ", qc_nm[iprm]));
        }
    } else {
        /*
         *  shift and write the ozone - it is a short in this usage
         */
        if ((float_SDS_oz = (int16 *) malloc(sizeof (int16) * npix * nlin))
                == NULL) pexit("malloc float_SDS_oz");

        l = 0;
        for (j = 0; j < nlin; j++) { /* lats */
            for (k = (npix / 2); k < npix; k++) { /* lons */
                if (*(datarr + k + j * npix) < 0 ||
                        *(datarr + k + j * npix) > SHRT_MAX)
                    float_SDS_oz[l] = 0;
                else
                    float_SDS_oz[l] = (int16)*(datarr + k + j * npix);
                l++; /* index counter */
            } /* for k */

            for (k = 0; k < (npix / 2); k++) { /* lons */
                if (*(datarr + k + j * npix) < 0 ||
                        *(datarr + k + j * npix) > SHRT_MAX)
                    float_SDS_oz[l] = 0;
                else
                    float_SDS_oz[l] = (int16)*(datarr + k + j * npix);
                l++; /* index counter */
            } /* for k */
        } /* for j */
        /*
         *  interpolate the array to the traditional TOVS ozone size
         */
        resize_oz((short *) float_SDS_oz, npix, nlin,
                npix_out, nlin_out, (short *) oz_newsiz);
        free(float_SDS_oz);

        sdsname = "ozone";
        longname = "Total ozone";
        units = "Dobson units";
        datafmt = "int16";
        datatype = DFNT_INT16;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                oz_newsiz, gridid)) != 0)
            pexit("SDSinFile ozone");

        free(oz_newsiz);
        /* ---- Write QC flag SDS, units attribute, and add to grid Vgroup -- */
        sdsname = "ozone_QC";
        longname = "Total ozone Q/C flag";
        units = "";
        datafmt = "int8";
        datatype = DFNT_INT8;

        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape,
                qc_out, gridid)) != 0)
            pexit("SDSinFile ozone_QC");

    }

    /*
     * free storage for QC array (used for all products)
     */
    free(qc_out);
    free(datarr);
    free(float_out);
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

void shift_180(float *in, int npix, int nlin, float fact, float *out)
/*****************************************************************
   shift_180  - shift a float array to have longitude range change 
       from 0 -> 360 to -180 -> 180

  Returns: none

  Parameters: (in calling order)
    Type              Name            I/O     Description
    ----              ----            ---     -----------
    float *           in               I      input array
    int               npix             I      array size in pixels
    int               nlin             I      array size in lines
    float             fact             I      conversion factor to apply
                                              to the data values
    float *           out              O      finished output array,
                                              already allocated storage

  W. Robinson, SAIC, 19 Mar 2015

 *****************************************************************/ {
    int iln, ipx, oloc;

    oloc = 0;
    for (iln = 0; iln < nlin; iln++) {
        for (ipx = npix / 2; ipx < npix; ipx++) {
            *(out + oloc) = *(in + ipx + npix * iln) * fact;
            oloc++;
        }

        for (ipx = 0; ipx < npix / 2; ipx++) {
            *(out + oloc) = *(in + ipx + npix * iln) * fact;
            oloc++;
        }
    }
}

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
            W. Robinson, SAIC, 20 Mar 2015  expand with ANC_METOZ type
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
         */
        if (argc - optind == 7) {
            *anctyp = ANC_MET;
            if (!opt_source_found) strcpy(source_name, "NCEP");
        } else if (argc - optind == 3) {
            *anctyp = ANC_OZONE;
            if (!opt_source_found) strcpy(source_name, "TOVS");
        } else if (argc - optind == 8) {
            *anctyp = ANC_METOZ;
            if (!opt_source_found) strcpy(source_name, "NCEP");
        } else if (argc - optind == 10) {
            *anctyp = ANC_METOZ2;
            if (!opt_source_found) strcpy(source_name, "NCEP");
        } else {
            errflg++;
        }
    }
    if (errflg) {
        printf("\n\nUsage:\n");
        printf("\t%s <metadata> <flist> <file> <outdir> [options]\n",
                argv[0]);
        printf("\nWhere:\n");
        printf("\tmetadata:   DAAC-style HDF metadata (ASCII file)\n");
        printf("\t            usually fillenv_met.ancnrt\n");
        printf("\tflist is a set of 1, 5, or 6 binary files created by wgrib\n");
        printf("\t    or wgrib2.  the 3 sets are:\n");
        printf("\t    1 file - for making an ozone dataset\n");
        printf("\t              file is just ozone (TOVS)\n");
        printf("\t    5 files - for making a MET dataset from the NCEP data,\n");
        printf("\t              files are z_wind, m_wind, press, p_water, and rel_hum\n");
        printf("\t    6 files - for making a MET and OZONE from the NCEP data,\n");
        printf("\t              files are MET set and ozone\n");
        printf("\t              This also has forecast file ingest ability\n");
        printf("\toutdir:     output directory ('./' can be used)\n");
        printf("\t   Options:\n");
        printf("\t-s <sname>   Optional source name to output file with:\n");
        printf("\t   for MET data: NYYYYDDD00_MET_<sname>_6h.hdf, default sname = NCEP\n");
        printf("\t   for MET and OZONE data:\n");
        printf("\t    NYYYYDDD00_MET_<sname>_<npix>x<nlin>_f<h_fcst>.hdf,\n");
        printf("\t    default sname = NCEP\n");
        printf("\t   for OZONE data: SYYYYDDD00DDD23_<sname>.OZONE, default sname = TOVS\n");
        printf("\t-2 <GRIB str>   Specify the -2 with the <GRIB str> to process\n");
        printf("\t        binary files derived from a GRIB 2 format file.  The\n");
        printf("\t        binary files are derived using wgrib2 with the -bin \n");
        printf("\t        option.  The <GRIB str> is a required date / time\n");
        printf("\t        identification string gotten using wgrib -t\n");
        printf("\t        (for 360 x 181 grid, 0 hour forecast)\n");
        printf("\t        or for MET and OZONE file: wgrib -t -nxny -ftime\n");
        printf("\t         (see $SWFAPP/scripts/met_ingest2.scr for full example)\n");
        printf("\n");
        exit(-1);
    }
    return (0);
}
