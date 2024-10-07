/*****************************************************************
 * File:        o3_toms.c
 *
 * Purpose:     create HDF ancillary real time datafile from
 *              OMI TOMS ascii input file
 *
 * Description: this program will create a HDF with the components:
 *
 *              HDF DAAC compliant global attributes
 *              Geometry Vdata for storage of LAT/LON endpoints and step
 *              2-D data parameter OZONE SDS
 *              2-D OZONE QC flag (padded all zeroes for TOMS and TOVS,
 *                with EARLIERVAL to EPTOMS missing points filled with
 *                prior day's data)
 *
 * Input parms:
 * char *infile    - input ozone file to process ("ga960903.ept")
 * char *directory - output file directory (".")
 * char *prevfile  - input ozone file for previous day ("OMITOMS" only)
 *
 * Returns:     false on error
 *
 * Modification history:
 *  W. Robinson, SAIC, Dec 2013  from o3nrt_omi
 *****************************************************************/
#include <mfhdf.h>
#include "ancil.h"
#include <time.h>
#include <unistd.h>
#include "ancnrt_proto.h"
#include "o3_toms.h"
#include "h5io.h"

/*
 * Ozone specific settings
 */

#define VGROUPNAME   "Geophysical Data"
#define EARLIERVAL    10  /* QC grid value for EPTOMS pts filled from prev day */
#define NEARESTVAL    11  /* QC grid value for EPTOMS pts filled from nearest
                             (file 1 or 2) non-zero point */

int main(int argc, char *argv[]) {
    char *datafmt;
    char *datalabel;
    char dataunit[MAXNAMELNG];
    char dataattr[MAXNAMELNG];
    char directory[MAXNAMELNG];
    char annotfile[MAXNAMELNG];
    char outfile[MAXNAMELNG];
    char outfilename[MAXNAMELNG];
    char infile[MAXNAMELNG];
    char prevfile[MAXNAMELNG];
    char start_time[17];
    char end_time[17];
    int i, j, k, l, kk;
    int rank;
    int result;
    int array_size = 0;
    int numannarr = 0;
    int modified_pts = 0; /* EPTOMS, track pts modded */
    int hour; /* data times */
    short startyear, startday; /* start data times */
    short endyear, endday; /* end data times */
    int startmsecday; /* millisec of day */;;
    int endmsecday; /* millisec of day */
    int pyear = 0; /* local proc times */
    int32 shape[2];
    time_t t; /* processing time */
    struct tm *localtm; /* processing time */
    struct annotation *xannot;
    toms_txt_info_struc toms_info, toms_info2;

    int32 sdfid, fid, gridid;
    int32 datatype;

    int8 *int8_SDSdataQC;
    int16 *int16_SDSdataO3;
    int16 *o3_prev, *o3;
    char *inst_name_frag[] = {"N7TOMS", "EPTOMS", "AURAOMI"};
    char *inst_data_src[] = {"Nimbus 7 / Total Ozone Mapping Spectrometer",
        "Earth Probe / Total Ozone Mapping Spectrometer",
        "Aura / Ozone Monitoring Instrument"};
    char *inst_platform[] = {"Nimbus 7", "Earth Probe", "Aura"};
    char gen_strs[24], *ocpath;
    int32_t bad_ct;

    /*
     * functions used
     */

    int wrtattr();
    int32 wrtsds();
    int startHDF();
    int setSDSref();
    int closeHDFstructs();
    void pexit();
    int8 check_usage();
    struct annotation * fillenv_annot();
    char *upcase();
    //void world_avg(int16 *, int32, int32, int, int, int16 *);
    int rd_oz_grid( char *, toms_txt_info_struc * );

    /*
     * ------- check command line arguments and set args  ------------------
     */

    if ((check_usage(argc, argv)) != 0) exit(-1);

    strcpy(annotfile, argv[1]);
    strcpy(infile, argv[1]);
    strcpy(directory, argv[2]);
    strcpy(prevfile, argv[3]);
    /*
     *  generate the annotation (fillenv...) file name internally
     */
    if ((ocpath = getenv("OCDATAROOT")) == NULL) {
        printf("%s, %d E: Cannot get OCDATAROOT definition\n", __FILE__,
                __LINE__);
        exit(-1);
    }
    strcpy(annotfile, ocpath);
    strcat(annotfile, "/common/fillenv.o3_toms");
    /*
     *  read the TOMS text file of data
     */
    
    if (rd_oz_grid( infile, &toms_info) != 0)
        pexit("rd_oz_grid");
    /*
     * ------- Create component date/time ----------
     *   TOVS input ASCII data provides start and end times in an additional
     *   header line.
     *   TOMS, EPTOMS and ADTOMS input ASCII does not provide start and end time,
     *   only one header line for year, Gregorian day, and Julian day.
     */

    hour = 00;
    sprintf(outfile, "%s/N%04d%03d%02d_O3_%s_24h.hdf",
            directory, toms_info.year, toms_info.doy, hour,
            inst_name_frag[ toms_info.toms_typ - 1 ]);

    sprintf(outfilename, "N%04d%03d%02d_O3_%s_24h.hdf",
            toms_info.year, toms_info.doy, hour,
            inst_name_frag[ toms_info.toms_typ - 1 ]);

    /* ending times based on 1st file in EPTOMS case */
    sprintf(end_time, "%04d%03d235959999", toms_info.year, toms_info.doy);

    endyear = toms_info.year;
    endday = toms_info.doy;
    endmsecday = 86399999;

    /* set the starting times as 0Z on day of primary file */
    sprintf(start_time, "%04d%03d000000000", toms_info.year, toms_info.doy);
    startyear = toms_info.year;
    startday = toms_info.doy;
    startmsecday = 0;

    /* Read previous day */

    if (rd_oz_grid(prevfile, &toms_info2) != 0)
        pexit("rd_oz_grid, #2");
    /*
     *  If the 2 files have different grid size, only use the primary day
     */
    if ((toms_info2.nlat != toms_info.nlat) || (toms_info2.nlon != toms_info.nlon)) {
        printf("%s, %d I: The 2 input files selected have different grid sizes\n",
                __FILE__, __LINE__);
        printf("1: %d X %d: %s\n", toms_info.nlon, toms_info.nlat, infile);
        printf("2: %d X %d: %s\n", toms_info2.nlon, toms_info2.nlat, prevfile);
        printf("    Using primary file only\n");
        toms_info2 = toms_info;
    }
    /*  WDR allow overwriting
     if (!access(outfile, F_OK))
       pexit("....output file exists.  Won't overwrite.");
     */

    /*
     * --- read annotation ASCII file for metadata placement in new HDF -----
     */

    if ((numannarr = count_annot(annotfile)) == 0)
        pexit("no annotations found");

    if ((xannot = (struct annotation *)
            malloc(sizeof (struct annotation) * numannarr))
            == NULL) pexit("malloc Annotation");

    xannot = fillenv_annot(annotfile);

    /*
     * Allocate space for 2D data arrays
     */

    rank = 2;
    shape[0] = toms_info.nlat; /* lat */
    shape[1] = toms_info.nlon; /* lon */
    array_size = shape[0] * shape[1];

    if ((int16_SDSdataO3 = (int16 *) malloc(sizeof (int16) * array_size))
            == NULL)
        pexit("malloc int16_SDSdataO3");

    if ((int8_SDSdataQC = (int8 *) malloc(sizeof (int8) * array_size)) == NULL)
        pexit("malloc int8_SDSdataQC");
    for (i = 0; i < array_size; i++) int8_SDSdataQC[i] = 0; /* init */

    /*
     * Create SDS array for writing
     */

    o3 = toms_info.datarr;
    o3_prev = toms_info2.datarr;
    l = 0;
    bad_ct = 0;
    for (j = 0; j < shape[0]; j++) /* lats */ {
        for (k = 0; k < shape[1]; k++) /* lons */ {
            /*
             *  Overwrite missings (0's) with matching
             *  position from previous day's data.  Mark QC array (as EARLIERVAL)
             */
            if (*(o3 + k + shape[1] * j) == 0) {
                if (*(o3_prev + k + shape[1] * j) != 0) {
                    *(o3 + k + shape[1] * j) = *(o3_prev + k + shape[1] * j);
                    int8_SDSdataQC[l] = EARLIERVAL;
                } else /* even 2nd file is missing for this point */ {
                    /* Scan lons in either longitudinal direction from current point */
                    for (kk = 1; kk <= 50; kk++) {
                        if ((k - kk >= 0) && (*(o3 + (k - kk) + shape[1] * j) > 0)) {
                            *(o3 + k + shape[1] * j) = *(o3 + (k - kk) + shape[1] * j);
                            int8_SDSdataQC[l] = NEARESTVAL;
                            break;
                        } else if ((k + kk < shape[1])
                                && (*(o3 + (k + kk) + shape[1] * j) > 0)) {
                            *(o3 + k + shape[1] * j) = *(o3 + (k + kk) + shape[1] * j);
                            int8_SDSdataQC[l] = NEARESTVAL;
                            break;
                        }
                    }
                }
                modified_pts++;
            }

            if (*(o3 + k + shape[1] * j) > 0) {
                int16_SDSdataO3[l] = *(o3 + k + shape[1] * j);
            } else {
                int16_SDSdataO3[l] = 0;
                int8_SDSdataQC[l] = 20; /* missing data bad ancillary */
                bad_ct++;
            }
            l++; /* index counter */
        } /* for k */
    } /* for j */
    /*
     *  if the file is all bad data, don't make it
     */
    if (shape[0] * shape[1] == bad_ct)
        pexit("Final ozone array is all bad");

    /*
     *  make the fill seamlessly join to the day's data
     */
    if (fill_smooth(int16_SDSdataO3, (char *) int8_SDSdataQC, shape[1], shape[0])
            != 0)
        pexit("fill_smooth problem");
    /*
     *  get both generation dates to gen_strs
     */
    strcpy(gen_strs, toms_info.gen_str);
    strcat(gen_strs, ",");
    strcat(gen_strs, toms_info2.gen_str);
    /*
     * -------  Assign metadata values to local descriptive variables ----
     * -------   insert dates and other values to metadata array  ---------
     */

    for (i = 0; i < numannarr; i++) {

        if (!strcmp(xannot[i].label, "Product Name"))
            sprintf(xannot[i].descr, "%s", outfilename);

        else if (!strcmp(xannot[i].label, "Processing Time")) {

            /* Determine processing time */
            (void) time(&t);
            localtm = localtime(&t);
            pyear = localtm->tm_year;
            if (pyear < 90) pyear = pyear + 2000; /* 2 digit to 4 */
            else pyear = pyear + 1900;

            sprintf(xannot[i].descr, "%04d%03d%02d%02d%02d%03d",
                    pyear, localtm->tm_yday + 1, localtm->tm_hour, localtm->tm_min,
                    localtm->tm_sec, 0);
        }
        else if (!strcmp(xannot[i].label, "Input Files"))
            sprintf(xannot[i].descr, "%s,%s", infile, prevfile);

        else if (!strcmp(xannot[i].label, "Processing Control"))
            sprintf(xannot[i].descr, "%s %s %s %s %s",
                argv[0], annotfile, infile, directory, prevfile);

        else if (!strcmp(xannot[i].label, "Start Time"))
            sprintf(xannot[i].descr, "%s", start_time);

        else if (!strcmp(xannot[i].label, "End Time"))
            sprintf(xannot[i].descr, "%s", end_time);

        else if (!strcmp(xannot[i].label, "Start Year"))
            sprintf(xannot[i].descr, "%04d", startyear);

        else if (!strcmp(xannot[i].label, "Start Day"))
            sprintf(xannot[i].descr, "%03d", startday);

        else if (!strcmp(xannot[i].label, "Start Millisec"))
            sprintf(xannot[i].descr, "%08d", startmsecday);

        else if (!strcmp(xannot[i].label, "End Year"))
            sprintf(xannot[i].descr, "%04d", endyear);

        else if (!strcmp(xannot[i].label, "End Day"))
            sprintf(xannot[i].descr, "%03d", endday);

        else if (!strcmp(xannot[i].label, "End Millisec"))
            sprintf(xannot[i].descr, "%08d", endmsecday);

        else if (!strcmp(xannot[i].label, "Node Crossing Time"))
            sprintf(xannot[i].descr, "%s", toms_info.node_time);

        else if (!strcmp(xannot[i].label, "Points Modified"))
            sprintf(xannot[i].descr, "%08d", modified_pts);

        else if (!strcmp(xannot[i].label, "Data Source"))
            sprintf(xannot[i].descr, "%s",
                inst_name_frag[ toms_info.toms_typ - 1 ]);

        else if (!strcmp(xannot[i].label, "Data Source Desc"))
            sprintf(xannot[i].descr, "%s",
                inst_data_src[ toms_info.toms_typ - 1 ]);

        else if (!strcmp(xannot[i].label, "Satellite Platform"))
            sprintf(xannot[i].descr, "%s",
                inst_platform[ toms_info.toms_typ - 1 ]);

        else if (!strcmp(xannot[i].label, "Latitude Step"))
            sprintf(xannot[i].descr, "%5.3f", toms_info.del_lat);

        else if (!strcmp(xannot[i].label, "Longitude Step"))
            sprintf(xannot[i].descr, "%5.3f", toms_info.del_lon);

        else if (!strcmp(xannot[i].label, "SW Point Latitude"))
            sprintf(xannot[i].descr, "%8.3f", toms_info.slat);

        else if (!strcmp(xannot[i].label, "SW Point Longitude"))
            sprintf(xannot[i].descr, "%8.3f", toms_info.slon);

        else if (!strcmp(xannot[i].label, "Number of Rows"))
            sprintf(xannot[i].descr, "%08d", toms_info.nlat);

        else if (!strcmp(xannot[i].label, "Number of Columns"))
            sprintf(xannot[i].descr, "%08d", toms_info.nlon);

        else if (!strcmp(xannot[i].label, "Raw Data Generation Date"))
            sprintf(xannot[i].descr, "%s", gen_strs);

    } /* for i */

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
     * --- Write OZONE SDS, units attribute, and add to grid Vgroup ----------
     */

    datalabel = "ozone";
    strcpy(dataunit, "Total ozone");
    strcpy(dataattr, "Dobson units");
    datafmt = "int16";
    datatype = DFNT_INT16;

    gridid = setupGrid(fid, VGROUPNAME);

    if ((SDSinFile(datalabel, dataunit, dataattr, datafmt,
            datatype, sdfid, rank, shape, int16_SDSdataO3, gridid)) != 0)
        pexit("SDSinFile ozone wrtsds");

    free(int16_SDSdataO3);

    /*
     * --- Write OZONE QC flag SDS, units attribute, and add to grid Vgroup ---
     */

    datalabel = "ozone_QC";
    strcpy(dataunit, "Total ozone Q/C flag");
    strcpy(dataattr, "");
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(datalabel, dataunit, dataattr, datafmt,
            datatype, sdfid, rank, shape, int8_SDSdataQC, gridid)) != 0)
        pexit("SDSinFile ozone QC wrtsds");

    free(int8_SDSdataQC);

    deattachHDFgrid(gridid); /* deattach HDF grid */

    /*
     * close HDF structures
     */

    result = 0;
    if ((result = closeHDFstructs(sdfid, fid)) != 0) pexit("closeHDFstructs");


    /*
     * Print info to STDOUT for SeaWiFS Project
     */

    printf("%s+%s+%s\n", outfile, start_time, end_time);

    exit(0);

} /* main */

int rd_oz_grid( char * infile, toms_txt_info_struc *toms_info )
/****************************************************************************
  rd_oz_grid - read the grid of ozone from a GES DISC HDF5 file of AURA
  OMI ozone

  Returns int, status 0 is good

  Parameters: (in calling order)
  Type              Name            I/O     Description
  ----              ----            ---     -----------
  char *            infile            I      file to read
  toms_txt_info_struc *  toms_info   O      All info about the grid

  W. Robinson, 3 Sep 2021
****************************************************************************/
  {
  static int npix = 360, nlin = 180;
  h5io_str fid, grp_id;
  int ret = 1, stat, iloc, oloc;
  int32_t year, doy;
  float *tmp_flt, fval;
  char dat_grp[] =
    "HDFEOS/GRIDS/OMI Column Amount O3/Data Fields/ColumnAmountO3";
  char attr_path[] = "HDFEOS/ADDITIONAL/FILE_ATTRIBUTES";
 /*
  *  at first, assume the file is hdf 5 and try reading it that way
  *  open file, get the ozone array
  */
  if( ( stat = h5io_openr( infile, 0, &fid ) ) != 0 ) {
    printf("%s, %d: Open error on file: %s\n", __FILE__, __LINE__, infile );
    }
  if( stat == 0 )
    {
   /*
    *  the file is h5.  assume it is ONI TOMS
    */
    toms_info->toms_typ = TOMS_TYP_OMITOMS; 
    // make a temporary float for the read data
    if( ( tmp_flt = (float *) malloc( 180 * 360 * sizeof(float) ) ) == NULL )
      {
      printf("%s, %d: Unable to malloc the oz storage space\n",
        __FILE__, __LINE__ );
      return ret;
      }
    if( ( stat = h5io_grab_ds( &fid, dat_grp, (float *)tmp_flt ) ) != 0 ) {
      printf("%s, %d: read error on file: %s\n", __FILE__, __LINE__, infile );
      return ret;
      }
   /*
    *  Set up the final int16 storage for the ozone data
    */
    if( ( toms_info->datarr = (int16 *) malloc( npix * nlin * sizeof(int16) ) )
      == NULL )
      {
      printf("%s, %d: Unable to allocate the storage for the ozone data\n", 
        __FILE__, __LINE__ );
      return ret;
      }
    for( iloc = 0; iloc < npix * nlin; iloc++ )
      {
     /* to match the old AURA O3 flip in latitude */
      oloc = ( iloc % npix ) + npix * ( nlin - 1 - iloc / npix );
      fval = *( tmp_flt + iloc );
      *( toms_info->datarr + oloc ) = ( ( fval < 50. ) || ( fval > 700. ) )
        ? 0 : (int16) ( fval + 0.5 );
      }
    free( tmp_flt );
   /*
    *  set group to HDFEOS/ADDITIONAL/FILE_ATTRIBUTES and get the GranuleYear,
    *  GranuleDayOfYear
    */
    if( ( stat = h5io_set_grp( &fid, attr_path, &grp_id ) ) != 0 ) {
      printf("%s, %d: path could not be set on file: %s\n", __FILE__,
        __LINE__, infile );
      return ret;
      }
  
    if( ( stat = h5io_rd_attr( &grp_id, "GranuleYear", 
      &year ) )!= 0 ) {
      printf("%s, %d: unable to read attr GranuleYear on file: %s\n", __FILE__,
        __LINE__, infile );
      return ret;
      }
    if( ( stat = h5io_rd_attr( &grp_id, "GranuleDayOfYear", 
      &doy ) )!= 0 ) {
      printf("%s, %d: unable to read attr GranuleDayOfYear on file: %s\n", __FILE__,
        __LINE__, infile );
      return ret;
      }

   /*
    *  set up various other info
    */
    toms_info->year = year;
    toms_info->doy = doy;
    sprintf( toms_info->node_time, "%04d%03d000000000", year, doy );
    toms_info->nlon = npix;
    toms_info->slon = -179.5;
    toms_info->elon = 179.5;
    toms_info->del_lon = 1.;
    toms_info->nlat = nlin;
    toms_info->slat = -89.5;
    toms_info->elat = 89.5;
    toms_info->del_lat = 1.;
    strcpy( toms_info->gen_str, "None" );
   /*
    *  Close it up
    */
    h5io_close( &fid );
    ret = 0;
    return ret;
    }
  else
    {
   /*
    *  if it fails, assume ascii and try to read
    */
    if (rd_toms_ascii(infile, toms_info) != 0)
        pexit("rd_toms_ascii");
    }
  ret = 0;
  return ret;
  }

/*****************************************************************
 * File:        check_usage
 *
 * Purpose:     check command line arguments
 *
 * Description: check command line arguements and report proper usage
 *              on argument count error or no argumnts.
 *
 * Input parms:
 * char *annotfile - input file of DAAC compliant annots ("fillenv.tomsnrt")
 * char *infile    - input ozone file to process ("nimoz100.y93")
 * char *directory - output file directory (".")
 *
 * Output parms:none
 *
 * Returns:     success (0) or failure(-1)
 *
 * Subs called: none
 *
 * History:     none
 *
 * Note:        This routine is custom for each program that uses it since
 *              the arguments are different of course.
 *
 * Author:      Brian D. Schieber, GSC, 9/93
 *
 * Modification history:
 *  W. Robinson, SAIC  Dec 2013  mod for use in o3_toms
 *****************************************************************/


int8 check_usage(int argc, char *argv[]) {
    if (argc < 4) {
        printf("\n\nUsage:\n");
        printf("\t%s <file> <directory> <prevfile>\n", argv[0]);
        printf("\nWhere:\n");
        printf("\tfile:       file to process\n");
        printf("\tdirectory:  output directory (no ending slash)\n");
        printf("\tprevfile:   previous day \n");
        printf("\n\n");
        printf("\tExample:\n");
        printf("\t\t o3_toms $SDSDEMO/fillenv.eptomsnrt $SDSDEMO/ga960903.ept\n");
        printf("\t\t ./ $SDSDEMO/ga960902.ept\n\n");
        return -1;
    }
    return 0;
}
