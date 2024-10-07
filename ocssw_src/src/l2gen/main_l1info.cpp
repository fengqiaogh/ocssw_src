/* ========================================================================
 *
 *   MSl1b2info input-l1b-filename
 *
 * Description:
 * 
 * Modification history:
 *
 *     Programmer     Organization      Date       Description of change
 *   --------------   ------------    --------     ---------------------
 *   Bryan A. Franz   GSC             28 March 1999 Original development
 *   Dan Knowles      SAIC            2004          Complete rewrite to create geoboxes
 *   W. Robinson      SAIC            10 dec 2004   add CZCS sensor
 *   Dan Knowles      SAIC            18 Aug 2005   Changed daynight logic to be based on solar zenith angle
 *   Sean Bailey      NASA            23 Oct 2015   Added logic to determine if sun in SD for VIIRS
 *   Gwyn Fireman     SAIC            2018-04-19    Determine if non-nominal pointing
 *   Steve Lockhart   SAIC	      2018-04-27    Add support for gring, changing from C to C++ 
 *
 * ======================================================================== */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>

#include "l12_proto.h"
#include <timeutils.h>
#include "version.h"
#include <stdbool.h>
#include "gringHelper.h"
#include "allocate2d.h"
#include <string>

#define NUM_SPACECRAFT_DIRECTION_MODES 3
#define DEFAULT_SPACECRAFT_DIRECTION 0   /* set to an invalid value */
#define ASCENDING_TRACK  1
#define DESCENDING_TRACK 2

#define NORTH 1
#define SOUTH 0

#ifdef TRUE
#undef TRUE
#endif

#define TRUE 1
#define FALSE 0

#define TMP_FILENAME_MAX 255

#define DEFAULT_COORD_VALUE -999.   /* set to an invalid value */

#define DEFAULT_DAYNIGHT_MODE 0
#define DAY_MODE (1<<0)
#define NIGHT_MODE (1<<1)

/* equatorial radius, didn't bother to use polar value also since precision not that important for current use */
#define EARTH_RADIUS_EQUATORIAL 6378   

#define MAX_ATTERR 1.0  /* degrees deviation considered as non-nominal attitude */

enum {
    ALL_NAVFAIL = 100,
    NO_VALID_NAV = 101,
    DIRECTION_NOT_SET = 102,
    DAY_NIGHT_NOT_SET = 103,
    GEOBOX_NOT_FOUND = 104,
    NON_NADIR = 105
};

typedef struct {
    float32 north_lat;
    float32 south_lat;
    float32 west_lon;
    float32 east_lon;
    unsigned char daynightflag;
} box_t;

typedef struct {
    int sensor_id;
    int day_node;
} day_node_t;


void set_north_south_boundaries(float32, float32 *, float32 *);
void set_west_east_boundaries(float32, float32 *, float32 *);
void set_west_east_boundaries2(float32, float32 *, float32 *);
int check_if_in_west_east_boundaries(float32, float32, float32);
int check_if_in_box(float32, float32, float32, float32, float32, float32);

#define USAGESTR \
"%s %d.%d.%d-%s (%s %s)\n\n\
Usage: %s [-n number-of-boxes] [-s] L1-filename [met-file]\n\n\
Where:\n\
\t-d degree limits to set on boxes (use instead of -n)\n\
\t-i set subsampling increment (default=1 gets all scan lines and pixels)\n\
\t-n generates up to this number of lat/lon boxes per orbit-side\n\
\t-s prints only data needed for database fields\n\
\t-v verify number of pixels that match in boxed region\n\
\t-o [output file]\n\n\
Return status codes:\n\
\t0 - good\n\
\t1 - general error\n\
\t100 - all pixels have NAVFAIL or NAVWARN flag set\n\
\t101 - no 'good' navigation was found\n\
\t102 - direction (ascending/descending) not set\n\
\t103 - day night flag not set\n\
\t104 - total box count = 0\n\
\t105 - non-nominal pointing\n\
(For multiple conditions, lowest value status code takes priority)\n\
\n\
"


#define PRINT_USAGE(x)  printf(USAGESTR,(x),VERSION_MAJOR,VERSION_MINOR,VERSION_PATCH,GITSHA,__DATE__,__TIME__,(x))

/* -------------------------------------------------------------------- *
 *                              main                                    *
 * -------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    l1str *l1rec; /* generic level-1b scan structure      */
    filehandle *l1file; /* input file handle                    */

    char ofile[FILENAME_MAX];
    FILE *fp;
    char isodatetime_first[30], isodatetime_last[30];

#if 0
    char geofile[FILENAME_MAX] = "";
    char geofile[TMP_FILENAME_MAX] = "";
#endif
    float32 mem_st_lat=DEFAULT_COORD_VALUE;
    float32 mem_st_lon=DEFAULT_COORD_VALUE;
    float32 mem_en_lat=DEFAULT_COORD_VALUE;
    float32 mem_en_lon=DEFAULT_COORD_VALUE;

    int32_t spix = 0; /* start pixel for subscene process   */
    int32_t cpix = -1; /* center pixel of subscene scan      */
    int32_t epix = -1; /* end pixel for subscene process     */
    int32_t sscan = 0; /* start scan for subscene process    */
    int32_t cscan = -1; /* center scan for subscene process   */
    int32_t escan = -1; /* end scan for subscene process      */
    int32_t iscan;
    int32_t ip;
    int32_t num_match = 0;
    int32_t num_no_match = 0;
    int32_t num_nav_fail = 0;
    int32_t num_pixels = 0;
    int32_t navwarn_all_set = 1;
    int32_t good_spix, good_epix, good_cpix;

    int box_index;
    int spacecraft_direction_index;
    int curr_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
    int prev_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
    int initial_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
    int final_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
    int match;
    int in_box;
    int num_boxes = 0;
    float num_degrees = 20;
    float max_num_degrees = 20;
    int show_standard_info = TRUE;
    int show_sdps_info = FALSE;
    int case_n = FALSE;
    int case_s = FALSE;
    int case_d = FALSE;
    int case_v = FALSE;
    int case_o = FALSE;
    int user_defined_increment = 1;
    int exitflag = EXIT_SUCCESS;
    int line_nav_ok;
    static int first_good_scan = -1;
    static int first_good_scantime = FALSE;
    double scantime_first, scantime_last;

    double percent_match;

    unsigned char curr_daynight_mode = DAY_MODE;
    unsigned char master_daynightflag = DEFAULT_DAYNIGHT_MODE;

    static char *node[3] = {"Unidentified", "Ascending", "Descending"};
    static char *daynight_string[4] = {"Unidentified", "Day", "Night", "Both"};

    float32 epix_northern_lat = DEFAULT_COORD_VALUE;
    float32 epix_southern_lat = DEFAULT_COORD_VALUE;
    float32 spix_northern_lat = DEFAULT_COORD_VALUE;
    float32 spix_southern_lat = DEFAULT_COORD_VALUE;
    float32 northern_lat = DEFAULT_COORD_VALUE;
    float32 southern_lat = DEFAULT_COORD_VALUE;
    float32 eastern_lon = DEFAULT_COORD_VALUE;
    float32 western_lon = DEFAULT_COORD_VALUE;
    float32 center_lon = DEFAULT_COORD_VALUE;
    float32 center_lat = DEFAULT_COORD_VALUE;
    float32 prev_lat_cpix = DEFAULT_COORD_VALUE;
    float32 lat_breakpoint;
    float32 increment;

    box_t ****box = NULL; /* will become 3 dimensional box array in [DNF][direction][box_num][lonFlag]  */
    int num_args;

    char *option_string = "d:hi:n:sv:o:";
    int options = 0;
    while ((options = getopt(argc, argv, option_string)) != -1) {
        switch (options) {
        case 'd':
            case_d = TRUE;
            num_degrees = atof(optarg);
            break;
        case 'i':
            user_defined_increment = atoi(optarg);
            break;
        case 'h':
            PRINT_USAGE(argv[0]);
            exit(EXIT_FAILURE);
            break;
        case 'n':
            case_n = TRUE;
            num_boxes = atoi(optarg);
            break;
        case 's':
            case_s = TRUE;
            break;
        case 'v':
            case_v = TRUE;
            break;
        case 'o':
            case_o = TRUE;
            snprintf(ofile, FILENAME_MAX, "%s", optarg);
            break;
        default:
            break;
        }
    }

    want_verbose = 0;

    /*******************************************************************************************
     *    Check options and set num_boxes
     *******************************************************************************************/
    if (case_d == TRUE) {
        if (num_degrees > 0 && num_degrees <= 180) {
            if (case_n != TRUE) {
                num_boxes = (int) (180 / num_degrees);
                printf("num_boxes=%d\n", num_boxes);
            } else {
                printf("INPUT ERROR: Cannot choose both -d and -n options\n");
                PRINT_USAGE(argv[0]);
                exit(EXIT_FAILURE);
            }
        } else {
            printf("INPUT ERROR: -d (degrees lat) option must be between 1 and 180\n");
            PRINT_USAGE(argv[0]);
            exit(EXIT_FAILURE);
        }
    } else if (case_n == TRUE) {
        if (num_boxes < 0) {
            printf("INPUT ERROR: -n (number of lat divisions) option must be greater than or equal to 0\n");
            PRINT_USAGE(argv[0]);
            exit(EXIT_FAILURE);
        }
    } else {
        num_boxes = 0;
    }

    /*******************************************************************************************
     *    Setup print options
     *******************************************************************************************/

    if (case_s == TRUE) {
        show_sdps_info = TRUE;
        show_standard_info = FALSE;
    }
    if (case_o == TRUE) {
        fp = fopen(ofile, "w");
    } else {
        fp = stdout;
    }

    num_args = argc - optind;

    // make the l1 readers quiet
    want_verbose = 0;

    cdata_();
    msl12_input_init();
    l1file = (filehandle*)allocateMemory(sizeof (filehandle), "filehandle structure");
    filehandle_init(l1file);
    l1rec = (l1str*)allocateMemory(sizeof (l1str), "l1rec structure");


    /*									*/
    /* Process required command-line arguments                          */
    /*									*/
    switch (num_args) {
    case 2:
        strcpy(l1file->name, argv[optind + 0]);
        /* strcpy(geofile,argv[optind+1]);*/
        l1file->geofile = argv[optind + 1];
        break;
    case 1:
        strcpy(l1file->name, argv[optind + 0]);
        break;
    case 0:
        PRINT_USAGE(argv[0]);
        exit(EXIT_SUCCESS);
        break;
    default:
        PRINT_USAGE(argv[0]);
        exit(EXIT_FAILURE);
        break;
    }

    /* Set default input parameters */
    if (msl12_input_defaults(l1file) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    input->proc_sst = 0;
    input->proc_ocean = 0;
    input->atmocor = 0;
    input->proc_land = 0;


    /*
     * The following tidbit calculates if the Sun is visible to the solar
     * diffuser for VIIRS so that the OBC files can be flagged to keep
     * Since this requires reading from the geolocation file, this bit is
     * done prior to opening the L1 file (which would also open the geo file)
     */
    int sun_in_sd = FALSE;
    int non_nadir = FALSE;
    int filledScans = 0;
    std::string instName = sensorId2InstrumentName(l1file->sensorID);

    // The bowtie effect...if an extract is small enough and on the scene
    // edge, the bowtie effect can confuse the geobox logic, so bump the
    // subsampling increment to make things happy

    //TODO: Check l1file->ndets ...should be a better way to do one per scan 
    if ((instName == "MODIS") &&
            (user_defined_increment < 10)) {
        user_defined_increment = 10;
    }
    if(instName == "VIIRS") {
        if (user_defined_increment < 16)
            user_defined_increment = 16;
        // Get the ACTUAL scans for printing the real number not the 
        // arbitrary dimension size
        int32_t ncid, status;
        if ((nc_open((const char*) input->ifile, NC_NOWRITE, &ncid)) == NC_NOERR) {
            nc_get_att_int(ncid, NC_GLOBAL, "number_of_filled_scans", &filledScans);
            status = nc_close(ncid);
            if (status != NC_NOERR) {
                printf("-E- %s: Error closing %s.\n", argv[0], (char*) input->ifile);
                exit(EXIT_FAILURE);
            }
        }
        int32_t geoid;
        //Only bother to do this if provided a geolocation file
        if (input->geofile[0] && (nc_open(input->geofile, NC_NOWRITE, &geoid)) == NC_NOERR) {

            int32_t grpid, varid, dimid;
            size_t nscans, nvectors;

            // don't try this for the SDR files...only the NASA L1 format files
            if ((nc_inq_ncid(geoid, "All_Data", &grpid)) != NC_NOERR) {
                status = nc_inq_dimid(geoid, "number_of_scans", &dimid);
                if (status != NC_NOERR) {
                    printf("-E- Error reading number_of_scans attribute.\n");
                    exit(EXIT_FAILURE);
                }
                nc_inq_dimlen(geoid, dimid, &nscans);
                status = nc_inq_dimid(geoid, "vector_elements", &dimid);
                if (status != NC_NOERR) {
                    printf("-E- Error reading vector_elements attribute.\n");
                    exit(EXIT_FAILURE);
                }
                nc_inq_dimlen(geoid, dimid, &nvectors);

                if ((nc_inq_ncid(geoid, "navigation_data", &grpid)) == NC_NOERR) {
                    int i = 0;
                    size_t start[] = {0, 0}; /* start at first value */
                    size_t cnt[] = {nscans, nvectors};

                    // sun in solar diffuser?
                    float **solar_inst = allocate2d_float(nscans, nvectors);
                    status = nc_inq_varid(grpid, "solar_inst", &varid);
                    if (status != NC_NOERR) {
                        printf("-E- Error finding solar_inst variable.\n");
                        exit(EXIT_FAILURE);
                    }
                    status = nc_get_vara_float(grpid, varid, start, cnt, &solar_inst[0][0]);
                    if (status != NC_NOERR) {
                        printf("-E- Error reading solar_inst variable.\n");
                        exit(EXIT_FAILURE);
                    }
                    for (i = 0; i < (int)nscans; i++) {
                        if (solar_inst[i][0] >= -1 && solar_inst[i][0] <= 1) {
                            float solar_inst_az = RADEG * atan2f(solar_inst[i][0], -solar_inst[i][2]);
                            if (solar_inst_az >= 100.0 && solar_inst_az <= 110.0) {
                                sun_in_sd = TRUE;
                                break;
                            }
                        }
                    }
                    free2d_float(solar_inst);

                    // attitude not nominal?
                    float **att_ang= allocate2d_float(nscans, nvectors);
                    status = nc_inq_varid(grpid, "att_ang_mid", &varid);
                    if (status != NC_NOERR) {  // try old name
                        status = nc_inq_varid(grpid, "att_ang", &varid);
                        if (status != NC_NOERR) {
                            printf("-E- Error finding attitude angles.\n");
                            exit(EXIT_FAILURE);
                        }
                    }
                    status = nc_get_vara_float(grpid, varid, start, cnt, &att_ang[0][0]);
                    if (status != NC_NOERR) {
                        printf("-E- Error reading attitude angles.\n");
                        exit(EXIT_FAILURE);
                    }
                    for (i = 0; (i < (int)nscans) && !non_nadir; i++) {
                        for (int j = 0; j < (int)nvectors; j++) {
                            if (fabs(att_ang[i][j]) > MAX_ATTERR) {
                                non_nadir = TRUE;
                                break;
                            }
                        }
                    }
                    free2d_float(att_ang);
                }
            }

            status = nc_close(geoid);
            if (status != NC_NOERR) {
                printf("-E- %s: Error closing %s.\n", argv[0], input->geofile);
                exit(EXIT_FAILURE);
            }
        }
    }

    /*									*/
    /* Open input file and get sensor and scan information from handle  */
    /*                                                                  */
    if (openl1(l1file) != 0) {
        printf("-E- %s: Error opening %s for reading.\n", argv[0], l1file->name);
        exit(EXIT_FAILURE);
    }

    /*									*/
    /* Allocate memory for L1 scan data			  	        */
    /*									*/
    if (alloc_l1(l1file, l1rec) == 0) {
        printf("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    spix = 0;
    epix = l1file->npix - 1;
    cpix = spix + (epix - spix + 1) / 2;

    sscan = 0;
    escan = l1file->nscan - 1;
    cscan = sscan + (escan - sscan + 1) / 2;

    l1file->spix = spix;
    l1file->epix = epix;

    if (show_standard_info == TRUE) {
        fprintf(fp, "Sensor=%s\n", sensorId2SensorName(l1file->sensorID));
    }

    if (show_sdps_info == TRUE || show_standard_info == TRUE) {
        fprintf(fp, "Orbit_Number=%d\n", l1file->orbit_number);
        if (l1file->format == FT_VIIRSL1A || l1file->format == FT_VIIRSL1BNC) {
            fprintf(fp, "Number_of_Scans=%d\n", filledScans);
        } else {
            fprintf(fp, "Number_of_Scans=%d\n", l1file->nscan);
        }
    }

    if (show_standard_info == TRUE) {
        fprintf(fp, "Number_of_Pixels_per_Scan=%d\n", l1file->npix);
    }


    if (num_boxes) {
        /***********************************************************************************************
         *    Allocate boxes
         ***********************************************************************************************/
        box = (box_t ****) calloc(2, sizeof (box_t***));
        if (box == NULL) {
            printf("can not allocate box\n");
            exit(EXIT_FAILURE);
        }
        for(int day_night = 0; day_night < 2; day_night++) {
            box[day_night] = (box_t ***) calloc(NUM_SPACECRAFT_DIRECTION_MODES, sizeof (box_t**));
            if (box[day_night] == NULL) {
                printf("can not allocate box[day_night]\n");
                exit(EXIT_FAILURE);
            }
            for(int direction = 0; direction < NUM_SPACECRAFT_DIRECTION_MODES; direction++) {
                box[day_night][direction] = (box_t **) calloc(num_boxes, sizeof (box_t*));
                if (box == NULL) {
                    printf("can not allocate box[day_night][direction]\n");
                    exit(EXIT_FAILURE);
                }
                for(box_index=0; box_index<num_boxes; box_index++) {
                    box[day_night][direction][box_index] = (box_t *) calloc(2, sizeof (box_t));
                    if (box[day_night][direction][box_index] == NULL) {
                        printf("can not allocate box[day_night][direction][box_index]\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }

        /***********************************************************************************************
         *    Initialize boxes
         ***********************************************************************************************/
        for(int day_night = 0; day_night < 2; day_night++) {
            for(int direction = 0; direction < NUM_SPACECRAFT_DIRECTION_MODES; direction++) {
                for (box_index = 0; box_index < num_boxes; box_index++) {
                    for(int lonFlag = 0; lonFlag < 2; lonFlag++) {
                        box[day_night][direction][box_index][lonFlag].north_lat = DEFAULT_COORD_VALUE;
                        box[day_night][direction][box_index][lonFlag].south_lat = DEFAULT_COORD_VALUE;
                        box[day_night][direction][box_index][lonFlag].west_lon = DEFAULT_COORD_VALUE;
                        box[day_night][direction][box_index][lonFlag].east_lon = DEFAULT_COORD_VALUE;
                        box[day_night][direction][box_index][lonFlag].daynightflag = DEFAULT_DAYNIGHT_MODE;
                    }
                }
            }
        }

    }

    /***************************************************************************
     *    Loop through all scan lines
     ***************************************************************************/

    /*
        for (iscan=sscan; iscan<=escan; iscan++) 
        {  */

    iscan = sscan;
    line_nav_ok = 0;

    // Before looping on scans, prepare to instantiate gringHelper
    gringConfig_t *gringConfig = new gringConfig_t;
    // Get scan direction for this sensor
    // Most instruments scan to the left (so, use default scan_direction = -1); however, some exceptions scan to the right. 
    // Set scan_direction = 1 for these exceptions.
    if (l1file->sensorID == CZCS || l1file->sensorID == OCI || l1file->sensorID == HARP2 || l1file->sensorID == SPEXONE) {
        gringConfig->scan_direction = 1;
    } else {
        gringConfig->scan_direction = -1;
    }
    
    // set the number of lines between direction test
    if (strcmp(sensorId2InstrumentName(l1file->sensorID), "MODIS") == 0) {
        gringConfig->direction_delta = 10;
    } else if(strcmp(sensorId2InstrumentName(l1file->sensorID), "VIIRS") == 0) {
        gringConfig->direction_delta = 16;
    } else {
        gringConfig->direction_delta = 5;
    }    
    
    if (num_degrees < max_num_degrees) {
        gringConfig->delta_degrees_lat = num_degrees;
    } else {
        gringConfig->delta_degrees_lat = max_num_degrees;
        //printf("Setting delta_degrees_lat = %d\n", (int)gringConfig->delta_degrees_lat);
    }
    // instantiate a gringHelper
    gringHelper *gring_helper = new gringHelper(gringConfig);

    // loop on scans
    while (iscan <= escan) {
        good_spix = -1;
        good_epix = -1;
        good_cpix = -1;

        readl1_lonlat(l1file, iscan, l1rec);
        /*
         *  get info for each line and remember the info if end 
         *  and center geolocation are OK
         */

        if (l1rec->scantime > 0) {
            scantime_last = l1rec->scantime;
            if (first_good_scantime == FALSE) {
                first_good_scantime = TRUE;
                scantime_first = l1rec->scantime;
            }
        }

        /*
         * establish good start, end, center pixels if any
         */
        for (ip = spix; ip <= epix; ip++) {
            if ((l1rec->flags[ip] & NAVFAIL) == 0) {
                if (good_spix == -1) {
                    good_spix = ip;
                }
                good_epix = ip;
                if ((ip <= cpix) || (good_cpix == -1))
                    good_cpix = ip;
                if ((l1rec->flags[ip] & NAVWARN) == 0)
                    navwarn_all_set = 0;
            }
        }
        if (good_spix != -1) {
            line_nav_ok = 1;
            if (first_good_scan == -1) {
                first_good_scan = iscan;
                if (show_standard_info == TRUE) {
                    fprintf(fp, "Upper_Left_Lon=%f\n", l1rec->lon[good_spix]);
                    fprintf(fp, "Upper_Left_Lat=%f\n", l1rec->lat[good_spix]);
                    fprintf(fp, "Upper_Right_Lon=%f\n", l1rec->lon[good_epix]);
                    fprintf(fp, "Upper_Right_Lat=%f\n", l1rec->lat[good_epix]);
                }
            }


            mem_st_lat = l1rec->lat[good_spix];
            mem_st_lon = l1rec->lon[good_spix];
            mem_en_lat = l1rec->lat[good_epix];
            mem_en_lon = l1rec->lon[good_epix];

        }

        if (iscan == escan) {
            // Last scan
            if (scantime_first > scantime_last) {
                strncpy(isodatetime_first, unix2isodate(scantime_last, 'G'), 30);
                strncpy(isodatetime_last, unix2isodate(scantime_first, 'G'), 30);
            } else {
                strncpy(isodatetime_last, unix2isodate(scantime_last, 'G'), 30);
                strncpy(isodatetime_first, unix2isodate(scantime_first, 'G'), 30);
            }
            if (show_standard_info == TRUE) {
                fprintf(fp, "\n");
                fprintf(fp, "INFO: FIRST SCAN LINE\n");
                fprintf(fp, "Start_Date=%s\n", isodatetime_first);
                fprintf(fp, "\n");
                fprintf(fp, "INFO: LAST SCAN LINE\n");
                fprintf(fp, "End_Date=%s\n", isodatetime_last);
                fprintf(fp, "Lower_Left_Lon=%f\n", mem_st_lon);
                fprintf(fp, "Lower_Left_Lat=%f\n", mem_st_lat);
                fprintf(fp, "Lower_Right_Lon=%f\n", mem_en_lon);
                fprintf(fp, "Lower_Right_Lat=%f\n", mem_en_lat);
            }

            if (show_sdps_info == TRUE) {
                fprintf(fp, "Start_Date=%.*s\n", 10, isodatetime_first);
                fprintf(fp, "Start_Time=%.*s\n", 12, isodatetime_first + 11);
                fprintf(fp, "End_Date=%.*s\n", 10, isodatetime_last);
                fprintf(fp, "End_Time=%.*s\n", 12, isodatetime_last + 11);
            }
        }


        /***************************************************************************************************
         *    Determine MAX and MIN values for epix_lat and spix_lat based on current and all preceeding scan lines
         ***************************************************************************************************/


        if (line_nav_ok) {
            set_north_south_boundaries(l1rec->lat[good_spix], &spix_northern_lat, &spix_southern_lat);
            set_north_south_boundaries(l1rec->lat[good_epix], &epix_northern_lat, &epix_southern_lat);
        }

        /**************************************************************
         *    Find lat and lon of the center pixel of the middle scan
         **************************************************************/

        if (center_lat == DEFAULT_COORD_VALUE || center_lon == DEFAULT_COORD_VALUE) {
            if (iscan >= cscan && (l1rec->flags[good_cpix] & NAVFAIL) == 0) {
                center_lat = l1rec->lat[good_cpix];
                center_lon = l1rec->lon[good_cpix];
            }
        }


        /***************************************************************************
         *    Determine spacecraft direction based on center pixel of the scan
         ***************************************************************************/

        if ((l1rec->flags[good_cpix] & NAVFAIL) == 0) {
            if (prev_lat_cpix != DEFAULT_COORD_VALUE) {
                if (l1rec->lat[good_cpix] > prev_lat_cpix) {
                    curr_spacecraft_direction = ASCENDING_TRACK;
                } else if (l1rec->lat[good_cpix] < prev_lat_cpix) {
                    curr_spacecraft_direction = DESCENDING_TRACK;
                } else {
                    curr_spacecraft_direction = prev_spacecraft_direction;
                }

                if (initial_spacecraft_direction == DEFAULT_SPACECRAFT_DIRECTION &&
                        curr_spacecraft_direction != DEFAULT_SPACECRAFT_DIRECTION) {
                    initial_spacecraft_direction = curr_spacecraft_direction;
                }

            }
            prev_lat_cpix = l1rec->lat[good_cpix];
            prev_spacecraft_direction = curr_spacecraft_direction;

        }


        /***************************************************************************
         *    Loop through all pixels 
         ***************************************************************************/

        if (good_spix != -1) {
            ip = good_spix;
            while (ip <= good_epix) {

                if ((l1rec->flags[ip] & NAVFAIL) == 0) {
                    /***************************************************************************
                     *    Determine day/night mode of pixel 
                     ***************************************************************************/

                    // needed to test for valid solz - invalid ones were incorrectly
                    // changing the daynight mode
                    if (l1rec->solz[ip] >= 0.0 && l1rec->solz[ip] <= 180.0) {
                        if (l1rec->solz[ip] < SOLZNIGHT) {
                            curr_daynight_mode = DAY_MODE;
                        } else {
                            curr_daynight_mode = NIGHT_MODE;
                        }
                        master_daynightflag |= curr_daynight_mode;
                    }

                    /***************************************************************************************************
                     *    Determine MAX and MIN values for both lat and lon based on current and all preceeding scan lines 
                     ***************************************************************************************************/
                    if (line_nav_ok) {
                        set_north_south_boundaries(l1rec->lat[ip], &northern_lat, &southern_lat);
                        set_west_east_boundaries2(l1rec->lon[ip], &western_lon, &eastern_lon);
                    }


                    if (num_boxes) {
                        /***************************************************************************
                         *    Determine and set boundaries for regular boxes
                         ***************************************************************************/

                        if (l1rec->lat[ip] >= -90. && l1rec->lat[ip] <= 90. && l1rec->lon[ip] >= -180. && l1rec->lon[ip] <= 180.) {
                            int i = 0;
                            box_index = num_boxes;
                            increment = 180. / num_boxes;

                            for (lat_breakpoint = 90. - increment; lat_breakpoint >= -90.; lat_breakpoint -= increment) {
                                if (l1rec->lat[ip] >= lat_breakpoint && l1rec->lat[ip] < (lat_breakpoint + increment)) {
                                    box_index = i;
                                    break;
                                }
                                i++;
                            }

                            if (box_index < num_boxes) {
                                int lonFlag;
                                if(l1rec->lon[ip] < 0)
                                    lonFlag = 0;
                                else
                                    lonFlag = 1; 
                                set_north_south_boundaries(l1rec->lat[ip],
                                        &box[curr_daynight_mode-1][curr_spacecraft_direction][box_index][lonFlag].north_lat,
                                        &box[curr_daynight_mode-1][curr_spacecraft_direction][box_index][lonFlag].south_lat);
                                set_west_east_boundaries(l1rec->lon[ip],
                                        &box[curr_daynight_mode-1][curr_spacecraft_direction][box_index][lonFlag].west_lon,
                                        &box[curr_daynight_mode-1][curr_spacecraft_direction][box_index][lonFlag].east_lon);
                                box[curr_daynight_mode-1][curr_spacecraft_direction][box_index][lonFlag].daynightflag = curr_daynight_mode;
                            }
                        } else {

                            if (l1rec->lat[ip] < -90 || l1rec->lat[ip] > 90 || l1rec->lon[ip] < -180 || l1rec->lon[ip] > 180) {
                                fprintf(stderr, "-W- lat/lon pixel out of range: lat=%f, lon=%f, pixel_index=%d\n", l1rec->lat[ip], l1rec->lon[ip], ip);
                            }

                        }
                    }
                } else {
                    num_nav_fail++; /* count failed pixels  */
                }

                num_pixels++; /* count total pixels */

                if (ip < epix - user_defined_increment) {
                    ip += user_defined_increment;
                } else {
                    ip++;
                }
            }
        }
        
        // gring-related processing (per scan)
        gring_helper->process_scan(l1rec, iscan, escan); 

        if (iscan < escan - user_defined_increment) {
            iscan += user_defined_increment;
        } else {
            iscan++;
        }
    }

    final_spacecraft_direction = curr_spacecraft_direction;



    /***************************************************************************************************
     *    Merge initial scans box with appropriate box now that we know the initial_spacecraft_direction 
     ***************************************************************************************************/

    if (num_boxes) {
        for(int day_night = 0; day_night < 2; day_night++) {
            for (box_index = 0; box_index < num_boxes; box_index++) {
                for(int lonFlag=0; lonFlag<2; lonFlag++) {
                    if (box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].north_lat != DEFAULT_COORD_VALUE) {
                        set_north_south_boundaries(box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].north_lat,
                                &box[day_night][initial_spacecraft_direction][box_index][lonFlag].north_lat,
                                &box[day_night][initial_spacecraft_direction][box_index][lonFlag].south_lat);
                        set_north_south_boundaries(box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].south_lat,
                                &box[day_night][initial_spacecraft_direction][box_index][lonFlag].north_lat,
                                &box[day_night][initial_spacecraft_direction][box_index][lonFlag].south_lat);
                        set_west_east_boundaries(box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].west_lon,
                                &box[day_night][initial_spacecraft_direction][box_index][lonFlag].west_lon,
                                &box[day_night][initial_spacecraft_direction][box_index][lonFlag].east_lon);
                        set_west_east_boundaries(box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].east_lon,
                                &box[day_night][initial_spacecraft_direction][box_index][lonFlag].west_lon,
                                &box[day_night][initial_spacecraft_direction][box_index][lonFlag].east_lon);
                        box[day_night][initial_spacecraft_direction][box_index][lonFlag].daynightflag = box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].daynightflag;

                    /* re initialize default spacecraft direction box */
                        box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].north_lat = DEFAULT_COORD_VALUE;
                        box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].south_lat = DEFAULT_COORD_VALUE;
                        box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].west_lon = DEFAULT_COORD_VALUE;
                        box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].east_lon = DEFAULT_COORD_VALUE;
                        box[day_night][DEFAULT_SPACECRAFT_DIRECTION][box_index][lonFlag].daynightflag = DEFAULT_DAYNIGHT_MODE;
                    }
                }
            }
        }

        if (num_pixels == num_nav_fail) {
            exitflag = ALL_NAVFAIL;
        } else if (northern_lat == DEFAULT_COORD_VALUE ||
                southern_lat == DEFAULT_COORD_VALUE ||
                eastern_lon == DEFAULT_COORD_VALUE ||
                western_lon == DEFAULT_COORD_VALUE) {
            exitflag = NO_VALID_NAV;
        } else if (case_n == TRUE && (initial_spacecraft_direction == DEFAULT_SPACECRAFT_DIRECTION ||
                final_spacecraft_direction == DEFAULT_SPACECRAFT_DIRECTION)
                ) {
            exitflag = DIRECTION_NOT_SET;
        }
    }

    if (navwarn_all_set == 1) {
        exitflag = ALL_NAVFAIL;
    }

    /***************************************************************************************************
     *    Print summary info 
     ***************************************************************************************************/

    if (show_standard_info == TRUE) {
        fprintf(fp, "\n");
        fprintf(fp, "SUMMARY STATS\n");
        fprintf(fp, "Northernmost_Lat=%f\n", northern_lat);
        fprintf(fp, "Southernmost_Lat=%f\n", southern_lat);
        fprintf(fp, "Easternmost_Lon=%f\n", eastern_lon);
        fprintf(fp, "Westernmost_Lon=%f\n", western_lon);
        fprintf(fp, "Center_Lat=%f\n", center_lat);
        fprintf(fp, "Center_Lon=%f\n", center_lon);
        fprintf(fp, "Start_Node=%s\n", node[initial_spacecraft_direction]);
        fprintf(fp, "End_Node=%s\n", node[final_spacecraft_direction]);
        fprintf(fp, "Daynight=%s\n", daynight_string[master_daynightflag]);
        fprintf(fp, "Moon_in_SV=%1d\n", l1file->sv_with_moon);
        if(instName == "VIIRS") {
            fprintf(fp, "Sun_in_SD=%1d\n", sun_in_sd);
            fprintf(fp, "Non_Nadir=%1d\n", non_nadir);
        }
    }

    if (show_sdps_info == TRUE) {
        fprintf(fp, "DayNightFlag=%d\n", master_daynightflag);
        fprintf(fp, "Moon_in_SV=%1d\n", l1file->sv_with_moon);
        if(instName == "VIIRS") {
            fprintf(fp, "Sun_in_SD=%1d\n", sun_in_sd);
            fprintf(fp, "Non_Nadir=%1d\n", non_nadir);
        }
    }

    if (exitflag == 0 && master_daynightflag == 0) {
        exitflag = DAY_NIGHT_NOT_SET;
    }


    /***************************************************************************************************
     *    Print box info 
     ***************************************************************************************************/

    if (show_sdps_info == TRUE && num_boxes) {
        int geobox_found = 0;
        int box_count = 1;
        for(int day_night = 0; day_night < 2; day_night++) {
            for (spacecraft_direction_index = 0; spacecraft_direction_index < NUM_SPACECRAFT_DIRECTION_MODES; spacecraft_direction_index++) {
                for (box_index = 0; box_index < num_boxes; box_index++) {
                    for(int lonFlag = 0; lonFlag < 2; lonFlag++) {
                        if (box[day_night][spacecraft_direction_index][box_index][lonFlag].north_lat != DEFAULT_COORD_VALUE) {
                        geobox_found = 1;
                        fprintf(fp, "GeoBox_%d=%f,%f,%f,%f,%d\n",
                                box_count,
                                    box[day_night][spacecraft_direction_index][box_index][lonFlag].north_lat,
                                    box[day_night][spacecraft_direction_index][box_index][lonFlag].south_lat,
                                    box[day_night][spacecraft_direction_index][box_index][lonFlag].west_lon,
                                    box[day_night][spacecraft_direction_index][box_index][lonFlag].east_lon,
                                    box[day_night][spacecraft_direction_index][box_index][lonFlag].daynightflag);
                        box_count++;
                        }
                    }
                }
            }
        }
        if (exitflag == 0 && case_n == TRUE && !geobox_found) {
            exitflag = GEOBOX_NOT_FOUND;
        }
    }

    if (case_v == TRUE) {

        for (iscan = sscan; iscan <= escan; iscan++) {
            readl1_lonlat(l1file, iscan, l1rec);

            for (ip = spix; ip <= epix; ip++) {
                if ((l1rec->flags[ip] & NAVFAIL) == 0) {
                    match = FALSE;

                    for(int day_night = 0; day_night < 2; day_night++) {
                        for (spacecraft_direction_index = 0; spacecraft_direction_index < NUM_SPACECRAFT_DIRECTION_MODES; spacecraft_direction_index++) {
                            for (box_index = 0; box_index < num_boxes; box_index++) {
                                for(int lonFlag = 0; lonFlag < 2; lonFlag++) {
                                in_box = check_if_in_box(l1rec->lat[ip],
                                        l1rec->lon[ip],
                                            box[day_night][spacecraft_direction_index][box_index][lonFlag].north_lat,
                                            box[day_night][spacecraft_direction_index][box_index][lonFlag].south_lat,
                                            box[day_night][spacecraft_direction_index][box_index][lonFlag].west_lon,
                                            box[day_night][spacecraft_direction_index][box_index][lonFlag].east_lon);
                                if (in_box == TRUE) {
                                    match = TRUE;
                                }
                            }
                        }
                    }
                    }
                    if (match == TRUE) {
                        num_match++;
                    } else {
                        num_no_match++;
                    }

                }
            }
        }

        if ((num_match + num_no_match) > 0) {
            percent_match = 100. * num_match / (num_match + num_no_match);
        } else {
            percent_match = 0;
        }

        fprintf(fp, "percent_match=%f, num_match=%d, num_no_match=%d\n", percent_match, num_match, num_no_match);
    }


    /***************************************************************************************************
     *    Print gring info 
     ***************************************************************************************************/ 
    int gring_geobox_cnt = gring_helper->get_geobox_cnt();
    if (gring_geobox_cnt > 1) {
        // Proceed to print gring info
        // Prep for call to  gring_helper->get_geobox_cnt()   
        std::string gring_lon_string;
        std::string gring_lat_string;
        std::string gring_seq_string;
        // Call gring_helper->get_gring_strings
        if (gring_helper->get_gring_strings(gring_lon_string, gring_lat_string, gring_seq_string) == SUCCESS) {
            // fprintf gring-related strings
            fprintf(fp, "gringpointlongitude=%s\n",gring_lon_string.c_str());
            fprintf(fp, "gringpointlatitude=%s\n",gring_lat_string.c_str());
            fprintf(fp, "gringpointsequence=%s\n",gring_seq_string.c_str());
        } else {
            fprintf(stderr, "-W- Could not generate gring\n");
        }
        
    } else {
        fprintf(stderr, "-W- Not enough scans to form a gring; number of scans to include in gring = %d\n",gring_geobox_cnt);
    }
    
    

    /* free the dynamically allocated memory */
    delete(gringConfig);
    delete(gring_helper);

    if (num_boxes > 0) {
        for(int day_night = 0; day_night < 2; day_night++) {
            for (spacecraft_direction_index = 0; spacecraft_direction_index < NUM_SPACECRAFT_DIRECTION_MODES; spacecraft_direction_index++) {
                for(box_index=0; box_index<num_boxes; box_index++) {
                    free(box[day_night][spacecraft_direction_index][box_index]);
                }
                free(box[day_night][spacecraft_direction_index]);
            }
            free(box[day_night]);
        }
        free(box);
    }
    if (case_o == TRUE)
        fclose(fp);

    if (exitflag == 0 && non_nadir) {
        exitflag = NON_NADIR;
    }

    exit(exitflag);
}

/* take first arg (lat) and adjust northern_boundary or southern_boundary if needed  */

void
set_north_south_boundaries(float32 lat, float32 *northern_boundary, float32 *southern_boundary) {
    if (lat > 90.) {
        lat = 90.;
    }

    if (lat < -90.) {
        lat = -90.;
    }

    if (*northern_boundary != DEFAULT_COORD_VALUE) {
        *northern_boundary = MAX(lat, *northern_boundary);
    } else {
        *northern_boundary = lat;
    }

    if (*southern_boundary != DEFAULT_COORD_VALUE) {
        *southern_boundary = MIN(lat, *southern_boundary);
    } else {
        *southern_boundary = lat;
    }

}

int
check_if_in_box(float32 lat, float32 lon, float32 northern_boundary, float32 southern_boundary,
        float32 western_boundary, float32 eastern_boundary) {
    int results;
    int east_west_results;
    int north_south_results;

    east_west_results = check_if_in_west_east_boundaries(lon, western_boundary, eastern_boundary);

    if (lat >= southern_boundary && lat <= northern_boundary) {
        north_south_results = TRUE;
    } else {
        north_south_results = FALSE;
    }

    if (north_south_results == TRUE && east_west_results == TRUE) {
        results = TRUE;
    } else {
        results = FALSE;
    }

    return results;
}

int
check_if_in_west_east_boundaries(float32 lon, float32 western_boundary, float32 eastern_boundary) {

    int results = FALSE;

    if (eastern_boundary >= western_boundary) {
        /* no date line crossing */

        if (lon >= western_boundary &&
                lon <= eastern_boundary) {
            results = TRUE;
        }
    } else {
        /* date line crossing */

        if (lon >= western_boundary ||
                lon <= eastern_boundary) {
            results = TRUE;
        }
    }

    return results;
}

void 
set_west_east_boundaries(float32 lon, float32 * western_boundary, float32 * eastern_boundary) {
    if (lon > 180.) {
        lon = 180.;
    }

    if (lon < -180.) {
        lon = -180.;
    }

    if(lon > *eastern_boundary || *eastern_boundary == DEFAULT_COORD_VALUE)
        *eastern_boundary = lon;
    if(lon < *western_boundary || *western_boundary == DEFAULT_COORD_VALUE)
        *western_boundary = lon;
}

void
set_west_east_boundaries2(float32 lon, float32 * western_boundary, float32 * eastern_boundary) {
    float32 boundary_width;
    float32 width_to_west;
    float32 width_to_east;

    if (*western_boundary != -180. || *eastern_boundary != 180.) {
        if (*eastern_boundary != DEFAULT_COORD_VALUE && *western_boundary != DEFAULT_COORD_VALUE) {
            if (check_if_in_west_east_boundaries(lon, *western_boundary, *eastern_boundary) != TRUE) {
                /************************************************************************
                 *    Determine longitude width to east of current boundary
                 ************************************************************************/

                if (lon > *eastern_boundary) {
                    /* dateline not crossed */
                    width_to_east = lon - *eastern_boundary;
                } else {
                    /* dateline crossed */
                    width_to_east = 360. + lon - *eastern_boundary;
                }

                /************************************************************************
                 *    Determine longitude width to west of current boundary
                 ************************************************************************/

                if (*western_boundary > lon) {
                    /* dateline not crossed */
                    width_to_west = *western_boundary - lon;
                } else {
                    /* dateline crossed */
                    width_to_west = *western_boundary + 360. - lon;
                }

                /************************************************************************
                 *    Set closest west-east boundary
                 ************************************************************************/

                if (fabs(width_to_west) <= fabs(width_to_east)) {
                    *western_boundary = lon;
                } else {
                    *eastern_boundary = lon;
                }

                /************************************************************************
                 *    Determine longitude width between western_boundary and eastern_boundary
                 ************************************************************************/

                if (*eastern_boundary >= *western_boundary) {
                    /* no date line crossing */
                    boundary_width = *eastern_boundary - *western_boundary;
                } else {
                    /* date line crossing */
                    boundary_width = 360. + *eastern_boundary - *western_boundary;
                }

                /************************************************************************
                 *    if west-to-east span > 355 then just set span to -180 to 180
                 ************************************************************************/

                if (boundary_width > 355.) {
                    *western_boundary = -180.;
                    *eastern_boundary = 180.;
                }
            }
        } else {
            *eastern_boundary = lon;
            *western_boundary = lon;
        }
    }
}
