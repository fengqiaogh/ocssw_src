/*-----------------------------------------------------------------------------
    Function: main (MSl2brsgen)

    Returns: int32 (status)
        Returns a status code of 0 when successful.  Otherwise returns
        -1 - to indicate get_l2_open/close error
        -2 - to indicate get_l2_record error
        -3 - to indicate insufficient memory error

    Description:
        The function MSl2brsgen calls level-2 data read routines to access
        chlorophil data from the level 2 file.  After subsampling the data
        calls level-2 browse data write routine to ouput the level-2 browse
        product.

    Arguments: (in calling order)
      Type       Name             I/O     Description
      ----       ----             ---     -----------
      char *     l2_path           I      directory path and filename for the
                                          level-2 product
      char *     l2brs_path        I      directory path and filename for the
                                          level-2 browse product
      char *     xsub              I      Pixel subsampling rate
      char *     ysub              I      Scan subsampling rate
      char *     browse_prod       I      Browse product name
      char *     rflag		   I      Replacement Flag

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Norman Kuring	 GSFC (902.3) 		     Original development
        Lakshmi Kumar    Hughes STX      07/14/94    Modified according I/O
                                                        specs.
        Lakshmi Kumar	 Hughes STX	 12/10/94    Removed l2_scale and l2_
                                                     offset fr get_l2_open call
                                                     (see I/O specs. v4.1)
        Lakshmi Kumar    HITC		 05/18/95    Added Tilt and Navigation
                                                     Vgroups and nav datasets
        Lakshmi Kumar    Hughes STX      07/25/95    Added code to subsample
                                                     l2_flags and pass it to
                                                     put_l2brs routine
        Lakshmi Kumar    Hughes STX      11/01/95    Added output arguments
                                                     "orbit" & "mskflg" to
                                                     get_l2_open call
                                                     Added input arguments
                                                     flag_names & mskflg to
                                                     put_l2brs call
        Lakshmi Kumar	 Hughes STX	 01/11/96    Missing data pixels are
                                                     marked by setting its
                                                     value = -9E6 .
                                                     Lat/Lon for pixels with
                                                     nav_flag set are reset to
                                                     NAVERR_VAL = -999.0
        Norman Kuring	 GSFC (902)	 02/10/98    nav_flag no longer affects
                                                     outline coordinates.
        Robert Lindsay	 SAIC GSC	 08/17/00    modified to exit if file
                                                     file has no chlor_a SDS

        Joel Gales	 Futuretech	 03/01/01    Added support for other
                                                     sensors

        Joel Gales	 Futuretech	 09/05/06    Added support for sst
                                                     qual check

        Joel Gales	 Futuretech	 09/20/06    Added support for sst4
                                                     qual check

        Joel Gales	 SAIC    	 04/22/21    Added support for
                                                     multiple bands per pixel
----------------------------------------------------------------------------*/

#include "l2brsgen.h"

#include <unistd.h>
#include <strings.h>
#include "setupflags.h"
#include <timeutils.h>
#include <genutils.h>
#include <get_product_table.h>

// the real bad lat lon value is -999, but we don't want to do an exact check on
// a float, so value bigger than -900 is a good check
#define BAD_LAT_LON -900.0

char ERR_MSG[1024];

static float initialLonVal = -180.0;

/* this function is used when comparing lons to see if the lon is.  The lons
 * are rotated such that the first lon in the range is zero then normalized
 * from 0 to 360 */
static float normalizeLon(float lon) {
    lon -= initialLonVal;
    while (lon < 0.0)
        lon += 360.0;
    while (lon >= 360.0)
        lon -= 360.0;
    return lon;
}

/**
 * set the max value if val is bigger and not -999
 * @param max pointer to max val
 * @param val current value
 */
static void setMaxLon(float *max, float val) {
    if (val > BAD_LAT_LON) {
        if (normalizeLon(val) > normalizeLon(*max))
            *max = val;
    }
}

/**
 * set the min value if the normalized val is smaller and val not -999
 * @param min pointer to min val
 * @param val current value
 */
static void setMinLon(float *min, float val) {
    if (val > BAD_LAT_LON) {
        if (normalizeLon(val) < normalizeLon(*min))
            *min = val;
    }
}

/**
 * set the max value if the normalized val is bigger and val not -999
 * @param max pointer to max val
 * @param val current value
 */
static void setMaxLat(float *max, float val) {
    if (val > BAD_LAT_LON) {
        if (val > *max)
            *max = val;
    }
}

/**
 * set the min value if val is smaller and not -999
 * @param min pointer to min val
 * @param val current value
 */
static void setMinLat(float *min, float val) {
    if (val > BAD_LAT_LON) {
        if (val < *min)
            *min = val;
    }
}

/**
 * set the value only if not -999
 * @param min pointer to min val
 * @param val current value
 */
static void setGoodLL(float *good, float val) {
    if (val > BAD_LAT_LON) {
        *good = val;
    }
}

int main(int32 argc, char *argv[]) {
    int status;
    char l2_path[FILENAME_MAX];
    char ofile[FILENAME_MAX];
    short syear, sday, eyear, eday;
    int smsec, emsec;
    char dtype[8], *flag_names = NULL;
    int nrec, nsamp;
    int32_t brs_nrec, brs_nsamp;
    int32_t brs_crec;
    int32_t max_rec_used, max_samp_used;
    int ntilts;
    short *tilt_ranges;

    static meta_l2Type meta_l2;
    int i;
    int32 *l2brs_flags, *flags;
    float32 *pix, *qual_pix;
    int rec, samp, brs_rec, brs_samp;
    int start_samp, start_rec, xsub = 2, ysub = 2, qual = 999;
    int end_samp, end_rec;
    int center_samp;
    int spixl, epixl, sline, eline;
    float32 *l2brs_data, *data;
    float *px_ll_first, *px_ll_last, *sc_ll_first, *sc_ll_last;
    float *top, *bottom, *left, *right;
    int16 *lat, *lon;
    char infiles[MAXVAL], *sp;
    char proc_con[4096];
    char replaces[MAXVAL], ptime[17];
    int32 apply_pal = 0;
    char palfile[MAXVAL];
    char palette_dir[FILENAME_MAX];
    char product_table[FILENAME_MAX];
    const char* oformat;
    int firstLoop;

    int16 *brs_lat = NULL, *brs_lon = NULL;
    char *flaguse, *chl_flags, *sst_flags;

    const float MISSING_DATAFLAG = -9.1E6;
    char browse_prod[128] = "chlor_a";
    char qual_prod1[128] = "qual_sst";
    char qual_prod2[128] = "qual_sst4";
    unsigned char palette[256 * 3];
    short r[256], g[256], b[256];

    static l2_prod l2_str;
    int32 prod_index;
    int32 qual_index;

    char* timePtr;

    product_table_t *ptable = NULL;
    int32_t ptable_size = 0;
    int32_t ptable_index = -1;
    product_table_t *ptable_rec;
    product_table_t tmp_ptable;

    setlinebuf(stdout);

    /* hold all of the command line options */
    clo_setEnableExtraOptions(0);
    clo_optionList_t* list = clo_createList();
    l2brsgen_init_options(list);

    /* make it print help if no options supplied */
    if (argc == 1) {
        clo_printUsage(list);
        exit(1);
    }

    /* make it print version if -version or --version */
    if (strcmp(argv[1], "-version") == 0 || strcmp(argv[1], "--version") == 0) {
        clo_printVersion(list);
        exit(1);
    }

    if (l2brsgen_read_options(list, argc, argv, &l2_str, &meta_l2) != 0) {
        clo_printUsage(list);
        exit(1);
    }

    proc_con[0] = 0;
    for (i = 0; i < argc; i++) {
        strcat(proc_con, argv[i]);
        strcat(proc_con, " ");
    }

    parse_file_name(clo_getString(list, "ifile"), l2_path);
    parse_file_name(clo_getString(list, "ofile"), ofile);
    strcpy(browse_prod, clo_getString(list, "prod"));
    qual = clo_getInt(list, "quality");
    strcpy(replaces, clo_getString(list, "rflag"));
    spixl = clo_getInt(list, "spixl");
    epixl = clo_getInt(list, "epixl");
    xsub = clo_getInt(list, "dpixl");
    sline = clo_getInt(list, "sline");
    eline = clo_getInt(list, "eline");
    ysub = clo_getInt(list, "dline");
    apply_pal = clo_getBool(list, "apply_pal");
    if (clo_isSet(list, "palfile")) {
        apply_pal = 1;
    }
    strcpy(palfile, clo_getString(list, "palfile"));
    parse_file_name(clo_getString(list, "palette_dir"), palette_dir);
    parse_file_name(clo_getString(list, "product_table"), product_table);

    /* Read product table */
    ptable = get_product_table(product_table, &ptable_size);
    if (ptable == NULL) {
        fprintf(stderr, "l2brsgen - error reading product table \"%s\"\n", product_table);
        exit(1);
    }

    ptable_index = search_product_table(ptable, ptable_size, browse_prod);
    if (ptable_index == -1) {
        if (clo_isSet(list, "datamin") &&
                clo_isSet(list, "datamax") &&
                clo_isSet(list, "stype")) {
            tmp_ptable.description = strdup(browse_prod);
            tmp_ptable.min = clo_getFloat(list, "datamin");
            tmp_ptable.max = clo_getFloat(list, "datamax");
            tmp_ptable.name = strdup(browse_prod);
            tmp_ptable.palette = strdup(palfile);
            tmp_ptable.precision = strdup("F");
            if (clo_getInt(list, "stype") == 2)
                tmp_ptable.scaling = strdup("logarithmic");
            else
                tmp_ptable.scaling = strdup("linear");

            tmp_ptable.units = strdup("unknown");
            ptable_rec = &tmp_ptable;
        } else {
            fprintf(stderr, "l2brsgen - error finding \"%s\" in product table \"%s\" and datamin, datamax, stype not set.\n",
                    browse_prod, product_table);
            exit(1);
        }
    } else {
        ptable_rec = ptable + ptable_index;
    }

    /* over ride values in the product table */
    if (strcasecmp(palfile, "default") != 0) {
        free(ptable_rec->palette);
        ptable_rec->palette = strdup(palfile);
    }
    if (clo_isSet(list, "datamin"))
        ptable_rec->min = clo_getFloat(list, "datamin");

    if (clo_isSet(list, "datamax"))
        ptable_rec->max = clo_getFloat(list, "datamax");

    if (clo_isSet(list, "stype")) {
        free(ptable_rec->scaling);
        if (clo_getInt(list, "stype") == 2)
            ptable_rec->scaling = strdup("logarithmic");
        else
            ptable_rec->scaling = strdup("linear");
    }

    flaguse = NULL;
    if (clo_isSet(list, "flaguse")) {
        flaguse = clo_getRawString(list, "flaguse");
    }
    chl_flags = clo_getRawString(list, "chl_flags");
    sst_flags = clo_getRawString(list, "sst_flags");

    oformat = getFileFormatName(clo_getString(list, "oformat"));

    /* print out some useful information */
    if (clo_isSet(list, "par")) {
        fprintf(stderr, "par=%s\n", clo_getString(list, "par"));
    }
    fprintf(stderr, "ifile=%s\n", l2_path);
    fprintf(stderr, "ofile=%s\n", ofile);
    fprintf(stderr, "prod=%s\n", browse_prod);
    fprintf(stderr, "oformat=%s\n", oformat);

    /* Read palette file */
    strcpy(palfile, palette_dir);
    strcat(palfile, "/");
    strcat(palfile, ptable_rec->palette);
    strcat(palfile, ".pal");

    if (getlut_file(palfile, r, g, b)) {
        fprintf(stderr, "l2brsgen: Error reading palette file %s\n", palfile);
        exit(1);
    }
    for (i = 0; i < 256; i++) {
        palette[i * 3] = r[i];
        palette[i * 3 + 1] = g[i];
        palette[i * 3 + 2] = b[i];
    }

    // already opened by l2brsgen_read_options
    //status = openL2(l2_path, 0x0, &l2_str);

    syear = l2_str.syear;
    sday = l2_str.sday;
    smsec = l2_str.smsec;
    eyear = l2_str.eyear;
    eday = l2_str.eday;
    emsec = l2_str.emsec;
    nrec = l2_str.nrec;
    nsamp = l2_str.nsamp;

    flag_names = (char *) calloc(strlen(l2_str.flagnames) + 1, sizeof (char));
    strcpy(flag_names, l2_str.flagnames);

    ntilts = l2_str.ntilts;
    tilt_ranges = (short *) calloc(20 * 2, sizeof (short));
    for (i = 0; i < ntilts; i++) {
        tilt_ranges[2 * i] = l2_str.tilt_ranges[0][i];
        tilt_ranges[2 * i + 1] = l2_str.tilt_ranges[1][i];
    }

    strcpy(dtype, l2_str.dtype);

    // already loaded by l2brsgen_read_options
    //status = readL2meta(&meta_l2, 0);

    /* Get the current time and represent it in the format yyyydddhhmmssfff */
    get_time(ptime);

    /*  Set starting (X,Y) position of browse data  */
    start_samp = spixl - 1;
    if (start_samp < 0)
        start_samp = 0;
    end_samp = epixl - 1;
    if (end_samp < 0)
        end_samp = nsamp - 1;
    if (start_samp >= end_samp) {
        fprintf(stderr, "l2brsgen: spixel (%d) needs to be less than epixl (%d)",
                start_samp + 1, end_samp + 1);
        fprintf(stderr, " cntl_pt_cols ");
        exit(1);
    }

    start_rec = sline - 1;
    if (start_rec < 0)
        start_rec = 0;
    end_rec = eline - 1;
    if (end_rec < 0)
        end_rec = nrec - 1;
    if (start_rec >= end_rec) {
        fprintf(stderr, "l2brsgen: sline (%d) needs to be less than eline (%d)",
                start_rec + 1, end_rec + 1);
        fprintf(stderr, " cntl_pt_cols ");
        exit(1);
    }

    /* Compute browse image dimensions */
    brs_nrec = (end_rec - start_rec) / ysub + 1;
    brs_nsamp = (end_samp - start_samp) / xsub + 1;

    /* Compute the maximum record and sample numbers used to make the subsampled
            browse output.  */
    max_rec_used = (brs_nrec - 1) * ysub + start_rec;
    max_samp_used = (brs_nsamp - 1) * xsub + start_samp;

    /* compute the center coodinates */
    brs_crec = brs_nrec / 2;
    center_samp = ((max_samp_used - start_samp) / 2) + start_samp;

    if ((brs_lat = (int16 *) calloc(brs_nrec * brs_nsamp, sizeof (int16))) == NULL) {
        fprintf(stderr, "\nl2brsgen: Calloc error while allocating memory for");
        fprintf(stderr, " brs_lat ");
        exit(MEMERR);
    }
    lat = brs_lat;

    if ((brs_lon = (int16 *) calloc(brs_nrec * brs_nsamp, sizeof (int16))) == NULL) {
        fprintf(stderr, "\nl2brsgen: Calloc error while allocating memory for");
        fprintf(stderr, " brs_lon ");
        exit(MEMERR);
    }
    lon = brs_lon;

    /* Allocate memory for accumulation of output data. */
    l2brs_data = (float32 *) malloc(brs_nrec * brs_nsamp * sizeof (float32));
    if (l2brs_data == NULL) {
        fprintf(stderr, "\nError: In allocating memory\n");
        exit(MEMERR);
    }
    data = l2brs_data;

    /* Allocate memory for accumulation of level2 flags */
    l2brs_flags = (int32 *) malloc(brs_nrec * brs_nsamp * sizeof (int32));
    if (l2brs_flags == NULL) {
        fprintf(stderr, "\nError: In allocating memory\n");
        exit(MEMERR);
    }
    flags = l2brs_flags;

    /* Allocate memory for accumulation of output coordinates. */
    px_ll_first = (float *) malloc(2 * brs_nsamp * sizeof (float));
    px_ll_last = (float *) malloc(2 * brs_nsamp * sizeof (float));
    sc_ll_first = (float *) malloc(2 * brs_nrec * sizeof (float));
    sc_ll_last = (float *) malloc(2 * brs_nrec * sizeof (float));
    if (!(px_ll_first && px_ll_last && sc_ll_first && sc_ll_last)) {
        fprintf(stderr, "\nError: In allocating memory\n");
        exit(MEMERR);
    }

    /* The following four names are somewhat arbitrary but intuitive. */
    top = px_ll_first;
    bottom = px_ll_last;
    left = sc_ll_first;
    right = sc_ll_last;

    prod_index = findprod(&l2_str, browse_prod);
    qual_index = findprod(&l2_str, qual_prod1);
    if (qual_index == -1) {
        qual_index = findprod(&l2_str, qual_prod2);
    }

    /* if no browse SDS in input file, quit and go home */

    if (prod_index == -1) {
        fprintf(stderr, "\n Error: Browse product (%s) not found in input file \n", browse_prod);

        /* Close the level 2 input file. */
        status = closeL2(&l2_str, 0);

        /* Deallocate some memory */
        free(l2brs_data);
        free(px_ll_first);
        free(px_ll_last);
        free(sc_ll_first);
        free(sc_ll_last);
        free(l2brs_flags);
        free(flag_names);
        free(brs_lat);
        free(brs_lon);
        free(tilt_ranges);

        status = freeL2meta(&meta_l2);

        exit(MEMERR);
    }

    /* set the initial lon before the loop and before resetting westlon.  Also
     * nudge the init point to the west a little bit for roundoff error. */
    setGoodLL(&initialLonVal, meta_l2.westlon);

    firstLoop = 1;

    /* Read the level 2 data into memory one line at a time
     * skipping lines according to the line subsampling factor.  */
    rec = start_rec;
    for (brs_rec = 0; brs_rec < brs_nrec; brs_rec++) {

        status = readL2(&l2_str, 0, rec, -1, NULL);
        pix = l2_str.l2_data[prod_index];
        if (strncmp(browse_prod, "sst", 3) == 0 && qual_index != -1)
            qual_pix = l2_str.l2_data[qual_index];
        else
            qual_pix = NULL;

        samp = start_samp;
        for (brs_samp = 0; brs_samp < brs_nsamp; brs_samp++) {

            if (qual_pix && qual_pix[samp] > qual)
                *data = MISSING_DATAFLAG;
            else
                *data = pix[samp];

            *flags = l2_str.l2_flags[samp];
            if (l2_str.latitude[samp] > BAD_LAT_LON)
                *lat = (int16) roundf(l2_str.latitude[samp] * 360.0);
            else
                *lat = -32768;
            if (l2_str.longitude[samp] > BAD_LAT_LON)
                *lon = (int16) roundf(l2_str.longitude[samp] * 180.0);
            else
                *lon = -32768;

            if (rec == start_rec) {
                *top++ = l2_str.latitude[samp];
                *top++ = l2_str.longitude[samp];
            }
            if (rec == max_rec_used) {
                *bottom++ = l2_str.latitude[samp];
                *bottom++ = l2_str.longitude[samp];
            }

            /* reset the lat and lon  limits first time through */
            if (firstLoop) {
                firstLoop = 0;
                meta_l2.northlat = l2_str.latitude[samp];
                meta_l2.southlat = l2_str.latitude[samp];
                meta_l2.westlon = l2_str.longitude[samp];
                meta_l2.eastlon = l2_str.longitude[samp];
            }

            /* set the min and max lat and lon*/
            setMaxLat(&meta_l2.northlat, l2_str.latitude[samp]);
            setMinLat(&meta_l2.southlat, l2_str.latitude[samp]);
            setMaxLon(&meta_l2.eastlon, l2_str.longitude[samp]);
            setMinLon(&meta_l2.westlon, l2_str.longitude[samp]);

            lat++;
            lon++;
            flags++;
            data++;
            samp += xsub;
        } // for brs_samp

        /* Save the ground coordinates of the two ends of this scan. */
        *left++ = l2_str.latitude[start_samp];
        *left++ = l2_str.longitude[start_samp];
        *right++ = l2_str.latitude[max_samp_used];
        *right++ = l2_str.longitude[max_samp_used];

        /* start record */
        if (brs_rec == 0) {
            syear = l2_str.year;
            sday = l2_str.day;
            smsec = l2_str.msec;

            // only change the metadata values if the lat lons are valid
            setGoodLL(&meta_l2.startclat, l2_str.latitude[center_samp]);
            setGoodLL(&meta_l2.startclon, l2_str.longitude[center_samp]);
        }

        /* center record */
        if (brs_rec == brs_crec) {
            meta_l2.ncrec = rec;
            if (l2_str.year != -999 && l2_str.day != -999 && l2_str.msec != -999) {
                timePtr = unix2ydhmsf(yds2unix(l2_str.year, l2_str.day, l2_str.msec / 1000.0), 'G');
                if (meta_l2.ctime)
                    free(meta_l2.ctime);
                meta_l2.ctime = strdup(timePtr);
            }
            // only change the metadata values if the lat lons are valid

            // How do I get the solar zenith ??
            // meta_l2.scsol_z = l2_str.

        } // if center record

        // store these for every line, that way the last valid value will
        // be saved.  Sometimes the last few lines are bad.
        setGoodLL(&meta_l2.endclat, l2_str.latitude[center_samp]);
        setGoodLL(&meta_l2.endclon, l2_str.longitude[center_samp]);
        if (l2_str.year != -999 && l2_str.day != -999 && l2_str.msec != -999) {
            eyear = l2_str.year;
            eday = l2_str.day;
            emsec = l2_str.msec;
        }

        rec += ysub;
    } // for brs_rec



    /* Write out a level 2 browse file */

    sp = strrchr(l2_path, '/'); /* write product name without path*/
    if (sp == NULL)
        strcpy(infiles, l2_path);
    else
        strcpy(infiles, ++sp);

    if (flaguse == NULL) {
        if (strncmp(browse_prod, "chl", 3) == 0) {
            flaguse = chl_flags;
        } else if (strncmp(browse_prod, "sst", 3) == 0) {
            flaguse = sst_flags;
        }
    }

    status = put_l2brs(ofile, replaces, ptime, infiles, start_samp + 1, end_samp + 1, xsub,
            brs_nsamp, start_rec + 1, end_rec + 1, ysub, brs_nrec, browse_prod, l2brs_data,
            l2brs_flags, flag_names,l2_str.l2_bits, flaguse, palette, px_ll_first, px_ll_last,
            sc_ll_first, sc_ll_last, proc_con, syear, sday, smsec, eyear, eday,
            emsec, dtype, nrec, nsamp, ntilts, &l2_str.tilt_flags[0],
            tilt_ranges, brs_lat, brs_lon, &meta_l2, ptable_rec, oformat, apply_pal);

    if (status < 0)
        fprintf(stderr, "\n Error: put_l2brs unsuccessful\n");

    /* Close the level 2 input file. */
    status = closeL2(&l2_str, 0);

    /* Deallocate some memory */

    free(l2brs_data);
    free(px_ll_first);
    free(px_ll_last);
    free(sc_ll_first);
    free(sc_ll_last);
    free(l2brs_flags);
    free(flag_names);
    free(brs_lat);
    free(brs_lon);
    free(tilt_ranges);

    status = freeL2meta(&meta_l2);

    return SUCCEED;
}

