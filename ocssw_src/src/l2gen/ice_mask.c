#include "l12_proto.h"

#define ICE_TYPE_OLD   0
#define ICE_TYPE_NSIDC 1
#define ICE_TYPE_OISST 2

static int ice_initalized = 0;
static int ice_file_type = ICE_TYPE_OLD;
static float ice_threshold = 0;

/***********************************************
 *
 * global variabls for the old ice mask
 *
 ***********************************************/

#define NX 4096
#define NY 2048

typedef char ice_t[NX];
static ice_t* ice;


/***********************************************
 *
 * global variables for NSIDC ice data
 *
 ***********************************************/

#define NAMING_CONVENTION_REFERENCE \
"http://nsidc.org/data/docs/daac/nsidc0051_gsfc_seaice.gd.html#namingconvention"

#define NROWS 448
#define NCOLS 304
#define SROWS 332
#define SCOLS 316

#define RE 6378.273 /* Earth's equatorial radius */
#define EC 0.081816153 /* Eccentricity */
#define E2 0.006693883 /* Eccentricity squared */
#define SLAT 70.0  /* Standard latitude */
#define CELL 25.0  /* Stereographic pixel dimension (km) */
#define NODATA -1.0

typedef unsigned char north_t[NCOLS];
static north_t* north;

typedef unsigned char south_t[SCOLS];
static south_t* south;

static double slat = 0.0, sinslat, tc, mc;


/***********************************************
 *
 * global variables for NSIDC the OISST ice data file.
 *
 ***********************************************/

static size_t OI4NX = 1440;
static size_t OI4NY = 720;

//typedef float iceref_t[OI4NX + 2];
static float **iceref;

/**
 * Initialize the ice data using the given old ice HDF file.  The data set
 * has a bit for each month of the year.
 *
 * @param file name of the ice file
 * @param day day of year
 * @param sd_id opened HDF file id
 * @return 0=good, else error
 */
int ice_mask_init_old(char *file, int day, int32 sd_id) {
    int32 sds_id;
    int32 status;
    int32 sds_index;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int32 start[2];
    int32 edges[2];
    char name[H4_MAX_NC_NAME];

    int16 ice_data[NX];
    int16 mon = (int) day / 31; /* month of year (no need for perfection) */
    int16 bit = pow(2, mon); /* bit mask associated with this month    */
    int16 i, j;

    // allocate ice global data
    ice = (ice_t*) allocateMemory(NY * sizeof (ice_t), "ice");

    /* Get the SDS index */
    sds_index = SDnametoindex(sd_id, "ice_mask");
    if (sds_index == -1) {
        printf("-E- %s:  Error seeking ice_mask SDS from %s.\n",
                __FILE__, file);
        return (1);
    }

    /* Select the SDS */
    sds_id = SDselect(sd_id, sds_index);

    /* Verify the characteristics of the array */
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS info for ice_mask from %s.\n",
                __FILE__, file);
        return (1);
    }
    if (dims[0] != NY || dims[1] != NX) {
        printf("-E- %s:  Dimension mis-match in ice_mask file %s.\n  Expecting %d x %d\n  Reading   %d x %d\n",
                __FILE__, file, NX, NY, (int) dims[1], (int) dims[0]);
        return (1);
    }

    for (j = 0; j < NY; j++) {

        /* Read one row from the array     */
        start[0] = j; /* row offset */
        start[1] = 0; /* col offset */
        edges[0] = 1; /* row count  */
        edges[1] = dims[1]; /* col count  */

        status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) ice_data);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS ice_mask from %s.\n",
                    __FILE__, file);
            return (1);
        }

        /* Set byte array for if monthly bit set in ice data */
        for (i = 0; i < NX; i++)
            ice[j][i] = (ice_data[i] & bit) > 0;
    }

    /* Terminate access to the array */
    status = SDendaccess(sds_id);

    printf("Loaded old ice climatology HDF file.\n");
    return (0);
}

char ice_mask_old(float lon, float lat) {
    int16 i, j;

    i = MAX(MIN((int16) ((lon + 180.0) / 360.0 * NX), NX - 1), 0);
    j = MAX(MIN((int16) ((90.0 - lat) / 180.0 * NY), NY - 1), 0);
    return (ice[j][i]);
}

float get_icefrac_old(float lon, float lat) {
    if (ice_mask_old(lon, lat))
        return 1.0;

    return 0.0;
}


#define COPYNAME(orig,copy){     \
(copy) = (char *)strdup( orig );    \
if( (copy) == NULL ){      \
  fprintf(stderr,"-E- %s line %d: Memory allocation error\n", \
  __FILE__,__LINE__);      \
  exit(EXIT_FAILURE);      \
}        \
}

#define READIT(file,handle,rows,cols,array){    \
if( ( (handle) = fopen((file),"rb")) == NULL ){    \
  fprintf(stderr,"-E- %s line %d: Could not open file, \"%s\" .\n", \
  __FILE__,__LINE__,(file));      \
  perror("");        \
  exit(EXIT_FAILURE);       \
}         \
if( fseek( (handle) ,300,SEEK_SET) == -1 ){    \
  fprintf(stderr,       \
  "-E- %s line %d: Could not skip header of file, \"%s\" .\n",  \
  __FILE__,__LINE__,(file));      \
  perror("");        \
  exit(EXIT_FAILURE);       \
}         \
if(         \
fread((array), sizeof(unsigned char), (rows)*(cols), handle)  \
!= (rows)*(cols)       \
){         \
  fprintf(stderr,"-E- %s line %d: Error reading file, \"%s\" .\n", \
  __FILE__,__LINE__,(file));      \
  perror("");        \
  exit(EXIT_FAILURE);       \
}         \
}

int ice_mask_init_nsidc_raw(char *file) {

    /*
      Read ice data from the north and south pole files the first time
      this function is called and whenever the icefile name changes.
      Subsequent calls just refer to ice values stored in static arrays.
     */

    char *nfile, *sfile, *p;
    FILE *fh;

    COPYNAME(file, nfile);
    COPYNAME(file, sfile);

    // allocate global arrays
    north = (north_t*) allocateMemory(NROWS * sizeof (north_t), "ice north");
    south = (south_t*) allocateMemory(SROWS * sizeof (south_t), "ice south");

    /*
      According to the file naming convention, the filename for the
      north polar file has an 'n' immediately befor the last '.' in
      the name; the south polar file has an 's' in the same position
      with all other characters being the same.  Both files are expected
      to be present in the same directory.
     */
    p = strrchr(nfile, '.');
    if (p == NULL || p <= nfile || (*(p - 1) != 'n' && *(p - 1) != 's')) {
        printf("-E- %s line %d: ", __FILE__, __LINE__);
        printf("File, \"%s\", doesn't follow naming convention described at %s .\n",
                file, NAMING_CONVENTION_REFERENCE);
        return 1;
    }
    --p;
    *p = 'n';
    p = strrchr(sfile, '.');
    --p;
    *p = 's';

    READIT(nfile, fh, NROWS, NCOLS, north);
    fclose(fh);
    free(nfile);

    READIT(sfile, fh, SROWS, SCOLS, south);
    fclose(fh);
    free(sfile);

    printf("Loaded raw NSIDC ice files.\n");
    return 0;
}

int ice_mask_init_monthly(char *file, int year, int day, int32 sd_id) {
    int32 sds_id;
    int32 status;
    int32 sds_index;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int32 start[2];
    int32 edges[2];
    char name[H4_MAX_NC_NAME];
    int32 startYear;
    int32 endYear;

    int16 month, dom;

    char northName[32];
    char southName[32];

    // allocate global arrays
    north = (north_t*) allocateMemory(NROWS * sizeof (north_t), "ice north");
    south = (south_t*) allocateMemory(SROWS * sizeof (south_t), "ice south");

    /* Find the file attribute named "start_year". */
    if (SDreadattr(sd_id, SDfindattr(sd_id, "start_year"), (VOIDP) (&startYear))) {
        printf("-E- %s:  Error reading start_year from monthly file %s.\n",
                __FILE__, file);
        return 1;
    }

    /* Find the file attribute named "end_year". */
    if (SDreadattr(sd_id, SDfindattr(sd_id, "end_year"), (VOIDP) (&endYear))) {
        printf("-E- %s:  Error reading end_year from monthly file %s.\n",
                __FILE__, file);
        return 1;
    }

    // if the requested year is not in the monthly file
    // use the closest year avaliable

    if (year < startYear)
        year = startYear;
    if (year > endYear)
        year = endYear;

    // create the data set name for the monthly data
    yd2md(year, day, &month, &dom);
    sprintf(northName, "ice_%d_%02d_north", year, month);
    sprintf(southName, "ice_%d_%02d_south", year, month);

    //
    // read north array
    //

    /* Get the SDS index */
    sds_index = SDnametoindex(sd_id, northName);
    if (sds_index == -1) {
        printf("-E- %s:  Error seeking %s SDS from %s.\n",
                __FILE__, northName, file);
        return (1);
    }

    /* Select the SDS */
    sds_id = SDselect(sd_id, sds_index);

    /* Verify the characteristics of the array */
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS info for %s from %s.\n",
                __FILE__, northName, file);
        return (1);
    }
    if (dims[0] != NROWS || dims[1] != NCOLS) {
        printf("-E- %s:  Dimension mis-match in %s file %s.\n",
                __FILE__, northName, file);
        printf("  Expecting %d x %d\n", NROWS, NCOLS);
        printf("  Reading   %d x %d\n", (int) dims[0], (int) dims[1]);
        return (1);
    }

    // read the north array
    start[0] = 0; /* row offset */
    start[1] = 0; /* col offset */
    edges[0] = dims[0]; /* row count  */
    edges[1] = dims[1]; /* col count  */

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) north);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS %s from %s.\n",
                __FILE__, northName, file);
        return (1);
    }

    /* Terminate access to the north array */
    status = SDendaccess(sds_id);

    //
    // read south array
    //

    /* Get the SDS index */
    sds_index = SDnametoindex(sd_id, southName);
    if (sds_index == -1) {
        printf("-E- %s:  Error seeking %s SDS from %s.\n",
                __FILE__, southName, file);
        return (1);
    }

    /* Select the SDS */
    sds_id = SDselect(sd_id, sds_index);

    /* Verify the characteristics of the array */
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS info for %s from %s.\n",
                __FILE__, southName, file);
        return (1);
    }
    if (dims[0] != SROWS || dims[1] != SCOLS) {
        printf("-E- %s:  Dimension mis-match in %s file %s.\n",
                __FILE__, southName, file);
        printf("  Expecting %d x %d\n", SROWS, SCOLS);
        printf("  Reading   %d x %d\n", (int) dims[0], (int) dims[1]);
        return (1);
    }

    // read the south array
    start[0] = 0; /* row offset */
    start[1] = 0; /* col offset */
    edges[0] = dims[0]; /* row count  */
    edges[1] = dims[1]; /* col count  */

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) south);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS %s from %s.\n",
                __FILE__, southName, file);
        return (1);
    }

    /* Terminate access to the south array */
    status = SDendaccess(sds_id);

    printf("Loaded monthly NSIDC ice climatology HDF file.\n");
    return 0;
}

int ice_mask_init_daily(char *file, int year, int day, int32 sd_id) {
    int32 sds_id;
    int32 status;
    int32 sds_index;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int32 start[2];
    int32 edges[2];
    char name[H4_MAX_NC_NAME];

    // allocate global arrays
    north = (north_t*) allocateMemory(NROWS * sizeof (north_t), "ice north");
    south = (south_t*) allocateMemory(SROWS * sizeof (south_t), "ice south");

    //
    // read north array
    //

    /* Get the SDS index */
    sds_index = SDnametoindex(sd_id, "north");
    if (sds_index == -1) {
        printf("-E- %s:  Error seeking north SDS from %s.\n", __FILE__, file);
        return (1);
    }

    /* Select the SDS */
    sds_id = SDselect(sd_id, sds_index);

    /* Verify the characteristics of the array */
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS info for north from %s.\n",
                __FILE__, file);
        return (1);
    }
    if (dims[0] != NROWS || dims[1] != NCOLS) {
        printf("-E- %s:  Dimension mis-match in north file %s.\n",
                __FILE__, file);
        printf("  Expecting %d x %d\n", NROWS, NCOLS);
        printf("  Reading   %d x %d\n", (int) dims[0], (int) dims[1]);
        return (1);
    }

    // read the north array
    start[0] = 0; /* row offset */
    start[1] = 0; /* col offset */
    edges[0] = dims[0]; /* row count  */
    edges[1] = dims[1]; /* col count  */

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) north);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS north from %s.\n",
                __FILE__, file);
        return (1);
    }

    /* Terminate access to the north array */
    status = SDendaccess(sds_id);

    //
    // read south array
    //

    /* Get the SDS index */
    sds_index = SDnametoindex(sd_id, "south");
    if (sds_index == -1) {
        printf("-E- %s:  Error seeking south SDS from %s.\n", __FILE__, file);
        return (1);
    }

    /* Select the SDS */
    sds_id = SDselect(sd_id, sds_index);

    /* Verify the characteristics of the array */
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS info for south from %s.\n",
                __FILE__, file);
        return (1);
    }
    if (dims[0] != SROWS || dims[1] != SCOLS) {
        printf("-E- %s:  Dimension mis-match in south file %s.\n",
                __FILE__, file);
        printf("  Expecting %d x %d\n", SROWS, SCOLS);
        printf("  Reading   %d x %d\n", (int) dims[0], (int) dims[1]);
        return (1);
    }

    // read the south array
    start[0] = 0; /* row offset */
    start[1] = 0; /* col offset */
    edges[0] = dims[0]; /* row count  */
    edges[1] = dims[1]; /* col count  */

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) south);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS south from %s.\n",
                __FILE__, file);
        return (1);
    }

    /* Terminate access to the south array */
    status = SDendaccess(sds_id);

    printf("Loaded near real time NSIDC ice HDF file.\n");
    return 0;
}

float get_icefrac_nsidc(float lon, float lat) {

    double sinlat, t, rho, x, y, ii, jj, i, j, u, v, d[4], sum, interp;
    float delta, xdist, ydist;
    int sgn, rows, cols, ix, jy, cnt, n;
    unsigned char *data;


    if (lat >= 90 || lat <= -90) {
        return 0.0;
    }

    /*
    The north polar stereographic projection has the 45 West meridian
    running vertically from pole to bottom center;  the south polar
    stereographic projection has the 0 meridian running vertically
    from pole to top center.
     */
    if (lat >= 0) {
        delta = 45;
        sgn = 1;
        xdist = 3850.0; /* distance from pole to left edge */
        ydist = 5350.0; /* distance from pole to bottom edge */
        rows = NROWS;
        cols = NCOLS;
        data = &north[0][0];
    } else {
        delta = 0;
        sgn = -1;
        lat = -lat;
        xdist = 3950.0; /* distance from pole to left edge */
        ydist = 3950.0; /* distance from pole to bottom edge */
        rows = SROWS;
        cols = SCOLS;
        data = &south[0][0];
    }

    while (lon < 0) lon += 360;
    while (lon > 360) lon -= 360;

    lon += delta;

    /* Convert input coordinates to radians. */
    lat *= PI / 180;
    lon *= PI / 180;

    sinlat = sin((double) lat);

    /*
    The FORTRAN routines from the sidads.colorado.edu site use the
    ellipsoidal form of the stereographic projection equations, so
    I use them here as well; however, I think this may be a bit of
    overkill at 25-kilometer resolution.
     */

    /*
    Equation 15-9 from page 161 of
    Map Projections -- A Working Manual
    John P. Snyder
    U.S. Geological Survey Professional Paper 1395
    Fourth Printing 1997
     */
    t = tan((double) (PI / 4 - lat / 2))
            / pow((double) ((1 - EC * sinlat) / (1 + EC * sinlat)), (double) (EC / 2));

    /*
    The variables, tc and mc, are derived from the latitude of true
    scale, slat, and do not depend on the input parameters, so I only
    need to compute them once.  Actually, I could just replace these
    by constants, but the following lines make it clearer where the
    constants come from.
     */
    if (slat == 0.0) {
        slat = SLAT * PI / 180;
        sinslat = sin(slat);

        /*
        Equation 15-9 from page 161 of aforementioned publication.
         */
        tc = tan((double) (PI / 4 - slat / 2))
                / pow((double) ((1 - EC * sinslat) / (1 + EC * sinslat)), (double) (EC / 2));

        /*
        Equation 14-15 from page 101 of aforementioned publication.
         */
        mc = cos((double) slat) / sqrt((double) (1 - E2 * sinslat * sinslat));
    }

    /*
    Equation 21-34 from page 161 of aforementioned publication.
     */
    rho = RE * mc * t / tc;

    /*
    Equations 21-30 and 21-31 from page 161 of aforementioned publication.
     */
    x = rho * sgn * sin((double) (sgn * lon));
    y = -rho * sgn * cos((double) (sgn * lon));

    ii = (x + xdist - CELL / 2) / (CELL);
    jj = rows - (y + ydist - CELL / 2) / (CELL);

    /*
    Interpolate data from ice-cover files to get ice fraction
    at the specified coordinate.  Make sure to handle cases
    where data is flagged as missing (value > 250) for one
    reason or another.  Also handle cases where input coordinate
    maps outside of the area covered by the ice data files.
     */

    if (ii < 0 || ii > cols || jj < 0 || jj > rows) {
        /*
        The input coordinate pair was outside of the areas
        represented by the northern and southern data files.
         */
        return 0;
    }

    d[0] = *(data + (int) jj * cols + (int) ii);
    if (d[0] > 250) {
        /*
        The input coordinate pair mapped to a pixel that
        is marked as something other than an ice fraction.
         */
        return 0.0;
    }

    /*
    Choose the 4 pixels whose centers mark the corners of a
    square that encloses the user's specified input coordinate.
    Put the values of the 4 pixels in the array, d.
     */
    u = modf(ii, &i);
    v = modf(jj, &j);

    if (u < 0.5) {
        i--;
        u += 0.5;
    } else {
        u -= 0.5;
    }
    if (v < 0.5) {
        j--;
        v += 0.5;
    } else {
        v -= 0.5;
    }
    sum = 0;
    cnt = 0;
    n = 0;
    for (jy = j; jy < j + 2; jy++) {
        for (ix = i; ix < i + 2; ix++) {
            if (ix >= 0 && ix < cols && jy >= 0 && jy < rows) {

                /*
                This pixel is within the bounds of the
                reference array, so store its value.
                 */
                d[n] = *(data + jy * cols + ix);

                if (d[n] <= 250) {
                    /*
                    This is a valid value so include it in the mean that
                    will be used to fill in flagged or out-of-bound pixels.
                     */
                    sum += d[n];
                    cnt++;
                }
            } else {
                /*
                Use a flag value of 255 to indicate
                that this pixel is out of bounds.
                 */
                d[n] = 255;
            }
            n++;
        }
    }
    for (n = 0; n < 4; n++) {
        /*
        Replace missing values with mean of other non-missing values.
        This only affects pixels adjacent to the one the actually
        contains the desired Earth coordinate.  If the desired coordinate
        did not fall within a valid ice pixel, then this function already
        returned a NODATA value higher up in the code.
         */
        if (d[n] > 250) {
            d[n] = sum / cnt;
        }
    }

    /* Bi-linear interpolation. */
    interp = (1 - u)*(1 - v) * d[0]
            + u * (1 - v) * d[1]
            + (1 - u) * v * d[2]
            + u * v * d[3];

    return (float) (interp / 250.0);
}

char ice_mask_nsidc(float lon, float lat) {
    if (get_icefrac_nsidc(lon, lat) > ice_threshold)
        return 1;
    else
        return 0;
}

/* ------------------------------------------------------------------ *
 * read and interpolate Reynolds 0.25-deg daily V2 netcdf OI files
 *
 * return 0 if OK or 1 if error                                                                                     *
 * ------------------------------------------------------------------ */
int ice_init_oisst(char *icefile) {

    int i, j, ii;
    char* varName = "ice";
    int ncid, ndims, nvars, ngatts, unlimdimid;
    int rotateLon = 1;
    int varid;
    size_t length;
    float slope;
    float offset;


    if (nc_open(icefile, NC_NOWRITE, &ncid) != NC_NOERR) {
        return 1;
    }

    if (nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid) != NC_NOERR) {
        nc_close(ncid);
        return 1;
    }

    if (nc_inq_attlen(ncid, NC_GLOBAL, "title", &length) != NC_NOERR) {
        nc_close(ncid);
        return 1;
    }
    char titleStr[length + 2];
    if (nc_get_att_text(ncid, NC_GLOBAL, "title", titleStr) != NC_NOERR) {
        nc_close(ncid);
        return 1;
    }
    titleStr[length] = 0;

    if (strstr(titleStr, "Daily-OI") ||
        strstr(titleStr, "NOAA/NCEI 1/4 Degree Daily Optimum Interpolation Sea Surface Temperature (OISST) Analysis")) {
        printf("Loading Daily V2 0.25-deg OISST Ice field from %s\n\n", icefile);
    } else if (strstr(titleStr, "CMC")) {
        varName = "sea_ice_fraction";
        rotateLon = 0;
        int XdimID, YdimID;
        nc_inq_dimid(ncid, "lat", &YdimID);
        nc_inq_dimlen(ncid, YdimID, &OI4NY);
        nc_inq_dimid(ncid, "lon", &XdimID);
        nc_inq_dimlen(ncid, XdimID, &OI4NX);
        printf("Loading Daily CMC Ice field from %s\n\n", icefile);
    } else {
        nc_close(ncid);
        return 1;
    }

    if (nc_inq_varid(ncid, varName, &varid) != NC_NOERR) {
        nc_close(ncid);
        return 1;
    }

    // Allocate space for the ice data array
    static short **icetmp;
    if (iceref == NULL){
        iceref = allocate2d_float(OI4NY+2, OI4NX+2);
        icetmp = allocate2d_short(OI4NY, OI4NX);
    }
//    typedef int16 icetmp_t[OI4NX];
//    icetmp_t *icetmp;
//    icetmp = (icetmp_t*) malloc(OI4NY * sizeof (icetmp_t));

    /* Read the data. */
    if (nc_get_var_short(ncid, varid, &icetmp[0][0]) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                __FILE__, __LINE__, varName, icefile);
        exit(EXIT_FAILURE);
    }

    if (nc_get_att_float(ncid, varid, "scale_factor", &slope) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: error reading scale factor.\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    if (nc_get_att_float(ncid, varid, "add_offset", &offset) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: error reading scale offset.\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    nc_close(ncid);

    /* rotate 180-deg and add wrapping border to simplify interpolation */
    /* new grid is -180.125,180.125 in i=0,1441 and -90.125,90.125 in j=0,721 */
    for (j = 0; j < OI4NY; j++) {
        for (i = 0; i < OI4NX; i++) {
            if (rotateLon) {
                ii = (i < OI4NX / 2) ? i + OI4NX / 2 : i - OI4NX / 2;
            } else {
                ii = i;
            }
            if (icetmp[j][i] > -999)
                iceref[j + 1][ii + 1] = icetmp[j][i] * slope + offset;
            else
                iceref[j + 1][ii + 1] = BAD_FLT;
        }
        iceref[j + 1][0] = iceref[j + 1][OI4NX];
        iceref[j + 1][OI4NX + 1] = iceref[j + 1][1];
    }
    for (i = 0; i < OI4NX + 2; i++) {
        iceref[0][i] = iceref[1][i];
        iceref[OI4NY + 1][i] = iceref[OI4NY][i];
    }

    free2d_short(icetmp);
    return 0;
}

float get_icefrac_oisst(float lon, float lat) {
    float dx = 360.0 / OI4NX;
    float dy = 180.0 / OI4NY;

    float ice = 0;
    int i, j;
    int ntmp;
    float xx, yy;
    float t, u;
    float reftmp[2][2];
    float ftmp;


    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), OI4NX + 1), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), OI4NY + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, replacing missing (land) values with average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    ftmp = 0.0;
    ntmp = 0;
    if (iceref[j ][i ] > BAD_FLT + 1) {
        ftmp += iceref[j ][i ];
        ntmp++;
    }
    if (iceref[j ][i + 1] > BAD_FLT + 1) {
        ftmp += iceref[j ][i + 1];
        ntmp++;
    }
    if (iceref[j + 1][i + 1] > BAD_FLT + 1) {
        ftmp += iceref[j + 1][i + 1];
        ntmp++;
    }
    if (iceref[j + 1][i ] > BAD_FLT + 1) {
        ftmp += iceref[j + 1][i ];
        ntmp++;
    }
    if (ntmp > 0) {
        ftmp /= ntmp;
        reftmp[0][0] = (iceref[j ][i ] > BAD_FLT + 1 ? iceref[j ][i ] : ftmp);
        reftmp[0][1] = (iceref[j ][i + 1] > BAD_FLT + 1 ? iceref[j ][i + 1] : ftmp);
        reftmp[1][1] = (iceref[j + 1][i + 1] > BAD_FLT + 1 ? iceref[j + 1][i + 1] : ftmp);
        reftmp[1][0] = (iceref[j + 1][i ] > BAD_FLT + 1 ? iceref[j + 1][i ] : ftmp);

        ice = (1 - t)*(1 - u) * reftmp[0][0] + t * (1 - u) * reftmp[0][1] + t * u * reftmp[1][1] + (1 - t) * u * reftmp[1][0];

    } else
        ice = 0;

    return (ice);
}

char ice_mask_oisst(float lon, float lat) {
    if (get_icefrac_oisst(lon, lat) > ice_threshold)
        return 1;
    else
        return 0;
}

/*************************************************************
 * init the ice mask functions.  
 *
 * file      - file name of the ice data file
 * year      - year of the ice file
 * day       - day of year
 * threshold - ice fraction above which the ice mask returns 1
 *
 * return 0 if OK or 1 if error
 *
 *************************************************************/
int ice_mask_init(char *file, int year, int day, float threshold) {
    int32 sd_id;
    int32 attr_index;
    int32 status;
    int result = 1;
    char titleStr[256] = "";

    ice_initalized = 0;
    ice_file_type = ICE_TYPE_OLD;

    // set the threshold global variable
    ice_threshold = threshold;

    if (Hishdf(file)) {
        /* Open the file and initiate the SD interface */
        sd_id = SDstart(file, DFACC_RDONLY);
        if (sd_id == -1) {
            fprintf(stderr, "-E- %s line %d: Could not open %s as an HDF file.\n",
                    __FILE__, __LINE__, file);
            return (HDF_FUNCTION_ERROR);
        }
        // look for the global attribute "Title"
        attr_index = SDfindattr(sd_id, "Title");
        if (attr_index == -1) {

            // Title not found, assume it is an old ice HDF climatolgy file
            result = ice_mask_init_old(file, day, sd_id);
            if (result == 0) {
                ice_initalized = 1;
                ice_file_type = ICE_TYPE_OLD;
            }
            SDend(sd_id);
            return (result);
        }

        // Read the attribute data.
        status = SDreadattr(sd_id, attr_index, titleStr);
        if (status == FAIL) {
            fprintf(stderr, "-E- %s line %d: Could not read the global attribute \"Title\" from, %s .\n",
                    __FILE__, __LINE__, file);
            SDend(sd_id);
            return (1);
        }

        // check title to which type of file it is
        if (strcmp(titleStr, "Sea Ice Concentration Daily") == 0) {
            result = ice_mask_init_daily(file, year, day, sd_id);
            ice_file_type = ICE_TYPE_NSIDC;
        } else if (strcmp(titleStr, "Sea Ice Concentration Monthly") == 0) {
            result = ice_mask_init_monthly(file, year, day, sd_id);
            ice_file_type = ICE_TYPE_NSIDC;
        } else {
            fprintf(stderr, "-E- %s line %d: Global attribute \"Title\" not valid from, %s .\n",
                    __FILE__, __LINE__, file);
            result = 1;
        }

        if (result == 0)
            ice_initalized = 1;
        SDend(sd_id);
        return result;
    } // is HDF

    /* try netCDF */
    int ncid;
    if (nc_open(file, NC_NOWRITE, &ncid) == NC_NOERR) {
        nc_close(ncid);
        result = ice_init_oisst(file);
        if (result == 0) {
            ice_initalized = 1;
            ice_file_type = ICE_TYPE_OISST;
        }
        return result;
    }

    // if the file is not an HDF or NetCDF file then try to open it
    // as a raw NSIDC file
    result = ice_mask_init_nsidc_raw(file);
    if (result == 0) {
        ice_initalized = 1;
        ice_file_type = ICE_TYPE_NSIDC;
    }
    return ( result);

}

/**
 * Get the ice mask at this location. 1=ice, 0=water
 * 
 * @param lon longitude of interest
 * @param lat latitude of interest
 * @return 1 if icefrac above ice threshold else return 0
 */
char ice_mask(float lon, float lat) {
    if (!ice_initalized)
        return 0;

    switch (ice_file_type) {
    case ICE_TYPE_OLD:
        return ice_mask_old(lon, lat);
    case ICE_TYPE_NSIDC:
        return ice_mask_nsidc(lon, lat);
    case ICE_TYPE_OISST:
        return ice_mask_oisst(lon, lat);
    default:
        fprintf(stderr, "-E- %s line %d: Ice file type=%d is invalid.\n",
                __FILE__, __LINE__, ice_file_type);
        exit(EXIT_FAILURE);
    }
}

/**
 * function to return the ice fraction
 *
 * @param lon longitude of interest
 * @param lat latitude of interest
 * @return ice fraction coverage at this position (range of 0 to 1)
 */
float ice_fraction(float lon, float lat) {
    if (!ice_initalized)
        return NODATA;

    switch (ice_file_type) {
    case ICE_TYPE_OLD:
        return get_icefrac_old(lon, lat);
    case ICE_TYPE_NSIDC:
        return get_icefrac_nsidc(lon, lat);
    case ICE_TYPE_OISST:
        return get_icefrac_oisst(lon, lat);
    default:
        fprintf(stderr, "-E- %s line %d: Ice file type=%d is invalid.\n",
                __FILE__, __LINE__, ice_file_type);
        exit(EXIT_FAILURE);
    }
}

