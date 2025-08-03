/* =================================================================== */
/* module swim.c - Shallow Water Inversion Model                       */
/*                                                                     */
/* This module contains the functions to optimize and evaluate the     */
/* SWIM Bio-optical model. The SWIM algorithm  was designed to improve */
/* retrievals in optically shallow waters by considering both water    */
/* column depth and benthic albedo                                     */
/*                                                                     */
/* References:                                                         */
/* McKinna et al (2015) A semianalytical ocean color                   */
/* inversion algorithm with explicit water-column depth and substrate  */
/* reflectance parameterization, JGR Oceans, 120,                      */
/* doi: 10.1002/2014JC010224                                           */
/*                                                                     */
/* Reichstetter et al (2014) Seafloor brightness map of the Great      */
/* Barrier Reef, Australia, derived from biodiversity data, PANGEA,    */
/* doi:10.1594/PANGAEA.835979                                          */
/*                                                                     */
/* Implementation:                                                     */
/* L. McKinna and D. Shea, NASA/OBPG/SAIC, March 2015                  */
/*                                                                     */
/* Notes:                                                              */
/* In lieu of user-supplied data:                                      */
/* - Uses ETOP01 as default bathymetry.                                */
/* - Uses single albedo spectrum as default.                           */
/* =================================================================== */

#include <stdio.h>
#include <math.h>
#include <levmar.h>
#include <netcdf.h>

#include "l12_proto.h"

#define NPAR  3           //Number of free model parameters
#define MAXITR 1500          //Maximum number of LM optimisation iterations
static int32_t swimScanNum = -1; //last scan num that the model was calculated for
static int32_t LastRecNum = -1;

static int32_t* lambda; //array holding wavelengths for this sensor
static int32_t nbandVis; //number of visible bands for a given sensor

//output product pointers
static float *adg; // adg value for the scan line adg[pix, band]
static float *bbp; // bbp value for the scan line bbp[pix, band]
static float *aph; // aph value for the scan line aph[pix, band]
static float *atot; // atot value for the scan line atot[pix, band]
static float *bbtot; // bbtot value for the scan line bbtot[pix, band]
static int16 *iter;

//Bio-optical model parameters
static float *aphstar; //spectral shape of phyto abs (aphstar[band])
static float *adgstar; //spectral shape of detritus + cdom abs (adgstar[band])
static float *bbpstar; //spectral shape of particle backscat (bbpstar[band])
static float *aw; //pure water abs spectra (aw[band])
static float *bbw; //pure water backscat spectra (bbw[band])

//static float grd1 = 0.0949;        //the "g0" coefficient in rrs function Lee et al 1998
//static float grd2 = 0.0794;        //the "g1" coefficient in rrs function Lee et al 1998

static const double grd1 = 0.08945; //the "g0" coefficient in rrs function Lee et al 2002
static const double grd2 = 0.1247; //the "g1" coefficient in rrs function Lee et al 2002

static double invCosSolz; // 1.0 / cos(solar_zenith_angle)
static double invCosSenz; // 1.0/ cos(sensor_viewing_angle)
static double depth; // depth at this pixel

// variables for reading the benthic netCDF file
static int ncfile; // netCDF file ID
static int benthicProportionId; // netCDF variable id for benthicProportion
static int benthicROutOfBounds; // binary flag, pixel not within supplied benthic reflectance lat/lon grid
static int benthicRFileExist; // binary flag, benthic reflectance file supplied to l2gen or not.
static size_t numLat; // number of latitudes in the benthicProportion array
static size_t numLon; // number of longitudes in the benthicProportion array
static size_t numBottomTypes; // number of bottom types in the benthicProportion array
static double upperLat; // latitude of benthicProportion[0][0][bottomType]
static double lowerLat; // latitude of benthicProportion[numLat-1][0][bottomType]
static double deltaLat; // latitude grid spacing
static double leftLon; // longitude of benthicProportion[0][0][bottomType]
static double rightLon; // longitude of benthicProportion[0][numLon-1][bottomType]
static double deltaLon; // longitude grid spacing

static double* reflectance; // bottom reflectance for each sensor wavelength (reflectance[bottomType][band])
static double* benthicRefl; // bottom percent for this pixel (benthicRefl[band])

//Hard code in a regional chlorophyll-specific phytoplankton spectra
//for tropical waters (collected in southern GBR).
static float taphw[10] = {412.0, 443.0, 469.0, 488.0, 531.0, 551.0, 555.0,
    645.0, 667.0, 678.0};
static float taphs[10] = {0.840353901, 1.0, 0.821449702, 0.689387045,
    0.217028284, 0.14646583, 0.134941382, 0.122233711, 0.291530253,
    0.356272348};

//Hard code in a default benthic albedo in the case there is no default global albedo file found,
//or data is out of lat/lon bounds.
static float reflw[10] = {412.0, 443.0, 469.0, 488.0, 531.0, 551.0, 555.0,
    645.0, 667.0, 678.0};
static float refls[10] = {0.167663599, 0.197465145, 0.224745351, 0.24005827,
    0.28106442, 0.297352629, 0.30027302, 0.337329977, 0.309064323, 0.310894269};

/* ------------------------------------------------------------------------*/
/* swim_ran() - determine if swim has been run for scanline                */
/* ------------------------------------------------------------------------*/
int swim_ran(int recnum)
{ 
    if ( recnum == LastRecNum )
        return 1;
    else
        return 0; 
}

/* ------------------------------------------------------------------------*/
/* getAttr() - get an attribute from the global ncfile.  Use NC_GLOBAL for */
/* varid to get global attributes                                          */

/* ------------------------------------------------------------------------*/
double getAttr(int varid, char* attrName) {
    double val;
    if (nc_get_att_double(ncfile, varid, attrName, &val) != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: could get attribute %s from netCDF File.\n",
                __FILE__, __LINE__, attrName);
        exit(1);
    }
    return val;
}

/* -----------------------------------------------------------------------*/
/* getVarId() - get variable ID from the global ncfile.                   */

/* -----------------------------------------------------------------------*/
int getVarId(char* name) {
    int varid;
    if (nc_inq_varid(ncfile, name, &varid) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: could not find %s in netCDF File.\n",
                __FILE__, __LINE__, name);
        exit(1);
    }
    return varid;
}

/* -----------------------------------------------------------------------*/
/* getDimensionId() - get the dimension IDs from a variable               */

/* -----------------------------------------------------------------------*/
void getDimensionIds(int varid, int* dimIds) {
    if (nc_inq_vardimid(ncfile, varid, dimIds) != NC_NOERR) {
        char name[NC_MAX_NAME];
        nc_inq_varname(ncfile, varid, name);
        fprintf(stderr,
                "-E- %s line %d: could not find dimension Ids of %s in netCDF File.\n",
                __FILE__, __LINE__, name);
        exit(1);
    }
}

/* -----------------------------------------------------------------------*/
/* getDimensionLength() - get the length of the dimension ID              */

/* -----------------------------------------------------------------------*/
size_t getDimensionLength(int dimId) {
    size_t length;
    if (nc_inq_dimlen(ncfile, dimId, &length) != NC_NOERR) {
        char name[NC_MAX_NAME];
        nc_inq_dim(ncfile, dimId, name, &length);
        fprintf(stderr,
                "-E- %s line %d: could not get size of dimension \"%s\" in netCDF File.\n",
                __FILE__, __LINE__, name);
        exit(1);
    }
    return length;
}

/* -----------------------------------------------------------------------*/
/* initBenthicFile() - init the global variables from the bottom          */
/* reflectance netCDF file                                                */

/* -----------------------------------------------------------------------*/
void initBenthicFile(char* fileName) {
    int varid;
    int dimIds[3];
    size_t bottom;
    int band;
    size_t numWavelengths; // number of wavelengths in reflectance array
    double firstWavelength; // wavelength of reflectance[0]
    double deltaWavelength; // delta between each wavelength of reflectance
    double reflectanceScale; // reflectance scale factor
    double reflectanceOffset; // reflectance offset

    //Was a benthic reflectance filename supplied?
    if (strlen(fileName) != 0) {
        printf("\n");
        printf("Reading benthic reflectance from: %s \n", fileName);
        printf("\n");
        benthicRFileExist = 1;
    } else {
        printf("\n");
        printf("-E- No benthic reflectance file supplied: use default benthic albedo spectra \n");
        printf("\n");
        benthicRFileExist = 0;
        //Return out of function
        return;
    }

    //If file exists and cannot be read, exit and print error
    if (benthicRFileExist) {
        if (nc_open(fileName, NC_NOWRITE, &ncfile) != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d: could not open netCDF File \"%s\".\n",
                    __FILE__, __LINE__, fileName);
            exit(1);
        }
    }

    lowerLat = getAttr(NC_GLOBAL, "lower_lat");
    upperLat = getAttr(NC_GLOBAL, "upper_lat");
    leftLon = getAttr(NC_GLOBAL, "left_lon");
    rightLon = getAttr(NC_GLOBAL, "right_lon");

    // normalize the lons
    while (leftLon > 180.0)
        leftLon -= 360.0;
    while (leftLon < -180.0)
        leftLon += 360.0;
    while (rightLon > 180.0)
        rightLon -= 360.0;
    while (rightLon < -180.0)
        rightLon += 360.0;
    while (rightLon <= leftLon)
        rightLon += 360.0;

    firstWavelength = getAttr(NC_GLOBAL, "firstWavelength");
    deltaWavelength = getAttr(NC_GLOBAL, "deltaWavelength");

    // make sure deltaWavelength is valid
    if (deltaWavelength <= 0) {
        fprintf(stderr,
                "-E- %s line %d: deltaWavelength must be > 0 in netCDF File \"%s\".\n",
                __FILE__, __LINE__, fileName);
        exit(1);
    }

    varid = getVarId("reflectance");

    // get reflectance scale and offset
    reflectanceScale = 1.0 / getAttr(varid, "scale_factor");
    reflectanceOffset = getAttr(varid, "add_offset");

    getDimensionIds(varid, dimIds);
    numBottomTypes = getDimensionLength(dimIds[0]);
    numWavelengths = getDimensionLength(dimIds[1]);

    // allocate global variable
    reflectance = (double*) allocateMemory(
            numBottomTypes * nbandVis * sizeof (double), "reflectance");

    size_t start[2];
    size_t count[2];
    count[0] = 1;

    // load up reflectance for each bottom type for each sensor wavelength
    for (bottom = 0; bottom < numBottomTypes; bottom++) {
        start[0] = bottom;
        for (band = 0; band < nbandVis; band++) {

            // average an 11nm band width
            if (lambda[band] < firstWavelength) {
                reflectance[bottom * nbandVis + band] = 0.0;
                continue;
            }
            int firstIndex = round(
                    (lambda[band] - 5 - firstWavelength) / deltaWavelength);
            if (firstIndex < 0)
                firstIndex = 0;
            if (firstIndex >= numWavelengths) {
                reflectance[bottom * nbandVis + band] = 0.0;
                continue;
            }
            int numIndex = round(11.0 / deltaWavelength);
            if (firstIndex + numIndex >= numWavelengths) {
                numIndex = numWavelengths - firstIndex;
            }
            if (numIndex < 1)
                numIndex = 1;

            start[1] = firstIndex;
            count[1] = numIndex;

            ushort tmpShort[numIndex];
            if (nc_get_vara_ushort(ncfile, varid, start, count,
                    tmpShort) != NC_NOERR) {
                fprintf(stderr,
                        "-E- %s line %d: could not read values of reflectance in netCDF File \"%s\".\n",
                        __FILE__, __LINE__, fileName);
                exit(1);
            }

            // average the 11nm wide readings and scale
            int i;
            double sum = 0.0;
            for (i = 0; i < numIndex; i++)
                sum += tmpShort[i];
            reflectance[bottom * nbandVis + band] = sum / numIndex
                    * reflectanceScale + reflectanceOffset;

        } // for bands
    } // for bottoms

    // get dimensions of the benthicProportion variable
    benthicProportionId = getVarId("benthicProportion");

    getDimensionIds(benthicProportionId, dimIds);
    numLat = getDimensionLength(dimIds[0]);
    numLon = getDimensionLength(dimIds[1]);
    size_t tmpSize = getDimensionLength(dimIds[2]);

    if (tmpSize != numBottomTypes) {
        fprintf(stderr,
                "-E- %s line %d: numBottomTypes in benthicProportion is not the same as in reflectance in netCDF File \"%s\".\n",
                __FILE__, __LINE__, fileName);
        exit(1);
    }

    deltaLat = (upperLat - lowerLat) / numLat;
    deltaLon = (rightLon - leftLon) / numLon;

}

/* -----------------------------------------------------------------------*/
/* getDefaultBenthicR() - get default benthic reflectance spectra and     */
/* interpolate to sensor wavelengths                                      */

/* -----------------------------------------------------------------------*/
void getDefaultBenthicR() {

    static int32_t refln = -1;
    int band;

    /* read default benthic reflectance table size  */
    refln = 0;
    for (refln = 0; refln < nbandVis; refln++) {
        if (reflw[refln] < 0) {
            break;
        }
    }

    /*Interploate the hard-coded benthic albedo to sensor wavelenghth*/
    for (band = 0; band < nbandVis; band++) {
        benthicRefl[band] = linterp(reflw, refls, refln, lambda[band]);
    }

}

/* -----------------------------------------------------------------------*/
/* getBottomReflectance() - get the benthic reflectance for this pixel    */
/*  into benthicRefl global array                                         */

/* -----------------------------------------------------------------------*/
void getBottomReflectance(float lat, float lon) {
    int latIndex;
    int lonIndex;
    int latFlag = 0;
    int lonFlag = 0;

    int band;
    int bottom;
    double sum;

    // see if coordinate is covered by the file
    if (lat < lowerLat || lat > upperLat) {
        latFlag = 1;
    }

    // see if coordinate is covered by the file
    if (lon < leftLon || lon > rightLon) {
        lonFlag = 1;
    }

    //If the pixel is outside the lat/lon range of the user input benthic albedo file,
    // a default hard-coded sand benthic albedo is used.
    if (latFlag || lonFlag) {
        //Set benthic out of bounds flag to 1;
        benthicROutOfBounds = 1;

        //Call function to interpolate default benthic reflectance
        getDefaultBenthicR();

        //Break out of getBottomReflectance function
        return;

    } else {
        /*Else the pixel is within the user-supplied benthic albedo file*/
        //Set benthic out of bounds flag to 0;
        benthicROutOfBounds = 0;
    }

    //printf("benthic out of bounds  = %d \n",benthicROutOfBounds);
    latIndex = (lat - lowerLat) / deltaLat;
    lonIndex = (lon - leftLon) / deltaLon;

    size_t start[3];
    start[0] = latIndex;
    start[1] = lonIndex;
    start[2] = 0;
    size_t count[3];
    count[0] = 1;
    count[1] = 1;
    count[2] = numBottomTypes;

    uint8_t benthicPercent[numBottomTypes];
    if (nc_get_vara_uchar(ncfile, benthicProportionId, start, count,
            benthicPercent) != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: could not read values from benthicProportion in netCDF File.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    for (band = 0; band < nbandVis; band++) {
        sum = 0.0;
        for (bottom = 0; bottom < numBottomTypes; bottom++) {
            sum += reflectance[bottom * nbandVis + band] * benthicPercent[bottom]
                    / 100.0;
        }

        benthicRefl[band] = sum;
    }
}

/* ------------------------------------------------------------------- */
/* swim_func() - optically shallow semianlytical reflectance model     */

/* ------------------------------------------------------------------- */
void swim_func(double *initialParams, double *rrsTotal, int numParams,
        int numBands, void* dataPtr) {

    int iw; //iterator
    double rrsDeep; //deep water rrs term
    double rrsColumn; //water column rrs term
    double rrsBenthos; //bottom contribution rrs term
    double bb; //total backscattering coefficient
    double ac; //total absorption coefficient
    double kappa; //kappa = ac + bb
    double u; //u = bb_total/(kappa)
    double duC; // Upward attenuation term
    double duB; // Downward attenuation term

    //-------------------------------------------//
    /* Shallow water sub-surface reflectance model*/
    //-------------------------------------------//
    // Iterate over wavelengths
    for (iw = 0; iw < numBands; iw++) {

        //Calculate the bulk IOPs
        ac = aw[iw] + initialParams[0] * aphstar[iw]
                + initialParams[1] * adgstar[iw];
        bb = bbw[iw] + initialParams[2] * bbpstar[iw];


        //Attenuation coefficients
        kappa = bb + ac;
        u = bb / kappa;

        //Pathlength elongation factors
        duC = 1.03 * pow((1 + 2.48 * u), 0.5);
        duB = 1.04 * pow((1.54 + 5.4 * u), 0.5);

        //Deep water reflectance
        rrsDeep = grd1 * u + grd2 * pow(u, 2);

        //Water column term
        rrsColumn = rrsDeep
                * (1.0
                + exp(
                -1.0 * kappa * depth
                * (invCosSolz + (duC * invCosSenz))));

        //Bottom contribution term
        rrsBenthos = (benthicRefl[iw] / OEL_PI)
                * exp(-1.0 * kappa * depth * (invCosSolz + (duB * invCosSenz)));

        rrsTotal[iw] = rrsColumn + rrsBenthos;
        
    }
}

/* ---------------------------------------------------------------------- */
/* run_swim() - calculates the SWIM model for full scanline               */

/* ---------------------------------------------------------------------- */
void run_swim(l2str *l2rec) {

    //Declare variables
    static int32_t taphn = -1;
    //static float *taphw = NULL;
    //static float *taphs = NULL;
    static float adg_s = -1;
    static float bbp_s = -1;

    //Iterators indexes
    int i, iw;
    int validr;
    int statusLM;
    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    double Rrs; //above-water observed remote sensing reflectance
    double rrs_s[l1rec->l1file->nbands]; //Array holding sub-surface observed remote sensing reflectance
    //    double rrs_a[l1rec->l1file->nbands];  //Array holding above-water observed remote sensing reflectance
    double fitParams[NPAR]; //Array for holding initial guesses for aph440, adg440, bbp550
    double opts[LM_OPTS_SZ]; //Atray for holding Levmar optionsfloat get_benthic_r(float wave, float chl)
    double info[LM_INFO_SZ]; //Array for holding Levmar info

    float aph443; //the value of phyto abs at 443 nm
    float bbp443; //the value of detritus + cdom abs at 443 nm
    float adg443; //the value of particle backscatter 443 nm

    int defaultRefl; //binary flag, was default reflectance used?

    // first check to see if this is the first time we have run
    // so we can allocate memory and init stuff
    if (swimScanNum == -1) {

        /* limit to visible wavelengths */
        nbandVis = rdsensorinfo(l1file->sensorID, l1_input->evalmask,
                "NbandsVIS", NULL);

        /*Get the wavelengths - needed for IOP spectral models*/
        rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Lambda",
                (void**) &lambda);

        /*Allocate memory*/
        adg = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "adg");
        bbp = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "bbp");
        aph = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "aph");
        atot = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "atot");
        bbtot = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "bbtot");
        iter = (int16*) allocateMemory(l1rec->npix * nbandVis * sizeof (int16),
                "iter");
        aphstar = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "aphstar");
        adgstar = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "adgstar");
        bbpstar = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "bbpstar");
        aw = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "aw");
        bbw = (float*) allocateMemory(l1rec->npix * nbandVis * sizeof (float),
                "bbw");
        benthicRefl = (double*) allocateMemory(
                l1rec->npix * nbandVis * sizeof (double), "benthicRefl");

        /*--------------------------------------------------------*/
        /*IOP models ar defined once here and allocated to memory */
        /*--------------------------------------------------------*/

        /*Need to load up the model coefficient shapes (aph*, adg*, bbp*)*/
        /* set adg slope and bbp power law static coefficients */
        adg_s = 0.017;
        bbp_s = 1.0;

        /* read aphstar table size  */
        taphn = 0;
        for (taphn = 0; taphn < nbandVis; taphn++)
            if (taphw[taphn] < 0)
                break;

        /*Loop over wavelengths */
        for (iw = 0; iw < nbandVis; iw++) {

            /*Phytoplankton absorption shape*/
            aphstar[iw] = linterp(taphw, taphs, taphn, lambda[iw]);

            //printf("aphi = %f @ lambda = %d \n", aphstar[iw], lambda[iw] );

            /*Coloured dissolved and detrital matter absorption shape*/
            adgstar[iw] = exp(-adg_s * (lambda[iw] - 443.0));

            /*Particulate matter backscatter shape*/
            bbpstar[iw] = pow((443.0 / lambda[iw]), bbp_s);

            /*Pure water absorption coefficient*/
            aw[iw] = aw_spectra(lambda[iw], BANDW);

            /*Pure water backscatter coefficient*/
            bbw[iw] = bbw_spectra(lambda[iw], BANDW);
        }

        /*Begin reading benthic reflectance data*/
        initBenthicFile(input->breflectfile);
    }

    /*record that we have run the model on this scan line*/
    swimScanNum = l1rec->iscan;

    /*----------------------------------------*/
    /*  Iterate through pixels            */
    /*----------------------------------------*/
    for (i = 0; i < l1rec->npix; i++) {

        // init the spectral outputs
        for (iw = 0; iw < nbandVis; iw++) {
            adg[i * nbandVis + iw] = BAD_FLT;
            bbp[i * nbandVis + iw] = BAD_FLT;
            aph[i * nbandVis + iw] = BAD_FLT;
            atot[i * nbandVis + iw] = BAD_FLT;
            bbtot[i * nbandVis + iw] = BAD_FLT;
            iter[i * nbandVis + iw] = 0;
        }

        //Set status to default value of 0.
        validr = 0;

        //Use the quality control from GSM to flag bad pixelsparamFit
        if (l1rec->mask[i] == 0) {

            /* Solar and Sensor Viewing Geometry */
            double solzRad = (l1rec->solz[i] * OEL_PI) / 180.0; /*Convert to radians*/
            double senzRad = (l1rec->senz[i] * OEL_PI) / 180.0; /*convert to radians*/
            invCosSolz = 1.0 / cos(solzRad); /*Check rad or deg */
            invCosSenz = 1.0 / cos(senzRad); /*Check rad or deg */

            depth = 0 - l1rec->dem[i];

            for (iw = 0; iw < nbandVis; iw++) {
                Rrs = l2rec->Rrs[i * l1file->nbands + iw];
                //                                rrs_a[iw] = Rrs;
                //Sub-surface remote sensing reflectance
                rrs_s[iw] = Rrs / (0.52 + 1.7 * Rrs);
                if (rrs_s[iw] >= 0.0)
                    validr = validr + 1;
            }

            /* if less than 3 valid radiances, fail pixel    */
            if (validr <= 3) {
                l1rec->flags[i] |= PRODFAIL;
                continue;
            }

            /*---------------------------------------------------*/
            /* fill up the benthicRefl global array */
            defaultRefl = -1;
            if (!benthicRFileExist) {

                if (defaultRefl < 0) {
                    //call function to interpolate default reflectance
                    getDefaultBenthicR();

                    //Assign this flag so the default reflectance does not
                    //need to be calculated multiple times.
                    defaultRefl = 1;
                }

            } else {
                getBottomReflectance(l1rec->lat[i], l1rec->lon[i]);
            }


            /*If the hard-coded reflectance is used, set l2flag to PRODWARN
            if ( benthicROutOfBounds ) {
                    l2rec->flags[i] |= PRODWARN;
            };*/

            /*---------------------------------------------------*/
            /* fit model for this pixel using levmar optimisation*/
            /*---------------------------------------------------*/

            /*Initial guess (LM starting points)*/
            fitParams[0] = 0.2;
            fitParams[1] = 0.01;
            fitParams[2] = 0.001;

            /*Optimisation control parameters*/
            opts[0] = 1E-3; //1E-3 LM_INIT_MU; /* scale factor for initial mu */
            opts[1] = 1E-10; //1E-10 LM_INIT_MU;  convergence in gradient
            opts[2] = 1E-5; //1E-5 relative error desired in the approximate solution
            opts[3] = 1E-7; //1E-7 relative error desired in the sum of squares
            opts[4] = 1E-6; //1E-6  LM_DIFF_DELTA;  /*step used in difference approximation to the Jacobian*/

            /*Set box constraints*/
            /*L/U Bounds inferred from Werdell et al. 2013*/
            double lowerBounds[NPAR] = {-0.0003755, -0.0003755, -0.0001195625}; //Lower bounds for three parameters
            double upperBounds[NPAR] = {5.0, 5.0, 5.0}; //Upper bounds for three parameters
            int m = NPAR;
            int n = nbandVis;
            int maxit = MAXITR;

            /*Run optimization*/
            /*levmar routine (dlevmar), box contstraints (bc), numerical partials (diff)*/
            //For further details, see source code in: lmbc_core.c //
            statusLM = dlevmar_bc_dif(swim_func, fitParams, rrs_s, m, n,
                    lowerBounds, upperBounds, NULL, maxit, opts, info, NULL,
                    NULL, NULL);

            /*---------------------------------------------------*/
            /*             Write products to arrays              */
            /*---------------------------------------------------*/

            /* save results and flag failure */
            //printf("num iteration %f \n",info[5]);
            //if (statusLM <= 0 || info[5] == maxit) {
            if (statusLM <= 0) {
                l1rec->flags[i] |= PRODFAIL;

            } else {
                /*Output results and calculate IOPs*/
                aph443 = fitParams[0]; //absorption of phytoplankton at 443 nm band
                adg443 = fitParams[1]; //absorption of detritus + cdom at 443 nm band
                bbp443 = fitParams[2]; //particle backscattering at 443 nm band

                /*calculate spectral IOPs for output*/
                for (iw = 0; iw < nbandVis; iw++) {

                    //absorption of phytoplankton
                    //Note: have normalised aphi440=1.0 - fix later with regional aphi
                    //Free parameter in model thus aph443 rather than Chl
                    aph[i * nbandVis + iw] = aph443 * aphstar[iw];
                    adg[i * nbandVis + iw] = adg443 * adgstar[iw]; //absorption detritus + cdom
                    bbp[i * nbandVis + iw] = bbp443 * bbpstar[iw]; //backscattering of particles
                    atot[i * nbandVis + iw] = aw[iw] + aph443 * aphstar[iw]
                            + adg443 * adgstar[iw]; //total absorption
                    bbtot[i * nbandVis + iw] = bbw[iw] + bbp443 * bbpstar[iw]; //total backscattering

                }
            }

        }
    }
    
    LastRecNum = l1rec->iscan;

}

/* -------------------------------------------------------------------- */
/* get_swim() - returns requested swim product for full scanline        */

/* ---------------------------------------------------------------------*/
void get_swim(l2str *l2rec, l2prodstr *p, float prod[]) {
    
    int32_t ip;
    l1str *l1rec = l2rec->l1rec;
    
    if (!swim_ran(l1rec->iscan)) {
        run_swim(l2rec);
    }

    int bandNum = p->prod_ix;
    if (bandNum < 0) {
        printf("-E- %s line %d : prod_ix must be positive, prod_ix=%d\n",
                __FILE__, __LINE__, p->prod_ix);
        exit(1);
    }
    if (bandNum >= nbandVis) {
        printf(
                "-E- %s line %d : prod_ix must be less than numBands=%d, prod_ix=%d\n",
                __FILE__, __LINE__, nbandVis, p->prod_ix);
        exit(1);
    }

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {

        switch (p->cat_ix) {

        case CAT_a_swim:
            prod[ip] = atot[ip * nbandVis + bandNum];
            break;

        case CAT_bb_swim:
            prod[ip] = bbtot[ip * nbandVis + bandNum];
            break;

        case CAT_adg_swim:
            prod[ip] = adg[ip * nbandVis + bandNum];
            break;

        case CAT_aph_swim:
            prod[ip] = aph[ip * nbandVis + bandNum];
            break;

        case CAT_bbp_swim:
            prod[ip] = bbp[ip * nbandVis + bandNum];
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to swim.\n",
                    __FILE__, __LINE__, p->cat_ix);
            exit(1);
        } 
    } 
}

/* ------------------------------------------------------------------- */
/* interface to convl12() to return SWIM bulk iops */

/* ------------------------------------------------------------------- */
void iops_swim(l2str *l2rec) {
    
    int ib, ip;
    l1str *l1rec = l2rec->l1rec;
    
    if (!swim_ran(l1rec->iscan)) {
        run_swim(l2rec);
    }

    for (ip = 0; ip < l1rec->npix; ip++) {
        for (ib = 0; ib < nbandVis; ib++) {
            l2rec->a[ip * l1rec->l1file->nbands + ib] = atot[ip * nbandVis + ib];
            l2rec->bb[ip * l1rec->l1file->nbands + ib] = bbtot[ip * nbandVis + ib];

        }
    }

    return;
}

/* ------------------------------------------------------------------- */
