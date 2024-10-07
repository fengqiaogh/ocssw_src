#include "l3mapgen.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <genutils.h>
#include <clo.h>
#include <dfutils.h>
#include <sensorInfo.h>
#include <unistd.h>

/** add all of the accepted command line options to list */
int l3mapgen_init_options(clo_optionList_t* list, const char* softwareVersion) {

    clo_setVersion2("l3mapgen", softwareVersion);

    clo_setHelpStr("Usage: l3mapgen argument-list"
            "\n"
            "\n  This program takes a product (or products if netCDF output) from an L3 bin"
            "\n  or SMI file, reprojects the data using proj.4 and writes a mapped file in"
            "\n  the requested output format."
            "\n"
            "\n  Return values"
            "\n    0 = All Good"
            "\n    1 = Error"
            "\n    110 = No valid data to map"
            "\n"
            "\n  The argument list is a set of keyword=value pairs.  Arguments can"
            "\n  be specified on the command line, or put into a parameter file, or the"
            "\n  two methods can be used together, with command line overriding."
            "\n"
            "\nThe list of valid keywords follows:"
            "\n");

    clo_addOption(list, "suite", CLO_TYPE_STRING, NULL, "suite for default parameters");

    // files
    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L3 bin filename"
            "\n");

    clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output filename");
    clo_addOption(list, "oformat", CLO_TYPE_STRING, "netcdf4",
            "output file format"
            "\n        netcdf4: netCDF4 file, can contain more than one product"
            "\n        hdf4:    HDF4 file (old SMI format)"
            "\n        png:     PNG image file"
            "\n        ppm:     PPM image file"
            "\n        tiff:    TIFF file with georeference tags");

    clo_addOption(list, "ofile_product_tag", CLO_TYPE_STRING, "PRODUCT",
            "sub-string in ofile name that will be substituted by the product name");
    clo_addOption(list, "ofile2", CLO_TYPE_OFILE, NULL,
            "second output filename");
    clo_addOption(list, "oformat2", CLO_TYPE_STRING, "png",
            "second output file format"
            "\n        same options as oformat");

    clo_addOption(list, "deflate", CLO_TYPE_INT, "4", "netCDF4 deflation level"
            "\n");

    // data, projection, coverage
    clo_addOption(list, "product", CLO_TYPE_STRING, NULL,
            "comma separated list of products."
            "\n        Each product can have an optional colon and modifier appended."
            "\n        For example, \"product=chlor_a,chlor_a:stdev,Kd_490:nobs\""
            "\n        Available modifiers:"
            "\n            avg       average value (default)"
            "\n            stdev     standard deviation"
            "\n            var       variance"
            "\n            nobs      number of observations in the bin"
            "\n            nscenes   number of contributing scenes"
            "\n            obs_time  average observation time (TAI93)"
            "\n            bin_num   bin ID number"
            "\n");
    clo_addOption(list,"wavelength_3d",CLO_TYPE_STRING,NULL, "comma separated list of wavelengths for 3D products");
    clo_addOption(list, "resolution", CLO_TYPE_STRING, NULL,
            "size of output pixel (default from input file)"
            "\n        in meters or SMI dimensions"
            "\n        90km: 432 x 216 image for full globe"
            "\n        36km: 1080 x 540"
            "\n        18km: 2160 x 1080"
            "\n         9km: 4320 x 2160"
            "\n         4km: 8640 x 4320"
            "\n         2km: 17280 x 8640"
            "\n         1km: 34560 x 17280"
            "\n         hkm: 69120 x 34560"
            "\n         qkm: 138240 x 69120"
            "\n         smi: 4096 x 2048"
            "\n        smi4: 8192 x 4096"
            "\n        land: 8640 x 4320"
            "\n         #.#:  width of a pixel in meters"
            "\n       #.#km:  width of a pixel in kilometers"
            "\n      #.#deg:  width of a pixel in degrees");
    clo_addOption(list, "width", CLO_TYPE_INT, NULL,
            "width of output image in pixels; supercedes resolution parameter.\n");

    clo_addOption(list, "projection", CLO_TYPE_STRING, "platecarree",
            "proj.4 projection string or one"
            "\n        of the following predefined projections:"
            "\n        smi:       Standard Mapped image, cylindrical projection,"
            "\n                   uses central_meridian.  NSEW defaults to whole globe."
            "\n                   projection=\"+proj=eqc +lon_0=<central_meridian>\""
            "\n        platecarree: Plate Carree image, cylindrical projection,"
            "\n                   uses central_meridian."
            "\n                   projection=\"+proj=eqc +lon_0=<central_meridian>\""
            "\n        mollweide: Mollweide projection"
            "\n                   projection=\"+proj=moll +lon_0=<central_meridian>\""
            "\n        lambert:   Lambert conformal conic (2SP) projection"
            "\n                   projection=\"+proj=lcc +lon_0=<central_meridian>"
            "\n                                +lat_0=<scene center latitude>"
            "\n                                +lat_1=<scene south latitude>"
            "\n                                +lat_2=<scene north latitude>\""
            "\n        albersconic: Albers Equal Area Conic projection"
            "\n                   projection=\"+proj=aea +lon_0=<central_meridian>"
            "\n                                +lat_0=<scene center latitude>"
            "\n                                +lat_1=<scene south latitude>"
            "\n                                +lat_2=<scene north latitude>\""
            "\n        aeqd:      Azimuthal Equidistant projection"
            "\n                   projection=\"+proj=aeqd +lon_0=<central_meridian>"
            "\n                                +lat_0=<scene center latitude>\""
            "\n        mercator:  Mercator cylindrical map projection"
            "\n                   projection=\"+proj=merc +lon_0=<central_meridian>\""
            "\n        transmerc:  Transverse Mercator cylindrical map projection"
            "\n                   projection=\"+proj=tmerc +lon_0=<central_meridian>"
            "\n                                +lat_0=<scene center latitude>\""
            "\n        utm:  Universal Transverse Mercator cylindrical map projection"
            "\n                   projection=\"+proj=utm +zone=<utm_zone> [+south]\""
            "\n        obliquemerc:  Oblique Mercator cylindrical map projection"
            "\n                   projection=\"+proj=omerc +gamma=0 +lat_0=<lat_0>"
            "\n                          +lonc=<central_meridian> +alpha=<azimuth>"
            "\n                          +k_0=1 +x_0=0 +y_0=0\""
            "\n        ease2:     EASE-Grid 2.0 projection"
            "\n                   projection=\"EPSG:6933\""
            "\n        ortho:     Orthographic projection"
            "\n                   projection=\"+proj=ortho +lat_0=<lat_0> +lon_0=<central_meridian>"
            "\n                          +ellps=GRS80 +units=m\""
            "\n        stere:     Stereographic projection"
            "\n                   projection=\"+proj=stere +lat_0=<lat_0> +lat_ts=<lat_ts>"
            "\n                          +lon_0=<central_meridian>"
            "\n                          +ellps=WGS84 +datum=WGS84 +units=m\""
            "\n        conus:     USA Contiguous Albers Equal Area Conic USGS version"
            "\n                   projection=\"+proj=aea +lat_1=29.5 +lat_2=45.5"
            "\n                         +lat_0=23.0 +lon_0=-96 +x_0=0 +y_0=0"
            "\n                         +ellps=GRS80 +datum=NAD83 +units=m\""
            "\n        alaska:    Alaskan Albers Equal Area Conic USGS version"
            "\n                   projection=\"EPSG:3338\""
            "\n        gibs:      latitudinally dependent projection"
            "\n                   Plate Carree between 60S and 60N"
            "\n                   else use Polar Sterographic"
            "\n                   North Polar: projection=\"EPSG:3413\""
            "\n                   South Polar: projection=\"EPSG:3031\""
            "\n         raw:       Raw dump of bin file contents."
            "\n");
    clo_addOption(list, "write_projtext", CLO_TYPE_BOOL, "no", "write projection information to a text file.");
    clo_addOption(list, "central_meridian", CLO_TYPE_FLOAT, "-999",
            "central meridian for projection in deg east."
            "\n        Used only for raw dump and predefined projections as above.");
    clo_addOption(list, "lat_ts", CLO_TYPE_FLOAT, NULL,
            "latitude of true scale for projection in deg north."
            "\n        Used only for predefined projections above as required.");
    clo_addOption(list, "lat_0", CLO_TYPE_FLOAT, NULL,
            "latitude of origin for projection in deg north."
            "\n        Used only for predefined projections above as required.");
    clo_addOption(list, "lat_1", CLO_TYPE_FLOAT, NULL,
            "latitude of first standard parallel (south)."
            "\n        Used only for predefined projections above as required.");
    clo_addOption(list, "lat_2", CLO_TYPE_FLOAT, NULL,
            "latitude of second standard parallel (north)."
            "\n        Used only for predefined projections above as required.");
    clo_addOption(list, "azimuth", CLO_TYPE_FLOAT, NULL,
            "projection rotation angle in deg north."
            "\n        Used only for predefined projections above as required.");
    clo_addOption(list, "utm_zone", CLO_TYPE_STRING, NULL,
            "UTM zone number."
            "\n        Used only for the UTM projection;"
            "\n        Append 'S' for southern hemisphere zones (e.g. 59S).");
    clo_addOption(list, "north", CLO_TYPE_FLOAT, "-999",
            "Northernmost Latitude (default: file north)");
    clo_addOption(list, "south", CLO_TYPE_FLOAT, "-999",
            "Southernmost Latitude (default: file south)");
    clo_addOption(list, "east", CLO_TYPE_FLOAT, "-999",
            "Easternmost Longitude (default: file east)");
    clo_addOption(list, "west", CLO_TYPE_FLOAT, "-999",
            "Westernmost Longitude (default: file west)");
    clo_addOption(list, "trimNSEW", CLO_TYPE_BOOL, "no",
            "should we trim output"
            "\n        to match input NSEW range"
            "\n");

    clo_addOption(list, "interp", CLO_TYPE_STRING, "nearest",
            "interpolation method:"
            "\n        nearest: use the value of the nearest bin for the pixel"
            "\n        bin:     bin all of the pixels that intersect the area of the"
            "\n                  output pixel"
            "\n        area:    bin weighted by area of all the pixels that intersect"
            "\n                  the area of the output pixel"
            "\n");

    // color table, scaling
    clo_addOption(list, "apply_pal", CLO_TYPE_BOOL, "yes",
            "apply color palette:"
            "\n        yes: color image"
            "\n         no: grayscale image");
    clo_addOption(list, "palfile", CLO_TYPE_IFILE, NULL,
            "palette filename (default from product.xml)");
    clo_addOption(list, "use_transparency", CLO_TYPE_BOOL, "no",
            "make missing data transparent (only valid for color PNG and TIFF)");
    clo_addOption(list, "datamin", CLO_TYPE_FLOAT, NULL,
            "minimum value for scaling (default from product.xml)");
    clo_addOption(list, "datamax", CLO_TYPE_FLOAT, NULL,
            "maximum value for scaling (default from product.xml)");
    clo_addOption(list, "scale_type", CLO_TYPE_STRING, NULL,
            "data scaling type (default from product.xml)"
            "\n        linear:  linear scaling"
            "\n        log:     logarithmic scaling"
            "\n        arctan:  arc tangent scaling"
            "\n");

    // other options
    clo_addOption(list, "quiet", CLO_TYPE_BOOL, "false",
            "stop the status printing");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, "Unspecified",
            "processing version string");
    clo_addOption(list, "use_quality", CLO_TYPE_BOOL, "yes",
            "should we do quality factor processing");
    clo_addOption(list, "quality_product", CLO_TYPE_STRING, NULL,
            "product to use for quality factor processing");
    clo_addOption(list, "use_rgb", CLO_TYPE_BOOL, "no",
            "should we use product_rgb to make a"
            "\n        pseudo-true color image");
    clo_addOption(list, "product_rgb", CLO_TYPE_STRING,
            "rhos_670,rhos_555,rhos_412",
            "\n        Three products to use for RGB.  Default is sensor-specific.");
    clo_addOption(list, "fudge", CLO_TYPE_FLOAT, "1.0",
            "fudge factor used to modify size of L3 pixels");
    clo_addOption(list, "threshold", CLO_TYPE_FLOAT, "0",
            "minimum percentage of filled pixels before"
            "\n        an image is generated");
    clo_addOption(list, "num_cache", CLO_TYPE_INT, "500",
            "number of rows to cache in memory.");
    clo_addOption(list, "mask_land", CLO_TYPE_BOOL, "no",
            "set land pixels to pixel value 254");
    clo_addOption(list, "rgb_land", CLO_TYPE_STRING, "160,82,45",
            "RGB value to use for land mask; comma separate string");
    clo_addOption(list, "land", CLO_TYPE_IFILE,
                  "$OCDATAROOT/common/landmask_GMT15ARC.nc", "land mask file");
    clo_addOption(list, "full_latlon", CLO_TYPE_BOOL, "yes",
                  "write full latitude and longitude arrays."
            "\n");
    clo_addOption(list, "doi", CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");

    return 0;
}

int getSensorId(const char* fileName) {
    idDS dsId;
    int sensorId = -1;

    dsId = openDS(fileName);
    if (dsId.fid == -1) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.  \n",
                __FILE__, fileName);
        exit(EXIT_FAILURE);
    }

    if (dsId.fftype == DS_NCDF) {
        char* instrumentStr = readAttrStr(dsId, "instrument");
        if (instrumentStr) {
            char* platformStr = readAttrStr(dsId, "platform");
            if (platformStr) {
                sensorId = instrumentPlatform2SensorId(instrumentStr,
                        platformStr);
            }
        } else {
            char* sensorStr = readAttrStr(dsId, "Sensor");
            if (sensorStr) {
                sensorId = sensorName2SensorId(sensorStr);
            }
        }
    } else {
        char* sensorNameStr = readAttrStr(dsId, "Sensor Name");
        if (sensorNameStr) {
            sensorId = sensorName2SensorId(sensorNameStr);
        }
    }

    if (sensorId == -1) {
        printf("Did not find a valid sensor ID - using OCRVC as the sensor ID.\n");
        sensorId = OCRVC;
    }

    endDS(dsId);
    return sensorId;
}

/*
 Read the command line option and all of the default parameter files.

 This is the order for loading the options:
 - load the main program defaults file
 - load the command line (including specified par files)
 - re-load the command line disabling file descending so command
 line arguments will over ride

 */
int l3mapgen_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.  \n");
        return (-1);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);

    // load program defaults
    sprintf(tmpStr, "%s/common/l3mapgen_defaults.par", dataRoot);
    clo_readFile(list, tmpStr);

    // read all arguments
    clo_readArgs(list, argc, argv);

    // get sensor directory
    int sensorId = getSensorId(clo_getString(list, "ifile"));
    int subsensorId = sensorId2SubsensorId(sensorId);
    const char* sensorDir = sensorId2SensorDir(sensorId);

    // load the sensor specific defaults file
    sprintf(tmpStr, "%s/%s/l3mapgen_defaults.par", dataRoot, sensorDir);
    int defaultLoaded = 0;
    if (access(tmpStr, R_OK) != -1) {
        if (want_verbose)
            printf("Loading default parameters from %s\n", tmpStr);
        clo_readFile(list, tmpStr);
        defaultLoaded = 1;
    }

    if (subsensorId != -1) {
        sprintf(tmpStr, "%s/%s/%s/l3mapgen_defaults.par", dataRoot, sensorDir,
                subsensorId2SubsensorDir(subsensorId));
        if (access(tmpStr, R_OK) != -1) {
            if (want_verbose)
                printf("Loading default parameters from %s\n", tmpStr);
            clo_readFile(list, tmpStr);
            defaultLoaded = 1;
        }
    }

    if (!defaultLoaded) {
        printf("-E- Failed to load sensor program defaults for %s \n",
                sensorId2SensorName(sensorId));
        exit(EXIT_FAILURE);
    }

    // load the suite specific defaults file
    clo_option_t *option = clo_findOption(list, "suite");
    if (clo_isOptionSet(option)) {
        int suiteLoaded = 0;
        sprintf(tmpStr, "%s/common/l3mapgen_defaults_%s.par", dataRoot, 
                clo_getOptionString(option));
        if (access(tmpStr, R_OK) != -1) {
            if (want_verbose)
                printf("Loading default parameters from %s\n", tmpStr);
            clo_readFile(list, tmpStr);
            suiteLoaded = 1;
        }

        sprintf(tmpStr, "%s/%s/l3mapgen_defaults_%s.par", dataRoot, sensorDir,
                clo_getOptionString(option));
        if (access(tmpStr, R_OK) != -1) {
            if (want_verbose)
                printf("Loading default parameters from %s\n", tmpStr);
            clo_readFile(list, tmpStr);
            suiteLoaded = 1;
        }

        if (subsensorId != -1) {
            sprintf(tmpStr, "%s/%s/%s/l3mapgen_defaults_%s.par", dataRoot,
                    sensorDir, subsensorId2SubsensorDir(subsensorId),
                    clo_getOptionString(option));
            if (access(tmpStr, R_OK) != -1) {
                if (want_verbose)
                    printf("Loading default parameters from %s\n", tmpStr);
                clo_readFile(list, tmpStr);
                suiteLoaded = 1;
            }
        }

        if (!suiteLoaded) {
            printf("-E- Failed to load parameters for suite %s for sensor %s\n", clo_getOptionString(option),
                    sensorId2SensorName(sensorId));
            exit(EXIT_FAILURE);
        }
    }
    // enable the dump option
    clo_setEnableDumpOptions(1);
    // make the command line over ride
    clo_readArgs(list, argc, argv);

    return 0;
}
