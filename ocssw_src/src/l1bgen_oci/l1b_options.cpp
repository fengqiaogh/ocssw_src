
#include "l1b_options.hpp"
#include <dirent.h>
#include <iostream>
#include <genutils.h>
#include <sys/stat.h>

std::string ifile = "ifile";
std::string ofile = "ofile";
std::string cal_lut = "cal_lut";
std::string geo_lut = "geo_lut";
std::string doi = "doi";
std::string crosstalk_lut = "crosstalk_lut";
std::string pversion = "pversion";
std::string demfile = "demfile";
std::string radiance = "radiance";
std::string disable_aggregation = "disable_aggregation";
std::string disable_geo = "disable_geo";
std::string floating_point = "floating_point";
std::string ephfile = "ephfile";
std::string sline = "sline";
std::string eline = "eline";
std::string enable_crosstalk = "enable_crosstalk";
std::string deflate = "deflate";

oel::L1bOptions::L1bOptions(int argc, char *argv[], const char *version) {
    using namespace std;

    clo_setVersion2("l1bgen_oci", version);
    optionList = clo_createList();
    clo_addOption(optionList, ifile.c_str(), CLO_TYPE_IFILE, NULL, "Input L1A file");
    clo_addOption(optionList, ofile.c_str(), CLO_TYPE_OFILE, NULL, "Output L1B file");
    clo_addOption(optionList, cal_lut.c_str(), CLO_TYPE_IFILE, NULL, "CAL LUT file");
    clo_addOption(optionList, geo_lut.c_str(), CLO_TYPE_IFILE, NULL, "GEO LUT file");
    clo_addOption(optionList, doi.c_str(), CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");
    clo_addOption(optionList, crosstalk_lut.c_str(), CLO_TYPE_IFILE,
                  "$OCDATAROOT/oci/cal/PACE_OCI_cross_coef_LUT_2025-05-02.nc", "Crosstalk LUT file");
    clo_addOption(optionList, pversion.c_str(), CLO_TYPE_STRING, "Unspecified", "processing version string");
    clo_addOption(optionList, demfile.c_str(), CLO_TYPE_STRING, "$OCDATAROOT/common/gebco_ocssw_v2025.nc",
                  "Digital elevation map file");
    clo_addOption(optionList, radiance.c_str(), CLO_TYPE_BOOL, "false",
                  "Generate radiances as opposed to reflectances");
    clo_addOption(optionList, disable_aggregation.c_str(), CLO_TYPE_BOOL, "false",
                  "Turn off aggregation for instrument bands");
    clo_addOption(optionList, disable_geo.c_str(), CLO_TYPE_BOOL, "false", "Disable geolocation");
    clo_addOption(optionList, floating_point.c_str(), CLO_TYPE_BOOL, "false",
                  "Whether to put out sensor/solar angles as floating point numbers");
    clo_addOption(optionList, ephfile.c_str(), CLO_TYPE_STRING, nullptr, "Definitive ephemeris file name");
    clo_addOption(optionList, sline.c_str(), CLO_TYPE_INT, "1", "Starting line number (1 based)");
    clo_addOption(optionList, eline.c_str(), CLO_TYPE_INT, "-1", "Ending line number (1 based, -1=last line)");
    clo_addOption(optionList, enable_crosstalk.c_str(), CLO_TYPE_BOOL, "false", "enable crosstalk correction");
    clo_addOption(optionList, deflate.c_str(), CLO_TYPE_INT, "5",
                  "Compression factor, 0-9. 0 is least compression, 9 is most");
    // Defaulting deflate to 5 results in no change to current version

    string help = "\nReturns:\n   0: All is well\n   100: Geolocation (or selenolocation) failed\n";
    clo_setHelpStr(help.c_str());

    clo_readArgs(optionList, argc, argv);

    if (argc == 1) {
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }

    // Grab the OPER LUTs
    char *ocVarRoot = getenv("OCVARROOT");
    string OCVARROOT;
    if (ocVarRoot) {
        OCVARROOT.assign(ocVarRoot);
    } else {
        cerr << "-E- Environment variable OCVARROOT not set" << endl;
        exit(EXIT_FAILURE);
    }
    const string CALDIR(OCVARROOT + "/oci/cal/OPER/");

    vector<string> lutFiles;  // Look up tables

    DIR *calibrationLutDir;
    struct dirent *caldirptr;
    if ((calibrationLutDir = opendir(CALDIR.c_str())) != NULL) {
        while ((caldirptr = readdir(calibrationLutDir)) != NULL) {
            lutFiles.push_back(string(caldirptr->d_name));
        }
        closedir(calibrationLutDir);
    }

    if (clo_isSet(optionList, ifile.c_str())) {
        l1aFilename = clo_getString(optionList, ifile.c_str());
    } else {
        cout << "-E- Input L1A file name was not provided" << endl;
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }
    if (clo_isSet(optionList, ofile.c_str())) {
        l1bFilename = clo_getString(optionList, ofile.c_str());
    } else {
        cout << "-E- Output L1B file name was not provided" << endl;
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }
    if (clo_isSet(optionList, cal_lut.c_str())) {
        calibrationLutFilename = clo_getString(optionList, cal_lut.c_str());
    } else {
        for (const string &lut : lutFiles) {
            if (lut.find("PACE_OCI_L1B_LUT") != std::string::npos) {
                calibrationLutFilename = CALDIR;
                calibrationLutFilename.append(lut);
                break;
            }
        }
    }
    {
        struct stat _;  // Just to fill a variable.
        if (calibrationLutFilename.empty() || (stat(calibrationLutFilename.c_str(), &_) != 0)) {
            cout << "Error: input CAL LUT file: " << calibrationLutFilename.c_str() << " does not exist\n";
            exit(EXIT_FAILURE);
        }
        if (clo_isSet(optionList, geo_lut.c_str())) {
            geolocationLutFilename = clo_getString(optionList, geo_lut.c_str());
        } else {
            for (const string &lut : lutFiles) {
                if (lut.find("PACE_OCI_GEO_LUT") != std::string::npos) {
                    geolocationLutFilename = CALDIR;
                    geolocationLutFilename.append(lut);
                    break;
                }
            }
        }
        if (geolocationLutFilename.empty() || (stat(geolocationLutFilename.c_str(), &_) != 0)) {
            cout << "Error: input GEO LUT file: " << geolocationLutFilename.c_str() << " does not exist\n";
            exit(EXIT_FAILURE);
        }

        char tmp_filename[FILENAME_MAX];
        parse_file_name(clo_getString(optionList, demfile.c_str()), tmp_filename);
        demFile = tmp_filename;
        if ((stat(demFile.c_str(), &_) != 0)) {
            cout << "Error: DEM file: " << demFile.c_str() << " does not exist\n";
            exit(EXIT_FAILURE);
        }
    }

    radianceGenerationEnabled = clo_getBool(optionList, radiance.c_str());

    if (clo_isSet(optionList, "doi")) {
        digitalObjectId = clo_getString(optionList, doi.c_str());
        if (digitalObjectId == "None") {
            digitalObjectId.clear();
        }
    }

    aggregationOff = clo_getBool(optionList, disable_aggregation.c_str());

    if (clo_isSet(optionList, ephfile.c_str())) {
        ephFile = clo_getString(optionList, ephfile.c_str());
    }

    disableGeolocation = clo_getBool(optionList, disable_geo.c_str());

    if (disableGeolocation && !radianceGenerationEnabled) {
        cout << " -E- Cannot generate reflectances without geolocation";
        exit(EXIT_FAILURE);
    }

    processingVersion = clo_getString(optionList, pversion.c_str());

    xtalkLutFilename = clo_getString(optionList, crosstalk_lut.c_str());
    char tmpFilename[FILENAME_MAX];
    parse_file_name(xtalkLutFilename.c_str(), tmpFilename);
    xtalkLutFilename = tmpFilename;

    enableCrosstalk = clo_getBool(optionList, enable_crosstalk.c_str());

    floatingAngles = clo_getBool(optionList, floating_point.c_str());

    // sline and eline come based on the first line=1,  First line=0 inside the program
    startingLine = clo_getInt(optionList, sline.c_str());
    if (startingLine <= 1)
        startingLine = 0;
    else
        startingLine--;

    endingLine = clo_getInt(optionList, eline.c_str());
    if (endingLine < 1)
        endingLine = -1;

    deflateLevel = clo_getInt(optionList, deflate.c_str());
    if (deflateLevel < 0 || 9 < deflateLevel) {
        cerr << " -E- Invalid deflate level: " << to_string(deflateLevel) << endl;
        exit(EXIT_FAILURE);
    }
}

oel::L1bOptions::~L1bOptions() {  // Would delete the list, but this will run at the end of the program anyway
}
