
#include "l1b_options.hpp"
#include <dirent.h>
#include <iostream>
#include <genutils.h>
#include <sys/stat.h>

oel::L1bOptions::L1bOptions(int argc, char *argv[]) {
    using namespace std;

    optionList = clo_createList();
    clo_addOption(optionList, "ifile", CLO_TYPE_IFILE, NULL, "Input L1A file");
    clo_addOption(optionList, "ofile", CLO_TYPE_OFILE, NULL, "Output L1B file");
    clo_addOption(optionList, "cal_lut", CLO_TYPE_IFILE, NULL, "CAL LUT file");
    clo_addOption(optionList, "geo_lut", CLO_TYPE_IFILE, NULL, "GEO LUT file");
    clo_addOption(optionList, "doi", CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");
    clo_addOption(optionList, "crosstalk_lut", CLO_TYPE_IFILE,
                  "$OCDATAROOT/oci/cal/PACE_OCI_cross_coef_LUT_2025-05-02.nc", "Crosstalk LUT file");
    clo_addOption(optionList, "pversion", CLO_TYPE_STRING, "Unspecified", "processing version string");
    clo_addOption(optionList, "demfile", CLO_TYPE_STRING, "$OCDATAROOT/common/gebco_ocssw_v2020.nc",
                  "Digital elevation map file");
    clo_addOption(optionList, "radiance", CLO_TYPE_BOOL, "false",
                  "Generate radiances as opposed to reflectances");
    clo_addOption(optionList, "disable_geo", CLO_TYPE_BOOL, "false", "Disable geolocation");
    clo_addOption(optionList, "ephfile", CLO_TYPE_STRING, nullptr, "Definitive ephemeris file name");
    clo_addOption(optionList, "sline", CLO_TYPE_INT, "1", "Starting line number (1 based)");
    clo_addOption(optionList, "eline", CLO_TYPE_INT, "-1", "Ending line number (1 based, -1=last line)");
    clo_addOption(optionList, "enable_crosstalk", CLO_TYPE_BOOL, "false", "enable crosstalk correction");
    clo_addOption(optionList, "deflate", CLO_TYPE_INT, "5",
                  "Compression factor, 0-9. 0 is least compression, 9 is most");
    // Defaulting deflate to 5 results in no change to current version

    clo_readArgs(optionList, argc, argv);

    if (argc == 1) {
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }

    // Grab the OPER LUTs
    const string OCVARROOT = getenv("OCVARROOT");
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

    if (clo_isSet(optionList, "ifile")) {
        l1aFilename = clo_getString(optionList, "ifile");
    } else {
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }
    if (clo_isSet(optionList, "ofile")) {
        l1bFilename = clo_getString(optionList, "ofile");
    } else {
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }
    if (clo_isSet(optionList, "cal_lut")) {
        calibrationLutFilename = clo_getString(optionList, "cal_lut");
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
        if (clo_isSet(optionList, "geo_lut")) {
            geolocationLutFilename = clo_getString(optionList, "geo_lut");
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
        parse_file_name(clo_getString(optionList, "demfile"), tmp_filename);
        demFile = tmp_filename;
        if ((stat(demFile.c_str(), &_) != 0)) {
            cout << "Error: DEM file: " << demFile.c_str() << " does not exist\n";
            exit(EXIT_FAILURE);
        }
    }

    radianceGenerationEnabled = clo_getBool(optionList, "radiance");
    if (clo_isSet(optionList, "doi")) {
        digitalObjectId = clo_getString(optionList, "doi");
        if (digitalObjectId == "None") {
            digitalObjectId.clear();
        }
    }

    if (clo_isSet(optionList, "ephfile")) {
        ephFile = clo_getString(optionList, "ephfile");
    }

    disableGeolocation = clo_getBool(optionList, "disable_geo");

    if (disableGeolocation && !radianceGenerationEnabled) {
        cout << " -E- Cannot generate reflectances without geolocation";
        exit(EXIT_FAILURE);
    }

    processingVersion = clo_getString(optionList, "pversion");

    xtalkLutFilename = clo_getString(optionList, "crosstalk_lut");
    char tmpFilename[FILENAME_MAX];
    parse_file_name(xtalkLutFilename.c_str(), tmpFilename);
    xtalkLutFilename = tmpFilename;

    enableCrosstalk = clo_getBool(optionList, "enable_crosstalk");

    // sline and eline come based on the first line=1,  First line=0 inside the program
    startingLine = clo_getInt(optionList, "sline");
    if (startingLine <= 1)
        startingLine = 0;
    else
        startingLine--;

    endingLine = clo_getInt(optionList, "eline");
    if (endingLine < 1)
        endingLine = -1;

    deflate = clo_getInt(optionList, "deflate");
    if (deflate < 0 || 9 < deflate) {
        cerr << " -E- Invalid deflate level: " << to_string(deflate) << endl;
        exit(EXIT_FAILURE);
    }
}

oel::L1bOptions::~L1bOptions() {
    if (optionList)
        clo_deleteList(optionList);
}
