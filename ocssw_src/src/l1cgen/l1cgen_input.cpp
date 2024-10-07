#include "l1cgen.h"

#include <l1.h>

/** add all of the accepted command line options to list */
int l1cgen_init_options(clo_optionList_t* list, const char* softwareVersion) {
    clo_setVersion2("l1cgen", softwareVersion);

    clo_setHelpStr(
        "Usage: l3mapgen argument-list"
        "\n"
        "\n  This program reads an L1 file"
        "\n"
        "\n  Return values"
        "\n    0 = All Good"
        "\n    1 = Error"
        "\n"
        "\n  The argument list is a set of keyword=value pairs.  Arguments can"
        "\n  be specified on the command line, or put into a parameter file, or the"
        "\n  two methods can be used together, with command line overriding."
        "\n"
        "\nThe list of valid keywords follows:"
        "\n");

    // files
    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L1 filename");

    clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output filename");
    clo_addOption(list, "oformat", CLO_TYPE_STRING, "netcdf4",
                  "output file format"
                  "\n        netcdf4: netCDF4 file, can contain more than one product"
                  "\n        hdf4:    HDF4 file (old SMI format)"
                  "\n        png:     PNG image file"
                  "\n        ppm:     PPM image file"
                  "\n        tiff:    TIFF file with georeference tags");

    clo_addOption(list, "spixl", CLO_TYPE_INT, "1", "start pixel number");
    clo_addOption(list, "epixl", CLO_TYPE_INT, "-1", "end pixel number (-1=the last pixel)");
    clo_addOption(list, "dpixl", CLO_TYPE_INT, "1", "pixel sub-sampling interval");
    clo_addOption(list, "sline", CLO_TYPE_INT, "1", "start line number");
    clo_addOption(list, "eline", CLO_TYPE_INT, "-1", "end line number (-1=the last line)");
    clo_addOption(list, "dline", CLO_TYPE_INT, "1", "line sub-sampling interval");
    clo_addOption(list, "doi", CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, "Unspecified", "processing version string");
    clo_addOption(list, "quiet", CLO_TYPE_BOOL, "false", "stop the status printing");

    l1_add_options(list);

    return 0;
}

/*
 Read the command line option and all of the default parameter files.
 */
int l1cgen_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    // read all arguments
    clo_readArgs(list, argc, argv);

    return 0;
}
