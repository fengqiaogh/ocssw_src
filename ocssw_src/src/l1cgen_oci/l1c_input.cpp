#include <string>
#include "l1c_processor.h"

static char* l1cgen_optionKeys[] = {"-help",           "-version",       "par",   "pversion", "doi",
                                    "verbose",         "ifile",          "ofile", "l1c_grid", "cloud_height",
                                    "cloud_anc_files", "area_weighting", "demfile", NULL};
/*
 * @brief Initialize command line options for L1C processing, using the CLO library
 * @param list - pointer to the command line option list to be initialized
 * @param prog - name of the program
 * @param version - version of the program
 * @return true if initialization is successful, false otherwise
 */
bool l1c_clo_ini(clo_optionList_t* list, const std::string& prog, const std::string& version) {
    std::string help_menu;
    clo_setSelectOptionKeys(l1cgen_optionKeys);
    help_menu = "Usage: " + prog + " argument-list\n\n";
    help_menu += "The argument list is a set of keyword=value pairs.  Arguments can\n";
    help_menu += "be specified on the command line, or put into a parameter file, or the\n";
    help_menu += "two methods can be used together, with command line overriding.\n\n";
    help_menu += "The list of valid keywords follows:\n";
    clo_setHelpStr(help_menu.c_str());
    std::string temp_str = "input L1B PACE OCI files";
    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, temp_str.c_str());
    temp_str = "input L1C grid file";
    clo_addOption(list, "l1c_grid", CLO_TYPE_IFILE, NULL, temp_str.c_str());
    temp_str = "input L2 PACE OCI files for cloud ancillary data";
    clo_addOption(list, "cloud_anc_files", CLO_TYPE_IFILE, NULL, temp_str.c_str());
    temp_str = "Cloud top height for cloud top correction, km";
    clo_addOption(list, "cloud_height", CLO_TYPE_FLOAT, NULL, temp_str.c_str());
    clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output for 1 file name");
    clo_addOption(list, "verbose", CLO_TYPE_BOOL, "off", "Allow more verbose screen messages");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, NULL, "processing version string");
    clo_addOption(list, "doi", CLO_TYPE_STRING, NULL, "Digital Object Identifier (DOI) string");
    clo_addOption(list, "area_weighting", CLO_TYPE_INT, "1",
                  "Area weighting for L1C binning, 0 for none, 1 for area weighting.");
    clo_addOption(list, "demfile", CLO_TYPE_STRING, "$OCDATAROOT/common/gebco_ocssw_v2025.nc", "Digital Elevation Model file");
    //*************************************************************************-

    clo_setVersion2(prog.c_str(), version.c_str());
    return true;
}

bool l1c_input(int argc, char** argv, L1CProcessor& input, const std::string& program_name,
               const std::string& version) {
    l1_input_init();
    clo_optionList_t* list = clo_createList();
    l1c_clo_ini(list, program_name, version);
    if (argc == 1) {
        clo_printUsage(list);
        exit(EXIT_SUCCESS);
    }
    clo_printVersion();
    clo_readArgs(list, argc, argv);
    return input.load_input(list);
}
