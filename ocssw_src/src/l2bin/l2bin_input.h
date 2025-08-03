#ifndef _INPUT_STR_H
#define _INPUT_STR_H

#include <stdio.h>
#include <vector>
#include <string>
#include "clo.h"

#define DEF_FLAG "ATMFAIL,LAND,HILT,HISATZEN,STRAYLIGHT,CLDICE,COCCOLITH,LOWLW,CHLFAIL,CHLWARN,NAVWARN,ABSAER,MAXAERITER,ATMWARN,HISOLZEN,NAVFAIL,FILTER"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @struct input_struct
 * @brief Structure containing input parameters and configuration for L2 binning
 */
typedef struct input_struct {
    char infile [FILENAME_MAX];        ///< Input file path
    char ofile [FILENAME_MAX];         ///< Output file path  
    char pfile [FILENAME_MAX];         ///< Parameter file path
    char fileuse [FILENAME_MAX];       ///< File usage specification
    char flaguse[2048];               ///< Flag usage specification
    char l3bprod[2048];               ///< L3 binned products list
    char prodtype[32];                ///< Product type
    char qual_prod[2048];             ///< Quality product specification
    char composite_prod[2048];         ///< Composite product specification
    char composite_scheme[2048];       ///< Composite scheme specification
    char pversion[16];                ///< Processing version
    char suite [32];                  ///< Processing suite
    char output_wavelengths[2048];     ///< Output wavelengths specification
    char output_product_names[2048];   ///< Output product names
    char parms [4096];                ///< Parameters string

    int32_t sday;                     ///< Start day
    int32_t eday;                     ///< End day
    char resolve[4];                  ///< Resolution specification
    int32_t rowgroup;                 ///< Row grouping
    int32_t meminfo;                  ///< Memory info flag
    int32_t dcinfo;                   ///< DC info flag
    int32_t night;                    ///< Night processing flag
    int32_t verbose;                  ///< Verbosity level
    int32_t minobs;                   ///< Minimum observations
    float deltaeqcross;               ///< Delta equator crossing time
    int32_t deflate;                  ///< Deflate flag

    float latsouth;                   ///< Southern latitude bound
    float latnorth;                   ///< Northern latitude bound
    float lonwest;                    ///< Western longitude bound
    float loneast;                    ///< Eastern longitude bound

    uint8_t qual_max;                 ///< Maximum quality level
    std::vector<std::string> files;   ///< List of input files

    int area_weighting;               ///< Area weighting flag
    char doi[1024];                   ///< Digital Object Identifier
    
} instr;

/**
 * @brief Processes command line arguments and initializes input structure
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 * @param input Pointer to input structure to populate
 * @param prog Program name
 * @param version Program version
 * @return 0 on success, non-zero on failure
 */
int l2bin_input(int argc, char **argv, instr *input, const char* prog, const char* version);

/**
 * @brief Initializes command line options
 * @param list Pointer to options list structure
 * @param prog Program name
 * @param version Program version
 * @return 0 on success, non-zero on failure
 */
int l2bin_init_options(clo_optionList_t* list, const char* prog, const char* version);

#ifdef __cplusplus
}
#endif
#endif
