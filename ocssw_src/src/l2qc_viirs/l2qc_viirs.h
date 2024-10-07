#ifndef L2QC_VIIRS_H /* avoid re-inclusion */
#define L2QC_VIIRS_H

#include <clo.h>

extern clo_optionList_t* optionList;

int l2qcviirs_init_options(clo_optionList_t* list, const char* softwareVersion);
int l2qcviirs_read_options(clo_optionList_t* list, int argc, char* argv[]);

#endif /* L2QC_VIIRS_H */
