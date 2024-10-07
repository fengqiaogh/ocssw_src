#ifndef L1CGEN_H /* avoid re-inclusion */
#define L1CGEN_H

#include <clo.h>

// extern clo_optionList_t* optionList;

int l1cgen_init_options(clo_optionList_t* list, const char* softwareVersion);
int l1cgen_read_options(clo_optionList_t* list, int argc, char* argv[]);

#endif /* L1CGEN_H */
