#include "common.h"

#define ANCSIZE 104

typedef struct {
    int32_t iyear;
    int32_t iday;
    double sec;
} time_struct;

int gen_time_string(time_struct &starttime, time_struct &endtime, std::string &startstring,
                    std::string &endstring);

extern "C" int isleap(int year);
extern "C" int32_t jday(int16_t i, int16_t j, int16_t k);
extern "C" int jdate(int32_t julian, int32_t *year, int32_t *doy);
