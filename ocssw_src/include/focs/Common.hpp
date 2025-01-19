#ifndef COMMON_UTIL_HPP
#define COMMON_UTIL_HPP
#include "genutils.h"
#include "timeutils.h"
#ifdef BUILD_ALL

#define ESDIST esdist_
#else
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <string>
#define ESDIST esdist
double esdist(int *iyr, int *iday, int *msec);
std::string call_sequence(int argc, char* argv[]);
#endif
#endif //COMMON_UTIL_HPP
