#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "version.h"
#define STRINGIFY(x) #x
#define TOSTRING(x)  STRINGIFY(x)
#define VERSION      TOSTRING(VERSION_MAJOR) "." TOSTRING(VERSION_MINOR) "." TOSTRING(VERSION_PATCH) "-" GITSHA
#define PROGRAM "l2bin"
#endif