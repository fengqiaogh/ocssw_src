#ifndef ObpgReadGrid_h
#define ObpgReadGrid_h

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#ifndef SUCCESS
#define SUCCESS 0
#endif
#ifndef FAIL
#define FAIL    1
#endif

typedef struct {
    size_t ny;
    double nyPerDeg;
    double latOrigin;
    size_t nx;
    double nxPerDeg;
    double lonOrigin;
    unsigned char **data;
} GridStruct;

void printGrid(GridStruct grid);
void deletenc2d(unsigned char**& data);
unsigned char** readnc2d(NcVar ncVar,
        size_t y0, size_t y1,
        size_t x0, size_t x1);

NcVar openGridFile(string filePath, string varName);
int getGridValue(GridStruct grid, float lat, float lon, unsigned char *value);
int loadGrid(NcVar var, float NSEW[4], GridStruct *grid);

#endif
