/*
 *  read a per pixel ancillary data file.
 */

#include <netcdf>
#include "read_pixel_anc_file.h"

#include <vector>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

// global variables
static bool firstRun = true;
NcFile *ncFile = nullptr;
NcVar windSpeedVar;
NcVar windDirectionVar;
NcVar windMeridionalVar;
NcVar windZonalVar;
NcVar surfacePressureVar;
NcVar ozoneVar;
NcVar waterVaporVar;
NcVar relativeHumidityVar;
//NcVar NO2Var;

extern "C" void read_pixel_anc_file(char* filename, l1str *l1rec) {
    
    if(firstRun) {
        firstRun = false;

        try {
            printf("Loading Pixel Ancillary Data from %s\n", filename);
            ncFile = new NcFile(filename, NcFile::read);
            string title;
            ncFile->getAtt("title").getValues(title);
            if(title != "Pixel Ancillary Data") {
                printf("-E- %s:%d - %s is not a per Pixel Ancillary Data file.\n", 
                        __FILE__, __LINE__, filename);
                exit(EXIT_FAILURE);
            }
            
            if(ncFile->getDim("number_of_lines").getSize() != (size_t)l1rec->l1file->nscan) {
                printf("-E- %s:%d - number_of_lines in %s, must be equal to number of line in the L1 file.\n", 
                        __FILE__, __LINE__, filename);
                exit(EXIT_FAILURE);
            }
            if(ncFile->getDim("pixels_per_line").getSize() != (size_t)l1rec->l1file->npix) {
                printf("-E- %s:%d - pixels_per_line in %s, must be equal to number of pixels in the L1 file.\n", 
                        __FILE__, __LINE__, filename);
                exit(EXIT_FAILURE);
            }

            windSpeedVar = ncFile->getVar("wind_speed");
            windDirectionVar = ncFile->getVar("wind_direction");
            windMeridionalVar = ncFile->getVar("wind_meridional");
            windZonalVar = ncFile->getVar("wind_zonal");
            surfacePressureVar = ncFile->getVar("surface_pressure");
            ozoneVar = ncFile->getVar("ozone");
            waterVaporVar = ncFile->getVar("water_vapor");
            relativeHumidityVar = ncFile->getVar("relative_humidity");
            //NO2Var = ncFile->getVar("NO2");
        
        } catch   (NcException& e) {
            printf("-E- %s:%d - Problem initializing pixel_anc_file %s\n", __FILE__, __LINE__, filename);
            printf("            %s\n", e.what());
            exit(EXIT_FAILURE);
        }  

     } // firstRun


    try {
        vector<size_t> start;
        vector<size_t> count;
        start.push_back(l1rec->iscan);
        start.push_back(l1rec->l1file->spix);
        count.push_back(1);
        count.push_back(l1rec->npix);
        
        if(!windSpeedVar.isNull()) {
            windSpeedVar.getVar(start, count, l1rec->ws);
        }
        if(!windDirectionVar.isNull()) {
            windDirectionVar.getVar(start, count, l1rec->wd);
        }
        if(!windMeridionalVar.isNull()) {
            windMeridionalVar.getVar(start, count, l1rec->mw);
        }
        if(!windZonalVar.isNull()) {
            windZonalVar.getVar(start, count, l1rec->zw);
        }
        if(!surfacePressureVar.isNull()) {
            surfacePressureVar.getVar(start, count, l1rec->pr);
        }
        if(!ozoneVar.isNull()) {
            ozoneVar.getVar(start, count, l1rec->oz);
        }
        if(!waterVaporVar.isNull()) {
            waterVaporVar.getVar(start, count, l1rec->wv);
        }
        if(!relativeHumidityVar.isNull()) {
            relativeHumidityVar.getVar(start, count, l1rec->rh);
        }
//        if(!NO2Var.isNull()) {
//            NO2Var.getVar(start, count, l1rec->no2_strat);
//        }

    } catch   (NcException& e) {
        printf("-E- %s:%d - Problem reading pixel_anc_file %s\n", __FILE__, __LINE__, filename);
        printf("            %s\n", e.what());
        exit(EXIT_FAILURE);
    }

}

