#include <netcdf>
#include "read_pixel_anc_file.h"
#include "scaled_nc_var.hpp"
#include <boost/algorithm/string.hpp>
#include <memory>
#include <vector>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

// global variables
static bool firstRun = true;
NcFile* ncFile;
ScaledNcVar windSpeedVar;
ScaledNcVar windDirectionVar;
ScaledNcVar windMeridionalVar;
ScaledNcVar windZonalVar;
ScaledNcVar surfacePressureVar;
ScaledNcVar ozoneVar;
ScaledNcVar waterVaporVar;
ScaledNcVar relativeHumidityVar;
ScaledNcVar NO2Var;
std::vector<float> no2_total;
std::vector<float> ozone;
template <typename T = NcFile>
ScaledNcVar find_var(T& nc_handle, const string& varname) {
    ScaledNcVar var;

    auto vars = nc_handle.getVars();
    if (vars.count(varname))
        return nc_handle.getVar(varname);
    auto grps = nc_handle.getGroups();
    for (const auto& key_value_pair : grps) {
        var = find_var(key_value_pair.second, varname);
        if (!var.isNull())
            return var;
    }
    return var;
}

extern "C" void read_pixel_anc_file(char* filename, l1str* l1rec, char * pixel_anc_vars) {
    if (firstRun) {
        firstRun = false;
        std::set<std::string> requested_vars;
        boost::split(requested_vars, pixel_anc_vars, boost::is_any_of(", "), boost::algorithm::token_compress_on);
        try {
            ncFile = new NcFile(std::string(filename), NcFile::read);
            string title;
            ncFile->getAtt("title").getValues(title);
            printf("Loading Pixel Ancillary Data from %s, file type %s\n", filename, title.c_str());
            if (ncFile->getDim("number_of_lines").getSize() != (size_t)l1rec->l1file->nscan) {
                printf("-E- %s:%d - number_of_lines in %s, must be equal to number of line in the L1 file.\n",
                       __FILE__, __LINE__, filename);
                ncFile->close();
                exit(EXIT_FAILURE);
            }
            if (ncFile->getDim("pixels_per_line").getSize() != (size_t)l1rec->l1file->npix) {
                printf(
                    "-E- %s:%d - pixels_per_line in %s, must be equal to number of pixels in the L1 file.\n",
                    __FILE__, __LINE__, filename);
                ncFile->close();
                exit(EXIT_FAILURE);
            }
            if (requested_vars.count("wind_speed") || std::string(pixel_anc_vars) == "ALL")
                windSpeedVar = find_var(*ncFile, "wind_speed");
            if (requested_vars.count("wind_direction") || std::string(pixel_anc_vars) == "ALL")
                windDirectionVar = find_var(*ncFile, "wind_direction");
            if (requested_vars.count("wind_meridional") || std::string(pixel_anc_vars) == "ALL")
                windMeridionalVar = find_var(*ncFile, "wind_meridional");
            if (requested_vars.count("wind_zonal") || std::string(pixel_anc_vars) == "ALL")
                windZonalVar = find_var(*ncFile, "wind_zonal");
            if (requested_vars.count("surface_pressure") || std::string(pixel_anc_vars) == "ALL")
                surfacePressureVar = find_var(*ncFile, "surface_pressure");
            if (requested_vars.count("ozone") || std::string(pixel_anc_vars) == "ALL")
                ozoneVar = find_var(*ncFile, "ozone");
            if (requested_vars.count("water_vapor") || std::string(pixel_anc_vars) == "ALL")
                waterVaporVar = find_var(*ncFile, "water_vapor");
            if (requested_vars.count("relative_humidity") || std::string(pixel_anc_vars) == "ALL")
                relativeHumidityVar = find_var(*ncFile, "relative_humidity");
            if (requested_vars.count("no2") || std::string(pixel_anc_vars) == "ALL")
                NO2Var = find_var(*ncFile, "no2");
            if (!NO2Var.isNull()) {
                no2_total = std::vector<float>(l1rec->npix, BAD_FLT);
            }
            if(!ozoneVar.isNull()){
                ozone = std::vector<float>(l1rec->npix, BAD_FLT);
            }
        } catch (NcException& e) {
            printf("-E- %s:%d - Problem initializing pixel_anc_file %s\n", __FILE__, __LINE__, filename);
            printf("            %s\n", e.what());
            ncFile->close();
            exit(EXIT_FAILURE);
        }

    }  // firstRun

    try {
        vector<size_t> start;
        vector<size_t> count;
        start.push_back(l1rec->iscan);
        start.push_back(l1rec->l1file->spix);
        count.push_back(1);
        count.push_back(l1rec->npix);

        if (!windSpeedVar.isNull()) {
            windSpeedVar.getVar(start, count, l1rec->ws);
        }
        if (!windDirectionVar.isNull()) {
            windDirectionVar.getVar(start, count, l1rec->wd);
        }
        if (!windMeridionalVar.isNull()) {
            windMeridionalVar.getVar(start, count, l1rec->mw);
        }
        if (!windZonalVar.isNull()) {
            windZonalVar.getVar(start, count, l1rec->zw);
        }
        if (!surfacePressureVar.isNull()) {
            surfacePressureVar.getVar(start, count, l1rec->pr);
        }
        if (!ozoneVar.isNull()) {
            ozoneVar.getVar(start, count, ozone.data());
            for (int i = 0; i < l1rec->npix; i++) {
                if (ozone[i] != BAD_FLT)
                    l1rec->oz[i] = ozone[i]/1e3; // l2gen does NOT store ozone in Dobson units, for some reason
            }
        }
        if (!waterVaporVar.isNull()) {
            waterVaporVar.getVar(start, count, l1rec->wv);
        }
        if (!relativeHumidityVar.isNull()) {
            relativeHumidityVar.getVar(start, count, l1rec->rh);
        }
        if (!NO2Var.isNull()) {
            NO2Var.getVar(start, count, no2_total.data());

            for (int i = 0; i < l1rec->npix; i++) {
                l1rec->no2_frac[i] = 1;
                if (no2_total[i] == BAD_FLT)
                    continue;
                if (l1rec->no2_strat[i] == BAD_FLT)
                    continue;
                float no2_tropo = no2_total[i] - l1rec->no2_strat[i];
                if (no2_tropo > 0)
                    l1rec->no2_tropo[i] = no2_total[i] - l1rec->no2_strat[i];
            }
        }

    } catch (NcException& e) {
        printf("-E- %s:%d - Problem reading pixel_anc_file %s\n", __FILE__, __LINE__, filename);
        printf("            %s\n", e.what());
        ncFile->close();
        exit(EXIT_FAILURE);
    }
}
