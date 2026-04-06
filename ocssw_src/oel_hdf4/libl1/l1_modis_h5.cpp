#include "libnav.h"
#include "l1.h"
#include <libgen.h>

#include <iostream>
#include <string>
#include <netcdf>
#include <algorithm>
#include <cassert>
#include "modis_bands.hpp"
#include "ncfile_utils.hpp"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

struct ScaledVar {
    NcVar var;
    vector<float> scale{1.0};
    vector<float> offset{0.0};
};

void copyvals(vector<float> inVec, size_t ibeg, size_t iend, float *outArr) {
    auto iterator = inVec.begin();
    std::copy(iterator + ibeg, iterator + iend, outArr);
}

////////////////////////////////////////////////////////////////////////////////

class Geo {
   private:
    netCDF::NcFile ncfile;
    size_t i;

   public:
    string filename;
    string shortname;   // "ShortName"
    string inputfiles;  // "InputDataProducts"
    short npixl;        // "Max Earth Frames"
    short nscan;        // "Number of Scans"
    short ndets = 10;
    short resolution = 1000;

    int orbitnumber;                 // "OrbitNumber"
    float equatorcrossinglongitude;  // "EquatorCrossingLongitude"
    string equatorcrossingdate;      // "EquatorCrossingDate"
    string equatorcrossingtime;      // "EquatorCrossingTime"
    string rangebeginningdate;       // "RangeBeginningDate"
    string rangebeginningtime;       // "RangeBeginningTime"
    vector<unsigned int> gflags;     // "Cumulated gflags"

    // preload per-scan variables
    vector<double> ev_start_time;        // "EV start time"
    vector<unsigned short> mirror_side;  // "Mirror side"
    vector<float> attitude_angles;       // "attitude_angles": [nscan,3]
    vector<float> T_inst2ECR;            // "T_inst2ECR": [nscan,3,3]

    // variables to be read one MODIS scan at a time
    vector<ScaledVar> scanvars;

    // constructor
    Geo(const string filepath) {
        filename.assign(filepath);

        try {
            ncfile.open(filepath, NcFile::read);
            NcGroup grp;
            scanvars.resize(NUM_GEO_VARS);
            gflags.resize(8);

            // read global attributes
            ncfile.getAtt("Number of Scans").getValues(&nscan);
            ncfile.getAtt("Max Earth Frames").getValues(&npixl);
            ncfile.getAtt("ShortName").getValues(shortname);
            ncfile.getAtt("InputDataProducts").getValues(inputfiles);

            ncfile.getAtt("OrbitNumber").getValues(&orbitnumber);
            ncfile.getAtt("EquatorCrossingLongitude").getValues(&equatorcrossinglongitude);
            ncfile.getAtt("EquatorCrossingDate").getValues(equatorcrossingdate);
            ncfile.getAtt("EquatorCrossingTime").getValues(equatorcrossingtime);
            ncfile.getAtt("RangeBeginningDate").getValues(rangebeginningdate);
            ncfile.getAtt("RangeBeginningTime").getValues(rangebeginningtime);
            ncfile.getAtt("Cumulated gflags").getValues(gflags.data());

            // check file type
            if (shortname.find("D03") == std::string::npos) {
                cerr << "Input file " << filename << " has type " << shortname
                     << "; does not contain MODIS geolocation." << endl;
                // TODO: throw an error
            }

            // read root-level variables
            readVar(ncfile, "EV start time", ev_start_time);
            readVar(ncfile, "Mirror side", mirror_side);
            readVar(ncfile, "attitude_angles", attitude_angles);
            readVar(ncfile, "T_inst2ECR", T_inst2ECR);

            // find vars in group Geolocation Fields
            // (stored as float)
            grp = ncfile.getGroup("Geolocation Fields", netCDF::NcGroup::AllChildrenGrps);
            scanvars[LONGITUDE].var = grp.getVar("Longitude");
            scanvars[LATITUDE].var = grp.getVar("Latitude");

            // find vars in group Data Fields
            // (stored as scaled short)
            grp = ncfile.getGroup("Data Fields", netCDF::NcGroup::AllChildrenGrps);
            scanvars[HEIGHT].var = grp.getVar("Height");
            scanvars[SOLAR_ZENITH].var = grp.getVar("SolarZenith");
            scanvars[SOLAR_AZIMUTH].var = grp.getVar("SolarAzimuth");
            scanvars[SENSOR_ZENITH].var = grp.getVar("SensorZenith");
            scanvars[SENSOR_AZIMUTH].var = grp.getVar("SensorAzimuth");

            // load scale and offset
            for (size_t i = 0; i < scanvars.size(); i++) {
                readAtt(scanvars[i].var, "scale_factor", &scanvars[i].scale);
                readAtt(scanvars[i].var, "add_offset", &scanvars[i].offset);
            }
        }

        catch (NcException &e) {
            cerr << "Error reading MODIS geolocation file " << filename << endl;
            cerr << e.what() << endl;
            exit(EXIT_FAILURE);
        }
    }  // end constructor

    // destructor
    ~Geo() {
        try {
            ncfile.close();
        } catch (NcException &e) {
            cerr << e.what() << endl;
        }
    }

    // member functions
    void print() {
        cout << endl;
        cout << "filename\t= " << filename << endl;
        cout << "shortname\t= " << shortname << endl;
        cout << "resolution\t= " << resolution << endl;
        cout << "ndets\t= " << ndets << endl;
        cout << "nscan\t= " << nscan << endl;
        cout << "npixl\t= " << npixl << endl;
        cout << "inputs:\t" << inputfiles << endl;
        for (ScaledVar v : scanvars)
            if (!v.var.isNull())
                cout << v.var.getName() << "\t" << v.var.getType().getName() << endl;
        cout << endl;
    }
};  // class Geo

static Geo *geo = nullptr;

////////////////////////////////////////////////////////////////////////////////

class L1b {
   private:
    netCDF::NcFile ncfile;

   public:
    string filename;
    string shortname;   // "ShortName"
    string inputfiles;  // "InputDataProducts"
    short npixl;        // "Max Earth View Frames"
    short nscan;        // "Number of Scans"
    short ndets;        // 10, 20, or 40
    short resolution;   // 1000, 500 or 250m
    short spixl = 1;    // "extract_pixel_start"

    // variables to be read one MODIS scan at a time
    vector<ScaledVar> scanvars;

    // constructor
    L1b(const string filepath) {
        filename.assign(filepath);

        try {
            ncfile.open(filepath, NcFile::read);
            NcGroup grp;
            scanvars.resize(NUM_L1B_VARS);

            // read scalar global attributes
            ncfile.getAtt("Number of Scans").getValues(&nscan);
            ncfile.getAtt("Max Earth View Frames").getValues(&npixl);
            ncfile.getAtt("ShortName").getValues(shortname);
            ncfile.getAtt("InputDataProducts").getValues(inputfiles);
            NcGroupAtt attr = ncfile.getAtt("extract_pixel_start");
            if (!attr.isNull())
                attr.getValues(&spixl);
            spixl--;  // make zero-based

            // find vars in group Data Fields
            // (stored as scaled short)
            grp = ncfile.getGroup("Data Fields", netCDF::NcGroup::AllChildrenGrps);

            if (shortname.find("D021KM") != std::string::npos) {
                resolution = 1000;
                ndets = 10;
                scanvars[RSB_250].var = grp.getVar("EV_250_Aggr1km_RefSB");
                scanvars[RSB_500].var = grp.getVar("EV_500_Aggr1km_RefSB");
                scanvars[RSB_1KM].var = grp.getVar("EV_1KM_RefSB");
                scanvars[CIR_1KM].var = grp.getVar("EV_Band26");
                scanvars[TEB_1KM].var = grp.getVar("EV_1KM_Emissive");

            } else if (shortname.find("D02HKM") != std::string::npos) {
                resolution = 500;
                ndets = 20;
                scanvars[RSB_250].var = grp.getVar("EV_250_Aggr500_RefSB");
                scanvars[RSB_500].var = grp.getVar("EV_500_RefSB");

            } else if (shortname.find("D02QKM") != std::string::npos) {
                resolution = 250;
                ndets = 40;
                scanvars[RSB_250].var = grp.getVar("EV_250_RefSB");

            } else {
                cerr << "Input file " << filename << " has type " << shortname
                     << "; does not contain MODIS calibrated radiances." << endl;
                // TODO: throw an error
            }

            // load scale and offset
            for (size_t i = 0; i < scanvars.size(); i++) {
                if (i == TEB_1KM) {  // read thermal bands as radiance
                    readAtt(scanvars[i].var, "radiance_scales", &scanvars[i].scale);
                    readAtt(scanvars[i].var, "radiance_offsets", &scanvars[i].offset);
                } else {  // read reflective bands as reflectance; convert to radiance later
                    readAtt(scanvars[i].var, "reflectance_scales", &scanvars[i].scale);
                    readAtt(scanvars[i].var, "reflectance_offsets", &scanvars[i].offset);
                }
            }
        }

        catch (NcException &e) {
            cerr << "Error reading MODIS L1b file " << filename << endl;
            cerr << e.what() << endl;
            exit(EXIT_FAILURE);
        }
    }

    // destructor
    ~L1b() {
        try {
            ncfile.close();
        } catch (NcException &e) {
            cerr << e.what() << endl;
        }
    }

    // member functions
    void print() {
        cout << endl;
        cout << "filename\t= " << filename << endl;
        cout << "shortname\t= " << shortname << endl;
        cout << "resolution\t= " << resolution << endl;
        cout << "ndets\t= " << ndets << endl;
        cout << "nscan\t= " << nscan << endl;
        cout << "npixl\t= " << npixl << endl;
        cout << "inputs:\t" << inputfiles << endl;
        for (ScaledVar v : scanvars)
            if (!v.var.isNull())
                cout << v.var.getName() << "\t" << v.var.getType().getName() << endl;
        cout << endl;
    }

};  // class L1b

static L1b *l1b = nullptr;

////////////////////////////////////////////////////////////////////////////////

class ScanVar {
   public:
    NcVar var;
    size_t scandim;
    size_t per_scan;
    size_t nvals;
    vector<size_t> start;
    vector<size_t> count;
    vector<float> data;

    // constructor
    ScanVar(const netCDF::NcVar var, const size_t scandim, const size_t per_scan) {
        this->var = var;
        this->scandim = scandim;
        this->per_scan = per_scan;

        // set up start and count indices
        size_t ndims = var.getDimCount();
        start.resize(ndims);
        count.resize(ndims);
        for (size_t i = 0; i < ndims; i++) {
            start[i] = 0;
            count[i] = var.getDim(i).getSize();
        }
        count[scandim] = per_scan;

        // total number of values in each scan
        nvals = 1;
        for (const size_t &n : count)
            nvals *= n;

        // allocate memory
        data.resize(nvals);
    }

    // member functions
    void load(const size_t iscan) {
        vector<size_t> startp = start;
        startp[scandim] = iscan * count[scandim];
        var.getVar(startp, count, data.data());
    }

};  // class ScanVar

class Scan {
   private:
    size_t i;

   public:
    // Data for a single MODIS scan
    size_t iscan;  // scan number (0-based)

    // constants
    size_t ndets;   // rows per scan
    size_t npixl;   // number of values per variable per line
    size_t nvals;   // number of values per variable per scan
    size_t nbands;  // number of bands
    size_t spixl;   // start pixel for extracted granlules

    // preloaded Geo vars; subset for each scan
    double taisec;     // EV start time
    int32_t mside;     // Mirror side
    double angles[3];  // attitude_angles
    double mnorm[3];   // T_inst2ECR (need only 1st 3 values)

    // geolocation data: [det, pixl]
    vector<ScanVar *> allGeo;

    // geophysical values for all bands: [band, det, pixl]
    vector<ScanVar *> allL1b;

    // constructor
    Scan() {
        ndets = l1b->ndets;
        npixl = l1b->npixl;
        nvals = ndets * npixl;
        nbands = bandList.size();
        spixl = l1b->spixl;
        mside = -1;
        size_t scandim;

        // initialize geo arrays
        allGeo.resize(NUM_GEO_VARS);
        for (i = 0; i < NUM_GEO_VARS; i++) {
            scandim = 0;
            allGeo[i] = new ScanVar(geo->scanvars[i].var, scandim, ndets);
        }

        // initialize l1b arrays
        allL1b.resize(NUM_L1B_VARS);
        for (i = 0; i < NUM_L1B_VARS; i++) {
            scandim = (i == CIR_1KM ? 0 : 1);
            allL1b[i] = new ScanVar(l1b->scanvars[i].var, scandim, ndets);
        }
    }

    // member functions
    void load(size_t iscan, int lonlat) {
        // Reads data for a scan, and scales from SI to physical value

        // ----- preloaded values for this scan -----
        this->iscan = iscan;
        taisec = geo->ev_start_time[iscan];
        mside = geo->mirror_side[iscan];
        for (i = 0; i < 3; i++) {
            size_t index = 3 * iscan + i;
            angles[i] = static_cast<double>(geo->attitude_angles[index]);
            mnorm[i] = static_cast<double>(geo->T_inst2ECR[3 * index]);
        }

        // ----- Geolocation values -----
        for (i = 0; i < NUM_GEO_VARS; i++) {
            // read data
            allGeo[i]->load(iscan);

            // apply traditional scaling: y = mx + b
            for (auto &val : allGeo[i]->data) {
                val *= geo->scanvars[i].scale[0];
                val += geo->scanvars[i].offset[0];
            }
        }

        // ----- Radiometric values -----
        if (!lonlat) {
            size_t iVar, iBand;
            for (iVar = 0; iVar < NUM_L1B_VARS; iVar++)
                allL1b[iVar]->load(iscan);

            // convert scaled integers to radiometric values
            for (i = 0; i < bandList.size(); i++) {
                iVar = bandList[i].varIndex;
                iBand = bandList[i].bandIndex;
                float scale = l1b->scanvars[iVar].scale[iBand];
                float offset = l1b->scanvars[iVar].offset[iBand];
                size_t ibeg = iBand * nvals;
                size_t iend = ibeg + nvals;  // 1 past end of range

                // apply MODIS-style "backwards" scaling: y = m(x - b)
                for (size_t j = ibeg; j < iend - 1; j++)
                    if (allL1b[iVar]->data[j] < MIN_BAD_SI) {
                        allL1b[iVar]->data[j] -= offset;
                        allL1b[iVar]->data[j] *= scale;
                    }

            }  // bands

        }  // not lonlat
    }  // load()

};  // class Scan

static Scan *scan = nullptr;

////////////////////////////////////////////////////////////////////////////////

int set_l1rec_scanvals(const int32_t iscan, l1str *l1rec) {
    // set all scan-based values in l1rec
    double esdist;
    static double last_taisec = -1;
    int16_t badatt, badmside;
    double taisec, dsec;
    int16_t year, day;
    int32_t i;

    // mirror side
    l1rec->mside = scan->mside;
    badmside = (l1rec->mside < 0);
    if (badmside)  // fix mirror side if set to -1
        l1rec->mside = 0;

    // check for non-nominal roll, pitch, or yaw
    badatt = (fabs(scan->angles[0]) > MAX_ATTERR) || (fabs(scan->angles[1]) > MAX_ATTERR) ||
             (fabs(scan->angles[2]) > MAX_ATTERR);

    // flag each pixel of scan for geolocation errors
    for (i = 0; i < l1rec->npix; i++) {
        l1rec->pixnum[i] = i + scan->spixl;
        l1rec->navfail[i] = badmside;
        l1rec->navwarn[i] = badatt;
    }

    // scan start time
    taisec = scan->taisec;
    if (taisec < 0) {  // workaround for scan time errors
        cerr << "-W- " << __FILE__ << ": Bad time for scan " << iscan << ";";
        if (iscan == 0) {
            cerr << "using granule start time." << endl;
            taisec = unix_to_tai93(
                isodate2unix((geo->rangebeginningdate + "T" + geo->rangebeginningtime).c_str()));
        } else {
            cerr << "extrapolating from last good scan." << endl;
            taisec = last_taisec + SCAN_TIME_INTERVAL;
        }
    }
    l1rec->scantime = tai93_to_unix(taisec);
    last_taisec = taisec;

    // Earth-Sun distance correction for this scan
    unix2yds(l1rec->scantime, &year, &day, &dsec);
    int32_t msec = (int32_t)(dsec * 1000.0);
    int32_t yr = (int32_t)year;
    int32_t dy = (int32_t)day;
    esdist = esdist_(&yr, &dy, &msec);
    l1rec->fsol = pow(1.0 / esdist, 2);

    return SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
static float *Fobar;

extern "C" int openl1_modis_h5(filehandle *l1file) {

    // make sure a geofile was supplied
    if ((l1file->geofile == nullptr) || (l1file->geofile[0] == '\0')) {
        cerr << "-E- No MODIS geolocation file specified." << endl;
        exit(EXIT_FAILURE);
    }

    // open files
    geo = new Geo(l1file->geofile);
    l1b = new L1b(l1file->name);

    assert(geo->npixl == l1b->npixl);
    assert(geo->nscan == l1b->nscan);
    assert(geo->ndets == l1b->ndets);
    assert(geo->resolution == l1b->resolution);

    scan = new Scan();

    // L1B file must have 1KM resolution
    int status = (l1b->resolution != 1000);
    if (status) {
        cerr << "Input file " << l1file->name << " has " << l1b->resolution << "m resolution; "
             << "please specify a 1km-resolution Level 1B file." << endl;
        return status;
    }

    // populate filehandle structure
    l1file->ndets = l1b->ndets;
    l1file->nscan = l1b->nscan * l1b->ndets;
    l1file->npix = l1b->npixl * l1b->ndets / 10;
    strcpy(l1file->spatialResolution, "1 km");  // TODO: revisit high resolution
    strcpy(l1_input->calfile, l1b->inputfiles.c_str());

    l1file->orbit_number = geo->orbitnumber;
    l1file->orbit_node_lon = geo->equatorcrossinglongitude;
    l1file->node_crossing_time =
        isodate2unix((geo->equatorcrossingdate + "T" + geo->equatorcrossingtime).c_str());
    l1file->terrain_corrected = ((int32_t)(geo->gflags[5]) == 0);

    // Read reflectance to radiance conversion factors
    rdsensorinfo(l1file->sensorID, l1_input->evalmask, "Fobar", (void **)&Fobar);

    return SUCCESS;
}

extern "C" int readl1_modis_h5(filehandle *l1file, const int32_t line, l1str *l1rec, const int lonlat) {
    static int32_t prev_scan = -1;
    int32_t ndets, iscan, idet;

    // Scan and detector numbers for this line
    ndets = l1b->ndets;
    iscan = line / ndets;
    idet = line % ndets;

    // Load all values for next scan as needed
    if (iscan != prev_scan) {
        prev_scan = iscan;
        scan->load(iscan, lonlat);
    }

    // The l1rec structure is reinitialized each line, so always reload basic info.
    l1rec->npix = l1file->npix;
    l1rec->detnum = (int32_t)idet;
    set_l1rec_scanvals(iscan, l1rec);

    // get pixel range within variable
    size_t ibeg = idet * scan->npixl;
    size_t iend = ibeg + scan->npixl;

    // ----- Geolocation values -----
    copyvals(scan->allGeo[LONGITUDE]->data, ibeg, iend, l1rec->lon);
    copyvals(scan->allGeo[LATITUDE]->data, ibeg, iend, l1rec->lat);
    copyvals(scan->allGeo[SOLAR_ZENITH]->data, ibeg, iend, l1rec->solz);

    // lonlat needs only time, lon, lat and solar zenith
    if (lonlat)
        return SUCCESS;

    copyvals(scan->allGeo[SOLAR_AZIMUTH]->data, ibeg, iend, l1rec->sola);
    copyvals(scan->allGeo[SENSOR_ZENITH]->data, ibeg, iend, l1rec->senz);
    copyvals(scan->allGeo[SENSOR_AZIMUTH]->data, ibeg, iend, l1rec->sena);
    copyvals(scan->allGeo[HEIGHT]->data, ibeg, iend, l1rec->height);

    // Compute polarization frame rotation angles
    compute_alpha(l1rec->lon, l1rec->lat, l1rec->senz, l1rec->sena, scan->mnorm, l1rec->npix, l1rec->alpha);

    // ----- Radiometric values -----
    size_t irsb = 0;
    size_t iteb = 0;
    float value;
    vector<float> thisBand;
    thisBand.resize(scan->nvals);
    size_t iVar, iBand;
    int32_t ip, ipb;

    for (size_t i = 0; i < bandList.size(); i++) {
        iVar = bandList[i].varIndex;
        iBand = bandList[i].bandIndex;

        // get pixel range within multi-band variable
        ibeg = (iBand * ndets + idet) * scan->npixl;
        iend = ibeg + scan->npixl;  // 1 past end of range

        // extract all data for this band & scan
        copyvals(scan->allL1b[iVar]->data, ibeg, iend, thisBand.data());

        // Thermal Emissive Bands
        if (iVar == TEB_1KM) {
            for (ip = 0; ip < l1rec->npix; ip++) {
                value = 0.0;

                if (thisBand[ip] < MIN_BAD_SI)
                    // convert radiance from Watts/m^2/micrometer/steradian to mW/cm2/um/sr
                    value /= 10.0;

                // TODO: Terra & Aqua adjustments?

                // populate l1rec->Ltir
                ipb = ip * l1rec->l1file->nbandsir + iteb;  // band varies fastest
                l1rec->Ltir[ipb] = value;

            }  // ip
            iteb++;
        }  // TEB

        // Cirrus Band
        else if (iVar == CIR_1KM) {
            for (ip = 0; ip < l1rec->npix; ip++) {
                value = BAD_FLT;

                if (thisBand[ip] < MIN_BAD_SI)
                    // normalize reflectance by solar zenith angle
                    value = thisBand[ip] / cos(l1rec->solz[ip] / OEL_RADEG);

                // populate l1rec->rho_cirrus
                l1rec->rho_cirrus[ip] = value;

            }  // ip
        }  // CIR

        // Reflective Solar Bands
        else {
            l1rec->Fo[irsb] = Fobar[irsb] * l1rec->fsol;

            for (ip = 0; ip < l1rec->npix; ip++) {
                value = BAD_FLT;

                // skip night data for visible bands
                if (SOLZNIGHT < l1rec->solz[ip])
                    continue;

                if (thisBand[ip] < MIN_BAD_SI)
                    // convert from reflectance to radiance
                    value = thisBand[ip] * l1rec->Fo[irsb] / OEL_PI;

                else
                    // flag any sentinel values
                    switch ((int32_t)thisBand[ip]) {
                        case SECTOR_ROTATION_SI:  // not in earth view
                            l1rec->navfail[ip] = 1;
                            break;
                        case TEB_OR_RSB_GT_MAX_SI:   // saturation
                        case SATURATED_DETECTOR_SI:  // saturation
                            value = 1000.0;
                            if (irsb < 13)  // set flag below 1000nm
                                l1rec->hilt[ip] = 1;
                            break;
                        default:  // all others default to BAD_FLT
                            break;
                    }

                // populate l1rec->Lt
                ipb = ip * l1rec->l1file->nbands + irsb;  // band varies fastest
                l1rec->Lt[ipb] = value;

            }  // ip
            irsb++;
        }  // RSB

    }  // for i

    return SUCCESS;
}

extern "C" int closel1_modis_h5() {
    return SUCCESS;
}

extern "C" int readl1_lonlat_modis_h5(filehandle *l1file, const int32_t line, l1str *l1rec) {
    int status = readl1_modis_h5(l1file, line, l1rec, 1);
    return status;
}
