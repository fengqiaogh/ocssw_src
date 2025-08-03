#include <stdio.h>
#include <math.h>
#include <libgen.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <netcdf>

#include <allocate2d.h>
#include <allocate3d.h>
#include <timeutils.h>
#include <clo.h>

#include "copy_var_utils.hpp"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;



int main(int argc, char *argv[]) {
    NcVar var, varin, varout;

    clo_optionList_t *list;
    clo_option_t *option;
    char *strVal;
    char keyword[50];
    char inAncillaryFilename[FILENAME_MAX];
    char outAncillaryFilename[FILENAME_MAX];

    list = clo_createList();

    option = clo_addOption(list, "inancfile", CLO_TYPE_IFILE, NULL, "Input ancillary file file");

    option = clo_addOption(list, "outancfile", CLO_TYPE_OFILE, NULL, "Output ancillary file file");

    clo_setVersion("0.1");

    if (argc == 1) {
        clo_printUsage(list);
        exit(1);
    }

    clo_readArgs(list, argc, argv);

    strVal = clo_getString(list, "inancfile");
    strcpy(inAncillaryFilename, strVal);

    strVal = clo_getString(list, "outancfile");
    strcpy(outAncillaryFilename, strVal);

    int numOptions, optionId;

    numOptions = clo_getNumOptions(list);
    for (optionId = 0; optionId < numOptions; optionId++) {
        option = clo_getOption(list, optionId);

        strcpy(keyword, option->key);
    }

    cout << endl << "Opening ancillary file" << endl;
    NcFile *inAncillaryfile = new NcFile(inAncillaryFilename, NcFile::read);

    NcDim acrossDim = inAncillaryfile->getDim("bins_across_track");
    uint32_t naCross = acrossDim.getSize();
    NcDim alongDim = inAncillaryfile->getDim("bins_along_track");
    uint32_t nalong = alongDim.getSize();

    NcDim levsDim = inAncillaryfile->getDim("levels");
    uint32_t nlev = levsDim.getSize();

    uint32_t naCrossOut = 25;

    float **out2d = allocate2d_float(nalong, naCrossOut);
    float ***out3d = allocate3d_float(nlev, nalong, naCrossOut);

    int binsAlongTrack = nalong;

    cout << "Creating ancillary file" << endl;
    NcFile *ncOutput;
    ncOutput = new NcFile(outAncillaryFilename, NcFile::replace);

    string dateCreated = string(unix2isodate(now(), 'G'));

    char buffer[1000];
    multimap<string, NcGroupAtt> AncillaryillaryGAtts = inAncillaryfile->getAtts();
    multimap<string, NcGroupAtt>::iterator itratt;
    for (itratt = AncillaryillaryGAtts.begin(); itratt != AncillaryillaryGAtts.end(); ++itratt) {
        string attrName = itratt->first;
        NcGroupAtt attin = itratt->second;

        NcType type = attin.getType();
        size_t attlen = attin.getAttLength();

        cout << itratt->first << " " << attlen << endl;

        if (attrName.compare("date_created") == 0) {
            ncOutput->putAtt("date_created", dateCreated);
        } else if (attrName.compare("product_name") == 0) {
            ncOutput->putAtt("product_name", basename(outAncillaryFilename));
        } else if (attrName.compare("title") == 0) {
            ncOutput->putAtt("title", "SPEXone L1C ancillary file");
        } else {
            attin.getValues((void *)buffer);
            ncOutput->putAtt(attrName, type, attlen, buffer);
        }
    }

    // Create netCDF dimensions
    NcDim rows = ncOutput->addDim("bins_along_track", binsAlongTrack);
    NcDim cols = ncOutput->addDim("bins_across_track", naCrossOut);
    NcDim levs = ncOutput->addDim("levels", nlev);

    vector<NcDim> dims;
    vector<size_t> start, count;

    multimap<string, NcVar> AncillaryillaryVars = inAncillaryfile->getVars();

    multimap<string, NcVar>::iterator itr;
    for (itr = AncillaryillaryVars.begin(); itr != AncillaryillaryVars.end(); ++itr) {
        string fieldName = itr->first;
        varin = itr->second;

        cout << itr->first << endl;

        int ndims = varin.getDimCount();

        if (ndims == 2) {
            start.clear();
            start.push_back(0);
            start.push_back((naCross - naCrossOut) / 2);

            count.clear();
            count.push_back(nalong);
            count.push_back(naCrossOut);

            varin.getVar(start, count, &out2d[0][0]);

            dims.clear();
            dims.push_back(rows);
            dims.push_back(cols);

            varout = ncOutput->addVar(fieldName, ncFloat, dims);
            copyVarAtts(&varin, &varout);
            varout.putVar(&out2d[0][0]);
        } else if (ndims == 3) {
            start.clear();
            start.push_back(0);
            start.push_back(0);
            start.push_back((naCross - naCrossOut) / 2);

            count.clear();
            count.push_back(nlev);
            count.push_back(nalong);
            count.push_back(naCrossOut);

            varin.getVar(start, count, &out3d[0][0][0]);

            dims.clear();
            dims.push_back(levs);
            dims.push_back(rows);
            dims.push_back(cols);

            varout = ncOutput->addVar(fieldName, ncFloat, dims);
            copyVarAtts(&varin, &varout);
            varout.putVar(&out3d[0][0][0]);
        }
    }

    free2d_float(out2d);
    free3d_float(out3d);

    return 0;
}
