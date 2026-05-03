#include "OutFileNetcdf4.h"
#include "l3mapgen.h"
#include <timeutils.h>
#include <nc4utils.h>
#include <L3FileSMI.h>
#include <vector>
#include <mtl_geometry.h>
#include "Wave3DParsing.hpp"
#include "CacheSize.h"
using namespace netCDF;
using namespace std;

std::string OutFileNetcdf4::removeWvFromLongName(const std::string& longNameWithWv) {
    vector<string> key_words;
    const string delim = " ";
    boost::split(key_words, longNameWithWv, boost::is_any_of(delim));
    const string start = "at";
    const string end = "nm";
    string out;
    string wv_prefix;
    stack<string> stack;
    for (const auto& word : key_words) {
        if (stack.empty()) {
            if (word != "at") {
                if (!out.empty())
                    out += delim;
                out += word;
            } else {
                stack.push(word);
            }
        } else {
            if (word != "nm") {
                stack.push(word);
            } else {
                while (!stack.empty()) {
                    if (!wv_prefix.empty())
                        wv_prefix += delim;
                    wv_prefix += stack.top();
                    stack.pop();
                }
            }
        }
    }
    boost::trim(out);
    return out;
}

OutFileNetcdf4::OutFileNetcdf4() : OutFile() {
    ncFile = nullptr;
    fileData = nullptr;
}

OutFileNetcdf4::~OutFileNetcdf4() {
    if (fileData) {
        free(fileData);
        fileData = nullptr;
    }
    delete ncFile;
    ncFile = nullptr;
}


void OutFileNetcdf4::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData) {
        free(fileData);
        fileData = nullptr;
    }
}

bool OutFileNetcdf4::open() {
    const char* tmpStr;
    currentLine = 0;

    try {
        // open file
        ncFile = new NcFile(fileName.c_str(), NcFile::replace);

        // write global metadata
        string prodName;
        std::size_t pos = fileName.find_last_of('/');
        if (pos == string::npos)
            prodName = fileName;
        else
            prodName = fileName.substr(pos + 1);
        ncFile->putAtt("product_name", prodName.c_str());

        ncFile->putAtt("instrument", metaData->sensor);
        string source = "satellite observations from ";
        source += metaData->sensor;
        ncFile->putAtt("title", metaData->title);
        ncFile->putAtt("project", PROJECT);

        if (strcmp(metaData->mission, "") != 0) {
            ncFile->putAtt("platform", metaData->mission);
            source += "-";
            source += metaData->mission;
        }
        ncFile->putAtt("source", source);

        ncFile->putAtt("temporal_range", metaData->prod_type);
        ncFile->putAtt("processing_version", metaData->pversion);

        strcpy(metaData->ptime, unix2isodate(time(nullptr), 'G'));
        ncFile->putAtt("date_created", metaData->ptime);
        ncFile->putAtt("history", metaData->proc_con);

        ncFile->putAtt("l2_flag_names", metaData->flag_names);

        ncFile->putAtt("time_coverage_start", unix2isodate(metaData->startTime, 'G'));
        ncFile->putAtt("time_coverage_end", unix2isodate(metaData->endTime, 'G'));
        ncFile->putAtt("start_orbit_number", ncInt, metaData->start_orb);
        ncFile->putAtt("end_orbit_number", ncInt, metaData->end_orb);
        ncFile->putAtt("map_projection", mapProjection.c_str());

        ncFile->putAtt("latitude_units", metaData->lat_units);
        ncFile->putAtt("longitude_units", metaData->lon_units);

        ncFile->putAtt("northernmost_latitude", ncFloat, metaData->north);
        ncFile->putAtt("southernmost_latitude", ncFloat, metaData->south);
        ncFile->putAtt("westernmost_longitude", ncFloat, metaData->west);
        ncFile->putAtt("easternmost_longitude", ncFloat, metaData->east);
        ncFile->putAtt("geospatial_lat_max", ncFloat, metaData->north);
        ncFile->putAtt("geospatial_lat_min", ncFloat, metaData->south);
        ncFile->putAtt("geospatial_lon_max", ncFloat, metaData->east);
        ncFile->putAtt("geospatial_lon_min", ncFloat, metaData->west);

        double latStep = (metaData->north - metaData->south) / (float)getHeight();
        double lonStep;
        if (metaData->east < metaData->west)
            lonStep = (360 + metaData->east - metaData->west) / (float)getWidth();
        else
            lonStep = (metaData->east - metaData->west) / (float)getWidth();
        ncFile->putAtt("latitude_step", ncFloat, latStep);
        ncFile->putAtt("longitude_step", ncFloat, lonStep);
        ncFile->putAtt("sw_point_latitude", ncFloat, metaData->south + latStep / 2.0);
        ncFile->putAtt("sw_point_longitude", ncFloat, metaData->west + lonStep / 2.0);

        string resolutionStr = resolution2string(resolution);
        ncFile->putAtt("spatialResolution", resolutionStr);
        ncFile->putAtt("geospatial_lon_resolution", resolutionStr);
        ncFile->putAtt("geospatial_lat_resolution", resolutionStr);
        ncFile->putAtt("geospatial_lat_units", metaData->lat_units);
        ncFile->putAtt("geospatial_lon_units", metaData->lon_units);

        ncFile->putAtt("number_of_lines", ncInt, height);
        ncFile->putAtt("number_of_columns", ncInt, width);
        ncFile->putAtt("measure", "Mean");  // FIXME

        ncFile->putAtt("suggested_image_scaling_minimum", ncFloat, productStuff[0]->productInfo->displayMin);
        ncFile->putAtt("suggested_image_scaling_maximum", ncFloat, productStuff[0]->productInfo->displayMax);
        if (!strcmp(productStuff[0]->productInfo->displayScale, "log"))
            tmpStr = "LOG";
        else if (!strcmp(productStuff[0]->productInfo->displayScale, "arctan"))
            tmpStr = "ATAN";
        else
            tmpStr = "LINEAR";
        ncFile->putAtt("suggested_image_scaling_type", tmpStr);
        ncFile->putAtt("suggested_image_scaling_applied", "No");

        ncFile->putAtt("_lastModified", metaData->ptime);
        ncFile->putAtt("Conventions", "CF-1.6 ACDD-1.3");
        ncFile->putAtt("institution", INSTITUTION);
        ncFile->putAtt("standard_name_vocabulary", STDNAME_VOCABULARY);
        ncFile->putAtt("naming_authority", NAMING_AUTHORITY);

        // create id
        string id = metaData->pversion;
        if (id == "Unspecified") {
            id = "L3/";
        } else {
            id += "/L3/";
        }
        id += metaData->product_name;
        ncFile->putAtt("id", id);

        ncFile->putAtt("license", LICENSE);
        ncFile->putAtt("creator_name", CREATOR_NAME);
        ncFile->putAtt("publisher_name", PUBLISHER_NAME);
        ncFile->putAtt("creator_email", CREATOR_EMAIL);
        ncFile->putAtt("publisher_email", PUBLISHER_EMAIL);
        ncFile->putAtt("creator_url", CREATOR_URL);
        ncFile->putAtt("publisher_url", PUBLISHER_URL);

        ncFile->putAtt("processing_level", "L3 Mapped");
        ncFile->putAtt("cdm_data_type", "grid");  // deprecated?
        if (proj4String.length() > 0)
            ncFile->putAtt("proj4_string", proj4String);

        // Some missions have DOIs
        if (clo_isSet(optionList, "doi")) {
            ncFile->putAtt("identifier_product_doi_authority", "http://dx.doi.org");
            ncFile->putAtt("identifier_product_doi", clo_getString(optionList, "doi"));
        }
        if (clo_isSet(optionList, "suite")) {
            const char* keywordStr = getGCMDKeywords(clo_getString(optionList, "suite"));
            if (keywordStr) {
                ncFile->putAtt("keywords", keywordStr);
                ncFile->putAtt("keywords_vocabulary", KEYWORDS_VOCABULARY);
            }
        }  // suite was set

        // create group: processing_control
        NcGroup grp1 = ncFile->addGroup("processing_control");
        grp1.putAtt("software_name", metaData->soft_name);
        grp1.putAtt("software_version", metaData->soft_ver);
        if ((tmpStr = strrchr(metaData->infiles, '/')) != nullptr)
            tmpStr++;
        else
            tmpStr = metaData->infiles;
        grp1.putAtt("input_sources", tmpStr);
        grp1.putAtt("l2_flag_names", metaData->flag_names);

        // create sub-group: input_parameters
        NcGroup grp2 = grp1.addGroup("input_parameters");
        char buf[2048];
        char* end_str;
        char* token = strtok_r(metaData->input_parms, "|", &end_str);
        while (token != nullptr) {
            char* end_token;
            strcpy(buf, token);
            char* name = strtok_r(token, "=", &end_token);
            for (uint32_t i = 0; i < strlen(name); i++) {
                if (name[i] == ' ') {
                    name[i] = 0;
                    break;
                }
            }
            char* val = strtok_r(nullptr, "|", &end_token);
            if (val == nullptr)
                val = "";
            strcpy(buf, val);
            grp2.putAtt(name, buf);
            token = strtok_r(nullptr, "|", &end_str);
        }

        // Define dimensions and coordinate variables
        vector<NcDim> dimIds;
        string coordinates;

        if (mapProjection == "Equidistant Cylindrical"  // SMI
            || mapProjection == "PlateCarree" || (proj4String.find("+proj=eqc") != string::npos)) {
            dimIds.push_back(ncFile->addDim("lat", height));
            dimIds.push_back(ncFile->addDim("lon", width));

            if (fullLatLon != LAT_LON_OFF)
                fullLatLon = LAT_LON_1D;
        } else {
            dimIds.push_back(ncFile->addDim("y", height));
            dimIds.push_back(ncFile->addDim("x", width));
            if (fullLatLon != LAT_LON_OFF) {
                fullLatLon = LAT_LON_2D;
                coordinates = "lat lon";
            }
        }

        // Define variables
        std::size_t dataSize = 1;
        vector<NcDim> dimIds_2D;
        vector<NcDim> dimIds_3D;
        copy(dimIds.begin(), dimIds.end(), back_inserter(dimIds_2D));

        if (!getWv3dName2dTo3dExpansion().empty()) {
            copy(dimIds.begin(), dimIds.end(), back_inserter(dimIds_3D));
            size_t number_of_wave = getLenWv3d();
            if (number_of_wave) {
                dimIds_3D.push_back(
                    ncFile->addDim("wavelength", number_of_wave));  // add the extra wavelength dimension

                vector<NcDim> dimIdswv3d = {dimIds_3D[2]};
                NcVar var = ncFile->addVar("wavelength", ncFloat, dimIdswv3d);
                var.putAtt("long_name", "wavelengths");
                var.putAtt("units", "nm");
                var.putAtt("_FillValue", ncFloat, -32767);
                var.putAtt("valid_min", ncFloat, 0);
                var.putAtt("valid_max", ncFloat, 20000);
                std::vector<float> data = (*getWv3dName2dTo3dExpansion().begin()).second;
                var.putVar(data.data());
            }
        }
        auto dimIds_to_set = &dimIds_2D;
        std::size_t count_3d_counter = 0;
        for (std::size_t prodNum = 0; prodNum < productStuff.size(); prodNum++) {
            // determine maximum product size
            NcType varDataType = getDataType(productStuff[prodNum]->dataStorage);
            dataSize = max(dataSize, varDataType.getSize());

            // create variable
            productInfo_t* pInfo = productStuff[prodNum]->productInfo;
            string temp_wave_name = pInfo->paramDesignator;
            string temp_prefix = pInfo->prefix;
            string temp_descp = pInfo->description;

            const auto wv3d_3d_name_to_2d_name = getWv3d3dNameTo2D();
            string full_prod_name = getProductNameFull(pInfo);
            const string suffix_prod = pInfo->suffix;
            if (!suffix_prod.empty()) {
                size_t pos = full_prod_name.find(suffix_prod);
                full_prod_name = full_prod_name.substr(0, pos);
            }
            if (wv3d_3d_name_to_2d_name.count(full_prod_name) > 0) {
                const auto prod_name_no_suffix = wv3d_3d_name_to_2d_name.at(full_prod_name);
                const auto prod_name = prod_name_no_suffix + suffix_prod;

                if (product3dAlreadySet.count(prod_name) > 0) {
                    product3dAlreadySet.at(prod_name)++;
                    slice2dInWv3d[prodNum] = product3dAlreadySet.at(prod_name);
                    index2d3d[prodNum] = index2d3d.at(product2dIndexesLastIndex.at(prod_name));
                    continue;
                } else {
                    product3dAlreadySet[prod_name] = 0;
                    slice2dInWv3d[prodNum] = 0;
                }
                product2dIndexesLastIndex[prod_name] = prodNum;
                dimIds_to_set = &dimIds_3D;
                if (pInfo->paramDesignator)
                    free(pInfo->paramDesignator);
                if (pInfo->prefix)
                    free(pInfo->prefix);
                const auto new_desc = removeWvFromLongName(pInfo->description);
                if (pInfo->description)
                    free(pInfo->description);
                pInfo->paramDesignator = strdup("none");
                pInfo->prefix = strdup(prod_name_no_suffix.c_str());
                pInfo->description = strdup(new_desc.c_str());
            } else {
                dimIds_to_set = &dimIds_2D;
            }
            NcVar var = createProduct(pInfo, varDataType, *dimIds_to_set);
            index2d3d[prodNum] = count_3d_counter;
            count_3d_counter++;
            // reset pInfo to default
            {
                if (wv3d_3d_name_to_2d_name.count(full_prod_name) > 0) {
                    free(pInfo->paramDesignator);
                    free(pInfo->prefix);
                    free(pInfo->description);
                    pInfo->paramDesignator = strdup(temp_wave_name.c_str());
                    pInfo->prefix = strdup(temp_prefix.c_str());
                    pInfo->description = strdup(temp_descp.c_str());
                }
            }
            initCompression(var);
            prodVars.push_back(var);

            // add display info
            var.putAtt("display_scale", pInfo->displayScale);
            var.putAtt("display_min", ncFloat, pInfo->displayMin);
            var.putAtt("display_max", ncFloat, pInfo->displayMax);

            // specify coordinates as needed
            if (coordinates.length() > 0)
                var.putAtt("coordinates", coordinates);

            // see if wavelength is defined
            auto& the2dMap = getWavelength2dMap();
            // const auto& item = the2dMap.find(full_prod_name);
            //if(item != the2dMap.end()) {
            if(the2dMap.count(full_prod_name)) {
                // float waveVal = the2dMap.at(full_prod_name);
                var.putAtt("wavelength", ncFloat, the2dMap.at(full_prod_name));
            }
        }

        if (fileData)
            free(fileData);
        fileData = allocateMemory(width * dataSize, "OutFileNetcdf4::open fileData");

        // add the quality variable
        if (qualityData) {
            productInfo_t* qInfo = allocateProductInfo();
            if (qualityName.empty()) {
                qualityName = (string) "qual_" + getProductNameFull(productStuff[0]->productInfo);
            }
            if (!findProductInfo(qualityName.c_str(), metaData->sensorID, qInfo)) {
                cerr << "-E- OutFileNetcdf4::open - cannot find product " << qualityName << " in product.xml"
                     << endl;
                exit(EXIT_FAILURE);
            }
            qInfo->fillValue = NC_FILL_UBYTE;
            qualVar = createProduct(qInfo, ncUbyte, dimIds);
            initCompression(qualVar);
            freeProductInfo(qInfo);
        }

        // add latitude and longitude variables

        productInfo_t* latInfo = allocateProductInfo();
        if (!findProductInfo("lat", metaData->sensorID, latInfo)) {
            cerr << "-E- OutFileNetcdf4::open - "
                    "cannot find \"lat\" in product.xml"
                 << endl;
            exit(EXIT_FAILURE);
        }
        productInfo_t* lonInfo = allocateProductInfo();
        if (!findProductInfo("lon", metaData->sensorID, lonInfo)) {
            cerr << "-E- OutFileNetcdf4::open - "
                    "cannot find \"lon\" in product.xml"
                 << endl;
            exit(EXIT_FAILURE);
        }

        if (fullLatLon == LAT_LON_1D) {
            float* latarray = (float*)allocateMemory(getHeight() * sizeof(float), "latarray");
            for (int i = 0; i < getHeight(); i++)
                latarray[i] = metaData->north - latStep * i - latStep / 2.0;
            vector<NcDim> latDim;
            latDim.push_back(dimIds[0]);
            NcVar lat = createProduct(latInfo, ncFloat, latDim);
            lat.putVar(latarray);
            free(latarray);

            float* lonarray = (float*)allocateMemory(getWidth() * sizeof(float), "lonarray");
            for (int i = 0; i < getWidth(); i++) {
                lonarray[i] = metaData->west + lonStep * i + lonStep / 2.0;
                if (lonarray[i] > 180.0)
                    lonarray[i] -= 360.0;
            }
            vector<NcDim> lonDim;
            lonDim.push_back(dimIds[1]);
            NcVar lon = createProduct(lonInfo, ncFloat, lonDim);
            lon.putVar(lonarray);
            free(lonarray);
        } else if (fullLatLon == LAT_LON_2D) {
            latVar = createProduct(latInfo, ncFloat, dimIds);
            lonVar = createProduct(lonInfo, ncFloat, dimIds);
            initCompression(latVar);
            initCompression(lonVar);
        }

        freeProductInfo(latInfo);
        freeProductInfo(lonInfo);

        // store palette
        if (red && green && blue) {
            uint8_t data[768];
            int j = 0;
            for (int i = 0; i < 256; i++) {
                data[j++] = red[i];
                data[j++] = green[i];
                data[j++] = blue[i];
            }
            vector<NcDim> dimIds;
            dimIds.push_back(ncFile->addDim("rgb", 3));
            dimIds.push_back(ncFile->addDim("eightbitcolor", 256));
            NcVar var = ncFile->addVar("palette", ncUbyte, dimIds);
            var.putVar(data);
        }

    } catch (exception const& e) {
        cerr << "Exception: " << e.what() << endl;
    }

    return true;
}

void OutFileNetcdf4::writeLine() {
    vector<std::size_t> start(2);
    vector<std::size_t> count(2);
    start[0] = currentLine;
    start[1] = 0;
    count[0] = 1;
    count[1] = width;
    vector<std::size_t> start_3D = {start[0], start[1], 0};
    vector<std::size_t> count_3D = {count[0], count[1], 1};
    // calculate and write one line of data for each product
    for (std::size_t prodNum = 0; prodNum < productStuff.size(); prodNum++) {
        productStuff[prodNum]->calcOutputLineVals(fileData);
        const auto index_in_products = index2d3d.at(prodNum);
        if (slice2dInWv3d.count(prodNum) == 0) {
            prodVars[index_in_products].putVar(start, count, fileData);
        } else {
            const auto slice_to_put = slice2dInWv3d.at(prodNum);
            start_3D[2] = slice_to_put;
            prodVars[index_in_products].putVar(start_3D, count_3D, fileData);
        }
    }

    // write quality data for this line
    if (qualityData) {
        qualVar.putVar(start, count, qualityData);
    }

    // write lat, lon for this line
    if (fullLatLon == LAT_LON_2D) {
        if (latData)
            latVar.putVar(start, count, latData);
        if (lonData)
            lonVar.putVar(start, count, lonData);
    }

    currentLine++;
}


bool OutFileNetcdf4::close() {
    if (metaData) {
        ncFile->putAtt("data_bins", ncInt64, metaData->data_bins);
        if (metaData->data_bins <= 0) {
            ncFile->putAtt("data_minimum", ncFloat, 0.0);
            ncFile->putAtt("data_maximum", ncFloat, 0.0);
        } else {
            ncFile->putAtt("data_minimum", ncFloat, (float)getFileMinVal());
            ncFile->putAtt("data_maximum", ncFloat, (float)getFileMaxVal());
        }
    }
    ncFile->close();
    prodVars.clear();

    if (fileData) {
        free(fileData);
        fileData = nullptr;
    }
    if (qualityData) {
        free(qualityData);
        qualityData = nullptr;
    }

    return true;
}

int32_t OutFileNetcdf4::addProduct(productInfo_t* productInfo, bool applyMask, const ProductL3Attributes & productAttr) {
    return addProductNonDisplay(productInfo, productAttr);
}

void OutFileNetcdf4::initCompression(NcVar var) {
    // if deflate level below 0 then don't deflate
    if (getDeflate() < 1)
        return;

    vector<NcDim> dims = var.getDims();
    vector<std::size_t> chunksize;
    // set chunk size along  first dimension, lines
    const size_t chunk_size_lines = get_chunk_along_lines();

    switch (dims.size()) {
        case 1:
            chunksize.push_back(MIN(dims[0].getSize(), 2048));
            break;
        case 2:
            chunksize.push_back(MIN(dims[0].getSize(), chunk_size_lines));
            chunksize.push_back(MIN(dims[1].getSize(), 1024));
            break;
        case 3:
            chunksize.push_back(MIN(dims[0].getSize(), chunk_size_lines));
            chunksize.push_back(MIN(dims[1].getSize(), 1024));
            chunksize.push_back(MIN(dims[2].getSize(), 8));
            break;
        default:
            printf("-E- OutFileNetcdf4::initCompression - bad rank for variable %s\n",
                   var.getName().c_str());
            exit(EXIT_FAILURE);
    }

    var.setChunking(NcVar::nc_CHUNKED, chunksize);
    var.setCompression(true, true, getDeflate());
}

NcVar OutFileNetcdf4::createProduct(productInfo_t* pInfo, const NcType& ncType,
                                     const vector<netCDF::NcDim> ncDim) {
    // create variable
    NcVar var = ncFile->addVar(getProductNameFull(pInfo), ncType, ncDim);

    // add standard metadata
    if (pInfo->description)
        var.putAtt("long_name", pInfo->description);
    if (pInfo->scaleFactor != 1.0 || pInfo->addOffset != 0.0) {
        var.putAtt("scale_factor", ncFloat, pInfo->scaleFactor);
        var.putAtt("add_offset", ncFloat, pInfo->addOffset);
    }
    if (pInfo->units != nullptr && strcmp(pInfo->units, "dimensionless") != 0)
        var.putAtt("units", pInfo->units);
    if (pInfo->standardName != nullptr && strcmp(pInfo->standardName, "") != 0)
        var.putAtt("standard_name", pInfo->standardName);
    try {
        if (pInfo->fillValue)
            var.putAtt("_FillValue", ncType, pInfo->fillValue);
    } catch (exception const& e) {
        cerr << "FillValue exception: " << e.what() << endl;
    }
    if (pInfo->validMin != PRODUCT_DEFAULT_validMin || pInfo->validMax != PRODUCT_DEFAULT_validMax) {
        var.putAtt("valid_min", ncType, (pInfo->validMin - pInfo->addOffset) / pInfo->scaleFactor);
        var.putAtt("valid_max", ncType, (pInfo->validMax - pInfo->addOffset) / pInfo->scaleFactor);
    }
    if (pInfo->reference != nullptr && strcmp(pInfo->reference, "") != 0)
        var.putAtt("reference", pInfo->reference);
    if (pInfo->comment != nullptr && strcmp(pInfo->comment, "") != 0)
        var.putAtt("comment", pInfo->comment);

    return var;
}


NcType OutFileNetcdf4::getDataType(DataStorage dataStorage) {
    switch (dataStorage) {
        case BYTE_DS:
            return NC_BYTE;  // ncByte
        case UBYTE_DS:
            return NC_UBYTE;  // ncUbyte
        case SHORT_DS:
            return NC_SHORT;  // ncShort
        case USHORT_DS:
            return NC_USHORT;  // ncUshort
        case INT_DS:
            return NC_INT;  // ncInt
        case UINT_DS:
            return NC_UINT;  // ncUint
        case FLOAT_DS:
            return NC_FLOAT;  // ncFloat
        case DOUBLE_DS:
            return NC_DOUBLE;  // ncDouble
        default:
            cerr << "-E- OutFileNetcdf4::getDataType - illegal data storage "
                    "type = "
                 << dataStorage << endl;
            exit(EXIT_FAILURE);
    }
}
