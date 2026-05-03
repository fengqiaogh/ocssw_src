#include "OutFileHdf4.h"
#include <hntdefs.h>
#include <mfhdf.h>
#include <genutils.h>
#include <timeutils.h>

using namespace std;

OutFileHdf4::OutFileHdf4() : OutFile() {
    fileData = nullptr;
    sdfid = -1;
    sdsid = -1;
    quality_sdsid = -1;
    hdfDataType = DFNT_FLOAT32;
}

OutFileHdf4::~OutFileHdf4() {
    if (fileData)
        free(fileData);
}

void OutFileHdf4::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = nullptr;
}

bool OutFileHdf4::open() {
    const char* tmpStr;
    float tmpFloat;
    int32_t tmpInt32;

    currentLine = 0;

    sdfid = SDstart(fileName.c_str(), DFACC_CREATE);

    if (sdfid < 0) {
        printf("-E- Could not create HDF4 file %s\n", fileName.c_str());
        exit(EXIT_FAILURE);
    }

    get_time(metaData->ptime);

    string prodName;
    size_t pos = fileName.find_last_of('/');
    if (pos == string::npos)
        prodName = fileName;
    else
        prodName = fileName.substr(pos + 1);
    DPTB(SDsetattr(sdfid, "Product Name", DFNT_CHAR, prodName.size() + 1, (VOIDP)prodName.c_str()));
    DPTB(SDsetattr(sdfid, "Sensor Name", DFNT_CHAR, strlen(metaData->sensor_name) + 1,
                   (VOIDP)metaData->sensor_name));
    DPTB(SDsetattr(sdfid, "Sensor", DFNT_CHAR, strlen(metaData->sensor) + 1, (VOIDP)metaData->sensor));
    DPTB(SDsetattr(sdfid, "Title", DFNT_CHAR, strlen(metaData->title) + 1, (VOIDP)metaData->title));
    DPTB(SDsetattr(sdfid, "Data Center", DFNT_CHAR, strlen(metaData->data_center) + 1,
                   (VOIDP)metaData->data_center));
    DPTB(
        SDsetattr(sdfid, "Station Name", DFNT_CHAR, strlen(metaData->station) + 1, (VOIDP)metaData->station));
    DPTB(SDsetattr(sdfid, "Station Latitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->station_lat));
    DPTB(SDsetattr(sdfid, "Station Longitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->station_lon));
    DPTB(SDsetattr(sdfid, "Mission", DFNT_CHAR, strlen(metaData->mission) + 1, (VOIDP)metaData->mission));
    DPTB(SDsetattr(sdfid, "Mission Characteristics", DFNT_CHAR, strlen(metaData->mission_char) + 1,
                   (VOIDP)metaData->mission_char));
    DPTB(SDsetattr(sdfid, "Sensor Characteristics", DFNT_CHAR, strlen(metaData->sensor_char) + 1,
                   (VOIDP)metaData->sensor_char));
    DPTB(SDsetattr(sdfid, "Product Type", DFNT_CHAR, strlen(metaData->prod_type) + 1,
                   (VOIDP)metaData->prod_type));
    DPTB(SDsetattr(sdfid, "Processing Version", DFNT_CHAR, strlen(metaData->pversion) + 1,
                   (VOIDP)metaData->pversion));
    DPTB(SDsetattr(sdfid, "Software Name", DFNT_CHAR, strlen(metaData->soft_name) + 1,
                   (VOIDP)metaData->soft_name));
    DPTB(SDsetattr(sdfid, "Software Version", DFNT_CHAR, strlen(metaData->soft_ver) + 1,
                   (VOIDP)metaData->soft_ver));
    DPTB(SDsetattr(sdfid, "Processing Time", DFNT_CHAR, strlen(metaData->ptime) + 1, (VOIDP)metaData->ptime));
    DPTB(SDsetattr(sdfid, "Input Files", DFNT_CHAR, strlen(metaData->infiles) + 1, (VOIDP)metaData->infiles));
    DPTB(SDsetattr(sdfid, "Processing Control", DFNT_CHAR, strlen(metaData->proc_con) + 1,
                   (VOIDP)metaData->proc_con));
    DPTB(SDsetattr(sdfid, "Input Parameters", DFNT_CHAR, strlen(metaData->input_parms) + 1,
                   (VOIDP)metaData->input_parms));
    DPTB(SDsetattr(sdfid, "L2 Flag Names", DFNT_CHAR, strlen(metaData->flag_names) + 1,
                   (VOIDP)metaData->flag_names));

    short syear, sday, eyear, eday;
    double ssec, esec;
    int32_t smsec, emsec;
    unix2yds(metaData->startTime, &syear, &sday, &ssec);
    smsec = (int32_t)(ssec * 1000.0);
    unix2yds(metaData->endTime, &eyear, &eday, &esec);
    emsec = (int32_t)(esec * 1000.0);
    DPTB(SDsetattr(sdfid, "Period Start Year", DFNT_INT16, 1, (VOIDP)&syear));
    DPTB(SDsetattr(sdfid, "Period Start Day", DFNT_INT16, 1, (VOIDP)&sday));
    DPTB(SDsetattr(sdfid, "Period End Year", DFNT_INT16, 1, (VOIDP)&eyear));
    DPTB(SDsetattr(sdfid, "Period End Day", DFNT_INT16, 1, (VOIDP)&eday));
    tmpStr = ydhmsf(metaData->startTime, 'G');
    DPTB(SDsetattr(sdfid, "Start Time", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    tmpStr = ydhmsf(metaData->endTime, 'G');
    DPTB(SDsetattr(sdfid, "End Time", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    DPTB(SDsetattr(sdfid, "Start Year", DFNT_INT16, 1, (VOIDP)&syear));
    DPTB(SDsetattr(sdfid, "Start Day", DFNT_INT16, 1, (VOIDP)&sday));
    DPTB(SDsetattr(sdfid, "Start Millisec", DFNT_INT32, 1, (VOIDP)&smsec));
    DPTB(SDsetattr(sdfid, "End Year", DFNT_INT16, 1, (VOIDP)&eyear));
    DPTB(SDsetattr(sdfid, "End Day", DFNT_INT16, 1, (VOIDP)&eday));
    DPTB(SDsetattr(sdfid, "End Millisec", DFNT_INT32, 1, (VOIDP)&emsec));

    DPTB(SDsetattr(sdfid, "Start Orbit", DFNT_INT32, 1, (VOIDP)&metaData->start_orb));
    DPTB(SDsetattr(sdfid, "End Orbit", DFNT_INT32, 1, (VOIDP)&metaData->end_orb));
    DPTB(SDsetattr(sdfid, "Orbit", DFNT_INT32, 1, (VOIDP)&metaData->orbit));
    DPTB(SDsetattr(sdfid, "Map Projection", DFNT_CHAR, mapProjection.size() + 1, mapProjection.c_str()));
    DPTB(SDsetattr(sdfid, "Latitude Units", DFNT_CHAR, 14, (VOIDP) "degrees North"));
    DPTB(SDsetattr(sdfid, "Longitude Units", DFNT_CHAR, 13, (VOIDP) "degrees East"));
    DPTB(SDsetattr(sdfid, "Northernmost Latitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->north));
    DPTB(SDsetattr(sdfid, "Southernmost Latitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->south));
    DPTB(SDsetattr(sdfid, "Westernmost Longitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->west));
    DPTB(SDsetattr(sdfid, "Easternmost Longitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->east));

    float latStep = (metaData->north - metaData->south) / (float)getHeight();
    DPTB(SDsetattr(sdfid, "Latitude Step", DFNT_FLOAT32, 1, (VOIDP)&latStep));
    float lonStep;
    if (metaData->east < metaData->west)
        lonStep = (360 + metaData->east - metaData->west) / (float)getWidth();
    else
        lonStep = (metaData->east - metaData->west) / (float)getWidth();
    DPTB(SDsetattr(sdfid, "Longitude Step", DFNT_FLOAT32, 1, (VOIDP)&lonStep));

    tmpFloat = metaData->south + latStep / 2.0;
    DPTB(SDsetattr(sdfid, "SW Point Latitude", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = metaData->west + lonStep / 2.0;
    DPTB(SDsetattr(sdfid, "SW Point Longitude", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));

    tmpInt32 = metaData->data_bins;
    DPTB(SDsetattr(sdfid, "Data Bins", DFNT_INT32, 1, (VOIDP)&tmpInt32));

    tmpInt32 = getHeight();
    DPTB(SDsetattr(sdfid, "Number of Lines", DFNT_INT32, 1, (VOIDP)&tmpInt32));
    tmpInt32 = getWidth();
    DPTB(SDsetattr(sdfid, "Number of Columns", DFNT_INT32, 1, (VOIDP)&tmpInt32));
    DPTB(SDsetattr(sdfid, "Parameter", DFNT_CHAR, strlen(productStuff[0]->productInfo->description) + 1,
                   (VOIDP)productStuff[0]->productInfo->description));
    DPTB(SDsetattr(sdfid, "Measure", DFNT_CHAR, 5, (VOIDP) "Mean"));
    DPTB(SDsetattr(sdfid, "Units", DFNT_CHAR, strlen(productStuff[0]->productInfo->units) + 1,
                   (VOIDP)productStuff[0]->productInfo->units));

    // we only use linear scaling for data storage
    tmpStr = "linear";
    DPTB(SDsetattr(sdfid, "Scaling", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    tmpStr = "(Slope*l3m_data) + Intercept = Parameter value";
    DPTB(SDsetattr(sdfid, "Scaling Equation", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));

    float tmpScale;
    float tmpOffset;
    const char* imageScalingApplied;

    switch (productStuff[0]->dataStorage) {
        case FLOAT_DS:
        case DOUBLE_DS:
            tmpScale = 1.0;
            tmpOffset = 0.0;
            imageScalingApplied = "No";
            break;
        default:
            tmpScale = productStuff[0]->scale;
            tmpOffset = productStuff[0]->offset;
            if (tmpScale == 1.0 && tmpOffset == 0.0)
                imageScalingApplied = "No";
            else
                imageScalingApplied = "Yes";
    }

    DPTB(SDsetattr(sdfid, "Slope", DFNT_FLOAT32, 1, (VOIDP)&tmpScale));
    DPTB(SDsetattr(sdfid, "Intercept", DFNT_FLOAT32, 1, (VOIDP)&tmpOffset));

    tmpFloat = productStuff[0]->productInfo->validMin;
    DPTB(SDsetattr(sdfid, "Data Minimum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = productStuff[0]->productInfo->validMax;
    DPTB(SDsetattr(sdfid, "Data Maximum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = productStuff[0]->productInfo->displayMin;
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Minimum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = productStuff[0]->productInfo->displayMax;
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Maximum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));

    tmpStr = strdup(productStuff[0]->productInfo->displayScale);
    upcase((char*)tmpStr);
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Type", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    free((void*)tmpStr);

    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Applied", DFNT_CHAR, strlen(imageScalingApplied) + 1,
                   (VOIDP)imageScalingApplied));
    DPTB(SDsetattr(sdfid, "_lastModified", DFNT_CHAR, strlen(metaData->ptime) + 1, (VOIDP)metaData->ptime));

    // delete file data
    if (fileData) {
        free(fileData);
        fileData = NULL;
    }

    int32_t dims[2];
    dims[0] = height;
    dims[1] = width;

    if (!strcmp(productStuff[0]->productInfo->dataType, "byte")) {
        hdfDataType = DFNT_INT8;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        int8_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(int8_t), "OutFileHdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "ubyte")) {
        hdfDataType = DFNT_UINT8;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        uint8_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(uint8_t), "OutFileHdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "short")) {
        hdfDataType = DFNT_INT16;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        int16_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(int16_t), "OutFileHdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "ushort")) {
        hdfDataType = DFNT_UINT16;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        uint16_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(uint16_t), "OutFileHdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "int")) {
        hdfDataType = DFNT_INT32;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        int32_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(int32_t), "OutFileHdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "uint")) {
        hdfDataType = DFNT_UINT32;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        uint32_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(uint32_t), "OutFileHdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "float")) {
        hdfDataType = DFNT_FLOAT32;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        float tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(float), "OutFileHdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "double")) {
        hdfDataType = DFNT_FLOAT64;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        double tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(double), "OutFileHdf4::open fileData");
    } else {
        printf("-E- Data type %s, not supported\n", productStuff[0]->productInfo->dataType);
        exit(EXIT_FAILURE);
    }

    tmpStr = "linear";
    DPTB(SDsetattr(sdsid, "Scaling", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    tmpStr = "(Slope*l3m_data) + Intercept = Parameter value";
    DPTB(SDsetattr(sdsid, "Scaling Equation", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    DPTB(SDsetattr(sdsid, "Slope", DFNT_FLOAT32, 1, (VOIDP)&tmpScale));
    DPTB(SDsetattr(sdsid, "Intercept", DFNT_FLOAT32, 1, (VOIDP)&tmpOffset));

    // create the SST quality data set
    if (qualityData) {
        quality_sdsid = SDcreate(sdfid, "l3m_qual", DFNT_UINT8, 2, dims);
        int32_t validRange[2];
        validRange[0] = 0;
        validRange[1] = 2;
        DPTB(SDsetattr(quality_sdsid, "valid_range", DFNT_INT32, 2, (VOIDP)validRange));
    }

    return true;
}

void OutFileHdf4::writeLine() {
    productStuff[0]->calcOutputLineVals(fileData);

    int32_t start[2];
    int32_t count[2];

    start[0] = currentLine;
    start[1] = 0;
    count[0] = 1;
    count[1] = width;

    if ((SDwritedata(sdsid, start, NULL, count, (VOIDP)fileData)) < 0) {
        printf("\n-E- OutFileHdf4::writeLine(): SDwritedata unsuccessful\n");
        exit(EXIT_FAILURE);
    }

    if (qualityData) {
        if ((SDwritedata(quality_sdsid, start, NULL, count, (VOIDP)qualityData)) < 0) {
            printf("\n-E- OutFileHdf4::writeLine(): SDwritedata unsuccessful\n");
            exit(EXIT_FAILURE);
        }
    }

    currentLine++;
}


bool OutFileHdf4::close() {
    if (fileData) {
        free(fileData);
        fileData = NULL;
    }

    if (metaData) {
        int32_t tmpInt = metaData->data_bins;
        DPTB(SDsetattr(sdfid, "Data Bins", DFNT_INT32, 1, (VOIDP)&tmpInt));
    }
    float tmpFloat = getFileMinVal();
    DPTB(SDsetattr(sdfid, "Data Minimum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = getFileMaxVal();
    DPTB(SDsetattr(sdfid, "Data Maximum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));

    if (qualityData) {
        SDendaccess(quality_sdsid);
        quality_sdsid = -1;
    }
    SDendaccess(sdsid);
    sdsid = -1;
    SDend(sdfid);
    sdfid = -1;

    /*-----------------------------------------------------------------------*/
    /*  write map_palette */

    // make sure a palette has been loaded
    if (red && green && blue) {
        uint8_t data[768];
        int j = 0;
        for (int i = 0; i < 256; i++) {
            data[j++] = red[i];
            data[j++] = green[i];
            data[j++] = blue[i];
        }

        if ((DFPaddpal(fileName.c_str(), (VOIDP)data)) < 0) {
            printf("-E- OutFileHdf4::close - Error writing map_palette.\n");
            return false;
        }

        uint16_t pal_ref;
        if ((pal_ref = DFPlastref()) > 0) {
            if ((DFANputlabel(fileName.c_str(), DFTAG_IP8, pal_ref, (char*)"palette")) < 0) {
                printf("-E- OutFileHdf4::close - Error writing palette label\n");
                return false;
            }
        }
    }
    return true;
}

int32_t OutFileHdf4::addProduct(productInfo_t* productInfo, bool applyMask, const ProductL3Attributes & productAttr) {
    return addProductNonDisplay(productInfo, productAttr);
}
