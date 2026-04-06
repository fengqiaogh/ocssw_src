#include "OutFileTiff.h"
#include "OutFile.h"
#include <stdexcept>
#include <string>
#include <geo_normalize.h>
#include <regex>
#include <genutils.h>

using namespace std;

OutFileTiff::~OutFileTiff() {
    close();
}

bool OutFileTiff::open() {
    currentLine = 0;

    // open TIFF file
    tiff = XTIFFOpen(fileName.c_str(), "w");
    if (tiff == nullptr) {
        cerr << "-E- Could not open outgoing TIFF image" << endl;
        exit(EXIT_FAILURE);
    }

    // extend the TIFF tags to write pversion
    ttag_t TIFFTAG_GDALMetadata = 42112;
    static const TIFFFieldInfo xtiffFieldInfo[] = {
        {TIFFTAG_GDALMetadata, -1, -1, TIFF_ASCII, FIELD_CUSTOM, true, false, "GDALMetadata"}};
    TIFFMergeFieldInfo(tiff, xtiffFieldInfo, 1);
    string tagVal = "<GDALMetadata>\n  <Item name=\"OBPG_version\">";
    tagVal += metaData->pversion;
    tagVal += "</Item>\n";
    tagVal += "</GDALMetadata>\n";
    TIFFSetField(tiff, TIFFTAG_GDALMetadata, tagVal.c_str());

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(tiff, TIFFTAG_PREDICTOR, PREDICTOR_NONE);
    TIFFSetField(tiff, TIFFTAG_ROWSPERSTRIP, height);

    // open GeoTIFF interface
    bool hasGeotiffInfo = false;
    gtif = GTIFNew(tiff);
    if (gtif == nullptr) {
        cerr << "-E- Could not create geoTIFF structure" << endl;
        exit(EXIT_FAILURE);
    }

    // define GeoTIFF keys for lat/lon projection
    if (mapProjection == "Equidistant Cylindrical" || mapProjection == "PlateCarree") {
        double north, south, east, west;
        if (metaData) {
            north = metaData->north;
            south = metaData->south;
            east = metaData->east;
            west = metaData->west;
        } else {
            north = 90.0;
            south = -90.0;
            east = 180.0;
            west = -180.0;
        }

        // pixel width, height in degrees
        pixscale[0] = (east - west) / width;
        pixscale[1] = (north - south) / height;

        // top left corner pixel lat, lon
        for (int i = 0; i < 6; i++)
            tiepoints[i] = 0;
        tiepoints[3] = west;
        tiepoints[4] = north;

        GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelGeographic);
        GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);
        hasGeotiffInfo = true;
    }
    // otherwise, parse the proj4 string
    else
        hasGeotiffInfo = (GTIFSetFromProj4(gtif, proj4String.c_str()));

    if (!hasGeotiffInfo) {
        if (proj4String.find("EPSG:") != string::npos)
            hasGeotiffInfo = true;
    }

    if (hasGeotiffInfo) {
        // write GeoTIFF keys
        GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
        static const string& epsg_regex_str{"EPSG:(\\d+)\\b"};
        static const regex epsg_regex{epsg_regex_str};
        smatch matches;
        if (regex_search(proj4String, matches, epsg_regex)) {
            unsigned short epsg = stoi(matches[1]);
            GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, epsg);
        }
        GTIFWriteKeys(gtif);

        // write GeoTIFF tags in map units
        TIFFSetField(tiff, GTIFF_PIXELSCALE, 3, pixscale);
        TIFFSetField(tiff, GTIFF_TIEPOINTS, 6, tiepoints);
    }

    // define colormap
    setTiffColor();
    return true;
}

bool OutFileTiff::close() {
    if (gtif) {
        GTIFFree(gtif);
        XTIFFClose(tiff);
        tiff = nullptr;
        gtif = nullptr;
    }
    return true;
}