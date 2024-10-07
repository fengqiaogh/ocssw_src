/************************************************************
*               Sentinel 2A MSI L1C PROCESSING              *
*                        By Jiaying He                      *
*               With modifications for seamless integration *
*                 in SeaDAS by Sudipta Sarkar SSAI          *
************************************************************/

/* Include head files */
#include "l1_msi.h"

#include "l1.h"
#include "jplaeriallib.h"

#include <libnav.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "l1_msi_private.h"
#include <cmath>
#include <algorithm>
#include <pugixml.hpp>

using namespace std;
using namespace pugi;

typedef enum msiBandIdx
{
    B01,B02,B03,B04,B05,B06,B07,B08,B8A,B09,B11,B12,B10
} msiBandIdx;

const static struct {
    msiBandIdx  val;
    const char *str;
} conversion [] = {
    {B01, "B01"},
    {B02, "B02"},
    {B03, "B03"},
    {B04, "B04"},
    {B05, "B05"},
    {B06, "B06"},
    {B07, "B07"},
    {B08, "B08"},
    {B8A, "B8A"},
    {B09, "B09"},
    {B11, "B11"},
    {B12, "B12"},
    {B10, "B10"}
};

void resample_msi(opj_image_t* image, filehandle* file, int recnum, int srcRes, int destRes);
int decodeMSI(filehandle *file, int32_t bandIdx, int32_t recnum);
void interpGPSpos(l1str *l1rec, double* pos, int detector, int band);
int inDetector (msi_t *data , float lat, float lon);
void interpViewAngles(l1str* l1rec, int pixel, int scan, int band, float *senz, float *sena);

void error_callback(const char *msg, void *client_data);
void warning_callback(const char *msg, void *client_data);
void info_callback(const char *msg, void *client_data);


msiBandIdx str2enum (const char *str)
{
     uint32_t j;
     for (j = 0;  j < sizeof (conversion) / sizeof (conversion[0]);  ++j)
         if (!strcmp (str, conversion[j].str))
             return conversion[j].val;

     fprintf(stderr,"-E- %s line %d: something is seriously wrong in Denmark...\n",
            __FILE__,__LINE__);
     exit(EXIT_FAILURE);
}

int inDetector (msi_t *data , float lat, float lon) {
    int detector = 0;
    PJ_COORD c, c_out;

    // set default z and t
    c.xyzt.x = 0.0;
    c.xyzt.t = HUGE_VAL;

    c.xy.x = lon;
    c.xy.y = lat;
    c_out = proj_trans(data->pj, PJ_INV, c);
    Point_t p(c_out.xy.x, c_out.xy.y);

    for (detector = 0; detector < numDetectors; detector++){
          if (boost::geometry::within(p, data->detectorFootprints[detector]))
              break;
    }
    if (detector == numDetectors)
        detector = -1;
    
    return detector;
}

void interpGPSpos(l1str* l1rec, double* pos, int detector, int band){
    msi_t *data = (msi_t*) l1rec->l1file->private_data;
    int i;
    int nelem = data->num_gps;
    // use GSL...why reinvent the wheel
    double pixelTime = l1rec->scantime;
    if (detector >= 0)
        pixelTime += data->detectorDeltaTime[detector][band];
           
    for(i = 0; i< 3; i++) {
        gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, nelem);
        gsl_interp_init(interpolation, data->gpstime, data->position[i], nelem);
        gsl_interp_accel *accelerator = gsl_interp_accel_alloc();
        pos[i]= gsl_interp_eval(interpolation, data->gpstime, data->position[i], l1rec->scantime, accelerator);
        gsl_interp_free(interpolation);
        gsl_interp_accel_free (accelerator);
    }
}

void interpViewAngles(msi_t* data, int pixel, int scan, int band, float *senz, float *sena) {

    float angles[2];

    int  nelem = 23;
    const gsl_interp2d_type *T = gsl_interp2d_bicubic;

    double *tieZenith = data->sensorZenith[band];
    double *tieAzimuth = data->sensorAzimuth[band];

    double tieX[nelem];
    double tieY[nelem];
    double incr = 5490. / (nelem -1);
    for (int i = 0; i < nelem; i++){
        tieX[i] = i*incr;
        tieY[i] = i*incr;
    }

    size_t nx = sizeof(tieX) / sizeof(tieX[0]);
    size_t ny = sizeof(tieY) / sizeof(tieY[0]);

    gsl_spline2d *splineZenith = gsl_spline2d_alloc(T, nx, ny);
    gsl_spline2d *splineAzimuth = gsl_spline2d_alloc(T, nx, ny);
    gsl_interp_accel *xacc = gsl_interp_accel_alloc();
    gsl_interp_accel *yacc = gsl_interp_accel_alloc();

    gsl_spline2d_init(splineZenith, tieX, tieY, tieZenith, nx, ny);
    gsl_spline2d_init(splineAzimuth, tieX, tieY, tieAzimuth, nx, ny);

    angles[0] = (float) gsl_spline2d_eval(splineZenith,pixel,scan,xacc, yacc);
    angles[1] = (float) gsl_spline2d_eval(splineAzimuth,pixel,scan,xacc, yacc);
    *senz = angles[0];
    *sena = angles[1];

    gsl_spline2d_free(splineZenith);
    gsl_spline2d_free(splineAzimuth);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);

}

static float *Fobar;

static geom_struc *gm_p_b = NULL;

//TODO: Replace the following callbacks with olog implementation

/**
  error callback expecting a FILE* client object
*/
void error_callback(const char *msg, void *client_data) {
    (void)client_data;
    fprintf(stdout, "[ERROR] %s", msg);
}
/**
 warning callback expecting a FILE* client object
*/
void warning_callback(const char *msg, void *client_data) {
    (void)client_data;
    fprintf(stdout, "[WARNING] %s", msg);
}
/**
  debug callback expecting no client object
*/
void info_callback(const char *msg, void *client_data) {
    (void)client_data;
    fprintf(stdout, "[INFO] %s", msg);
}
/************************************************/
/*  function createPrivateData:                 */
/*  create private data for Sentinel 2 MSI data */
/************************************************/
msi_t* createPrivateData_msi(int numBands){
    // Allocate memory for msi_t data struct
    msi_t* data = (msi_t*)calloc(1, sizeof(msi_t));
    if(data == NULL){
        fprintf(stderr,"-E- %s line %d: unable to allocate private data for MSI\n",
            __FILE__,__LINE__);
        exit(1);
    }

    // Allocate memory for storing upper left coordinates
    data->ULCoord = (int32_t*) malloc(2*sizeof(int32_t));

    return data;   
}


/*************************************************/
/* function resample_msi:                        */
/* Resamples 10m and 60m resolution bands to 20m */

/*************************************************/
void resample_msi(opj_image_t* image, filehandle* file, int recnum, int srcRes, int destRes) {

    int width;
    int i, i0;
    msi_t* data = (msi_t*) file->private_data;

    if (!data->buf) {
        data->buf = (uint32_t*) malloc(file->npix * sizeof (uint32_t));
    }
    // Set resizing scale
    float scale = (float) srcRes / (float) destRes;
    if (scale <= 0) {
        fprintf(stderr, "-E- %f scale is not calculated correctly\n", scale);
        exit(EXIT_FAILURE);
    }

    width = image->comps[0].w;

    // Resample the image
    for (i = 0; i < scale * width; i++) {
        if (scale > 1) {
            i0 = floor(i / scale);
            data->buf[i] = image->comps[0].data[recnum * width + i0];
        } else if (scale == 1) {
            data->buf[i] = image->comps[0].data[recnum * width + i];
        } else {
            i0 = floor(i / scale);
            data->buf[i] = (image->comps[0].data[recnum * width + i0]
                    + image->comps[0].data[recnum * width + (i0 + 1)]
                    + image->comps[0].data[(recnum + 1) * width + i0]
                    + image->comps[0].data[(recnum + 1) * width + (i0 + 1)]) / 4;
        }
    }
}

int32_t readTileMeta_msi(filehandle *file) {
    xml_document rootNode;
    xml_node dataNode;
    xml_node zenithNode;
    xml_node azimuthNode;
    xml_node valuesNode;
    char *delim = " "; // input separated by spaces
    char *token = NULL;
    msi_t *data = (msi_t*) file->private_data;
    char* tmpBuff;
    const char *xmlBuffer;
    if(!rootNode.load_file(data->granuleMetadataFile)) {
       return 0;
    }
    
    // fix orbit altitude
    data->alt = 786.;
    // Get geocoding information
    dataNode = rootNode.first_element_by_path("n1:Level-1C_Tile_ID/n1:Geometric_Info/Tile_Geocoding");
    // Get coordinate system EPSG code
    xmlBuffer = dataNode.first_element_by_path("HORIZONTAL_CS_CODE").first_child().value();
    std::string s = xmlBuffer;
    std::string delimiter = ":";
    s.erase(0, s.find(delimiter) + delimiter.length());
    data->CSCode = atoi(s.c_str());
    
    // Get coordinate system UTM zone
    xmlBuffer = dataNode.first_element_by_path("HORIZONTAL_CS_NAME").first_child().value();
    s = xmlBuffer;
    delimiter = "zone ";
    s.erase(0, s.find(delimiter) + delimiter.length());
    if(s.back() == 'S') {
        s.pop_back();
        s.append(" +south");
    } else if(s.back() == 'N') {
        s.pop_back();
    }
    data->UTMZone = strdup(s.c_str());
       
    xmlBuffer = dataNode.find_child_by_attribute("Size","resolution","20").child_value("NCOLS");
    file->npix = atoi(xmlBuffer);
    xmlBuffer = dataNode.find_child_by_attribute("Size","resolution","20").child_value("NROWS");
    file->nscan = atoi(xmlBuffer);
    xmlBuffer = dataNode.find_child_by_attribute("Geoposition","resolution","20").child_value("ULX");
    data->ULCoord[0] = atoi(xmlBuffer);    
    xmlBuffer = dataNode.find_child_by_attribute("Geoposition","resolution","20").child_value("ULY");
    data->ULCoord[1] = atoi(xmlBuffer);    
    
    // Set projection string
    char pjStr[FILENAME_MAX];
    sprintf(pjStr, "+proj=utm +ellps=WGS84 +datum=WGS84 +zone=%s +units=m", data->UTMZone);

    // init the proj4 projections
    PJ *pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                               pjStr,
                               "+proj=longlat +ellps=WGS84 +datum=WGS84",
                               NULL);
    if(pj == NULL) {
        printf("Error - MSI first PROJ projection failed to init\n");
        exit(1);
    }
    data->pj = proj_normalize_for_visualization(PJ_DEFAULT_CTX, pj);
    if(data->pj == NULL) {
        printf("Error - MSI visualization PROJ projection failed to init\n");
        exit(1);
    }
    proj_destroy(pj);
    
    // Read tie-point sensor view angles
    dataNode = rootNode.first_element_by_path("n1:Level-1C_Tile_ID/n1:Geometric_Info/Tile_Angles/Viewing_Incidence_Angles_Grids");
      // Set up tie point array sizes ([numDetectors][maxBands][HEIGHT][WIDTH])
    int tiepoint_height = 23;
    int tiepoint_width = 23;
    data->sensorZenith = (double **) malloc(maxBands * sizeof (double *));
    data->sensorAzimuth = (double **) malloc(maxBands * sizeof (double *));
    for (int i = 0; i < maxBands; i++) {
        data->sensorZenith[i] = (double *) calloc(tiepoint_width * tiepoint_height, sizeof (double ));
        data->sensorAzimuth[i] = (double *) calloc(tiepoint_width * tiepoint_height, sizeof (double ));
    }
    while (dataNode) {
        xmlBuffer = dataNode.attribute("bandId").value();
        int bandIdx = atoi(xmlBuffer);
        xmlBuffer = dataNode.attribute("detectorId").value();
        zenithNode = dataNode.first_element_by_path("Zenith/Values_List");
        valuesNode = zenithNode.first_element_by_path("VALUES");
        int i = 0;
        while (valuesNode){
            tmpBuff = strdup(valuesNode.first_child().value());
            int j = 0;
            for (token = strtok((char *) tmpBuff, delim); token != NULL; token = strtok(NULL, delim)) {
                char *unconverted;
                double value = strtof(token, &unconverted);
                if (!std::isnan(value))
                    data->sensorZenith[bandIdx][i*tiepoint_width + j] = value;
                j++;
            }
            free(tmpBuff);
            i++;
            valuesNode = valuesNode.next_sibling("VALUES");
        }
        azimuthNode = dataNode.first_element_by_path("Azimuth/Values_List");
        valuesNode = azimuthNode.first_element_by_path("VALUES");
        i = 0;
        while (valuesNode){
            tmpBuff = strdup(valuesNode.first_child().value());
            int j = 0;
            for (token = strtok((char *) tmpBuff, delim); token != NULL; token = strtok(NULL, delim)) {
                char *unconverted;
                double value = strtod(token, &unconverted);
                if (!std::isnan(value))
                    data->sensorAzimuth[bandIdx][i*tiepoint_width + j] = value;
                j++;
            }
            free(tmpBuff);
            i++;
            valuesNode = valuesNode.next_sibling("VALUES");
        }
        dataNode = dataNode.next_sibling("Viewing_Incidence_Angles_Grids");
    }
    return 1;
}

int32_t readDatastripMeta_msi(filehandle *file) {
    xml_document rootNode;
    xml_node dataNode;
    xml_node subDataNode;
    xml_node_iterator it;
    double unixseconds;
    char *delim = " "; // input separated by spaces
    char *token = NULL;
    msi_t *data = (msi_t*) file->private_data;
    char* tmpBuff;
    const char *xmlBuffer;
    if (!rootNode.load_file(data->datastripMetadataFile)) {
        return 0;
    }
    // Get start time
    dataNode = rootNode.first_element_by_path("n1:Level-1C_DataStrip_ID/n1:General_Info/Datastrip_Time_Info");
    xmlBuffer = rootNode.first_element_by_path("n1:Level-1C_DataStrip_ID/n1:General_Info/Datastrip_Time_Info/DATASTRIP_SENSING_START").first_child().value();
    data->scene_start_time = isodate2unix((char*) xmlBuffer);
    // Get detector line period
    dataNode = rootNode.first_element_by_path("n1:Level-1C_DataStrip_ID/n1:Image_Data_Info/Sensor_Configuration/Time_Stamp");
    xmlBuffer = dataNode.first_element_by_path("LINE_PERIOD").first_child().value();
    data->lineTimeDelta = atof(xmlBuffer) * 2e-3; //convert to seconds, for 20m band
    
    // Get detector relative startimes
    dataNode = rootNode.first_element_by_path("n1:Level-1C_DataStrip_ID/n1:Image_Data_Info/Sensor_Configuration/Time_Stamp/Band_Time_Stamp");
    while (dataNode) {
        xmlBuffer = dataNode.attribute("bandId").value();
        int bandIdx = atoi(xmlBuffer);
        subDataNode = dataNode.first_element_by_path("Detector");
        while (subDataNode) {
            xmlBuffer = subDataNode.attribute("detectorId").value();
            int detectorIdx = atoi(xmlBuffer) - 1;
            xmlBuffer = subDataNode.child_value("GPS_TIME");
            unixseconds = isodate2unix((char*) xmlBuffer);
            data->detectorDeltaTime[detectorIdx][bandIdx] = data->scene_start_time - unixseconds;
            subDataNode = subDataNode.next_sibling("Detector");
        }
        
        dataNode = dataNode.next_sibling("Band_Time_Stamp");
    }
    
    // Get number of GPS entries
    dataNode = rootNode.first_element_by_path("n1:Level-1C_DataStrip_ID/n1:Satellite_Ancillary_Data_Info/Ephemeris/GPS_Points_List");
    data->num_gps = 0;
    for (it = dataNode.begin(); it != dataNode.end(); it++) {
        data->num_gps++;
    }
    data->gpstime = (double *) malloc(data->num_gps * sizeof(double));
    for (int i = 0; i<3; i++)
        data->position[i] = (double *) malloc(data->num_gps * sizeof(double));
    // read GPS
    int i = 0;
    dataNode = rootNode.first_element_by_path("n1:Level-1C_DataStrip_ID/n1:Satellite_Ancillary_Data_Info/Ephemeris/GPS_Points_List/GPS_Point");
    while (dataNode) {
        tmpBuff = strdup(dataNode.child_value("POSITION_VALUES"));
        int j = 0;
        for (token = strtok((char *)tmpBuff, delim); token != NULL; token = strtok(NULL, delim)) {
            char *unconverted;
            data->position[j][i] = strtod(token, &unconverted);
            j++;
        }
        free(tmpBuff);
        xmlBuffer = dataNode.child_value("GPS_TIME");
        data->gpstime[i] = isodate2unix((char*) xmlBuffer);
        i++;
        dataNode = dataNode.next_sibling("GPS_Point");
    }
    return 1;
}

int32_t readDetectorFootprint_msi(filehandle *file, int band) {
    xml_document rootNode;
    xml_node dataNode;
    xml_node polyNode;
    msi_t *data = (msi_t*) file->private_data;
    std::vector<std::string> detstrs;
    std::vector<std::string> polypointstrs;
    std::vector<std::string>::iterator sit;
    const char *detectorName, *polygon;
    if (!rootNode.load_file(data->detectorFootprintFiles[band])) {
        return 0;
    }
    dataNode = rootNode.first_element_by_path("eop:Mask/eop:maskMembers/eop:MaskFeature");
    while (dataNode) {
        detectorName = dataNode.attribute("gml:id").value();
        //detector_footprint-B05-03-0
        boost::split(detstrs, detectorName, boost::is_any_of("-"));
        int detidx = std::atoi(detstrs[2].c_str()) - 1;
        polyNode = dataNode.first_element_by_path("eop:extentOf/gml:Polygon/gml:exterior/gml:LinearRing/gml:posList");
        polygon = polyNode.first_child().value();
        boost::split(polypointstrs, polygon, boost::is_any_of(" "));
        std::string polyWKT = "POLYGON((";
        int i = 1;
        for(sit=polypointstrs.begin() ; sit < polypointstrs.end(); sit++,i++ ) {
        // skip every third element
            if(i % 3) 
                polyWKT = polyWKT + *sit + " ";
            else
                polyWKT = polyWKT + ",";
        }
        polyWKT = polyWKT + "))";
        boost::geometry::read_wkt(polyWKT, data->detectorFootprints[detidx]);
        dataNode = dataNode.next_sibling("eop:MaskFeature");
    }
    
    return 1;
}


/************************************************/
/*  function: openl1_msi                       */
/*  Open msi l1c data with metadata file        */
/************************************************/
extern "C" int openl1_msi(filehandle *file){

    int i;
    xml_document rootNode;
    xml_node dataNode, metaNode;
    const char *productName, *dataName;
    msi_t *data;

    // fill up the private data
    file->private_data = data = createPrivateData_msi(maxBands); 
    
    if(!rootNode.load_file(file->name)) {
        printf("Could not open MSI file: %s\n", file->name);
        exit(EXIT_FAILURE);
    }

    if(want_verbose)
        printf("Input file: MSI Level-1C %s\n", file->name);
    
    metaNode = rootNode.first_element_by_path("xfdu:XFDU/metadataSection");
        
    // orbit direction
    dataNode = metaNode.find_child_by_attribute("metadataObject", "ID", "measurementOrbitReference");
    dataName = dataNode.first_element_by_path("metadataWrap/xmlData/safe:orbitReference/safe:orbitNumber").attribute("groundTrackDirection").value();
    data->orbit_direction = strdup(dataName);

    // Band image paths
    int nbandsImg = 0;
    int nbandsDetFoot = 0;
    int pickMe = 0;
    
    dataNode = rootNode.first_element_by_path("xfdu:XFDU/dataObjectSection/dataObject") ;
    string indir = file->name;
    size_t index = indir.find_last_of('/');
    if(index != string::npos)
        indir.erase(index);
    else
        indir.clear();
    string fileName;
    
    while (dataNode) {
        dataName = dataNode.attribute("ID").value();
        productName = dataNode.first_element_by_path("byteStream/fileLocation").attribute("href").value();
        if(indir.empty())
            fileName = productName;
        else
            fileName = indir + "/" + productName;
        index = fileName.find("./");
        if(index != string::npos)
            fileName.erase(index, 2);
        if (strstr(dataName, "S2_Level-1C_Tile1_Metadata")){
            data->granuleMetadataFile = strdup(fileName.c_str());
        }
        if (strstr(dataName, "S2_Level-1C_Datastrip1_Metadata")){
            data->datastripMetadataFile = strdup(fileName.c_str());
        }
        if (strstr(dataName, "IMG_DATA_Band") && !(strstr(dataName, "TCI"))) {
            if (nbandsImg > maxBands) {
                printf("%s, %d - E - Maximum number of radiance values (%d) reached\n",
                        __FILE__, __LINE__, maxBands);
                exit(EXIT_FAILURE);
            }
            data->jp2[nbandsImg] = strdup(fileName.c_str());
            if (!data->jp2[nbandsImg]) {
                printf("%s, %d - E - unable to set path for band %d\n",
                        __FILE__, __LINE__, nbandsImg);
                exit(EXIT_FAILURE);
            }
            nbandsImg++;
        }
        if (strstr(dataName, "DetectorFootprint")) {
            if (nbandsDetFoot > maxBands) {
                printf("%s, %d - E - Maximum number of radiance values (%d) reached\n",
                        __FILE__, __LINE__, maxBands);
                exit(EXIT_FAILURE);
            }
            data->detectorFootprintFiles[nbandsDetFoot] = strdup(fileName.c_str());
            if (!data->detectorFootprintFiles[nbandsDetFoot]) {
                printf("%s, %d - E - unable to set path for band %d\n",
                        __FILE__, __LINE__, nbandsDetFoot);
                exit(EXIT_FAILURE);
            }
            if (strstr(data->detectorFootprintFiles[nbandsDetFoot], "B07")) {
                pickMe = nbandsDetFoot;
            }
            nbandsDetFoot++;    
        }
        dataNode = dataNode.next_sibling("dataObject");
    }
    
    /* Read tile metadata file */
    if(!readTileMeta_msi(file))
        fprintf(stderr, "-E- %s line %d: unable read tile metadata file for MSI dataset %s\n",
            __FILE__,__LINE__,data->granuleMetadataFile);
    
    /* Read datastrip metadata file */
    if(!readDatastripMeta_msi(file))
        fprintf(stderr, "-E- %s line %d: unable read datastrip metadata file for MSI dataset %s\n",
            __FILE__,__LINE__,data->datastripMetadataFile);
    
    /* Read detector footprint file 
     for simplicity, read just one for a 20m band
     */
    if(!readDetectorFootprint_msi(file, pickMe))
        fprintf(stderr, "-E- %s line %d: unable read detector footprint file for MSI dataset %s\n",
            __FILE__,__LINE__,data->detectorFootprintFiles[pickMe]);

    strcpy(file->spatialResolution, "20 m");
    if(want_verbose){

        // Print out all MSI jp2 images 
        for (i = 0; i<maxBands; i++){
            printf("MSI file %d: %s\n",i, data->jp2[i]);
        }
    }
      
    /*
     *  get the Fobar here to set up Fo
     */
    rdsensorinfo(file->sensorID, l1_input->evalmask, "Fobar", (void **) &Fobar);
    file->terrain_corrected = 1;

    return 0;
}


/************************************************/
/*  function: readl1_msi                       */
/*  Get coordinates and convert to lon/lat      */
/************************************************/
extern "C" int readl1_msi_lonlat(filehandle *file, int recnum, l1str *l1rec)
{

    int ip;
    msi_t* data = (msi_t*) file->private_data; 
    
    // Convert lon and lon results of current scanline from radian to degree values
    PJ_COORD c, c_out;
    
    // set default z and t
    c.xyzt.z = 0.0;
    c.xyzt.t = HUGE_VAL;
    for (ip=0; ip<file->npix; ip++) {
        c.xy.x = data->ULCoord[0] + 10 + ip * 20;
        c.xy.y = data->ULCoord[1] - 10 - recnum * 20;
        c_out = proj_trans(data->pj, PJ_FWD, c);
        l1rec->lon[ip] = c_out.xy.x;
        l1rec->lat[ip] = c_out.xy.y;
    }

    return 0;
}

/************************************************/
/*  function: readl1_msi                       */
/*  Read MSI l1c image data line by line        */
/************************************************/
extern "C" int readl1_msi(filehandle *file, int recnum, l1str *l1rec, int lonlat)
{
    int i, ip, ib, ipb;
    msi_t* data = (msi_t*) file->private_data; 
    int16_t year, doy;
    float  sunpos[3];
    double secondOfDay;
    float sunDist;
    char bandStrBuffer[4];
    l1rec->scantime = data->scene_start_time + recnum * data->lineTimeDelta; //may want to to some math with recnum here

    // Get lat lon 
    if(readl1_msi_lonlat(file,recnum,l1rec)) {
        fprintf(stderr,"-E- %s line %d: unable to allocate lat/lon data for MSI\n",
            __FILE__,__LINE__);
        exit(1);
    }

    // Set information about data
    l1rec->npix = file->npix;
    l1rec->l1file->sensorID = file->sensorID;
    
    unix2yds(l1rec->scantime, &year, &doy, &secondOfDay);

    int32_t iyear, idoy, msec;
    iyear = (int32_t) year;
    idoy = (int32_t) doy;
    msec = (int32_t) secondOfDay*1e3;
    double esdist = esdist_(&iyear, &idoy, &msec);
    double fsol = pow(1.0 / esdist, 2);

    /*
    *  if required for that record, set up the geom_per_band storage
    */
    if (!lonlat) {
        if((l1_input->geom_per_band == 1 ) && ( l1rec->geom_per_band == NULL ) ) {
            init_geom_per_band( l1rec );
            gm_p_b = l1rec->geom_per_band; // store this address so that it can be later destroyed in close()
        }
    }
    
    l_sun_(&iyear, &idoy, &secondOfDay, sunpos, &sunDist); // get position vector for the sun

    for (i=0; i < 3; i++) {
        sunpos[i] *= 1.496e8; //convert to km for call to get_zenaz
    }
    // Assign viewing angles to l1rec struct
    for (ib = 0; ib<maxBands; ib++) {
        int len = strlen(data->jp2[ib])-7;
        strncpy(bandStrBuffer, data->jp2[ib]+len, 3);
        bandStrBuffer[3] = '\0';
        int bandIdx = str2enum(bandStrBuffer);
        if (bandIdx == B10)
            continue;
//        int detector = -1;
        for (ip=0; ip<file->npix; ip++) {
            // use boost.within to determine detector number
            // but only need to do this once - per band differences not significant
            if (ib == 0) {
//                detector = inDetector (data , l1rec->lat[ip], l1rec->lon[ip]);

                // interpolate GPS position vectors to scantime
                //  If we can get more accurate times, we can use this
                // interpGPSpos method for sensor angles
                // NOTE: the solar angles are also affected by the incorrect scantime
                //       but less egregiously so...

//                interpGPSpos(l1rec,pos,detector,ib);

//                for (i=0; i < 3; i++) {
//                    epos[i]    = pos[i] * 1e-6; //values are in mm, convert to km
//                }
//                get_zenaz(epos, l1rec->lon[ip], l1rec->lat[ip], &l1rec->senz[ip], &l1rec->sena[ip]);

                // Assign band independent solar angles to l1rec struct
                get_zenaz(sunpos, l1rec->lon[ip], l1rec->lat[ip], &l1rec->solz[ip], &l1rec->sola[ip]);
                if (!lonlat) {
                    // interpolate tiePoint sensor view angles
                    interpViewAngles(data, ip, recnum, ib, &l1rec->senz[ip], &l1rec->sena[ip]);
                }
            }
 
            // if getting only lonlat stuff, skip everything below and go to the next iteration
            if (lonlat)
                continue;
          
            if (l1_input->geom_per_band) {
                int ipb = ip*file->nbands + bandIdx;

                // re-interpolate GPS position vectors to per band scantime (see note above)
//                interpGPSpos(l1rec,pos,detector,ib);
//
//                for (i=0; i < 3; i++)
//                    epos[i]    = pos[i] * 1e-6; //values are in mm, convert to km
//
//                get_zenaz(epos, l1rec->lon[ip], l1rec->lat[ip], &l1rec->geom_per_band->senz[ipb], &l1rec->geom_per_band->sena[ipb]);
                // interpolate tiePoint sensor view angles
                interpViewAngles(data, ip, recnum, ib, &l1rec->geom_per_band->senz[ipb], &l1rec->geom_per_band->sena[ipb]);
                // per band solar angles are just repetitions of the nominal values
                l1rec->geom_per_band->solz[ipb] = l1rec->solz[ip];
                l1rec->geom_per_band->sola[ipb] = l1rec->sola[ip];
            }
        }

        // lonlat mode, only care about when ib = 0, which gets the sensor view angles
        if (lonlat)
            break;
    }

    // Calculate surface reflectance values for each band 
    for(ib = 0; ib < maxBands; ib++) {        
        int len = strlen(data->jp2[ib])-7;
        strncpy(bandStrBuffer, data->jp2[ib]+len, 3);
        bandStrBuffer[3] = '\0';
        int bandIdx = str2enum(bandStrBuffer);
        if (bandIdx == B10)
            continue;

        l1rec->Fo[bandIdx] = Fobar[bandIdx] * fsol;
   
        if(decodeMSI(file, bandIdx, recnum)!=0) {
            printf("-E-: Error decoding MSI jp2 files.\n");
            fprintf(stderr, "-E- %s line %d: Failed to read Lt, band %d, recnum %d\n",
                    __FILE__,__LINE__, bandIdx, recnum );
            exit(1);
        }

        for (ip=0; ip<file->npix; ip++) {
            ipb = ip*file->nbands+bandIdx;
            if(data->buf[ip] == 0) {
                l1rec->Lt[ipb] = BAD_FLT;   
                l1rec->navfail[ip] = 1;
            } else{
                // 10000 is the QUANTIFICATION_VALUE of MSI data to convert DN value to reflectance
                // This value is listed in the MTD_MSIL1C.xml file, but we're just hardcoding it...bad?
                float quant = 10000.;
                float rToa = (float) (data->buf[ip] / quant);

                l1rec->Lt[ipb] = (rToa * l1rec->Fo[bandIdx] * cos(l1rec->solz[ip]/RADEG)) / PI ;
            }
        }
    }
    
    // Skip everything else if lonlat
    // Surface reflectance should give a good start pixel for l1info
    if (lonlat)
        return 0;
  
    // Calculate rho_cirrus from cirrus band 10
    if(decodeMSI(file, B10, recnum)!=0) {
        fprintf(stderr, "-E- %s line %d: Failed to read cirrus band, recnum %d\n",
                __FILE__,__LINE__, recnum );
        exit(1);
    }
    for (ip=0;ip<file->npix; ip++) {
        if(data->buf[ip] == 0)
            l1rec->rho_cirrus[ip] = BAD_FLT;
        else
            l1rec->rho_cirrus[ip] = data->buf[ip] / (PI * 10000.0);
    }

    // Check lat and lon values
    for (ip=0; ip<file->npix; ip++) {
        l1rec->pixnum[ip] = ip;
        if ( std::isnan(l1rec->lat[ip]) )
            l1rec->lat[ip] = -999.0;
        if ( std::isnan(l1rec->lon[ip]) )
            l1rec->lon[ip] = -999.0;

        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
            l1rec->lat[ip] <  -91.0 || l1rec->lat[ip] >  91.0 )
        l1rec->navfail[ip] = 1;

    }
    return 0;
}

uint32_t scale_recnum( int32_t bandIdx, int32_t recnum){
    switch(bandIdx) {
        case B01:
        case B09:
        case B10:
            return floor(recnum/3);
        case B02:
        case B03:
        case B04:
        case B08:
            return recnum*2;
        default:
            return recnum;
    }
}
/************************************************/
/*  function: decodeMSI                         */
/*  Decode jp2 data                             */
/************************************************/
int decodeMSI(filehandle *file, int32_t bandIdx, int32_t recnum){

    msi_t* data = (msi_t*) file->private_data; 
    static int32_t initFile[13];

    int32_t fileIdx, tileIdx;
    int32_t scaledrecnum;

    for (fileIdx=0; fileIdx<maxBands; fileIdx++){
        if (strstr(data->jp2[fileIdx], conversion[bandIdx].str) != NULL)
            break;
    }
    char bandPath[FILENAME_MAX];
    strcpy(bandPath, data->imgDir);
    strcat(bandPath, data->jp2[fileIdx]);

    // Set decoding parameters to default values
    if (!initFile[bandIdx]){
//        memset(&(data->parameters[bandIdx]), 0, sizeof(opj_decompress_parameters));
        memset(&(data->parameters[bandIdx]), 0, sizeof(opj_dparameters_t));

        // Specify default decoding parameters
//        opj_set_default_decoder_parameters(&(data->parameters[bandIdx].core));
        opj_set_default_decoder_parameters(&(data->parameters[bandIdx]));

        data->image[bandIdx] = NULL;
        data->l_stream[bandIdx] = NULL;
        data->l_codec[bandIdx] = opj_create_decompress(OPJ_CODEC_JP2);
        data->cstr_info[bandIdx] = NULL;

//        opj_set_info_handler(data->l_codec[bandIdx], info_callback,00);
        opj_set_warning_handler(data->l_codec[bandIdx], warning_callback,00);
        opj_set_error_handler(data->l_codec[bandIdx], error_callback,00);

    } else if (scale_recnum(bandIdx,recnum) >= data->parameters[bandIdx].DA_y1) {
        opj_stream_destroy(data->l_stream[bandIdx]);
        opj_destroy_codec(data->l_codec[bandIdx]);
        opj_image_destroy(data->image[bandIdx]);
        data->l_codec[bandIdx] = opj_create_decompress(OPJ_CODEC_JP2);
        initFile[bandIdx] = 0;
    }

    // Setup the decoder decoding parameters using user parameters
    if(!initFile[bandIdx]){
        // Read input file and put it in memory
        data->l_stream[bandIdx] = opj_stream_create_default_file_stream(bandPath, 1);
        if (!data->l_stream[bandIdx]){
            fprintf(stderr, "ERROR -> failed to create the stream from the file %s\n", bandPath);
            free(&(data->parameters[bandIdx]));
            return EXIT_FAILURE;
        }
        
//        if ( !opj_setup_decoder(data->l_codec[bandIdx], &(data->parameters[bandIdx].core))){
        if ( !opj_setup_decoder(data->l_codec[bandIdx], &(data->parameters[bandIdx]))){
            fprintf(stderr, "ERROR -> opj_compress: failed to setup the decoder\n");
            opj_stream_destroy(data->l_stream[bandIdx]);
            opj_destroy_codec(data->l_codec[bandIdx]);
            return EXIT_FAILURE;
        }

        // Read the main header of the codestream and if necessary the JP2 boxes
        if(! opj_read_header(data->l_stream[bandIdx], data->l_codec[bandIdx], &(data->image[bandIdx]))){
            fprintf(stderr, "ERROR -> opj_decompress: failed to read the header\n");
            opj_stream_destroy(data->l_stream[bandIdx]);
            opj_destroy_codec(data->l_codec[bandIdx]);
            opj_image_destroy(data->image[bandIdx]);
            return EXIT_FAILURE;
        }
        data->cstr_info[bandIdx] = opj_get_cstr_info(data->l_codec[bandIdx]);
        initFile[bandIdx] = 1;
    }

    // Set boundaries for decoding MSI data based on different resolutions
    switch(bandIdx) {
        // For 60m resolution data: band 1, band 9 and band 10,
        case B01:
        case B09:
        case B10:
            scaledrecnum = scale_recnum(bandIdx,recnum);
            tileIdx = floor(scaledrecnum / (int32_t)data->cstr_info[bandIdx]->tdy);
            data->parameters[bandIdx].DA_x0 = 0;
            data->parameters[bandIdx].DA_x1 = floor(file->npix/3);
            data->parameters[bandIdx].DA_y0 = tileIdx*data->cstr_info[bandIdx]->tdy;
            data->parameters[bandIdx].DA_y1 = std::min((tileIdx+1)*(int32_t)data->cstr_info[bandIdx]->tdy,file->nscan/3);
            break;
        // For 10m resolution data: band 2, band 3, band 4 and band 8,
        case B02:
        case B03:
        case B04:
        case B08:
            scaledrecnum = scale_recnum(bandIdx,recnum);
            tileIdx = floor(scaledrecnum / (int32_t)data->cstr_info[bandIdx]->tdy);
            data->parameters[bandIdx].DA_x0 = 0;
            data->parameters[bandIdx].DA_x1 = 2 * file->npix;
            data->parameters[bandIdx].DA_y0 = tileIdx*data->cstr_info[bandIdx]->tdy;
            data->parameters[bandIdx].DA_y1 = std::min((tileIdx+1)*(int32_t)data->cstr_info[bandIdx]->tdy,file->nscan*2);
            break;
        // For 20m resolution data: band 5, band 6, band 7, band 8a, band 11 and band 12
        default:
            scaledrecnum = scale_recnum(bandIdx,recnum);
            tileIdx = floor(scaledrecnum / (int32_t)data->cstr_info[bandIdx]->tdy);
            data->parameters[bandIdx].DA_x0 = 0;
            data->parameters[bandIdx].DA_x1 = file->npix;
            data->parameters[bandIdx].DA_y0 = tileIdx*data->cstr_info[bandIdx]->tdy;
            data->parameters[bandIdx].DA_y1 = std::min((tileIdx+1)*(int32_t)data->cstr_info[bandIdx]->tdy,file->nscan);
            break;
    }


    // Decode the JPEG2000 stream
    if (data->image[bandIdx]->comps->data == NULL) {
        // Decode the image based on boundaries
        if (!opj_set_decode_area(data->l_codec[bandIdx], data->image[bandIdx], data->parameters[bandIdx].DA_x0,
                data->parameters[bandIdx].DA_y0, data->parameters[bandIdx].DA_x1, data->parameters[bandIdx].DA_y1)){
            fprintf(stderr, "ERROR -> opj_decompress: failed to set the decoded area\n");
            opj_stream_destroy(data->l_stream[bandIdx]);
            opj_destroy_codec(data->l_codec[bandIdx]);
            opj_image_destroy(data->image[bandIdx]);
            return EXIT_FAILURE;
        }
         if (!opj_decode(data->l_codec[bandIdx], data->l_stream[bandIdx], data->image[bandIdx]) &&
                opj_end_decompress(data->l_codec[bandIdx], data->l_stream[bandIdx])){
            fprintf(stderr, "ERROR -> opj_decompress: failed to set the decoded area\n");
            opj_stream_destroy(data->l_stream[bandIdx]);
            opj_destroy_codec(data->l_codec[bandIdx]);
            opj_image_destroy(data->image[bandIdx]);
            return EXIT_FAILURE;
        }
    }

    // Resampling data of each band to 20m
    int32_t relative_recnum = scaledrecnum - data->parameters[bandIdx].DA_y0;

    switch(bandIdx)
    {
        case B01:
        case B09:
        case B10:
            resample_msi(data->image[bandIdx], file, relative_recnum, 60, 20);
            break;
        case B02:
        case B03:
        case B04:
        case B08:
            resample_msi(data->image[bandIdx], file, relative_recnum, 10, 20);
            break;
        default:
           resample_msi(data->image[bandIdx], file, relative_recnum, 20, 20);
            break;
    }

    return 0;
};


/************************************************/
/*  function: freeMSIData                       */
/*  Free memory allocated in createPrivateData  */
/************************************************/
void freeMSIData(msi_t* data) {
    int i;
    free(data->ULCoord);
    for (i = 0; i< 3; i++)
        free(data->position[i]);

    free(data->gpstime);
    free(data->buf);

    for (int i = 0; i < maxBands; i++) {
        free(data->sensorZenith[i]);
        free(data->sensorAzimuth[i]);
    }
    free(data->sensorZenith);
    free(data->sensorAzimuth);

    free(data);
}


/************************************************/
/*  function: closel1_msi                      */
/*  Close opened files and free mempry          */
/************************************************/
extern "C" int closel1_msi(filehandle *file){

    msi_t* data = (msi_t*) file->private_data;
        
    freeMSIData(data);
    file->private_data = NULL;
    
    return 0;
}
