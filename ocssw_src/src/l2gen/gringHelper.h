#ifndef _GRING_HELPER_H
#define _GRING_HELPER_H

#ifndef MAX_GEOBOX_CNT
#define MAX_GEOBOX_CNT 100
#endif

#include "l1.h"
#include "allocate2d.h"

#include <string>

struct gringConfig_t {
    float delta_degrees_lat;                // To be derived from user input (e.g. l1info command line)
    int scan_direction;                     // CZCS is +1 (to the right); others are -1 (to the left)
    size_t direction_delta;                 // number of lines to wait before checking the direction
};

// The scanBounds structure is used to mark the end-points (and mid-point) of a scan that can contribute to a gring
struct scanBounds_t {
   float slat;		// Starting lat
   float slon;          // Starting lon
   float clat;		// Center lat
   float clon;          // Center lon
   float elat;		// End lat
   float elon;          // End lon
};

// gring_t struct
struct gring_t{
   float gring_lat[2*MAX_GEOBOX_CNT];
   float gring_lon[2*MAX_GEOBOX_CNT];
   size_t gring_seq[2*MAX_GEOBOX_CNT];
   int num_gring_pts;
};

//  The member functions of this class can be used by the CALLER to create a gring 
//  Call sequence:
//      CALLER (instantiation)  > gringHelper    
//      CALLER (per scan)       > process_scan      > includeScanInGring    > setScanBounds
//      CALLER                  > get_geobox_cnt
//      CALLER                  > get_gring_strings > createGring           > orderGeobox
//                                                  > xlat_floatArray2String
//                                                  > xlat_intArray2String


class gringHelper {
    
public:
    
    gringHelper(gringConfig_t *gringConfig);                                                            // Constructor
    virtual ~gringHelper();                                                                             // Destructor
    void process_scan(l1str *l1rec, size_t recnum, size_t escan);                                       // Process an individual scan line                                                                              // Print (to screen) gring-related information
    int get_geobox_cnt();
    int get_gring_strings(std::string& gring_lon_string, std::string& gring_lat_string, std::string& gring_seq_string);
    
    
private:
    
    // State information
    bool first_good_scan;
    int last_direction;
    float last_lat;
    float **geobox;                                                                 // Array 4 x MAX_GEOBOX_CNT
    float *centerLat;
    size_t geobox_cnt;
    size_t direction_line;  // last line where we checked the direction
    size_t direction_delta; // number of line to wait before checking the direction
    
    // Config
    float delta_degrees_lat;
    int scan_direction;
    
    // Structs
    //scanBounds_t *scanBounds;           // If scan should be included in gring, stuff the scanBoundsOut into a geobox.
    
    bool valid_lat(float lat);
    bool valid_lon(float lon);

    bool includeScanInGring(float *lat, float *lon, int32_t *flags, size_t num_pixels, size_t recnum, size_t escan,
                            scanBounds_t *scanBounds );

    void delete_extra_line();
    
    void createGring(gring_t *this_gring);

    void orderGeobox(float **geoboxOrdered);

    void xlat_floatArray2String(std::string& myString, float *myArray, int array_length);

    void xlat_intArray2String(std::string& myString, size_t *myArray, int array_length);

};

#endif

