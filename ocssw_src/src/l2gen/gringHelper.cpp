// Includes from main_l1info.c
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdbool.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include <unistd.h>

// Additional includes
#include "gringHelper.h"
#include "allocate2d.h"
#include "l2_flags.h"   // needed for NAVFAIL

/**
 * Constructor
 * @param gringConfig
 */
gringHelper::gringHelper(gringConfig_t *gringConfig) {
    
    // Allocate memory for geobox[4][MAX_GEOBOX_CNT]
    geobox = allocate2d_float(4, MAX_GEOBOX_CNT);
    centerLat = (float*) malloc(MAX_GEOBOX_CNT * sizeof(float));
    
    // State vars
    first_good_scan = false;
    last_direction = 0;
    last_lat = -999;
    geobox_cnt = 0;
    direction_line = -1;
    direction_delta = 5;
    
    // Unpack gringConfig
    delta_degrees_lat = gringConfig->delta_degrees_lat;
    scan_direction = gringConfig->scan_direction;
    direction_delta = gringConfig->direction_delta;
    
}


gringHelper::~gringHelper() {
    
    // Clean up
    free2d_float(geobox);

}

bool gringHelper::valid_lat(float lat) {
    if(lat > 90.0)
        return false;
    if(lat < -90.0)
        return false;
    return true;
}

bool gringHelper::valid_lon(float lon) {
    if(lon > 180.0)
        return false;
    if(lon < -180.0)
        return false;
    return true;
}

void gringHelper::process_scan(l1str *l1rec, size_t recnum, size_t escan) {
    bool gring_test_val;
    static scanBounds_t *scanBounds;
    if (scanBounds == NULL){
        scanBounds = new scanBounds_t;
    }
     
    gring_test_val = includeScanInGring(l1rec->lat,l1rec->lon,l1rec->flags,l1rec->l1file->npix,recnum,escan,scanBounds);
    if(gring_test_val) {
        if(geobox_cnt >= (MAX_GEOBOX_CNT - 1)) {
            printf("Error - Max number of gring points exceeded.\n");
            exit(EXIT_FAILURE);
        }
        //printf("Include last good scan %d where starting lat = %f, last direction = %d\n",(int)recnum,l1rec->lat[recnum],last_direction);
        //printf("        scanBoundsOut: slat = %f, slon = %f, elat = %f, elon = %f\n", scanBounds->slat, scanBounds->slon, scanBounds->elat, scanBounds->elon);
        // This scan should be included in the gring, so caller stuffs scanBoundsOut into a geobox array.
        // Note that a geobox (borrowed from l2_generic.c) is just the endpoints (slon,slat,elon,elat) of a scan line.
        geobox[0][geobox_cnt] = scanBounds->slon;          // Starting lon of this scan
        geobox[1][geobox_cnt] = scanBounds->slat;          // Starting lat of this scan
        geobox[2][geobox_cnt] = scanBounds->elon;          // Ending lon of this scan
        geobox[3][geobox_cnt] = scanBounds->elat;          // Ending lat of this scan
        centerLat[geobox_cnt] = scanBounds->clat;          // save the center lat also
        geobox_cnt++;
        
    }
    
}


/**
 * This function is called per scan. It returns true if this (recnum) scan should be included in a gring.
 * It borrows code from main_l1info.c (e.g. for finding good spix,epix,cpix) and from l2_generic.c (e.g. 
 * for building the gring from a geobox). 
 * 
 * NB: Allocation (and freeing) of memory for the output struct scanBounds is outside the scope of this function.
 *     If the return value is false, the output scanBounds should be ignored.
 * 
 * @param lat
 * @param lon
 * @param flags
 * @param num_pixels
 * @param recnum                    Index of this scan
 * @param escan                     Index of last scan
 * @param scanBounds                Output, a structure containing the good endpoints--if any--and center points of the scan line
 * @return 
 */
bool gringHelper::includeScanInGring(float *lat, float *lon, int32_t *flags, size_t num_pixels, size_t recnum, size_t escan,
            scanBounds_t *scanBounds ) { 
    
    int spix =  0;  // start pixel 
    int cpix = -1;  // center pixel     
    int epix = -1;  // end pixel 
    int ip, good_spix, good_epix, good_cpix;
    int this_direction = 0;
  
    // Initialize 
    spix = 0;
    epix = num_pixels - 1;
    cpix = num_pixels/2;
    good_spix = -1;
    good_epix = -1;
    good_cpix = -1;

    // Find the good start, center, and end pixels for this scan
    for (ip = spix; ip <= epix; ip++) {
        if ((flags[ip] & NAVFAIL) == 0 && valid_lat(lat[ip]) && valid_lon(lon[ip])) {
            // Store the first good pixel in good_spix
            if (good_spix == -1) {
                good_spix = ip;
            }
            // Store this good pixel in good_epix. Eventually (at or near the end), it will be the last good one.
            good_epix = ip;
            // Store this good pixel in cpix. Eventually (somewhere in the middle), it will be last good one at or before cpix.
            if ((ip <= cpix) || (good_cpix == -1)) {
                good_cpix = ip;
            }
        }
    }

    // If this is a good scan, set the scan bounds (i.e. endpoints), and test to see if it should be included in the gring
    if ((good_spix != -1) && (good_epix != -1) && (good_cpix != -1)) {
        //printf("good_spix = %d, good_cpix = %d, good_epix = %d\n", good_spix, good_cpix, good_epix);
        scanBounds->slat = lat[good_spix];
        scanBounds->slon = lon[good_spix];
        scanBounds->clat = lat[good_cpix];
        scanBounds->clon = lon[good_cpix];
        scanBounds->elat = lat[good_epix];
        scanBounds->elon = lon[good_epix];

        // CONDITION 1: If this is the FIRST good scan, return scanBounds, after resetting "state information".
        if (!first_good_scan) {
            first_good_scan = true;
            last_lat = lat[good_cpix];
            direction_line = recnum;
            //printf("includeScanInGring: First good scan returned, recnum = %d, lat[good_spix] = %f\n", (int)recnum, lat[good_spix]);
            return true;
        } else {
            // This is a good scan, but not the first good one. So, we can start tracking direction; i.e. no sense in 
            // tracking direction prior to the first good scan.
            // wait for a few lines before checking direction
            if((recnum - direction_line) > direction_delta) {
                direction_line = recnum;
                if ((lat[good_cpix] - last_lat) > 0) {
                    this_direction = 1;
                } else {
                    this_direction = -1;
                }
                if(last_direction == 0)
                    last_direction = this_direction;
            }
        }
        // CONDITION 2: If the latitude of this scan is more than delta_degrees_lat from the last_lat, return scanBounds, 
        // after resetting "state information". Also CONDITION 1 must have already been satisfied; i.e. last_lat
        // has a non-default value.
        if (first_good_scan && (fabs(last_lat - lat[good_cpix]) > delta_degrees_lat) ) {
            last_direction = this_direction;
            last_lat = lat[good_cpix];
            //printf("includeScanInGring: Scan returned due to delta lat, recnum = %d, lat[good_cpix] = %f\n", (int)recnum, lat[good_cpix]);
            return true;
        }
        // CONDITION 3: If the satellite passes over a pole, detect the change from e.g. ascending to descending (or vice versa)
        // and return scanBounds after resetting "state information". Also CONDITION 1 must have already been satisfied; i.e. last_lat
        // has a non-default value.  Only look in high latitudes.
        if (first_good_scan && (last_direction == -this_direction) && (this_direction != 0) && (fabs(lat[good_cpix]) > 75)) {
            last_direction = this_direction;
            last_lat = lat[good_cpix];
            //printf("includeScanInGring: Scan returned due to change in direction, recnum = %d, lat[good_cpix] = %f\n", (int)recnum, lat[good_cpix]);
            return true;
        }
    } 

   // If we get this far, this is not a scan to include in the gring. So, just return 0. However, if this is the last scan, we want to use
   // the last GOOD scan.
   if (recnum == escan) {
       // If we got this far (without returning) AND this is the last scan, use the last GOOD scan (still stored in scanBounds).
       return true;      
   } else {
       return false;
   }
}


void gringHelper::delete_extra_line() {
    
    if(geobox_cnt > 2) {
        if(fabs(centerLat[geobox_cnt-1] - centerLat[geobox_cnt-2]) < delta_degrees_lat / 2) {
            geobox_cnt--;
            geobox[0][geobox_cnt-1] = geobox[0][geobox_cnt];
            geobox[1][geobox_cnt-1] = geobox[1][geobox_cnt];
            geobox[2][geobox_cnt-1] = geobox[2][geobox_cnt];
            geobox[3][geobox_cnt-1] = geobox[3][geobox_cnt];
            centerLat[geobox_cnt-1] = centerLat[geobox_cnt];
        }
    }
}

/**
 * After accumulating the desired number of scanBounds into a geobox array (of dimension 4 x N), the caller 
 * (of includeScanInGring) can then call creatGring  to return a gring_t struct. This code is (mostly) borrowed 
 * from l2_generic.c. 
 * 
 * To create the clockwise gring, it goes forward through the ending lats/lons and backwards through the
 * starting lats/lons. (This approach works if the geobox is organized such that the scans are always in the -x 
 * direction, where the path is in the +y direction; i.e. the instrument scans "to the left". This function calls 
 * orderGeobox to guarantee that this is true for the geoboxOrdered.) 
 * 
 * @param geobox                geobox[4][MAX_GEOBOX_CNT] is a 2D array (borrowed from l2_generic.c), where
 *                                  geobox[0][i] is the starting longitude (slon) of the ith scan line to be included in the gring
 *                                  geobox[1][i] is the starting latitude  (slat) of the ith scan line to be included in the gring
 *                                  geobox[2][i] is the ending   longitude (elon) of the ith scan line to be included in the gring
 *                                  geobox[3][i] is the ending   latitude  (elat) of the ith scan line to be included in the gring
 * @param scan_direction        Indicates whether the instrument scans to the left (-1) or to the right (+1).
 * @param gringOut              Output, a struct having arrays gring_lat,gring_lon,gring_seq.
 */
void gringHelper::createGring(gring_t *gringOut) {
    delete_extra_line();
    
    // Declarations for packing gring struct
    float **geoboxOrdered;
    size_t i,j;
    
    // Allocate memory for geoboxOrdered[4][MAX_GEOBOX_CNT]
    geoboxOrdered = allocate2d_float(4, MAX_GEOBOX_CNT);

    // Order the geobox
    orderGeobox(geoboxOrdered);

    // Create the gring from the ordered geobox
    j = 1;
    // Start at first slon
    gringOut->gring_lon[0] = geoboxOrdered[0][0];
    // Get elons, going forward
    for (i = 0; i < geobox_cnt; i++) {
      gringOut->gring_lon[j++] = geoboxOrdered[2][i];
    }
    // Get slons, going backwards
    for (i = 0; i < geobox_cnt - 1; i++) {
      gringOut->gring_lon[j++] = geoboxOrdered[0][geobox_cnt - 1 - i];
    }
    // gring latitides and sequence numbers
    j = 1;
    // Start at first slat
    gringOut->gring_lat[0] = geoboxOrdered[1][0];
    gringOut->gring_seq[0] = j;
    // Get elats, going forward
    for (i = 0; i < geobox_cnt; i++) {
      gringOut->gring_seq[j] = j + 1;
      gringOut->gring_lat[j++] = geoboxOrdered[3][i];
    }
    // Get slats, going backwards
    for (i = 0; i < geobox_cnt - 1; i++) {
      gringOut->gring_seq[j] = j + 1;
      gringOut->gring_lat[j++] = geoboxOrdered[1][geobox_cnt - 1 - i];
    }
    // Verify
    for (i = 0; i < j; i++) {
      //printf("INFO: createGring: lat=%f lon=%f\n",gringOut->gring_lat[i],gringOut->gring_lon[i]);
    }
    // num_gring_pts
    gringOut->num_gring_pts = j;
    
    // Clean up
    free2d_float(geoboxOrdered);
    //delete[] geoboxOrdered[0];
    //delete[] geoboxOrdered;

}

/**
 * This function creates an ordered geobox (geoboxOrdered) from the possibly unordered geobox.
 * By ordered, we mean that the endpoints of each scan are ordered in a way that will make it
 * easy for createGring to create a gring in a clockwise order; i.e. geoboxOrdered is organized such that the
 * scans are always in the -x direction, where the path is in the +y direction; i.e. the instrument scans
 * "to the left".
 * 
 * @param geobox                    geobox[4][MAX_GEOBOX_CNT] is a 2D array (borrowed from l2_generic.c), where
 *                                      geobox[0][i] is the starting longitude (slon) of the ith scan line to be included in the gring
 *                                      geobox[1][i] is the starting latitude  (slat) of the ith scan line to be included in the gring
 *                                      geobox[2][i] is the ending   longitude (elon) of the ith scan line to be included in the gring
 *                                      geobox[3][i] is the ending   latitude  (elat) of the ith scan line to be included in the gring
 * @param scan_direction            Indicates whether the instrument scans to the left (-1) or to the right (+1).
 * @param geoboxOrdered             geoboxOrdered[4][MAX_GEOBOX_CNT] is like geobox except that endpoints may have been swapped to 
 *                                  make it easy for the caller to create the gring in clockwise order.
 */
void gringHelper::orderGeobox(float **geoboxOrdered) {
    //printf("INFO: orderGeobox: scan_direction = %d\n", scan_direction);

    // Declarations
    size_t scan_num,j;
    
    
    // Initialize geoboxOrdered = geobox
    for (scan_num=0; scan_num<geobox_cnt; scan_num++) {
        for (j=0; j<4; j++) {
            geoboxOrdered[j][scan_num] = geobox[j][scan_num];
        }
    }
    
    // For certain missions, the instrument scans to the right, so flip the endpoints
    if (scan_direction > 0) {
        for (scan_num=0; scan_num<geobox_cnt; scan_num++) {
            //printf("INFO: orderGeobox: Swapping endpoints for scan %d\n",(int)scan_num);
            // Lons
            geoboxOrdered[0][scan_num] = geobox[2][scan_num];
            geoboxOrdered[2][scan_num] = geobox[0][scan_num];
            // Lats
            geoboxOrdered[1][scan_num] = geobox[3][scan_num];
            geoboxOrdered[3][scan_num] = geobox[1][scan_num];
        }
    }
    
}



int gringHelper::get_geobox_cnt() {
    delete_extra_line();
    return (int)geobox_cnt;
    
}


/**
 * 
 * @param gring_lon_string: csv of gring longitudes, formatted 
 * @param gring_lat_string: csv of gring latitudes
 * @param gring_seq_string: csv of gring sequence numbers
 * @param max_string_length: Caller specifies max length for any of the above strings. So, caller must know the formats used
 *                           in xlat_floatArray2String, xlat_intArray2String
 * @return 0 on success, -1 on failure
 */
int gringHelper::get_gring_strings(std::string& gring_lon_string, std::string& gring_lat_string, std::string& gring_seq_string) {
    
    int num_gring_pts, gring_pt_num;
    
    if (geobox_cnt > 1) {
        // Enough scans for a gring, so keep going
        gring_t *this_gring = new gring_t;
        createGring(this_gring);
        // Declare arrays
        float *gring_lat, *gring_lon;
        size_t *gring_seq;
        // Allocate memory for arrays
        num_gring_pts = this_gring->num_gring_pts;
        gring_lat = (float *)calloc(num_gring_pts,sizeof(float));
        gring_lon = (float *)calloc(num_gring_pts,sizeof(float));
        gring_seq = (size_t *)calloc(num_gring_pts,sizeof(size_t));
        
        // Unpack the gring_t struct 
        for (gring_pt_num = 0; gring_pt_num < num_gring_pts; gring_pt_num++) {
            gring_lon[gring_pt_num] = this_gring->gring_lon[gring_pt_num];
            gring_lat[gring_pt_num] = this_gring->gring_lat[gring_pt_num];
            gring_seq[gring_pt_num] = this_gring->gring_seq[gring_pt_num];
            //printf("INFO: get_gring_strings: lat=%f lon=%f\n",gring_lat[gring_pt_num],gring_lon[gring_pt_num]);
        }
        
        // Convert arrays into strings and print
        // Lon
        xlat_floatArray2String(gring_lon_string, gring_lon, num_gring_pts); 
        // Lat
        xlat_floatArray2String(gring_lat_string, gring_lat, num_gring_pts);    
        // Seq
        xlat_intArray2String(gring_seq_string, gring_seq, num_gring_pts); 
        
        // Clean up
        delete(this_gring);
        free(gring_lat);
        free(gring_lon);
        free(gring_seq);
        // return success
        return(SUCCESS);
    } else {
        // NOT enough scans for a gring, so return -1
        return(-1);
    }
    
}
    
    
 
/**
 * To help display gring, convert an array of floats into a csv string, without whitespace, formatted like %6.2f
 * @param myString
 * @param myArray
 * @param array_length
 */
void gringHelper::xlat_floatArray2String(std::string& myString, float *myArray, int array_length) {

    int j;
    char tmp[64];                            // One value of myArray, formatted like %6.2f

    myString.clear();
    // Iterate through all but the last, adding a delimiter at the end
    for (j=0; j<array_length; j++) {
        if(j!=0)
            myString.append(",");
        sprintf(tmp, "%.5f",myArray[j]);
        myString.append(tmp);
    }
}


/**
 * To help display gring, convert an array of ints into a csv string, formatted like %d
 * @param myString
 * @param myArray
 * @param array_length
 */
void gringHelper::xlat_intArray2String(std::string& myString, size_t *myArray, int array_length) {

    int j;
    char tmp[64];

    myString.clear();
    // Iterate through all 
    for (j=0; j<array_length; j++) {
        if(j!=0)
            myString.append(",");
        sprintf(tmp, "%d",(int)myArray[j]);
        myString.append(tmp);
    }
}
