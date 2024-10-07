/******************************************************************************
 *
 * NAME: mds.h
 *
 * DESCRIPTION: description of the structure that constitutes a Map Data Set.
 * Defines prototypes for users of the map functions.
 *
 * PARAMETER        TYPE    DESCRIPTION
 * -------------    ------  --------------------------------
 * grid_type        short   identifies which map projection is being used
 *                           11=standard Mercator, cylindrical, coaxial
 *                           21=Polar Stereographic, Northern Hemisphere
 *                           25=Polar Stereographic, Southern Hemisphere
 *                           31=Lambert Conformal Conic, Tangent cone,
 *                              Northern Hemisphere
 *                           35=Lambert Conformal Conic, Tangent cone,
 *                              Southern Hemisphere
 *                           41=Lambert Conformal Conic, Secant cone,
 *                              Northern Hemisphere
 *                           45=Lambert Conformal Conic, Secant cone,
 *                              Southern Hemisphere
 *                           51=Cylindrical Equidistant
 *                           61=Northern Polar Azmuthal Equidistant
 *                           65=Southern Polar Azmuthal Equidistant
 *
 * wedge_rotation   short   used for Lambert Conformal Conic, position of the
 *                          empty wedge: 1=up,2=left,3=down,4=right
 *
 * mds_num          int     identification number of mds, optional
 *
 * stan_lat1        double  first standard latitude
 * stan_lat2        double  second standard latitude, needed for Lambert
 *                           Conformal Conics, secant cone projections
 * base_lon         double  base longitude of the X-Y coordinate system, not
 *                           always the same meaning, see notes for projection
 * grid_inc_constant double primarily determines scaling and grid sizes
 * grid_inc_c2      double  second constant of same purpose, needed in
 *                             cylindrical equidistant
 * grid_exponent    double  needed for Lamberts
 * grid_constant_A  double  needed for Lamberts,
 *                            -on mercators (grid_type=11) it is the width of
 *                            the MDS in radians
 * max_row          double  maximum row number in the grid
 * max_col          double  maximum column number in the grid
 * upr_left_lat     double  latitude of upper left corner of grid
 * upr_left_lon     double  longitude of upper left corner of grid
 * lwr_right_lat    double  latitude of lower right corner of grid
 * lwr_right_lon    double  longitude of lower right corner of grid
 * upr_left_x       double  x coordinate of upper left corner
 * upr_left_y       double  y coordinate of upper left corner
 * lwr_right_x      double  x coordinate of lower right corner
 * lwr_right_y      double  y coordinate of lower right corner
 * split_lon        double  for all lamberts, the longitude on both
 *                             sides of the empty sector
 * latsml           double  smallest latitude in the MDS
 * latbig           double  largest latitude in the MDS
 * lonsml[2]        double  smallest longitude in the MDS, need two because
 *                             there might be two longitude ranges in the MDS,
 *                             which happens when the MDS spans 180 degrees
 *                             longitude
 * lonbig[2]        double  biggest longitude in the MDS, need two for the
 *                             same reason as lonsml
 * num_lonrange     short   number of longitude ranges in the MDS, always
 *                             either 1 or 2
 * junk[3]          short   filler to make the length of an MDS an even
 *                          multiple of 8 bytes.
 * -------------    ------  --------------------------------
 *
 * REFERENCES: none
 *
 * LIMITATIONS: none
 *
 * NOTES (MISCELLANEOUS) SECTION:
 *
 * MDS = Map Data Set
 *
 ******************************************************************************/

#ifndef  MDS_TYPE_H
#define  MDS_TYPE_H

#ifdef __cplusplus
extern "C" {
#endif


/* Map Data Set Grid types */
const int NORTHERN_POLAR_STEREOGRAPHIC = 21;
const int SOUTHERN_POLAR_STEREOGRAPHIC = 25;

/***************************************************************************
Define the Map Data Set
 **********************/

typedef struct {
    short grid_type;
    short wedge_rotation;
    int mds_num;
    double stan_lat1;
    double stan_lat2;
    double base_lon;
    double grid_inc_constant;
    double grid_inc_c2;
    double grid_exponent;
    double grid_constant_A;
    double max_row;
    double max_col;
    double upr_left_lat;
    double upr_left_lon;
    double lwr_right_lat;
    double lwr_right_lon;
    double upr_left_x;
    double upr_left_y;
    double lwr_right_x;
    double lwr_right_y;
    double split_lon;
    double costanlat;
    double latsml;
    double latbig;
    double lonsml[2];
    double lonbig[2];
    short num_lonrange;
    char explicit_pad[6];
} mds_type;

/***************************************************************************
NOTE: All Map Data Sets (MDS) are rectangular, so the functions that
create them have to know how many pixels (grids) there are in the left-right
and up-down directions.  The upper left corner of the grid is row=0 and
column=0, rows increase down the grid, and columns increase toward the
right.

NOTE: The latitude or longitude in all of these functions refers to the
center of the pixel (intersection of the grid).

NOTE: The North-at-Top Method (nortop) - The functions whose names begin
with "bldmds" build map data sets for arbitrary grids.  In the "nortop"
method the first three arguments provide geodetic latitude and longitude
of the center of the map set, and the approximate number of kilometers
the map set will span in the north south direction, along the center
longitude.  By definition in the "nortop" method, the center longitude
is straight up and down the center of the map set, and north is towards
the top of that line.
 
NOTE: The Upper-Left-Lower-Right method (ullr) - Other functions that begin
with "bldmds" build MDS by specifying the upper left and lower right
corners of the MDS.  These functions have a middle part to their names
that is "ullr"

NOTE: Geodetic vs Geocentric and the Oblate Spheroid.  The functions in
this library do not account for the oblatness of the Earth because it is
assumed that the latitudes passed to this library are all Geodetic.
Once you have geodetic latitude, you have accounted for the oblatness of
the Earth.  Therefore, all the calculations in this library are based on
a spherical Earth model.
 **********************/

/***************************************************************************
These are the meanings of the variable names used in the parameter lists
of the functions:


   int   function_name - returns an integer error flag
                 0 means no errors
                 <100 means a warning, err_string contains message
                 >100 means an important error, err_string contains message
   double  clat  - center geodetic latitude, in radians, positive north,
                    range is -pi/2 to +pi/2
   double  clon  - center longitude, in radians, positive east,
                    range is -pi to +pi
   double  ul_lat  - upper left corner geodetic latitude, in radians,
                    positive north, range is -pi/2 to +pi/2,
                    must be > radian equivalent of 30 degrees south
   double  ul_lon  - upper left corner longitude, in radians,
                    positive east, range is -pi to +pi
   double  lr_lat  - lower right corner geodetic latitude, in radians,
                    positive north, range is -pi/2 to +pi/2,
                    must be > radian equivalent of 30 degrees south
   double  lr_lon  - lower right corner longitude, in radians,
                    positive east, range is -pi to +pi
   double  ns_dist  - approximate north-south distance, spanned by the
                    MDS, along the center longitude, in kilometers
   double row  - row coordinate of the input point
   double col  - column coordinate of the input point
   double* rlat  - input point geodetic latitude pointer, in radians,
                    positive north, range is -pi/2 to +pi/2,
   double* rlon  - input point longitude pointer, in radians,
                    positive east, range is -pi to +pi
   int     lftrgt_pix  - number of pixels (grids) in the left to right
                    direction
   int     updwn_pix  - number of pixels (grids) in the up to down
                    direction of the rectangle.
   mds_type* dmds  - Output MDS data structure pointer, memory allocated
                    by the calling function
   char*  err_string  - string pointer (256 bytes) containing any
                    error messages created, memory allocated by the calling
                    function
 **********************/


/***************************************************************************
bldmds_nortop_polen : build an MDS, by the the nortop method, on a
northern polar stereographic map projection.

lerr = bldmds_nortop_polen(clat,clon,ns_dist,lftrgt_pix,updwn_pix,
                           dmds,err_string);
 **********************/

int bldmds_nortop_polen(double clat,
        double clon,
        double ns_dist,
        int lftrgt_pix,
        int updwn_pix,
        mds_type *dmds,
        char *err_string);

/***************************************************************************
bldmds_nortop_poles : build an MDS, by the the nortop method, on a
southern polar stereographic map projection.

lerr = bldmds_nortop_poles(clat,clon,ns_dist,lftrgt_pix,updwn_pix,
                           dmds,err_string);
 **********************/

int bldmds_nortop_poles(double clat,
        double clon,
        double ns_dist,
        int lftrgt_pix,
        int updwn_pix,
        mds_type *dmds,
        char *err_string);

/***************************************************************************
bldmds_nortop_merc: build an MDS, by the the nortop method, on a
mercator map projection.

lerr = bldmds_nortop_merc(clat,clon,ns_dist,lftrgt_pix,updwn_pix,
                          dmds,err_string);
 **********************/

int bldmds_nortop_merc(double clat,
        double clon,
        double ns_dist,
        int lftrgt_pix,
        int updwn_pix,
        mds_type *dmds,
        char *err_string);

/***************************************************************************
bldmds_ullr_polen : build an MDS, by the "ullr" method, on a northern
polar stereographic map projection.

lerr = bldmds_ullr_polen(ul_lat,ul_lon,lr_lat,lr_lon,lftrgt_pix,updwn_pix,
                         dmds, err_string);
 **********************/

int bldmds_ullr_polen(double ul_lat,
        double ul_lon,
        double lr_lat,
        double lr_lon,
        int lftrgt_pix,
        int updwn_pix,
        mds_type *dmds,
        char *err_string);

/***************************************************************************
bldmds_ullr_poles : build an MDS, by the "ullr" method, on a southern
polar stereographic map projection.

lerr = bldmds_ullr_poles(ul_lat,ul_lon,lr_lat,lr_lon,lftrgt_pix,updwn_pix,
                         dmds, err_string);
 **********************/

int bldmds_ullr_poles(double ul_lat,
        double ul_lon,
        double lr_lat,
        double lr_lon,
        int lftrgt_pix,
        int updwn_pix,
        mds_type *dmds,
        char *err_string);

/***************************************************************************
bldmds_ullr_merc : build an MDS, by the "ullr" method, on a mercator map
projection.

NOTE: The lower right longitude is deliberately not present.  Once the
lower right latitude is defined, the lower right longitude is no longer
arbitrary.  It is determined by the upper left corner and the
dimensions of the grid.  We could have allowed specification of the
lower right longitude instead of the latitude, but we arbitrarily chose
the latitude.

lerr = bldmds_ullr_merc(ul_lat,ul_lon,lr_lat,lftrgt_pix,updwn_pix,
                        dmds, err_string);
 **********************/

int bldmds_ullr_merc(double ul_lat,
        double ul_lon,
        double lr_lat,
        int lftrgt_pix,
        int updwn_pix,
        mds_type* dmds,
        char* err_string);

/***************************************************************************
grid_to_latlon : convert an input grid coordinate row and column to latitude
and longitude.

lerr = grid_to_latlon(row,col,imds,&rlat,&rlon,err_string);
 **********************/

int grid_to_latlon(double row,
        double col,
        mds_type *imds,
        double *rlat,
        double *rlon,
        char *err_string);

/***************************************************************************
latlon_to_grid : convert an input latitude and longitude to grid coordinate
row and column.  The output coordinates are doubles, but when an integer
grid position is needed, the calling function would simply apply the rint
or round function to the coordinates.

lerr = latlon_to_grid(rlat,rlon,imds,&row,&col,err_string);
 **********************/

int latlon_to_grid(double rlat,
        double rlon,
        mds_type *imds,
        double *row,
        double *col,
        char *err_string);

/***************************************************************************
ll2g_polen : converts a latitude and longitude to a row and column 
position on a northern polar stereographic MDS.

lerr = ll2g_polen(rlat,rlon,imds,&row,&col,err_string);
 **********************/

int ll2g_polen(double rlat,
        double rlon,
        mds_type *imds,
        double *row,
        double *col,
        char *err_string);

/***************************************************************************
ll2g_poles : converts a latitude and longitude to a row and column 
position on a southern polar stereographic MDS.

lerr = ll2g_poles(rlat,rlon,imds,&row,&col,err_string);
 **********************/

int ll2g_poles(double rlat,
        double rlon,
        mds_type *imds,
        double *row,
        double *col,
        char *err_string);

/***************************************************************************
ll2g_merc : converts a latitude and longitude to a row and column 
position on a mercator MDS.

lerr = ll2g_merc(rlat,rlon,imds,&row,&col,err_string);
 **********************/

int ll2g_merc(double rlat,
        double rlon,
        mds_type *imds,
        double *row,
        double *col,
        char *err_string);

/***************************************************************************
g2ll_merc : converts a grid position on a mercator MDS to latitude and 
longitude.

lerr = g2ll_merc(row,col,imds,&rlat,&rlon,err_string);
 **********************/

int g2ll_merc(double row,
        double col,
        mds_type *imds,
        double *rlat,
        double *rlon,
        char *err_string);

/***************************************************************************
g2ll_polen : converts a grid position on a northern polar stereographic MDS
to latitude and longitude.

lerr = g2ll_polen(row,col,imds,&rlat,&rlon,err_string);
 **********************/

int g2ll_polen(double row,
        double col,
        mds_type *imds,
        double *rlat,
        double *rlon,
        char *err_string);

/***************************************************************************
g2ll_poles : converts a grid position on a southern polar stereographic MDS
to latitude and longitude.

lerr = g2ll_poles(row,col,imds,&rlat,&rlon,err_string);
 **********************/

int g2ll_poles(double row,
        double col,
        mds_type *imds,
        double *rlat,
        double *rlon,
        char *err_string);

/***************************************************************************
fllr_merc : finds and sets the ranges of latitude and longitude in the input
mercator MDS.

lerr = fllr_merc(imds,err_string);
 **********************/

int fllr_merc(mds_type *imds,
        char *err_string);

/***************************************************************************
fllr_polar : finds and sets the ranges of latitude and longitude in the
input polar stereographic MDS.
 
lerr = fllr_polar(imds,err_string);
 **********************/

int fllr_polar(mds_type *imds,
        char *err_string);

/***************************************************************************
NOTE: "ced" is cylindrical equidistant.
 **********************/

/***************************************************************************
bldmds_nortop_ced: build an MDS, by the the nortop method, on a
cylindrical equidistant map projection.

lerr = bldmds_nortop_ced(clat,clon,ns_dist,lftrgt_pix,updwn_pix,
                         dmds,err_string);
 **********************/

int bldmds_nortop_ced(double clat,
        double clon,
        double ns_dist,
        int lftrgt_pix,
        int updwn_pix,
        mds_type *dmds,
        char *err_string);

/***************************************************************************
bldmds_ullr_ced: build an MDS, by the "ullr" method, on a cylindrical
equidistant map projection.

lerr = bldmds_ullr_ced(ul_lat,ul_lon,lr_lat,lr_lon,lftrgt_pix,updwn_pix,
                       dmds, err_string);
 **********************/

int bldmds_ullr_ced(double ul_lat,
        double ul_lon,
        double lr_lat,
        double lr_lon,
        int lftrgt_pix,
        int updwn_pix,
        mds_type* dmds,
        char* err_string);

/***************************************************************************
fllr_ced: finds and sets the ranges of latitude and longitude in the input
cylindrical equidistant MDS.

lerr = fllr_ced(imds,err_string);
 **********************/

int fllr_ced(mds_type *imds,
        char *err_string);

/***************************************************************************
g2ll_ced : converts a grid position on a cylindrical equidistant MDS to 
latitude and longitude.

lerr = g2ll_ced(row,col,imds,&rlat,&rlon,err_string);
 **********************/

int g2ll_ced(double row,
        double col,
        mds_type *imds,
        double *rlat,
        double *rlon,
        char *err_string);

/***************************************************************************
ll2g_ced : converts a latitude and longitude to a row and column 
position on a cylindrical equidistant MDS.

lerr = ll2g_ced(rlat,rlon,imds,&row,&col,err_string);
 **********************/

int ll2g_ced(double rlat,
        double rlon,
        mds_type *imds,
        double *row,
        double *col,
        char *err_string);

/***************************************************************************
NOTE: "polazn" or "polazs" are related to Polar Azimuthal Equidistant
 **********************/

/***************************************************************************
bldmds_ullr_polazn: build an MDS, by the "ullr" method, on a northern polar
azimuthal equidistant map projection.

lerr = bldmds_ullr_polazn(ul_lat,ul_lon,lr_lat,lr_lon,lftrgt_pix,updwn_pix,
                          dmds, err_string);
 **********************/

int bldmds_ullr_polazn(double ul_lat,
        double ul_lon,
        double lr_lat,
        double lr_lon,
        int lftrgt_pix,
        int updwn_pix,
        mds_type* dmds,
        char* err_string);

/***************************************************************************
bldmds_ullr_polazs: build an MDS, by the "ullr" method, on a southern polar
azimuthal equidistant map projection.

lerr = bldmds_ullr_polazs(ul_lat,ul_lon,lr_lat,lr_lon,lftrgt_pix,updwn_pix,
                          dmds, err_string);
 **********************/

int bldmds_ullr_polazs(double ul_lat,
        double ul_lon,
        double lr_lat,
        double lr_lon,
        int lftrgt_pix,
        int updwn_pix,
        mds_type* dmds,
        char* err_string);

/***************************************************************************
bldmds_nortop_polazn: build an MDS, by the the nortop method, on a
northern polar azimuthal equidistant map projection.

lerr = bldmds_nortop_polazn(clat,clon,ns_dist,lftrgt_pix,updwn_pix,
                            dmds,err_string);
 **********************/

int bldmds_nortop_polazn(double clat,
        double clon,
        double ns_dist,
        int lftrgt_pix,
        int updwn_pix,
        mds_type* dmds,
        char* err_string);

/***************************************************************************
bldmds_nortop_polazs: build an MDS, by the the nortop method, on a
southern polar azimuthal equidistant map projection.

lerr = bldmds_nortop_polazs(clat,clon,ns_dist,lftrgt_pix,updwn_pix,
                            dmds,err_string);
 **********************/

int bldmds_nortop_polazs(double clat,
        double clon,
        double ns_dist,
        int lftrgt_pix,
        int updwn_pix,
        mds_type* dmds,
        char* err_string);

/***************************************************************************
g2ll_polazn : converts a grid position on a northern polar azimuthal
equidistant MDS to latitude and longitude.

lerr = g2ll_polazn(row,col,imds,&rlat,&rlon,err_string);
 **********************/

int g2ll_polazn(double row,
        double col,
        mds_type *imds,
        double *rlat,
        double *rlon,
        char *err_string);

/***************************************************************************
g2ll_polazs : converts a grid position on a southern polar azimuthal
equidistant MDS to latitude and longitude.

lerr = g2ll_polazs(row,col,imds,&rlat,&rlon,err_string);
 **********************/

int g2ll_polazs(double row,
        double col,
        mds_type *imds,
        double *rlat,
        double *rlon,
        char *err_string);

/***************************************************************************
ll2g_polazn : converts a latitude and longitude to a row and column 
position on a northern polar azimuthal equidistant MDS.

lerr = ll2g_polazn(rlat,rlon,imds,&row,&col,err_string);
 **********************/

int ll2g_polazn(double rlat,
        double rlon,
        mds_type *imds,
        double *row,
        double *col,
        char *err_string);

/***************************************************************************
ll2g_polazs : converts a latitude and longitude to a row and column 
position on a southern polar azimuthal equidistant MDS.

lerr = ll2g_polazs(rlat,rlon,imds,&row,&col,err_string);
 **********************/

int ll2g_polazs(double rlat,
        double rlon,
        mds_type *imds,
        double *row,
        double *col,
        char *err_string);

/***************************************************************************
end of include file
 **********************/

#ifdef __cplusplus
}
#endif

#endif
