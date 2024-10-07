#include "l12_proto.h"
#include "nc_gridutils.h"

/* global variables */
static const char* landnames[] = {"watermask", "landmask", "z", NULL};
static const char* dem_names[] = {"height", "z", "depth", NULL};
static const char* watersurface_names[] = {"water_surface_height", NULL};


static grid_info_t* landmask_grid = {0};
static grid_info_t* dem_grid={0};
static grid_info_t* dem_aux_grid = {0};

static grid_info_t* watersurface_grid={0};

// support for legacy b128 land/shallow water mask
static int landindex = 0;
static int bathindex = 1;
static int use_b128_land = FALSE;
static int use_b128_bath = FALSE;

/* Land Mask tools */
int land_mask_init() {
    int status;

    /* load structure from NetCDF file*/
    landmask_grid = allocate_gridinfo();
    status = init_gridinfo(input->land, landnames, landmask_grid);

    /* if not NetCDF, free structure and read binary file */
    if (status != NC_NOERR) {
        free(landmask_grid);
        landmask_grid = NULL;
        status = b128_msk_init(input->land, landindex);
        if (status == 0)
            use_b128_land = TRUE;
    }
        
    return status;

}

int dem_init() {
    int status;
    dem_grid = allocate_gridinfo();
    status = init_gridinfo(input->demfile, dem_names, dem_grid);
    if (status == NC_NOERR) {
        watersurface_grid = allocate_gridinfo();
        status = init_gridinfo(input->demfile, watersurface_names, watersurface_grid);
        if (status != NC_NOERR) {
            free(watersurface_grid);
            watersurface_grid = NULL;
            printf("-E-  Unable to initialize grid for water surface height \n");
            return EXIT_FAILURE;
        }
        if (input->dem_auxfile != NULL && input->dem_auxfile[0] != 0) {
            dem_aux_grid = allocate_gridinfo();
            status = init_gridinfo(input->dem_auxfile, dem_names, dem_aux_grid);

            if (status) {
                printf("-E- Unable to initialize grid for auxiliary digital elevation\n");
                exit(1);
            }
        }
    } else {
        free(dem_grid);
        dem_grid = NULL;
        status = b128_msk_init(input->water, bathindex);
        use_b128_bath = TRUE;
    }

    if (status != 0){
        printf("-E- %s : Unable to initialize digital elevation map\n", __FILE__);
        return EXIT_FAILURE;
    }

    return status;
}

int land_mask(float lat, float lon) {
    double value;
    int status;
    static int messagePrinted = 0;
    /* get value from NetCDF file*/
    if (landmask_grid != NULL) {
        status = get_bylatlon(landmask_grid, lat, lon, &value);
        if (status) {
            if (!messagePrinted) {
                printf("-W- file contains locations not contained in the landmask: %s\n...assuming those are water.\n",
                        landmask_grid->file);
                messagePrinted = 1;
            }
            return EXIT_SUCCESS;
        }
        return ( (short) value != 1); // convert from water=1 to land=1
    }
        /* otherwise get from binary file */
    else return b128_msk_get(lat, lon, landindex);
}

/* Bathymetry placeholders */
int bath_mask_init(char *file) {
    int status;
    status = b128_msk_init(file, bathindex);
    return status;
}

int bath_mask(float lat, float lon) {
    return b128_msk_get(lat, lon, bathindex);
}

float get_dem(float lat, float lon) {
    int status;
    double dem;
    status = get_bylatlon(dem_aux_grid, lat, lon, &dem);
    if ((status != 0) || (dem == BAD_FLT)) {
        status = get_bylatlon(dem_grid, lat, lon, &dem);
    }
    return (float) dem;
}

/*
 * Set land and shallow water flags
 */
int land_bath_mask(l1str *l1rec,int32_t ip)
{
    double value;
    int status=0;
    static int landMessagePrinted = 0;
    static int bathMessagePrinted = 0;

    double height_floor, height_surface;
    height_floor = (double) l1rec->dem[ip];

    float lat=l1rec->lat[ip];
    float lon=l1rec->lon[ip];

    if (use_b128_land) {
        l1rec->land[ip] = b128_msk_get(lat, lon, landindex);
    } else {
        /* get value from NetCDF file*/
        if (landmask_grid != NULL) {
            status = get_bylatlon(landmask_grid, lat, lon, &value);
            if (status) {
                if (!landMessagePrinted) {
                    fprintf(stderr, "-W- file contains locations not contained in the landmask: %s\n...assuming those are water.\n",
                            landmask_grid->file);
                    landMessagePrinted = 1;
                }
                return EXIT_SUCCESS;
            }
            l1rec->land[ip]=( (short) value != 1); // convert from water=1 to land=1
        }
    }
    //if land mask is off, set bathymetry flags
    if(!l1rec->land[ip]){
        if (use_b128_bath) {
            l1rec->swater[ip] = l1rec->land[ip] = b128_msk_get(lat, lon, bathindex);
            if (input->shallow_water_depth != 30) {
                if (!bathMessagePrinted) {
                    fprintf(stderr, "-W- shallow water value (%6.2f) is not consistent with the bathymetry file used: shallow water flag set to 30m",
                            input->shallow_water_depth);
                    bathMessagePrinted = 1;
                }
            } 
        } else {
            //status = get_bylatlon(dem_grid, lat, lon, &height_floor);
            status = get_bylatlon(watersurface_grid, lat, lon, &height_surface);

            if(fabs(height_floor-dem_grid->FillValue)>DBL_EPSILON && fabs(height_surface-watersurface_grid->FillValue)>DBL_EPSILON){
                if((height_surface-height_floor)<input->shallow_water_depth)
                    l1rec->swater[ip]=ON;
            }
        }
    }

    return EXIT_SUCCESS;
}

