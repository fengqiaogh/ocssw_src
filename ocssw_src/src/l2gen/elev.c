#include "l12_proto.h"
#include "nc_gridutils.h"

/* global variables */
static const char* elevNames[] = {"height", "z", "depth", NULL};
static int elev_initialized = 0;
static grid_info_t* elev_global_gridinfo = {0};
static grid_info_t* elev_aux_gridinfo = {0};
static int elev_global_depthmode = -1;
static int elev_aux_depthmode = -1;

int initElevFile(char* elevFilename, grid_info_t* elevGrid, int* depth_mode) {
    int status;
    nc_type type;
    size_t len;

    /* initial file tests */
    if (elevFilename == NULL) {
        fprintf(stderr, "-E- %s line %d: Elevation file is NULL.\n",
                __FILE__, __LINE__);
        return 1;
    }
    status = nc_open(elevFilename, NC_NOWRITE, &elevGrid->ncid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: Could not open netCDF File \"%s\".\n",
                __FILE__, __LINE__, elevFilename);
        return 1;
    }

    /* allocate and load grid info */
    status = init_gridinfo(elevFilename, elevNames, elevGrid);
    if (status != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: Could not read grid info from \"%s\".\n",
                __FILE__, __LINE__, elevFilename);
        return status;
    }

    /* read depth_mode attribute */
    /*   positive = "up"    means pos values for mountains   */
    /*   positive = "down"  means pos values for ocean depth */
    status = nc_inq_att(elevGrid->ncid, elevGrid->varid,
            "positive", &type, &len);
    if ((status == NC_NOERR) && (type == NC_CHAR) && (len > 0)) {
        char* value = (char*) malloc(len + 1);
        status = nc_get_att_text(elevGrid->ncid, elevGrid->varid,
                "positive", value);
        value[len] = '\0'; // null terminate
        lowcase(value);
        depth_mode[0] = (strcmp(value, "down") == 0); // 0=up, 1=down
        free(value);
    }        /* otherwise get depth_mode from variable name */
    else {
        /* (commented out: not valid for SWIM input)
            char name[NC_MAX_NAME];
            status = nc_inq_varname(elevGrid->ncid, elevGrid->varid, name);
            printf("varname = %s\n",name);
            depth_mode[0] = (strcmp(name, "depth") == 0);
         */
        status = 0;
        depth_mode[0] = 0;
    }
    return status;
}

void elev_init(char* elevGlobalFilename, char* elevAuxFilename) {
    int status;
    if (elevGlobalFilename != NULL && elevGlobalFilename[0] != 0) {
        elev_global_gridinfo = allocate_gridinfo();
        status = initElevFile(elevGlobalFilename,
                elev_global_gridinfo,
                &elev_global_depthmode);
        if (status) {
            fprintf(stderr, "-E- %s line %d: "
                    "Could not read Global Elevation file \"%s\".\n",
                    __FILE__, __LINE__, elevGlobalFilename);
            exit(1);
        }
    }
    if (elevAuxFilename != NULL && elevAuxFilename[0] != 0) {
        elev_aux_gridinfo = allocate_gridinfo();
        status = initElevFile(elevAuxFilename,
                elev_aux_gridinfo,
                &elev_aux_depthmode);
        if (status) {
            fprintf(stderr, "-E- %s line %d: "
                    "Could not read Auxiliary Elevation file \"%s\".\n",
                    __FILE__, __LINE__, elevAuxFilename);
            exit(1);
        }
    }
    elev_initialized = 1;
}

float readElev(grid_info_t* elevGrid, int depth_mode, float lat, float lon) {
    int status;
    double elev = BAD_FLT;

    if (elevGrid != NULL) {
        status = get_bylatlon(elevGrid, lat, lon, &elev);
        if ((!status) && (depth_mode))
            elev *= -1;
    }

    return (float) elev;
}

float get_elev(float lat, float lon) {
    float elev;
    elev = readElev(elev_aux_gridinfo, elev_aux_depthmode, lat, lon);
    if (elev == BAD_FLT) {
        elev = readElev(elev_global_gridinfo, elev_global_depthmode, lat, lon);
    }
    /* printf(" = %d (get_elev)\n",elev); */
    return elev;
}
