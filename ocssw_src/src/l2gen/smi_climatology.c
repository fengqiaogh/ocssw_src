#include "l12_proto.h"
#include "smi_climatology.h"
#include "allocate2d.h"

#define NX 4096
#define NY 2048

static uint8 **map [NPROD];
static float slope [NPROD];
static float offset[NPROD];
static int loaded[NPROD];

int smi_climatology_init(char *file, int day, int prodID) {
    int32 sd_id;
    int32 sds_id;
    int32 status;
    int32 sds_index;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int32 start[2];
    int32 edges[2];
    char name[H4_MAX_NC_NAME];
    char sdsname[256];

    int16 mon = (int) day / 31; /* month of year (no need for perfection) */

    /* Make sure this product was not already loaded */
    if (loaded[prodID])
        return (0);

    printf("Loading climatology file %s\n", file);

    if (prodID < 0 || prodID >= NPROD) {
        printf("-E- %s:  Invalid SMI product ID %d.\n",
                __FILE__, prodID);
        return (1);
    }


    /* Allocate space for requested product */
    if (map[prodID] == NULL) {
        map[prodID] = allocate2d_uchar(NY, NX);
        if (map[prodID] == NULL) {
            printf("-E- %s:  Error allocating space for %s.\n",
                    __FILE__, file);
            return (1);
        }
    }


    /* Open the file and initiate the SD interface */
    sd_id = SDstart(file, DFACC_RDONLY);
    if (sd_id == -1) {
        printf("-E- %s:  Error opening file %s.\n",
                __FILE__, file);
        return (1);
    }


    /* Get the SDS index */
    sprintf(sdsname, "data%2.2i", mon + 1);
    sds_index = SDnametoindex(sd_id, sdsname);
    if (sds_index == -1) {
        printf("-E- %s:  Error seeking %s SDS from %s.\n",
                __FILE__, sdsname, file);
        return (1);
    }

    /* Select the SDS */
    sds_id = SDselect(sd_id, sds_index);

    /* Verify the characteristics of the array */
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS info for %s from %s.\n",
                __FILE__, sdsname, file);
        return (1);
    }
    if (dims[0] != NY || dims[1] != NX) {
        printf("-E- %s:  Dimension mis-match on %s array from %s.\n",
                __FILE__, sdsname, file);
        printf("  Expecting %d x %d\n", NX, NY);
        printf("  Reading   %d x %d\n", dims[1], dims[0]);
        return (1);
    }

    start[0] = 0; /* row offset */
    start[1] = 0; /* col offset */
    edges[0] = NY; /* row count  */
    edges[1] = NX; /* col count  */

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) * map[prodID]);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS %s from %s.\n",
                __FILE__, sdsname, file);
        return (1);
    }

    /* Read the scaling info */
    status = SDreadattr(sds_id, SDfindattr(sds_id, "slope"), (VOIDP) & slope [prodID]);
    if (status != 0) {
        printf("-E- %s:  Error reading slope attribute for SDS %s from %s.\n",
                __FILE__, sdsname, file);
        return (1);
    }
    status = SDreadattr(sds_id, SDfindattr(sds_id, "intercept"), (VOIDP) & offset[prodID]);
    if (status != 0) {
        printf("-E- %s:  Error reading intercept attribute for SDS %s from %s.\n",
                __FILE__, sdsname, file);
        return (1);
    }

    /* Terminate access to the array */
    status = SDendaccess(sds_id);

    /* Terminate access to the SD interface and close the file */
    status = SDend(sd_id);

    loaded[prodID] = 1;

    return (0);
}

float smi_climatology(float lon, float lat, int prodID) {
    int16 i, j;

    i = MAX(MIN((int16) ((lon + 180.0) / 360.0 * NX), NX - 1), 0);
    j = MAX(MIN((int16) ((90.0 - lat) / 180.0 * NY), NY - 1), 0);

    return (map[prodID][j][i] * slope[prodID] + offset[prodID]);
}


