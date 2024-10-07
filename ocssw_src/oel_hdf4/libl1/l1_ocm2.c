#include "l1_ocm2.h"
#include "l1.h"
#include <mfhdf.h>
#include <allocate2d.h>

static char dname[8][6] = {"L_412", "L_443", "L_490", "L_510",
    "L_555", "L_620", "L_740", "L_865"};

static float *data = NULL;
static float32 **solz = NULL;
static float32 **sola = NULL;
static float32 **senz = NULL;
static float32 **sena = NULL;


#define READ_GLBL_ATTR(nam,ptr) {                                           \
  if(SDreadattr(sd_id,SDfindattr(sd_id,(nam)),(VOIDP)(ptr))){               \
    fprintf(stderr,                                                         \
    "-E- %s line %d: Could not get global attribute, %s.\n", \
    __FILE__,__LINE__,(nam));        \
  }                                                                         \
}

#define READ_SDS(nam,ptr,s0,s1,s2,e0,e1,e2) {                           \
  int32 start[3];                                                       \
  int32 edge[3];                                                        \
  edge[0]=(e0); edge[1]=(e1); edge[2]=(e2);                             \
  start[0]=(s0); start[1]=(s1); start[2]=(s2);                          \
  if(SDreaddata(SDselect(sd_id, SDnametoindex(sd_id, (nam))),           \
  start, NULL, edge, (VOIDP)(ptr)) == FAIL){                            \
    fprintf(stderr,"-E- %s line %d: Could not read SDS, %s.\n",         \
    __FILE__,__LINE__,(nam));                                           \
  }                                                                     \
}

void interp_ocm2_geo(int32 sd_id, int32 nx, int32 ny, char *sdsname, float32 **geo) {
    int32 sds_id;
    int32 rank;
    int32 dims[3];
    int32 type;
    int32 numattr;
    int32 i, j, jout;
    int32 mx, my;
    int32 nx_ctl;
    int32 ny_ctl;

    float32 *x_ctl, *dx_ctl;
    float32 *y_ctl, *dy_ctl;
    float32 *v_ctl;

    // determine size of path geometry control points        
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, (sdsname)));
    if (SDgetinfo(sds_id, NULL, &rank, dims, &type, &numattr) == -1) {
        fprintf(stderr, "-E- %s line %d: error getting dimension info.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    nx_ctl = dims[1];
    ny_ctl = dims[0];

    // set-up interpolation grid
    mx = nx / nx_ctl;
    my = ny / ny_ctl;
    x_ctl = (float32 *) calloc(nx_ctl, sizeof (float32));
    y_ctl = (float32 *) calloc(ny_ctl, sizeof (float32));
    for (i = 0; i < nx_ctl; i++) x_ctl[i] = (float32) i * mx + 1;
    for (j = 0; j < ny_ctl; j++) y_ctl[j] = (float32) j * my + 1;

    // interpolate each line to full pixel 
    v_ctl = (float32 *) calloc(nx_ctl, sizeof (float32));
    dx_ctl = (float32 *) calloc(nx_ctl, sizeof (float32));
    for (j = 0; j < ny_ctl; j++) {
        jout = j * my + 1;
        READ_SDS(sdsname, v_ctl, j, 0, 0, 1, nx_ctl, 1);
        spline(x_ctl, v_ctl, nx_ctl, 1e30, 1e30, dx_ctl);
        for (i = 0; i < nx; i++) {
            splint(x_ctl, v_ctl, dx_ctl, nx_ctl, (float) i, &geo[jout][i]);
            //printf("%d %d %f\n",jout,i,geo[jout][i]);
        }
    }
    free(v_ctl);
    free(dx_ctl);

    // interpolate each line to full pixel 
    v_ctl = (float32 *) calloc(ny_ctl, sizeof (float32));
    dy_ctl = (float32 *) calloc(ny_ctl, sizeof (float32));
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny_ctl; j++) {
            jout = j * my + 1;
            v_ctl[j] = geo[jout][i];
        }
        spline(y_ctl, v_ctl, ny_ctl, 1e30, 1e30, dy_ctl);
        for (j = 0; j < ny; j++) {
            splint(y_ctl, v_ctl, dy_ctl, ny_ctl, (float) j, &geo[j][i]);
        }
    }
    free(v_ctl);
    free(dy_ctl);

    free(x_ctl);
    free(y_ctl);
}

int openl1_ocm2(filehandle *file) {
    int32 npix;
    int32 nscan;
    int32 subsamp;
    int32 sd_id;

    /* Open the HDF input file */
    sd_id = SDstart(file->name, DFACC_RDONLY);
    if (sd_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, file->name, DFACC_RDONLY);
        return (HDF_FUNCTION_ERROR);
    }

    /* Read some of the level-1A global attributes. */
    READ_GLBL_ATTR("Pixels per Scan Line", &npix);
    READ_GLBL_ATTR("Number of Scan Lines", &nscan);
    READ_GLBL_ATTR("LAC Pixel Subsampling", &subsamp);

    file->npix = npix;
    file->nscan = nscan;
    file->nbands = 8;
    file->sd_id = sd_id;
    file->sensorID = OCM2;
    if (subsamp == 4)
        strcpy(file->spatialResolution, "4 km");
    else
        strcpy(file->spatialResolution, "360 m");

    return (LIFE_IS_GOOD);
}

int readl1_ocm2(filehandle *file, int32_t scan, l1str *l1rec) {
    static int firstCall = 1;

    int32 sd_id = file->sd_id;
    int32 npix = file->npix;
    int32 nscan = file->nscan;

    int32 ip, ib;

    if (firstCall) {

        firstCall = 0;

        if ((data = (float *) calloc(file->npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        // allocate space for radiant path geometries

        if ((solz = allocate2d_float(nscan, npix)) == NULL) {
            printf("-E- %s line %d: Error allocating geolocation space.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if ((sola = allocate2d_float(nscan, npix)) == NULL) {
            printf("-E- %s line %d: Error allocating geolocation space.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if ((senz = allocate2d_float(nscan, npix)) == NULL) {
            printf("-E- %s line %d: Error allocating geolocation space.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if ((sena = allocate2d_float(nscan, npix)) == NULL) {
            printf("-E- %s line %d: Error allocating geolocation space.\n",
                    __FILE__, __LINE__);
            return (1);
        }


        printf("Interpolating radiant path geometries\n");
        interp_ocm2_geo(sd_id, npix, nscan, "solz", solz);
        interp_ocm2_geo(sd_id, npix, nscan, "sola", sola);
        interp_ocm2_geo(sd_id, npix, nscan, "senz", senz);
        interp_ocm2_geo(sd_id, npix, nscan, "sena", sena);
    }

    int32_t year, day, msec;
    READ_SDS("year", &year, scan, 0, 0, 1, 1, 1);
    READ_SDS("day", &day, scan, 0, 0, 1, 1, 1);
    READ_SDS("msec", &msec, scan, 0, 0, 1, 1, 1);
    READ_SDS("longitude", l1rec->lon, scan, 0, 0, 1, npix, 1);
    READ_SDS("latitude", l1rec->lat, scan, 0, 0, 1, npix, 1);
    l1rec->scantime = yds2unix(year, day, (double) (msec / 1.e3));

    for (ip = 0; ip < npix; ip++) {
        l1rec->solz[ip] = solz[scan][ip];
        l1rec->sola[ip] = sola[scan][ip];
        l1rec->senz[ip] = senz[scan][ip];
        l1rec->sena[ip] = sena[scan][ip];
    }

    for (ib = 0; ib < file->nbands; ib++) {
        READ_SDS(dname[ib], data, scan, 0, 0, 1, npix, 1);
        for (ip = 0; ip < npix; ip++) {
            l1rec->Lt[ip * file->nbands + ib] = data[ip];
            //            if (data[ip] > 32767) l1rec->hilt[ip] = 1;
        }
    }

    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

int closel1_ocm2(filehandle *file) {
    if (SDend(file->sd_id)) {
        fprintf(stderr, "-E- %s line %d: SDend(%d) failed for file, %s.\n",
                __FILE__, __LINE__, file->sd_id, file->name);
        return (HDF_FUNCTION_ERROR);
    }

    return (LIFE_IS_GOOD);
}





