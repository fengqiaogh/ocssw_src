
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <hdf4utils.h>
#include "l1.h"
#include "filehdr_struc.h"
extern "C" {
    
    #include "l1_octs.h"
}
#include "l1_octs_netcdf.h"
#include <netcdf>
#include <iostream>


#include <hdf.h>
#include <mfhdf.h>

#define RADIUS 6378.137    /* Earth radius in km */
#define FL 1.0/298.257     /* Earth flatening factor */
#define NR 560             /* # rows */      
#define NC 30              /* max # columns */     

#define MAXOCLIN 6700      /* max # lines */
#define MAXOCPIX 2218      /* max # pixels */
#define MAXOCARR 10000
#define NOCBANDS 8


using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

static int32_t msec_start;
static int32_t msec[MAXOCLIN];
static float tilt[MAXOCLIN];
static int16_t year, day, nline, npix, sline;
static int32_t spix;
static int16_t maxday;
static float solz1, sola1, senz1, sena1;
static float *lon, *lat, *senz, *sena, *solz, *sola;
static int32_t ncol;
static int32_t scansPerScene;
static int32_t linesPerScan;
static uint16_t *tilt_flag;
static int16_t *gain;
static float *inst_temp;
static int16_t *miss_qual;
static int16_t num_tab;
static int32_t extr_line_offset = 0;
static int16_t *start_scan;
static int16_t last_table = 0;
static int16_t sample_table[3][8][2][400][2];

/* ------------------------------------------------------ */
/* CalcViewAngle() - calculates the solar and sensor view */
/*             angles from one geodetic longitude and     */
/*             latitude, position vector, and solar unit  */
/*             vector.                                    */
/*                                                        */
/* INPUT       DESCRIPTION                                */
/* -----       -----------                                */
/*  lon        longitude of pixel in degrees              */
/*  lat        latitude of pixel in degrees               */
/*  pos[3]     vector from earth center to spacecraft in  */
/*             ECR, km                                    */
/*  usun[3]    unit vector from earth center to sun in    */
/*             ECR                                        */
/*                                                        */
/* OUTPUT      DESCRIPTION                                */
/* ------      -----------                                */
/*  solz       solar zenith angle of pixel in degrees     */
/*  sola       solar azimuth angle of pixel in degrees    */
/*  senz       sensor zenith angle of pixel in degrees    */
/*  sena       sensor azimuth angle of pixel in degrees   */
/*                                                        */

/* ------------------------------------------------------ */
static int CalcViewAngle_nc(float lon1, float lat1, float pos[3],
        float usun[3]) {
    float rmat[3][3]; /* rotation matrix */
    float rlon, rlat;
    float up[3], upxy, ea[3], no[3];
    float phi, R;
    float gvec[3], scvec[3], senl[3], sunl[3];
    int i, j;

    rlon = lon1 * PI / 180.0;
    rlat = lat1 * PI / 180.0;

    /* First, we must define the axes (in ECR space) of a
       pixel-local coordinate system, with the z-axis along
       the geodetic pixel normal, x-axis pointing east, and
       y-axis pointing north.  */
    up[0] = cos(rlat) * cos(rlon);
    up[1] = cos(rlat) * sin(rlon);
    up[2] = sin(rlat);
    upxy = sqrt(up[0] * up[0 ] + up[1] * up[1]);
    ea[0] = -up[1] / upxy;
    ea[1] = up[0] / upxy;
    ea[2] = 0.0;
    /* calculate cross product of up and ea */
    no[0] = up[1] * ea[2] - ea[1] * up[2];
    no[1] = up[2] * ea[0] - ea[2] * up[0];
    no[2] = up[0] * ea[1] - ea[0] * up[1];


    /* Now we need the pixel-to-spacecraft vector in ECR,
       Compute geocentric pixel location vector in km and subtract
       from spacecraft position. */
    /* geocentric latitude */
    phi = atan(tan(rlat)*(1 - FL)*(1 - FL));
    /* dist to Earth surface */
    R = RADIUS * (1 - FL) / sqrt(1 - (2 - FL) * FL * (cos(phi) * cos(phi)));
    gvec[0] = R * cos(phi) * cos(rlon);
    gvec[1] = R * cos(phi) * sin(rlon);
    gvec[2] = R * sin(phi);
    for (i = 0; i < 3; i++) {
        scvec[i] = pos[i] - gvec[i];
    }

    /* Now we can transform the pixel-to-spacecraft and Sun
       vectors into the local frame. */
    for (i = 0; i < 3; i++) {
        rmat[0][i] = ea[i];
        rmat[1][i] = no[i];
        rmat[2][i] = up[i];
    }
    for (i = 0; i < 3; i++) {
        senl[i] = 0;
        sunl[i] = 0;
        for (j = 0; j < 3; j++) {
            senl[i] = senl[i] + rmat[i][j] * scvec[j];
            sunl[i] = sunl[i] + rmat[i][j] * usun[j];
        }
    }

    /* Compute the solar zenith and azimuth */
    solz1 = RADEG * atan(sqrt(sunl[0] * sunl[0] + sunl[1] * sunl[1]) / sunl[2]);
    if (solz1 < 0.0) solz1 += 180.0;

    if (solz1 > 0.05) {
        sola1 = RADEG * (atan2(sunl[0], sunl[1]));
    } else {
        sola1 = 0.0;
    }
    if (sola1 < 0.0) {
        sola1 = sola1 + 360.0;
    }

    /* Compute the sensor zenith and azimuth  */
    senz1 = RADEG * atan(sqrt(senl[0] * senl[0] + senl[1] * senl[1]) / senl[2]);
    if (senz1 < 0.0) senz1 += 180.0;
    if (senz1 > 0.05) {
        sena1 = RADEG * (atan2(senl[0], senl[1]));
    } else {
        sena1 = 0.0;
    }
    if (sena1 < 0.0) {
        sena1 = sena1 + 360.0;
    }
    return (0);
}



/* ------------------------------------------------------ */
/* navigation() - generates navigation info interpolated  */
/*              for each pixel, using info from a NASDA   */
/*              OCTS HDF file                             */
/*                                                        */

/* ------------------------------------------------------ */
extern "C" void navigation_netcdf(NcFile &dataFile) {
    int16_t det;

    int16_t *pxl;
    float *xctl;
    float *inlon, *inlat;
    float *insolz, *insola;
    float *insenz, *insena;
    float *usun, *pos;
    int32_t *row;
    int32_t *indx;
    float *yctl;
    float *in1, *in2;

    float iusun[3], ipos[3];
    float usun_sum;
    int32_t i, j, k, l;
    int32_t eline, nlon, nlat;
    int32_t ilat[MAXOCLIN], ilon[MAXOCLIN];
    float out1[MAXOCLIN], out2[MAXOCLIN];
    int32_t jout;
    float *spl_aux;
    float *x_ctl_ll, *y_ctl_ll, *z_ctl_ll;
    float *x_ctl_sol, *y_ctl_sol, *z_ctl_sol;
    float *x_ctl_sen, *y_ctl_sen, *z_ctl_sen;
    float * in_ptr[3][3], *out_ptr[3][2];

    inlon = (float *) calloc(ncol*scansPerScene, sizeof (float));
    inlat = (float *) calloc(ncol*scansPerScene, sizeof (float));
    insolz = (float *) calloc(ncol*scansPerScene, sizeof (float));
    insola = (float *) calloc(ncol*scansPerScene, sizeof (float));
    insenz = (float *) calloc(ncol*scansPerScene, sizeof (float));
    insena = (float *) calloc(ncol*scansPerScene, sizeof (float));
    usun = (float *) calloc(3 * scansPerScene, sizeof (float));
    pos = (float *) calloc(3 * scansPerScene, sizeof (float));
    row = (int32_t *) calloc(scansPerScene, sizeof (int));
    pxl = (int16_t *) calloc(ncol, sizeof (int16_t));
    xctl = (float *) calloc(ncol, sizeof (float));
    indx = (int32_t *) calloc(scansPerScene, sizeof (int));
    yctl = (float *) calloc(scansPerScene, sizeof (float));
    in1 = (float *) calloc(scansPerScene, sizeof (float));
    in2 = (float *) calloc(scansPerScene, sizeof (float));

    /* read the data sets */
    NcVar tempVar = dataFile.getVar("det");
    tempVar.getVar(&det);
    tempVar = dataFile.getVar("pxl"); /* control pt. pixels */
    tempVar.getVar(pxl);
    tempVar = dataFile.getVar("lon"); /* geodetic longitude */
    tempVar.getVar(inlon);   
    tempVar = dataFile.getVar("lat"); /* geodetic latitude  */
    tempVar.getVar(inlat);  
    tempVar = dataFile.getVar("sun_ref"); /* geocentric ECR solar refercence vector  */
    tempVar.getVar(usun);  
    tempVar = dataFile.getVar("orb_vec"); /* geocentric ECR S/C position */
    tempVar.getVar(pos); 

    det = det - 1;
    for (i = 0; i < ncol; i++) {
        pxl[i] = pxl[i] - 1;
    }

    /* convert sun vector to unit vector */
    for (i = 0; i < scansPerScene; i++) {
        usun_sum = 0;
        for (j = 0; j < 3; j++) {
            usun_sum = usun_sum + usun[3 * i + j] * usun[3 * i + j];
        }
        for (j = 0; j < 3; j++) {
            usun[3 * i + j] = usun[3 * i + j] / sqrt(usun_sum);
        }
    }

    /* define control point grid */
    if (want_verbose)
        printf("scansPerScene,ncol = %d %d \n", scansPerScene, ncol);
    for (i = 0; i < ncol; i++) {
        xctl[i] = (float) pxl[i];
    }
    for (i = 0; i < scansPerScene; i++) {
        row[i] = i * linesPerScan + det;
    }

    /* define output grid */
    eline = sline + nline - 1;
    nlon = npix; /* # of output pixs */
    if (nline < linesPerScan * scansPerScene) {
        nlat = nline; /* # of output lines */
    } else nlat = linesPerScan*scansPerScene;

    for (i = 0; i < nlon; i++) {
        ilon[i] = i + spix; /* pix #'s of output grid */
    }
    for (i = 0; i < nlat; i++) {
        ilat[i] = i + sline; /* line #'s of output grid */
    }

    /* create index of output lines to control point lines */
    k = 0;
    for (i = 0; i < nlat; i++) {
        for (j = 0; j < scansPerScene; j++) {
            if (ilat[i] == row[j]) {
                indx[k] = i;
                k++;
            }
        }
    }

    for (i = 0; i < nlat / linesPerScan; i++) {
        yctl[i] = (float) indx[i];
    }


    /* compute solar and sensor zenith and azimuth at input
    control points from info provided */
    if (want_verbose)
        printf("Computing view angles at input control points\n");
    for (i = 0; i < scansPerScene; i++) {
        for (j = 0; j < ncol; j++) {
            for (k = 0; k < 3; k++) {
                ipos[k] = pos[3 * i + k];
                iusun[k] = usun[3 * i + k];
            }
            CalcViewAngle_nc(inlon[i * ncol + j], inlat[i * ncol + j], ipos, iusun);
            insolz[i * ncol + j] = solz1;
            insola[i * ncol + j] = sola1;
            insenz[i * ncol + j] = senz1;
            insena[i * ncol + j] = sena1;
        }
    }


    /* Compute unit vectors from lon/lat of control points */
    x_ctl_ll = (float *) calloc(ncol*scansPerScene, sizeof (float));
    y_ctl_ll = (float *) calloc(ncol*scansPerScene, sizeof (float));
    z_ctl_ll = (float *) calloc(ncol*scansPerScene, sizeof (float));

    x_ctl_sol = (float *) calloc(ncol*scansPerScene, sizeof (float));
    y_ctl_sol = (float *) calloc(ncol*scansPerScene, sizeof (float));
    z_ctl_sol = (float *) calloc(ncol*scansPerScene, sizeof (float));

    x_ctl_sen = (float *) calloc(ncol*scansPerScene, sizeof (float));
    y_ctl_sen = (float *) calloc(ncol*scansPerScene, sizeof (float));
    z_ctl_sen = (float *) calloc(ncol*scansPerScene, sizeof (float));


    for (i = 0; i < scansPerScene; i++) {
        for (j = 0; j < ncol; j++) {
            inlon[i * ncol + j] = inlon[i * ncol + j] / RADEG;
            inlat[i * ncol + j] = inlat[i * ncol + j] / RADEG;

            x_ctl_ll[i * ncol + j] = cos(inlat[i * ncol + j]) * cos(inlon[i * ncol + j]);
            y_ctl_ll[i * ncol + j] = cos(inlat[i * ncol + j]) * sin(inlon[i * ncol + j]);
            z_ctl_ll[i * ncol + j] = sin(inlat[i * ncol + j]);


            insola[i * ncol + j] = insola[i * ncol + j] / RADEG;
            insolz[i * ncol + j] = insolz[i * ncol + j] / RADEG;

            x_ctl_sol[i * ncol + j] = cos(insolz[i * ncol + j]) * cos(insola[i * ncol + j]);
            y_ctl_sol[i * ncol + j] = cos(insolz[i * ncol + j]) * sin(insola[i * ncol + j]);
            z_ctl_sol[i * ncol + j] = sin(insolz[i * ncol + j]);


            insena[i * ncol + j] = insena[i * ncol + j] / RADEG;
            insenz[i * ncol + j] = insenz[i * ncol + j] / RADEG;

            x_ctl_sen[i * ncol + j] = cos(insenz[i * ncol + j]) * cos(insena[i * ncol + j]);
            y_ctl_sen[i * ncol + j] = cos(insenz[i * ncol + j]) * sin(insena[i * ncol + j]);
            z_ctl_sen[i * ncol + j] = sin(insenz[i * ncol + j]);
        }
    }

    in_ptr[0][0] = x_ctl_ll;
    in_ptr[0][1] = y_ctl_ll;
    in_ptr[0][2] = z_ctl_ll;
    in_ptr[1][0] = x_ctl_sol;
    in_ptr[1][1] = y_ctl_sol;
    in_ptr[1][2] = z_ctl_sol;
    in_ptr[2][0] = x_ctl_sen;
    in_ptr[2][1] = y_ctl_sen;
    in_ptr[2][2] = z_ctl_sen;

    out_ptr[0][0] = lon;
    out_ptr[0][1] = lat;
    out_ptr[1][0] = sola;
    out_ptr[1][1] = solz;
    out_ptr[2][0] = sena;
    out_ptr[2][1] = senz;


    /* we now have all the info at each control point, so we
    can interpolate to all pixels, all lines */
    /* interpolate angles across each control point line */
    if (want_verbose)
        printf("Interpolating rows for longitude/azimuth\n");
    spl_aux = (float *) calloc(ncol, sizeof (float));

    for (i = 0; i < scansPerScene; i++) {
        jout = row[i] - sline;
        if ((row[i] >= sline) && (row[i] <= eline)) {
            for (l = 0; l < 3; l++) {
                spline(xctl, in_ptr[l][0] + i*ncol, ncol, 1e30, 1e30, spl_aux);
                for (j = 0; j < nlon; j++)
                    splint(xctl, in_ptr[l][0] + i * ncol, spl_aux, ncol,
                        (float) ilon[j], out_ptr[l][0] + jout * npix + j);

                spline(xctl, in_ptr[l][1] + i*ncol, ncol, 1e30, 1e30, spl_aux);
                for (j = 0; j < nlon; j++)
                    splint(xctl, in_ptr[l][1] + i * ncol, spl_aux, ncol,
                        (float) ilon[j], out_ptr[l][1] + jout * npix + j);
            }
        }
    }
    free(spl_aux);


    /* fill missing lines by interpolating columns */
    if (want_verbose)
        printf("Interpolating columns for longitude/azimuth\n");
    spl_aux = (float *) calloc(nlat / linesPerScan, sizeof (float));

    for (i = 0; i < nlon; i++) {
        for (l = 0; l < 3; l++) {
            for (k = 0; k < nlat / linesPerScan; k++) {
                in1[k] = *(out_ptr[l][0] + indx[k] * npix + i);
                in2[k] = *(out_ptr[l][1] + indx[k] * npix + i);
            }
            spline(yctl, in1, nlat / linesPerScan, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in1, spl_aux, nlat / linesPerScan, (float) j, (float *) &out1[j]);

            spline(yctl, in2, nlat / linesPerScan, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in2, spl_aux, nlat / linesPerScan, (float) j, (float *) &out2[j]);

            for (j = 0; j < nlat; j++) {
                *(out_ptr[l][0] + j * npix + i) = atan2(out2[j], out1[j]) * RADEG;
                if (l >= 1 && *(out_ptr[l][0] + j * npix + i) < 0) {
                    *(out_ptr[l][0] + j * npix + i) += 360;
                }
            }
        }
    }
    free(spl_aux);


    if (want_verbose)
        printf("Interpolating rows for latitude/zenith\n");
    spl_aux = (float *) calloc(ncol, sizeof (float));

    for (i = 0; i < scansPerScene; i++) {
        jout = row[i] - sline;
        if ((row[i] >= sline) && (row[i] <= eline)) {
            for (l = 0; l < 3; l++) {
                spline(xctl, in_ptr[l][2] + i*ncol, ncol, 1e30, 1e30, spl_aux);
                for (j = 0; j < nlon; j++)
                    splint(xctl, in_ptr[l][2] + i * ncol, spl_aux, ncol,
                        (float) ilon[j], out_ptr[l][1] + jout * npix + j);
            }
        }
    }
    free(spl_aux);


    /* fill missing lines by interpolating columns */
    if (want_verbose)
        printf("Interpolating columns for latitude/zenith\n");
    spl_aux = (float *) calloc(nlat / linesPerScan, sizeof (float));

    for (i = 0; i < nlon; i++) {
        for (l = 0; l < 3; l++) {
            for (k = 0; k < nlat / linesPerScan; k++) {
                in1[k] = *(out_ptr[l][1] + indx[k] * npix + i);
            }
            spline(yctl, in1, nlat / linesPerScan, 1e30, 1e30, spl_aux);
            for (j = 0; j < nlat; j++)
                splint(yctl, in1, spl_aux, nlat / linesPerScan, (float) j, (float *) &out1[j]);

            for (j = 0; j < nlat; j++) {
                *(out_ptr[l][1] + j * npix + i) = asin(out1[j]) * RADEG;
            }
        }
    }
    free(spl_aux);


    free(x_ctl_ll);
    free(y_ctl_ll);
    free(z_ctl_ll);
    free(x_ctl_sol);
    free(y_ctl_sol);
    free(z_ctl_sol);
    free(x_ctl_sen);
    free(y_ctl_sen);
    free(z_ctl_sen);

    free(inlon);
    free(inlat);
    free(insolz);
    free(insola);
    free(insenz);
    free(insena);
    free(usun);
    free(pos);
    free(row);
    free(pxl);
    free(xctl);
    free(indx);
    free(yctl);
    free(in1);
    free(in2);
}

/* ------------------------------------------------------ */
/* openl1_octs() - opens a OCTS HDF level-1 file */
/*                    for reading, if not already opened. */
/*                    Reads the global attributes.        */
/*                    Calculates the subscene info, if    */
/*                    necessary, ane gets the navigation  */
/*                    info.                               */
/*                                                        */
/* CALLS: getHDFattr()                                    */
/*        navigation()                                    */
/*        rdSDS()                                         */
/*                                                        */

/* ------------------------------------------------------ */
extern "C" int32_t openl1_octs_netcdf(filehandle *l1file) {
    int32_t i, j;
    int32_t pixPerScan;
    int32_t linesPerScene;
    int32_t msec_temp1[MAXOCARR], msec_temp2[MAXOCARR];
    char buf[32];

    // Reading some Global Attributes and Dimentions
    try {
        NcFile dataFile(l1file->name, NcFile::read);
        char tempTimeCoverage[27];

        // Start Time
        dataFile.getAtt("time_coverage_start").getValues((void*) &tempTimeCoverage);

        // temp to hold the value and then reassign it to the 16 bit int
        int32_t tempDay, tempYear;
        isodate2ydmsec(tempTimeCoverage, &tempYear, &tempDay, &msec_start);
        day = tempDay;
        year = tempYear;

        /* get attributes pertaining to the full scene */
        pixPerScan = dataFile.getDim("nsamp").getSize();
        scansPerScene = dataFile.getDim("rec").getSize();
        ncol = dataFile.getDim("pxls").getSize();
        linesPerScan = 2;
        dataFile.getAtt("orbit_node_longitude").getValues((void*) &l1file->orbit_node_lon);
        dataFile.getAtt("orbit_number").getValues((void*) &l1file->orbit_number);
        dataFile.getAtt("node_crossing_time").getValues((char*) &buf);
        reform_octs_time(buf);
        l1file->node_crossing_time = zulu2unix(buf);

        linesPerScene = scansPerScene * linesPerScan;

        /* Adjust spix if this is an extract file */
        if (!dataFile.getAtt("extract_pixel_start").isNull()) {
            dataFile.getAtt("extract_pixel_start").getValues((void*) &spix);
            dataFile.getAtt("extract_line_start").getValues((void*) &extr_line_offset);
            printf("File is Level-1A extract starting on line %d.\n", extr_line_offset + 1);
        } else {
            extr_line_offset = 0;
            spix = 0;
        }
        sline = 0;
        nline = linesPerScene;
        npix = pixPerScan;


        /* Check that number of scene lines not greater than allocation */
        if (nline > MAXOCLIN) {
            printf("-E- %s: Number of scene lines: %d, greater than array allocation: %d.\n",
                    "openl1_octs", nline, MAXOCLIN);
            dataFile.close();
            exit(FATAL_ERROR);
        }


        /* Check that number of pixels per scan not greater than allocation */
        if (npix > MAXOCPIX) {
            printf("-E- %s: Number of scan pixels: %d, greater than array allocation: %d.\n",
                    "openl1_octs", npix, MAXOCPIX);
            dataFile.close();
            exit(FATAL_ERROR);
        }


        lon = (float *) calloc(npix*nline, sizeof (float));
        lat = (float *) calloc(npix*nline, sizeof (float));
        senz = (float *) calloc(npix*nline, sizeof (float));
        sena = (float *) calloc(npix*nline, sizeof (float));
        solz = (float *) calloc(npix*nline, sizeof (float));
        sola = (float *) calloc(npix*nline, sizeof (float));

        l1file->npix = npix;
        l1file->nscan = nline;
        strcpy(l1file->spatialResolution, "3.5 km");

        /* compute and store the view angles for the scene
        or subscene */
        navigation_netcdf(dataFile);

        /* get the time for each line of full scene */
        // Get a variable
        NcVar tempVar = dataFile.getVar("msec");
        tempVar.getVar(msec_temp1);

        for (i = 0; i < scansPerScene; i++) {
            for (j = 0; j < linesPerScan; j++) {
                msec_temp2[i * linesPerScan + j] = msec_temp1[i];
            }
        }
        msec_start = msec_temp2[0];
        j = 0;
        for (i = sline; i < (sline + nline); i++) {
            msec[j] = msec_temp2[i];
            j++;
        }
        // Get the tilt variable
        tempVar = dataFile.getVar("tilt");
        tempVar.getVar(tilt);

        maxday = 365 + abs(LeapCheck(year));

        /* Get tilt flag, instrument temperature, and gain index */
        tilt_flag = (uint16_t *) calloc(scansPerScene, sizeof (uint16_t));
        tempVar = dataFile.getVar("tilt_flag");
        tempVar.getVar(tilt_flag);
        inst_temp = (float *) calloc(4 * scansPerScene, sizeof (float));
        tempVar = dataFile.getVar("inst_temp");
        tempVar.getVar(inst_temp);
        gain = (int16_t *) calloc(8 * scansPerScene, sizeof (int16_t));
        tempVar = dataFile.getVar("gain");
        tempVar.getVar(gain);
        miss_qual = (int16_t *) calloc(scansPerScene, sizeof (int16_t));
        tempVar = dataFile.getVar("miss_qual");
        tempVar.getVar(miss_qual);

        /* Read GAC subsampling table */
        tempVar = dataFile.getVar("num_tables");
        tempVar.getVar(&num_tab);
        start_scan = (int16_t *) calloc(num_tab, sizeof (int16_t));
        tempVar = dataFile.getVar("start_line");
        tempVar.getVar(start_scan);        

        for (i = 0; i < num_tab; i++) {
            if ((extr_line_offset + sline) / linesPerScan < start_scan[i]) {
                last_table = i - 1;
                break;
            }
        }


        tempVar = dataFile.getVar("samp_table");

        vector<size_t> samp_start(tempVar.getDimCount());
        vector<size_t> samp_edges(tempVar.getDimCount());

        samp_start[0] = last_table;
        samp_edges[0] = 1;
        samp_start[1] = 0;
        samp_edges[1] = 3;
        samp_start[2] = 0;
        samp_edges[2] = 8;
        samp_start[3] = 0;
        samp_edges[3] = 2;
        samp_start[4] = 0;
        samp_edges[4] = 400;
        samp_start[5] = 0;
        samp_edges[5] = 2;

        tempVar.getVar(samp_start, samp_edges, &sample_table);


        return (0);
        }
    catch (NcException& e) {
        cout << "-E- Error in reading input NetCDF: " << e.what() << endl;
        exit(1);
    }
}

/* ------------------------------------------------------ */
/* readl1_octs() - reads a OCTS HDF level-1 record.   */
/*                                                        */

/* ------------------------------------------------------ */
extern "C" int32_t readl1_octs_netcdf(filehandle *l1file, int32_t recnum, l1str *l1rec) {
    int32_t i, status;
    uint16_t dataarr[NOCBANDS][MAXOCPIX];
    float scanarr[MAXOCPIX][NOCBANDS];
    int32_t tilt_deg;
    int32_t scan;
    int32_t ip, ib, iw;
    int32_t nwave = l1rec->l1file->nbands;
    int32_t *bindx = l1rec->l1file->bindx;
    char *cal_path = l1_input->calfile;
    static int32_t FirstCall = 1;
    int32_t current_scan = recnum / linesPerScan;

    /* load scan times */
    /* if the day changed between the start of the scene
    and the start of the subscene, adjust the day/year
    accordingly */
    static int32_t dayIncrimented = 0;
    if (msec[recnum] < msec_start) {
        if(!dayIncrimented) {
            dayIncrimented = 1;
            day = day + 1;
            if (day > maxday) {
                year = year + 1;
                day = day - maxday;
            }
        }
    }
    if (msec[recnum] >= 86400000L) {
        msec[recnum] = msec[recnum] - 86400000L;
        if(!dayIncrimented) {
            dayIncrimented = 1;
            day = day + 1;
            if (day > maxday) {
                year = year + 1;
                day = day - maxday;
            }
        }
    }
    l1rec->scantime = yds2unix(year, day, (double) (msec[recnum] / 1.e3));

    /* load standard navigation */
    memcpy(l1rec->lon, &lon [recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->lat, &lat [recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->solz, &solz[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->sola, &sola[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->senz, &senz[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->sena, &sena[recnum * npix], sizeof (float)*npix);

    /* Read the l1a data */
    try {
        NcFile dataFile(l1file->name, NcFile::read);

        NcVar l1aDataVar = dataFile.getVar("l1a_data");

        vector<size_t> start(l1aDataVar.getDimCount());
        vector<size_t> edges(l1aDataVar.getDimCount());
        start[1] = recnum;
        start[2] = 0;
        edges[0] = 1;
        edges[1] = 1;
        edges[2] = npix;

        l1aDataVar.getVar(start, edges, &dataarr);


        for (i = 0; i < NOCBANDS; i++) {
            start[0] = i;
            l1aDataVar.getVar(start, edges, &dataarr[i][0]);
        }

        scan = (recnum + extr_line_offset) / linesPerScan;

        l1rec->tilt = tilt[current_scan];

        tilt_deg = 0;
        if (tilt_flag[current_scan] == 1) tilt_deg = +20;
        if (tilt_flag[current_scan] == 2) tilt_deg = 0;
        if (tilt_flag[current_scan] == 4) tilt_deg = -20;

        /* Read the next subsample table if necessary */
        if (FirstCall == 1 && recnum > 0) {
            for (i = 0; i < num_tab; i++) {
                if (scan / linesPerScan < start_scan[i]) {
                    last_table = i;
                    FirstCall = 0;
                    break;
                }
            }
        }
        NcVar samp_tableVar = dataFile.getVar("samp_table");

        vector<size_t> samp_start(samp_tableVar.getDimCount());
        vector<size_t> samp_edges(samp_tableVar.getDimCount());

        samp_start[0] = last_table;
        samp_edges[0] = 1;
        samp_start[1] = 0;
        samp_edges[1] = 3;
        samp_start[2] = 0;
        samp_edges[2] = 8;
        samp_start[3] = 0;
        samp_edges[3] = 2;
        samp_start[4] = 0;
        samp_edges[4] = 400;
        samp_start[5] = 0;
        samp_edges[5] = 2;

        if (scan >= start_scan[last_table + 1]) {
            last_table++;
            // This is to keep this version bug compatible with the HDF4 reader
            if (last_table >= num_tab)
                last_table = num_tab - 1;
            samp_start[0] = last_table;

            samp_tableVar.getVar(samp_start, samp_edges, &sample_table);
        }

    }
    catch (NcException& e) {
        cout << "-E- Error in reading in readl1a_octs_netcdf: " << e.what() << endl;
        exit(1);
    }

    /* Calibrate the l1a data */
    status = get_octs_cal(cal_path, year, day, msec, recnum, npix, spix, tilt_deg,
            gain, inst_temp, sample_table, scansPerScene,
            dataarr, scanarr);

    if (status < 0) {
        fprintf(stderr,
                "-E- %s line %d: Error applying calibration table \"%s\".\n",
                __FILE__, __LINE__, cal_path);
        exit(status);
    }

    /* copy the scan array over all bands, and npix pixels */
    for (ip = 0; ip < npix; ip++) {

        // if solz NaN set navfail
        if (isnan(l1rec->solz[ip]))
            l1rec->navfail[ip] = 1;

        for (iw = 0; iw < nwave; iw++) {
            ib = bindx[iw];
            l1rec->Lt [ip * nwave + ib] = scanarr[ip][iw];

            /* Also set flags */
            if (gain[iw * scansPerScene + current_scan] > 0)
                l1rec->navwarn[ip] = 1;
            if (miss_qual[current_scan] > 0)
                l1rec->navwarn[ip] = 1;
            if (dataarr[iw][ip] > 1022)
                l1rec->hilt[ip] = 1;
            if (scanarr[ip][iw] <= 0.0)
                l1rec->navfail[ip] = 1;
        }
    }

    l1rec->npix = l1file->npix;

    return (0);
}

extern "C" int32_t readl1_lonlat_octs_netcdf(filehandle *l1file, int32_t recnum, l1str *l1rec) {
    int32_t ip;

    /* load scan times */
    /* if the day changed between the start of the scene
    and the start of the subscene, adjust the day/year
    accordingly */
    static int32_t dayIncrimented = 0;
    if (msec[recnum] < msec_start) {
        if(!dayIncrimented) {
            dayIncrimented = 1;
            day = day + 1;
            if (day > maxday) {
                year = year + 1;
                day = day - maxday;
            }
        }
    }
    if (msec[recnum] >= 86400000L) {
        msec[recnum] = msec[recnum] - 86400000L;
        if(!dayIncrimented) {
            dayIncrimented = 1;
            day = day + 1;
            if (day > maxday) {
                year = year + 1;
                day = day - maxday;
            }
        }
    }
    l1rec->scantime = yds2unix(year, day, (double) (msec[recnum] / 1.e3));

    /* load standard navigation */
    memcpy(l1rec->lon, &lon [recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->lat, &lat [recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->solz, &solz[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->sola, &sola[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->senz, &senz[recnum * npix], sizeof (float)*npix);
    memcpy(l1rec->sena, &sena[recnum * npix], sizeof (float)*npix);

    /* copy the scan array over all bands, and npix pixels */
    for (ip = 0; ip < npix; ip++) {

        // if solz NaN set navfail
        if (isnan(l1rec->solz[ip]))
            l1rec->navfail[ip] = 1;

    }

    l1rec->npix = l1file->npix;

    return (0);
}

/* ------------------------------------------------------ */
/* closel1_octs() - closes the level 1 OCTS HDF file  */
/*                                                        */

/* ------------------------------------------------------ */
extern "C" int32_t closel1_octs_netcdf(filehandle *l1file) {
    free(lon);
    free(lat);
    free(senz);
    free(sena);
    free(solz);
    free(sola);

    if (l1file->format == FT_OCTSL1A) {
        free(start_scan);
        free(tilt_flag);
        free(inst_temp);
    }

    /* End access to the HDF file */
    SDend(l1file->sd_id);

    return (0);
}
