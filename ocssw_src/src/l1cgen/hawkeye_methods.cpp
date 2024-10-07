#include <fstream>
#include <iostream>
#include <sstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <string.h>
#include <string>
#include "hawkeye_methods.h"
#include <netcdf>
#include "netcdf.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <timeutils.h>
#include <allocate2d.h>
#include "l1c_latlongrid.h"
#define RADEG 57.29577951
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std;

// function to calculate cross product of two vectors, 3 components each vector
void cross_product_double2(double vector_a[], double vector_b[], double temp[]) {
    temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
    temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
    temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
}

double cross_product_norm_double2(double vector_a[], double vector_b[]) {
    double temp[3], nvec;
    temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
    temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
    temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

    nvec = sqrt(temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2]);
    return nvec;
}

int orb_to_latlon(size_t ix_swt_ini,size_t ix_swt_end,size_t num_gridlines, int nbinx, double *orb_time_tot,
                  orb_array2 *p, orb_array2 *v, double mgv1, double *tmgv1, double *tmgvf, float **lat_gd,
                  float **lon_gd, float **alt,int FirsTerrain) {
    double rl2, pos_norm, clatg2, fe = 1 / 298.257;
    double v1[3], v2[3], vecout[3], orbnorm[3], nvec, vi, toff;
    float Re = 6378.137;  // Re earth radius in km at equator
    double oangle, G[3], glat, glon, gnorm, rem = 6371, omf2, omf2p, pxy, temp;

    size_t number_orecords_corr=ix_swt_end-ix_swt_ini+1;
 


 
    orb_array2 *posgrid = new orb_array2[num_gridlines]();  // these are doubles
    orb_array2 *velgrid = new orb_array2[num_gridlines]();
    orb_array2 *posgrid2 = new orb_array2[num_gridlines]();  // these are doubles
    orb_array2 *velgrid2 = new orb_array2[num_gridlines]();

    orb_array2 *pos2 = new orb_array2[number_orecords_corr]();  // these are doubles
    orb_array2 *vel2 = new orb_array2[number_orecords_corr]();
    double* orb_time_tot2 = (double*)calloc(number_orecords_corr, sizeof(double));
    int cc=0;
    int flag_twodays=0;

    for (size_t k=ix_swt_ini;k<=ix_swt_end;k++)
    {
        if(orb_time_tot[k+1]<orb_time_tot[k])
        {
        flag_twodays=1;
        }           
        if(flag_twodays) orb_time_tot[k+1]+=24*3600; 
           
    pos2[cc][0]=p[k][0];
    pos2[cc][1]=p[k][1];
    pos2[cc][2]=p[k][2];
    vel2[cc][0]=v[k][0];
    vel2[cc][1]=v[k][1];
    vel2[cc][2]=v[k][2];
    orb_time_tot2[cc]=orb_time_tot[k];
  
    cc++;
    }
  
    orb_interp2(number_orecords_corr, num_gridlines, orb_time_tot2, pos2, vel2, tmgv1, posgrid, velgrid);

    for (size_t i = 0; i < num_gridlines - 1; i++) {
        pos_norm = sqrt(posgrid[i][0] * posgrid[i][0] + posgrid[i][1] * posgrid[i][1] +
                        posgrid[i][2] * posgrid[i][2]);
        clatg2 = sqrt(posgrid[i][0] * posgrid[i][0] + posgrid[i][1] * posgrid[i][1]) / pos_norm;
        rl2 = Re * (1 - fe) / (sqrt(1 - (2 - fe) * fe * clatg2 * clatg2));

        v2[0] = velgrid[i][0] * rl2 / pos_norm;
        v2[1] = velgrid[i][1] * rl2 / pos_norm;
        v2[2] = velgrid[i][2] * rl2 / pos_norm;

        vi = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]) * 1000;
        toff = vi / mgv1;
        tmgvf[i] = tmgv1[i] + toff;       
    }
    orb_interp2(number_orecords_corr, num_gridlines, orb_time_tot2, pos2, vel2, tmgvf, posgrid2, velgrid2);
    
    // angle subsat track----
    omf2 = (1 - fe) * (1 - fe);

    for (size_t i = 0; i < num_gridlines - 1; i++) {
        pos_norm = sqrt(posgrid2[i][0] * posgrid2[i][0] + posgrid2[i][1] * posgrid2[i][1] +
                        posgrid2[i][2] * posgrid2[i][2]);
        clatg2 = sqrt(posgrid2[i][0] * posgrid2[i][0] + posgrid2[i][1] * posgrid2[i][1]) / pos_norm;
        rl2 = Re * (1 - fe) / (sqrt(1 - (2 - fe) * fe * clatg2 * clatg2));
        // ground pos
        v1[0] = (posgrid2[i][0]) * rl2 / pos_norm;
        v1[1] = (posgrid2[i][1]) * rl2 / pos_norm;
        v1[2] = (posgrid2[i][2]) * rl2 / pos_norm;
        // ground vel
        v2[0] = (velgrid2[i][0]) * rl2 / pos_norm;
        v2[1] = (velgrid2[i][1]) * rl2 / pos_norm;
        v2[2] = (velgrid2[i][2]) * rl2 / pos_norm;

        cross_product_double2(v1, v2, vecout);
        nvec = cross_product_norm_double2(v1, v2);  // length of orb norm vect

        orbnorm[0] = vecout[0] / nvec;
        orbnorm[1] = vecout[1] / nvec;
        orbnorm[2] = vecout[2] / nvec;

        for (int j = 0; j < nbinx; j++) {
            pos_norm = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);  // first gridline
            oangle = asin((j - (nbinx - 1) / 2) * 5.2 / pos_norm);
            // Geocentric vector
            G[0] = v1[0] * cos(oangle) - orbnorm[0] * pos_norm * sin(oangle);
            G[1] = v1[1] * cos(oangle) - orbnorm[1] * pos_norm * sin(oangle);
            G[2] = v1[2] * cos(oangle) - orbnorm[2] * pos_norm * sin(oangle);

            glon = atan2(G[1], G[0]) * 180 / M_PI;
            gnorm = sqrt(G[0] * G[0] + G[1] * G[1] + G[2] * G[2]);
            omf2p = (omf2 * rem + gnorm - rem) / gnorm;
            pxy = G[0] * G[0] + G[1] * G[1];
            temp = sqrt(G[2] * G[2] + omf2p * omf2p * pxy);
            glat = asin(G[2] / temp) * 180 / M_PI;
     
            lat_gd[i][j] = glat;
            lon_gd[i][j] = glon;    
            alt[i][j]=BAD_FLT;            
        }
    } 
  

    for (size_t i = 0; i < num_gridlines-1; i++)
    {
    if(flag_twodays)
    {
    tmgvf[i]-=24*3600;
    }
    }


    //================================================
//DEM info
 if(FirsTerrain)//terrain flag =0 so DEM  not added to L1C grid
   { 
   NcFile *nc_terrain;
   const char* dem_str="$OCDATAROOT/common/gebco_ocssw_v2020.nc";
   char demfile[FILENAME_MAX];
   parse_file_name(dem_str, demfile);

   try {       
                nc_terrain = new NcFile(demfile, NcFile::read);
            } catch (NcException& e) {
                cerr << e.what() << "l1cgen orb_to_latlon: Failure to open terrain height file: " << demfile
                     << endl;
                exit(1);
           }


   NcVar terrain,var1,var2,var3;
   vector<size_t> start,count,count2,start2,start3;
    //lon/lat dimensions
   NcDim londim = nc_terrain->getDim("lon");
   NcDim latdim = nc_terrain->getDim("lat");

   start.push_back(0);
   start.push_back(0);
   count.push_back(1);
   count.push_back(1);

   short oneheight=BAD_FLT;

   double onelat=-90.,onelon=-180.,dcoor=0.00416667;//dcoor in degrees , gridded gebco
   size_t ix_grid,ix_grid2;

   var1=nc_terrain->getVar("height");
   num_gridlines-=1;

   //COMPUTE L1C GRID MIN MAX LIMITS
   for(size_t i=0;i<num_gridlines;i++)
   {
    for(int j=0;j<nbinx;j++)
    {
    ix_grid=(-1*onelat+lat_gd[i][j])/dcoor;
    ix_grid2=(-1*onelon+lon_gd[i][j])/dcoor;
    start[0]=ix_grid;
    start[1]=ix_grid2; 
    var1.getVar(start,count,&oneheight);
    alt[i][j]=oneheight;
    if(oneheight==BAD_FLT) cout<<"WARNING: oneheight is a fillvalue"<<oneheight<<"demfile "<<"ybin# "<<i+1<<"xbin# "<<j+1<<endl;
    }
    if(i%500==0) cout<<"#gridline with DEM "<<i+1<<endl;
   }
  
    nc_terrain->close();
   }//if first terrain
 
    if (posgrid != nullptr)
        delete[] (posgrid);
    posgrid = nullptr;
    if (velgrid != nullptr)
        delete[] (velgrid);
    velgrid = nullptr;
    if (posgrid2 != nullptr)
        delete[] (posgrid2);
    posgrid2 = nullptr;
    if (velgrid2 != nullptr)
        delete[] (velgrid2);
    velgrid2 = nullptr;

    if (pos2 != nullptr)
        delete[] (pos2);
    pos2 = nullptr;
    if (vel2 != nullptr)
        delete[] (vel2);
    vel2 = nullptr;
    if ( orb_time_tot2 != nullptr)
        delete[] ( orb_time_tot2);
     orb_time_tot2 = nullptr;

    return 0;
}


int orb_interp2(size_t n_orb_rec, size_t sdim, double *torb, orb_array2 *p, orb_array2 *v, double *time,
                orb_array2 *posi, orb_array2 *veli) {
    double a0[3], a1[3], a2[3], a3[3];

    size_t ind = 0;
    for (size_t i = 0; i < sdim; i++) {
        //  Find input orbit vectors bracketing scan
        for (size_t j = ind; j < n_orb_rec; j++) {
            if (time[i] > torb[j] && time[i] <= torb[j + 1]) {
                ind = j;
                break;
            }
        }

        //  Set up cubic interpolation
        double dt = torb[ind + 1] - torb[ind];
        for (size_t j = 0; j < 3; j++) {
            a0[j] = p[ind][j];
            a1[j] = v[ind][j] * dt;
            a2[j] = 3 * p[ind + 1][j] - 3 * p[ind][j] - 2 * v[ind][j] * dt - v[ind + 1][j] * dt;
            a3[j] = 2 * p[ind][j] - 2 * p[ind + 1][j] + v[ind][j] * dt + v[ind + 1][j] * dt;
        }

        //  Interpolate orbit position and velocity components to scan line time
        double x = (time[i] - torb[ind]) / dt;
        double x2 = x * x;
        double x3 = x2 * x;
        for (size_t j = 0; j < 3; j++) {
            posi[i][j] = a0[j] + a1[j] * x + a2[j] * x2 + a3[j] * x3;
            veli[i][j] = (a1[j] + 2 * a2[j] * x + 3 * a3[j] * x2) / dt;
        }
    }  // i-loop

    return 0;
}

int orb_interp(size_t n_orb_rec, size_t sdim, double *torb, orb_array *p, orb_array *v, double *time,
               orb_array *posi, orb_array *veli) {
    double a0[3], a1[3], a2[3], a3[3];

    size_t ind = 0;
    for (size_t i = 0; i < sdim; i++) {
        //  Find input orbit vectors bracketing scan
        for (size_t j = ind; j < n_orb_rec; j++) {
            if (time[i] > torb[j] && time[i] <= torb[j + 1]) {
                ind = j;
                break;
            }
        }

        //  Set up cubic interpolation
        double dt = torb[ind + 1] - torb[ind];
        for (size_t j = 0; j < 3; j++) {
            a0[j] = p[ind][j];
            a1[j] = v[ind][j] * dt;
            a2[j] = 3 * p[ind + 1][j] - 3 * p[ind][j] - 2 * v[ind][j] * dt - v[ind + 1][j] * dt;
            a3[j] = 2 * p[ind][j] - 2 * p[ind + 1][j] + v[ind][j] * dt + v[ind + 1][j] * dt;
        }

        //  Interpolate orbit position and velocity components to scan line time
        double x = (time[i] - torb[ind]) / dt;
        double x2 = x * x;
        double x3 = x2 * x;
        for (size_t j = 0; j < 3; j++) {
            posi[i][j] = a0[j] + a1[j] * x + a2[j] * x2 + a3[j] * x3;
            veli[i][j] = (a1[j] + 2 * a2[j] * x + 3 * a3[j] * x2) / dt;
        }
    }  // i-loop

    return 0;
}

int j2000_to_ecr(int32_t iyr, int32_t idy, double sec, double ecmat[3][3]) {
    // Get J2000 to ECEF transformation matrix

    // Arguments:

    // Name               Type    I/O     Description
    // --------------------------------------------------------
    // iyr         I       I      Year, four digits
    // idy         I       I      Day of year
    // sec        R*8      I      Seconds of day
    // ecmat(3,3)  R       O      J2000 to ECEF matrix

    // Get transformation from J2000 to mean-of-date inertial
    double j2mod[3][3];
    j2000_to_mod(iyr, idy, sec, j2mod);

    // Get nutation and UT1-UTC (once per run)
    double xnut[3][3], ut1utc;
    get_nut(iyr, idy, xnut);
    get_ut1(iyr, idy, ut1utc);

    // Compute Greenwich hour angle for time of day
    double day = idy + (sec + ut1utc) / 86400;
    double gha, gham[3][3];
    gha2000(iyr, day, gha);

    gham[0][0] = cos(gha / RADEG);
    gham[1][1] = cos(gha / RADEG);
    gham[2][2] = 1;
    gham[0][1] = sin(gha / RADEG);
    gham[1][0] = -sin(gha / RADEG);

    gham[0][2] = 0;
    gham[2][0] = 0;
    gham[1][2] = 0;
    gham[2][1] = 0;

    // Combine all transformations
    gsl_matrix_view A = gsl_matrix_view_array(&gham[0][0], 3, 3);  // gham
    gsl_matrix_view B = gsl_matrix_view_array(&xnut[0][0], 3, 3);  // xnut
    gsl_matrix *C = gsl_matrix_alloc(3, 3);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &A.matrix, &B.matrix, 0.0, C);

    gsl_matrix_view D = gsl_matrix_view_array(&j2mod[0][0], 3, 3);  // j2mod
    gsl_matrix *E = gsl_matrix_alloc(3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C, &D.matrix, 0.0, E);
    double *ptr_E = gsl_matrix_ptr(E, 0, 0);

    memcpy(ecmat, ptr_E, 9 * sizeof(double));

    gsl_matrix_free(C);
    gsl_matrix_free(E);

    return 0;
}

int j2000_to_mod(int32_t iyr, int32_t idy, double sec, double j2mod[3][3]) {
    // Get J2000 to MOD (precession) transformation

    // Arguments:

    // Name               Type    I/O     Description
    // --------------------------------------------------------
    // iyr         I       I      Year, four digits
    // idy         I       I      Day of year
    // sec        R*8      I      Seconds of day
    // j2mod(3,3)  R       O      J2000 to MOD matrix

    int16_t iyear = iyr;
    int16_t iday = idy;

    double t = jday(iyear, 1, iday) - (double)2451545.5 + sec / 86400;
    t /= 36525;

    double zeta0 = t * (2306.2181 + 0.302 * t + 0.018 * t * t) / RADEG / 3600;
    double thetap = t * (2004.3109 - 0.4266 * t - 0.04160 * t * t) / RADEG / 3600;
    double xip = t * (2306.2181 + 1.095 * t + 0.018 * t * t) / RADEG / 3600;

    j2mod[0][0] = -sin(zeta0) * sin(xip) + cos(zeta0) * cos(xip) * cos(thetap);
    j2mod[0][1] = -cos(zeta0) * sin(xip) - sin(zeta0) * cos(xip) * cos(thetap);
    j2mod[0][2] = -cos(xip) * sin(thetap);
    j2mod[1][0] = sin(zeta0) * cos(xip) + cos(zeta0) * sin(xip) * cos(thetap);
    j2mod[1][1] = cos(zeta0) * cos(xip) - sin(zeta0) * sin(xip) * cos(thetap);
    j2mod[1][2] = -sin(xip) * sin(thetap);
    j2mod[2][0] = cos(zeta0) * sin(thetap);
    j2mod[2][1] = -sin(zeta0) * sin(thetap);
    j2mod[2][2] = cos(thetap);

    return 0;
}

int get_nut(int32_t iyr, int32_t idy, double xnut[3][3]) {
    int16_t iyear = iyr;
    int16_t iday = idy;

    double t = jday(iyear, 1, iday) - (double)2451545.5;

    double xls, gs, xlm, omega;
    double dpsi, eps, epsm;
    ephparms(t, xls, gs, xlm, omega);
    nutate(t, xls, gs, xlm, omega, dpsi, eps, epsm);

    xnut[0][0] = cos(dpsi / RADEG);
    xnut[1][0] = -sin(dpsi / RADEG) * cos(epsm / RADEG);
    xnut[2][0] = -sin(dpsi / RADEG) * sin(epsm / RADEG);
    xnut[0][1] = sin(dpsi / RADEG) * cos(eps / RADEG);
    xnut[1][1] =
        cos(dpsi / RADEG) * cos(eps / RADEG) * cos(epsm / RADEG) + sin(eps / RADEG) * sin(epsm / RADEG);
    xnut[2][1] =
        cos(dpsi / RADEG) * cos(eps / RADEG) * sin(epsm / RADEG) - sin(eps / RADEG) * cos(epsm / RADEG);
    xnut[0][2] = sin(dpsi / RADEG) * sin(eps / RADEG);
    xnut[1][2] =
        cos(dpsi / RADEG) * sin(eps / RADEG) * cos(epsm / RADEG) - cos(eps / RADEG) * sin(epsm / RADEG);
    xnut[2][2] =
        cos(dpsi / RADEG) * sin(eps / RADEG) * sin(epsm / RADEG) + cos(eps / RADEG) * cos(epsm / RADEG);

    return 0;
}

int ephparms(double t, double &xls, double &gs, double &xlm, double &omega) {
    //  This subroutine computes ephemeris parameters used by other Mission
    //  Operations routines:  the solar mean longitude and mean anomaly, and
    //  the lunar mean longitude and mean ascending node.  It uses the model
    //  referenced in The Astronomical Almanac for 1984, Section S
    //  (Supplement) and documented for the SeaWiFS Project in "Constants
    //  and Parameters for SeaWiFS Mission Operations", in TBD.  These
    //  parameters are used to compute the solar longitude and the nutation
    //  in longitude and obliquity.

    //  Sun Mean Longitude
    xls = (double)280.46592 + ((double)0.9856473516) * t;
    xls = fmod(xls, (double)360);

    //  Sun Mean Anomaly
    gs = (double)357.52772 + ((double)0.9856002831) * t;
    gs = fmod(gs, (double)360);

    //  Moon Mean Longitude
    xlm = (double)218.31643 + ((double)13.17639648) * t;
    xlm = fmod(xlm, (double)360);

    //  Ascending Node of Moon's Mean Orbit
    omega = (double)125.04452 - ((double)0.0529537648) * t;
    omega = fmod(omega, (double)360);

    return 0;
}

int nutate(double t, double xls, double gs, double xlm, double omega, double &dpsi, double &eps,
           double &epsm) {
    //  This subroutine computes the nutation in longitude and the obliquity
    //  of the ecliptic corrected for nutation.  It uses the model referenced
    //  in The Astronomical Almanac for 1984, Section S (Supplement) and
    //  documented for the SeaWiFS Project in "Constants and Parameters for
    //  SeaWiFS Mission Operations", in TBD.  These parameters are used to
    //  compute the apparent time correction to the Greenwich Hour Angle and
    //  for the calculation of the geocentric Sun vector.  The input
    //  ephemeris parameters are computed using subroutine ephparms.  Terms
    //  are included to 0.1 arcsecond.

    //  Nutation in Longitude
    dpsi = -17.1996 * sin(omega / RADEG) + 0.2062 * sin(2. * omega / RADEG) - 1.3187 * sin(2. * xls / RADEG) +
           0.1426 * sin(gs / RADEG) - 0.2274 * sin(2. * xlm / RADEG);

    //  Mean Obliquity of the Ecliptic
    epsm = (double)23.439291 - ((double)3.560e-7) * t;

    //  Nutation in Obliquity
    double deps = 9.2025 * cos(omega / RADEG) + 0.5736 * cos(2. * xls / RADEG);

    //  True Obliquity of the Ecliptic
    eps = epsm + deps / 3600;

    dpsi = dpsi / 3600;

    return 0;
}

int get_ut1(int32_t iyr, int32_t idy, double &ut1utc) {
    int16_t iyear = iyr;
    int16_t iday = idy;

    static int32_t ijd[25000];
    static double ut1[25000];
    string utcpole = "$OCVARROOT/var/modis/utcpole.dat";
    static bool first = true;
    int k = 0;
    if (first) {
        string line;
        expandEnvVar(&utcpole);
        istringstream istr;

        ifstream utcpfile(utcpole.c_str());
        if (utcpfile.is_open()) {
            getline(utcpfile, line);
            getline(utcpfile, line);
            while (getline(utcpfile, line)) {
                istr.clear();
                istr.str(line.substr(0, 5));
                istr >> ijd[k];
                istr.clear();
                istr.str(line.substr(44, 9));
                istr >> ut1[k];
                k++;
            }
            ijd[k] = 0;
            utcpfile.close();
            first = false;
        } else {
            cout << utcpole.c_str() << " not found" << endl;
            exit(1);
        }
    }  // if (first)

    k = 0;
    int mjd = jday(iyear, 1, iday) - 2400000;
    while (ijd[k] > 0) {
        if (mjd == ijd[k]) {
            ut1utc = ut1[k];
            break;
        }
        k++;
    }

    return 0;
}

int gha2000(int32_t iyr, double day, double &gha) {
    //  This subroutine computes the Greenwich hour angle in degrees for the
    //  input time.  It uses the model referenced in The Astronomical Almanac
    //  for 1984, Section S (Supplement) and documented for the SeaWiFS
    //  Project in "Constants and Parameters for SeaWiFS Mission Operations",
    //  in TBD.  It includes the correction to mean sidereal time for nutation
    //  as well as precession.
    //

    //  Compute days since J2000
    int16_t iday = day;
    double fday = day - iday;
    int jd = jday(iyr, 1, iday);
    double t = jd - (double)2451545.5 + fday;

    //  Compute Greenwich Mean Sidereal Time      (degrees)
    double gmst = (double)100.4606184 + (double)0.9856473663 * t + (double)2.908e-13 * t * t;

    //  Check if need to compute nutation correction for this day
    double xls, gs, xlm, omega;
    double dpsi, eps, epsm;
    ephparms(t, xls, gs, xlm, omega);
    nutate(t, xls, gs, xlm, omega, dpsi, eps, epsm);

    //  Include apparent time correction and time-of-day
    gha = gmst + dpsi * cos(eps / RADEG) + fday * 360;
    gha = fmod(gha, (double)360);

    return 0;
}

int expandEnvVar(std::string *sValue) {
    if ((*sValue).find_first_of("$") == string::npos)
        return 0;
    string::size_type posEndIdx = (*sValue).find_first_of("/");
    if (posEndIdx == string::npos)
        return 0;
    char *envVar_str = getenv((*sValue).substr(1, posEndIdx - 1).c_str());
    if (envVar_str == 0x0) {
        printf("Environment variable: %s not defined.\n", envVar_str);
        exit(1);
    }
    *sValue = envVar_str + (*sValue).substr(posEndIdx);

    return 0;
}
