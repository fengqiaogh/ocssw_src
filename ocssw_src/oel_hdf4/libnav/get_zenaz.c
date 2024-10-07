#include <libnav.h>

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define RAD_TO_DEG    57.295779513082321
#define DEG_TO_RAD   .017453292519943296

void euler(float *a,float *xm[]);

/*
 * get_zenaz
 *
 *  Created on: Aug 23, 2013
 *      Author: dshea
 *
 *  Calculate the zenith and  azimuth for input lat/lon values
 *  given the satellite position.
 *
 *   pos(3)   Orbit Position Vector (km)
 *   lat      Pixel geodetic latitudes
 *   lon      Pixel geodetic longitudes
 *   zenith   Pixel sensor zenith angle
 *   azimuth  Pixel sensor azimuth angle
 *
 */
void get_zenaz(float *pos, float lon, float lat, float *senz, float *sena) {
    int i;
    double gv[3];
    float ea[3];
    float no[3];
    float up[3];
    double rl[3];

    double re = 6378.137;
    double f = 1 / 298.257;
    double omf2 = (1.0 - f) * (1.0 - f);
    double xlat = lat * DEG_TO_RAD;
    double xlon = lon * DEG_TO_RAD;

    // Compute the local vertical, East and North unit vectors
    up[0] = cos(xlat) * cos(xlon);
    up[1] = cos(xlat) * sin(xlon);
    up[2] = sin(xlat);
    double upxy = sqrt(up[0] * up[0] + up[1] * up[1]);
    ea[0] = -up[1] / upxy;
    ea[1] = up[0] / upxy;
    ea[2] = 0.0;
    crossp_(up, ea, no);

    // Compute geocentric position vector
    double xlatg = atan(tan(xlat) * omf2);
    gv[0] = cos(xlatg) * cos(xlon);
    gv[1] = cos(xlatg) * sin(xlon);
    gv[2] = sin(xlatg);
    double r = re * (1.0 - f) / sqrt(1.0 - (2.0 - f) * f * pow(cos(xlatg), 2));

    //Transform the pixel-to-spacecraft and Sun vectors into the local frame
    gsl_matrix* xlMatrix = gsl_matrix_alloc(3, 3);
    gsl_matrix* rhMatrix = gsl_matrix_alloc(3, 1);
    gsl_matrix* rlMatrix = gsl_matrix_alloc(3, 1);

    for (i = 0; i < 3; i++) {
        gsl_matrix_set(xlMatrix, 0, i, ea[i]);
        gsl_matrix_set(xlMatrix, 1, i, no[i]);
        gsl_matrix_set(xlMatrix, 2, i, up[i]);

        gsl_matrix_set(rhMatrix, i, 0, pos[i] - r * gv[i]);
        gsl_matrix_set(rlMatrix, i, 0, 0.0);
    }
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, xlMatrix, rhMatrix, 0.0,
            rlMatrix);

    for (i = 0; i < 3; i++) {
        rl[i] = gsl_matrix_get(rlMatrix, i, 0);
    }

    gsl_matrix_free(xlMatrix);
    gsl_matrix_free(rhMatrix);
    gsl_matrix_free(rlMatrix);

    // Compute the sensor zenith and azimuth
    *senz = RAD_TO_DEG * atan2(sqrt(rl[0] * rl[0] + rl[1] * rl[1]), rl[2]);

    // Check for zenith close to zero
    if (*senz > 0.05)
        *sena = RAD_TO_DEG * atan2(rl[0], rl[1]);
    else
        *sena = 0.0;

    if (*sena < 0.0)
        *sena += 360.0;

}

void transpose(float *in[], float *out[]) {
	int i,j;

	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			out[j][i] = in[i][j];
}

void nav_get_vel(float ilat, float mlat, float ilon, float mlon, float *vel) {

	/*
	 * Determine the velocity vectors
	 *
	 * ilat, ilon (in) - latitude and longitude intercept of the lsq fitted lat lon positions
	 * mlat, mlon (in) - latitude and longitude slope     "   "   "     "    "   "      "
	 *
	 * float vel[3] - (out) velocity vector
	 *
	 * R. Healy 11/2/2016
	 */
	float heading, chead, shead;
	float clat,slat,clon,slon;

	heading = atan2(mlon*cos(ilat*DEG_TO_RAD),mlat);
	chead = cos(heading);
	shead = sin(heading);
	clat  = cos(ilat*DEG_TO_RAD);
	clon  = cos(ilon*DEG_TO_RAD);
	slat  = sin(ilat*DEG_TO_RAD);
	slon  = sin(ilon*DEG_TO_RAD);
	// V is computed from a 3-2-1 Euler rotation
	// lon is phi, lat is -theta, head is -psi
	vel[0] = -clon*slat*chead - slon*shead;
	vel[1] = -slon*slat*chead + clon*shead;
	vel[2] = clat*chead;

}
void nav_get_pos(float lat, float lon, float alt, float *p) {

	/*
	 * Determine the position vectors
	 *
	 * lat, lon (in)- latitude and longitude of the pixel
	 * latr, lonr (in) - latitude and longitude of position at time > time of lat, lon
	 *
	 * float pos[3] - (out) position and velocity vectors
	 *
	 * R. Healy 11/2/2016
	 */
	double re = 6378.137;
	double f = 1 / 298.257;
	double omf2 = (1.0 - f) * (1.0 - f);
	double glat,cglat,rl;
	double olon=lon*DEG_TO_RAD,olat=lat*DEG_TO_RAD;

	glat = atan(tan(olat)*omf2);
	cglat = cos(glat);
	rl = re*(1.-f)/sqrt(1.-(2.-f)*f*cglat*cglat);
	p[0] = cos(olon)*(rl*cglat + alt*cos(olat));
	p[1] = sin(olon)*(rl*cglat + alt*cos(olat));
	p[2] = rl*sin(glat) + alt*sin(olat);

}
void nav_gd_orient(float *pos, float *vel, float *att, float *smat[], float *coef) {
	/*
	  This subroutine performs a simple calculation of the sensor
	  orientation from the orbit position vector and input values of the
	  attitude offset angles.  The calculations assume that the angles
	  represent the roll, pitch and yaw offsets between the orbit
	  reference frame (at the spacecraft position) and the sensor frame.
	  Sensor tilt angles are assumed to be included in the pitch angle.
	  The outputs are the matrix which represents the transformation from
	  the geocentric rotating to sensor frame, and the coefficients which
	  represent the Earth scan track in the sensor frame.

	  The reference ellipsoid uses an equatorial radius of 6378.137 km and
	  a flattening factor of 1/298.257 (WGS 1984).

	  Calling Arguments

	  Name     Type  I/O   Description

	  pos(3)       R*4      I ECR Orbit Position Vector (km)
	  vel(3)       R*4      I ECR Orbit Velocity Vector (km/sec)
	  att(3)       R*4      I Attitude angles roll, pitch, yaw
	  smat(3,3)    R*4      O Sensor Orientation Matrix
	  coef(10)     R*4      O Scan path coefficients

	Translated from IDL function
	R. Healy 11/2/2016


	*/
    double re = 6378.137,rem = 6371.0;
    double f = 1 / 298.257;
    double omf2 = (1.0 - f) * (1.0 - f);
	double omegae = 7.29211585494e-5;
    double rd = 1.0/omf2;

    double velm[3],pm,omf2p,pxy,temp,x[3],y[3],z[3];
    double on[3],onm;
    float *sm1[3],sm2[3][3];
    int i,j;

    //  Determine local vertical reference axes
    velm[0] = vel[0] - pos[1]*omegae;
    velm[1] = vel[1] + pos[0]*omegae;
    velm[2] = vel[2];

    for (i=0;i<3;i++)
    	sm1[i] = (float *)malloc(3*sizeof(float));

    //Compute Z axis as local nadir vector
     pm = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
     omf2p = (omf2*rem + pm - rem)/pm;
     pxy = pos[0]*pos[0]+pos[1]*pos[1];
     temp = sqrt(pos[2]*pos[2] + omf2p*omf2p*pxy);
     z[0] = -omf2p*pos[0]/temp;
     z[1] = -omf2p*pos[1]/temp;
     z[2] = -pos[2]/temp;

     //Compute Y axis along negative orbit normal
     // Do the cross product

     on[0] = velm[1]*z[2]-velm[2]*z[1];
     on[1] = velm[2]*z[0]-velm[0]*z[2];
     on[2] = velm[0]*z[1]-velm[1]*z[0];

     onm=sqrt(on[0]*on[0]+on[1]*on[1]+on[2]*on[2]);

     for (i=0;i<3;i++) y[i] = -on[i]/onm;

     // Compute X axis to complete orthonormal triad (velocity direction)

     x[0] = y[1]*z[2]-y[2]*z[1];
     x[1] = y[2]*z[0]-y[0]*z[2];
     x[2] = y[0]*z[1]-y[1]*z[0];

     euler(att,sm1);

     //  Compute attitude matrix in geocentric frame
     //   Store local vertical reference vectors in matrix
     for (i=0;i<3;i++) {
            sm2[0][i]=x[i];
            sm2[1][i]=y[i];
            sm2[2][i]=z[i];
     }

 	for (i=0;i<3;i++)
 		for (j=0;j<3;j++)
 			smat[i][j]=sm1[i][0]*sm2[0][j]+sm1[i][1]*sm2[1][j]+sm1[i][2]*sm2[2][j];

 	//  Compute coefficients of intersection ellipse in scan plane
 	   coef[0] = 1.0+(rd-1.0)*smat[0][2]*smat[0][2];
 	   coef[1] = 1.0+(rd-1.0)*smat[1][2]*smat[1][2];
 	   coef[2] = 1.0+(rd-1.0)*smat[2][2]*smat[2][2];
 	   coef[3] = (rd-1.0)*smat[0][2]*smat[1][2]*2.0;
 	   coef[4] = (rd-1.0)*smat[0][2]*smat[2][2]*2.0;
 	   coef[5] = (rd-1.0)*smat[1][2]*smat[2][2]*2.0;
 	   coef[6] = (smat[0][0]*pos[0]+smat[0][1]*pos[1]+smat[0][2]*pos[2]*rd)*2.0;
 	   coef[7] = (smat[1][0]*pos[0]+smat[1][1]*pos[1]+smat[1][2]*pos[2]*rd)*2.0;
 	   coef[8] = (smat[2][0]*pos[0]+smat[2][1]*pos[1]+smat[2][2]*pos[2]*rd)*2.0;
 	   coef[9] = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]*rd-re*re;

}

void euler(float *a,float *xm[]) {
	//  Computes coordinate transformation matrix corresponding to Euler
	//  sequence; assumes order of rotations is X, Y, Z
	//
	//  Reference:  Wertz, Appendix E
	//	Translated from IDL function
	//	R. Healy 11/2/2016
	//

	double xm1[3][3],xm2[3][3],xm3[3][3],xmm[3][3];
	double c1,c2,c3,s1,s2,s3;
	int i,j;
	c1=cos(a[0]*DEG_TO_RAD);
	s1=sin(a[0]*DEG_TO_RAD);
	c2=cos(a[1]*DEG_TO_RAD);
	s2=sin(a[1]*DEG_TO_RAD);
	c3=cos(a[2]*DEG_TO_RAD);
	s3=sin(a[2]*DEG_TO_RAD);
	//  Convert individual rotations to matrices
	xm1[0][0]=1.0;
	xm1[1][1]= c1;
	xm1[2][2]= c1;
	xm1[1][2]= s1;
	xm1[2][1]=-s1;
	xm2[1][1]= 1.0;
	xm2[0][0]= c2;
	xm2[2][2]= c2;
	xm2[2][0]= s2;
	xm2[0][2]=-s2;
	xm3[2][2]= 1.0;
	xm3[1][1]= c3;
	xm3[0][0]= c3;
	xm3[0][1]= s3;
	xm3[1][0]=-s3;

	//  Compute total rotation as xm3*xm2*xm1
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			xmm[i][j]=xm2[i][0]*xm1[0][j]+xm2[i][1]*xm1[1][j]+xm2[i][2]*xm1[2][j];

	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			xm[i][j]=xm3[i][0]*xmm[0][j]+xm3[i][1]*xmm[1][j]+xm3[i][2]*xmm[2][j];

	return;

}
void nav_get_geonav(float *sunr, float *pos, float *pview, float *coef, float *smat[], float *xlon, float *xlat, float *solz, float *sola, float *senz, float *sena) {
//	;  This subroutine performs navigation of a scanning sensor on the
//	;  surface of an ellipsoid based on an input orbit position vector and
//	;  spacecraft orientation matrix.  It uses a closed-form algorithm for
//	;  determining points on the ellipsoidal surface which involves
//	;  determining the intersection of the scan plan with the ellipsoid.
//	;  The sensor view vectors in the sensor frame are passed in as a 3xN array.
//	;
//	;  The reference ellipsoid is set according to the scan
//	;  intersection coefficients in the calling sequence; an equatorial
//	;  radius of 6378.137 km. and a flattening factor of 1/298.257 are
//	;  used by both the Geodetic Reference System (GRS) 1980 and the
//	;  World Geodetic System (WGS) 1984.
//	;
//	;  It then computes geometric parameters using the pixel locations on
//	;  the Earth, the spaecraft position vector.  The outputs are arrays of
//	;  geodetic latitude and longitude, sensor
//	;  zenith and azimuth.  The azimuth angles are measured from local
//	;  North toward East.  Flag values of 999. are returned for any pixels
//	;  whose scan angle is past the Earth's horizon.
//	;
//	;  Reference: "Exact closed-form geolocation algorithm for
//	;  Earth survey sensors", F. S. Patt and W. W. Gregg, IJRS, Vol. 15
//	;  No. 18, 1994.
//
//	Translated from IDL function uni_geonav - solar angles removed
//	R. Healy 11/2/2016
//
    int i;
    float ea[3];
    float no[3];
    double up[3];
    double rh[3];
    double rl[3];
    double sl[3];

    double f = 1 / 298.257;
    double omf2 = (1.0 - f) * (1.0 - f);
    double temp,d,o,p,q,r,uxy,x1[3],vv[3],geovec[3];
    float *smatt[3];
    for (i=0;i<3;i++)
    	smatt[i] = (float *) malloc(3*sizeof(float));

    // Move scan view vectors to local array
    vv[0] = pview[0];
    vv[1] = pview[1];
    vv[2] = pview[2];

    //Compute sensor-to-surface vectors for all scan angles
    // Compute terms for quadratic equation
    o =   coef[0]*vv[0]*vv[0] + coef[1]*vv[1]*vv[1]
		+ coef[2]*vv[2]*vv[2] + coef[3]*vv[0]*vv[1]
		+ coef[4]*vv[0]*vv[2] + coef[5]*vv[1]*vv[2];

    p = coef[6]*vv[0] + coef[7]*vv[1] + coef[8]*vv[2];
    q = coef[9];
    r = p*p-4.0*q*o;

    //  Solve for magnitude of sensor-to-pixel vector and compute components

    d = (-p-sqrt(r))/(2.0*o);
    x1[0]=d*vv[0];
    x1[1]=d*vv[1];
    x1[2]=d*vv[2];
    //  Transform vector from sensor to geocentric frame

    transpose(smat, smatt);
	for (i=0;i<3;i++)
		rh[i]=smatt[i][0]*x1[0]+smatt[i][1]*x1[1]+smatt[i][2]*x1[2];
	for (i=0;i<3;i++)
		geovec[i] = pos[i] + rh[i];
	// Compute the local vertical, East and North unit vectors
	uxy = geovec[0]*geovec[0]+geovec[1]*geovec[1];
	temp = sqrt(geovec[2]*geovec[2] + omf2*omf2*uxy);

    up[0] = omf2*geovec[0]/temp;
    up[1] = omf2*geovec[1]/temp;
    up[2] = geovec[2]/temp;
    double upxy = sqrt(up[0] * up[0] + up[1] * up[1]);
    ea[0] = -up[1] / upxy;
    ea[1] =  up[0] / upxy;
    ea[2] =  0.0;
    no[0] = -up[2]*ea[1];
    no[1] =  up[2]*ea[0];
    no[2] =  up[0]*ea[1] - up[1]*ea[0];

    //    Compute geodetic latitude and longitude
    *xlat = asin(up[2])*RAD_TO_DEG;
    *xlon = atan2(up[1],up[0])*RAD_TO_DEG;
    //range = d

    //     Transform the pixel-to-spacecraft and Sun vectors into the local
    //     frame
    rl[0] = -ea[0]*rh[0] - ea[1]*rh[1] - ea[2]*rh[2];
    rl[1] = -no[0]*rh[0] - no[1]*rh[1] - no[2]*rh[2];
    rl[2] = -up[0]*rh[0] - up[1]*rh[1] - up[2]*rh[2];

    sl[0] = ea[0]*sunr[0] + ea[1]*sunr[1] + ea[2]*sunr[2];
    sl[1] = no[0]*sunr[0] + no[1]*sunr[1] + no[2]*sunr[2];
    sl[2] = up[0]*sunr[0] + up[1]*sunr[1] + up[2]*sunr[2];

    *solz = RAD_TO_DEG * atan2(sqrt(sl[0] * sl[0] + sl[1] * sl[1]), sl[2]);
    *sola = RAD_TO_DEG * atan2(sl[0], sl[1]);

    // Compute the sensor zenith and azimuth
    *senz = RAD_TO_DEG * atan2(sqrt(rl[0] * rl[0] + rl[1] * rl[1]), rl[2]);

    //      Check for zenith close to zero
    if (*senz > 0.02)
        *sena = RAD_TO_DEG * atan2(rl[0], rl[1]);
    else
	   *sena = 0.0;


}
