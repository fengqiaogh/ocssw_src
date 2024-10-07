/*
 * get_npp.h
 *
 *  Created on: Jan 9, 2015
 *      Author: rhealy
 */

float npp_vgpm(float chl, float irr, float sst, float dayL);
float npp_eppley(float chl, float irr, float sst, float dayL);
double npp_cbpm2(double chl, double bbp, double irr, double k490, double mld, double zno3, double daylength);
float get_zno3(float lon, float lat, int day);
double Kd490_obpg(l2str *l2rec, float *kd);
void get_par_clim(float *parin, float *lat, float *lon, int nlat, int nlon, float *latp, float *lonp, int32_t npix, float *par);
