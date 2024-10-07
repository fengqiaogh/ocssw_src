/*******************************************************************

   met_cvt.h

   purpose: include file for the use of the meteorological conversion
     routines

   Parameters: 
      Type              Name            I/O     Description
      ----              ----            ---     -----------

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 12-Mar-2009     Original development
      W. Robinson, SAIC 5 Aug 2013      rename for met_cvt

 *******************************************************************/

#include <math.h>
/*
 *  set up some definitions met data units
 */
#define MET_UNITS__P_PA  0  /* Pressure in Pascals */
#define MET_UNITS__P_HPA  1 /* Pressure in Hectopascals */
#define MET_UNITS__T_K  10   /* Temperature in deg K */
#define MET_UNITS__T_C  11   /* Temperature in deg C */
#define MET_UNITS__Q_KG_KG 20 /* Specific humidity in Kg / Kg */
#define MET_UNITS__Q_G_KG 21  /* Specific humidity in g / Kg */

/*
 *  common physical constants
 */
#define M_DRY 28.9644   /* molar mass of dry air */
#define M_WET 18.01534   /* molar mass of water */
#define C_IN_K 273.15    /* the value of 0 C in K */
#define R_W 461.5       /* gas constanr for water vapor in K^-1 kg^-1 */
#define L_ENTHALPY 2.38e6  /* enthalpy of vaporization (varies from
             2.501 at 273.15K to 2.257 at 373.15)  */
#define MAGNUS_A1 17.625  /* the A1 value for the Magnus equation for es */
#define MAGNUS_B1 243.04  /* the B1 value for the Magnus equation for es */
/*  coefficients for computing saturated vapor pressure */
static const double es_coef_wet[] = {6.107799961, 4.436518521e-1, 1.428945805e-2,
    2.650648471e-4, 3.031240396e-6, 2.034080948e-8, 6.136820929e-11};
static const double es_coef_ice[] = {6.109177956, 5.034698970e-1, 1.886013408e-2, 
    4.176223716e-4, 5.824720280e-6, 4.838803174e-8, 1.838826904e-10 };
/*
 *  standard pressure level set
 */
static const float std_p_lvls_42[] = {1000., 975., 950., 925., 900., 875., 850.,
    825., 800., 775., 750., 725., 700., 650., 600., 550., 500., 450., 400.,
    350., 300., 250., 200., 150., 100., 70., 50., 40., 30., 20., 10., 7., 5.,
    4., 3., 2., 1., 0.7, 0.5, 0.4, 0.3, 0.1};
/*
 *  prototypes
 */
int met_cvt_q_to_rh(int, float *, int, float *, int, float *, int, float *);
int met_cvt_rh_to_q( int, float *, int, float *, int, float *,  float *, int);
double met_cvt_p_cvt(double, int, int);
double met_cvt_t_cvt(double, int, int);
double met_cvt_q_cvt(double, int, int);
int met_cvt_ttd_to_rh(int, float *, int, float *, int, float *);
