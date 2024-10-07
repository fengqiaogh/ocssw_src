#include "met_cvt.h"
#include <genutils.h>

/*

  met_cvt.c is a collection of meteorological conversion routines

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       17-Aug-2006     Original development
      W. Robinson       5 Aug 2013      change to met_cvt
 */

int met_cvt_q_to_rh(int nval, float *pres, int p_type, float *temp,
        int t_type, float *q, int q_type, float *rh)
/*******************************************************************

   met_cvt_q_to_rh

   purpose: convert q (specific humidity) into relative humidity

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               nval             I      # of values in following 
                                                arrays to convert
      float *           pres             I      pressure array
      int               p_type           I      pressure type see met_convert.h
                                                for the types
      float *           temp             I      temperature array
      int               t_type           I      temperature type
      float *           q                I      specific humidity array
      int               q_type           I      specific humidity type
      float *           rh               O      relative humidity in %

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       12-May-2009     Original development

 *******************************************************************/ {
    int ival;
    double x_h20; /* volume mixing ratio */
    double es, es_wet, es_ice; /* saturation vapor pressure  */
    double p_lcl, t_lcl, q_lcl, rh_lcl, t_es;
    /*
     *  loop through the inputs
     */
    for (ival = 0; ival < nval; ival++) {
        /*
         *  get local parameter values and convert them to 
         *  hPa (p), C (t), and kg kg^-1 (q)
         */
        p_lcl = *(pres + ival);
        t_lcl = *(temp + ival);
        q_lcl = *(q + ival);
        if ((p_lcl != BAD_FLT) && (t_lcl != BAD_FLT) &&
                (q_lcl != BAD_FLT)) {

            p_lcl = met_cvt_p_cvt(p_lcl, p_type, MET_UNITS__P_HPA);
            t_lcl = met_cvt_t_cvt(t_lcl, t_type, MET_UNITS__T_C);
            q_lcl = met_cvt_q_cvt(q_lcl, q_type, MET_UNITS__Q_KG_KG);
            /*
             *  then, create the RH
             *  Note that for T < -50 C, the es is not good.  To have something
             *  compute es at -50 if t is lower
             */
            t_es = (t_lcl < -50.) ? -50. : t_lcl;
            es_wet = es_coef_wet[0] + t_es * (es_coef_wet[1] + t_es *
                    (es_coef_wet[2] + t_es *
                    (es_coef_wet[3] + t_es * (es_coef_wet[4] + t_es *
                    (es_coef_wet[5] + t_es * es_coef_wet[6])))));
            es_ice = es_coef_ice[0] + t_es * (es_coef_ice[1] + t_es *
                    (es_coef_ice[2] + t_es *
                    (es_coef_ice[3] + t_es * (es_coef_ice[4] + t_es *
                    (es_coef_ice[5] + t_es * es_coef_ice[6])))));
            es = (es_wet < es_ice) ? es_wet : es_ice;

            x_h20 = q_lcl * M_DRY / (q_lcl * M_DRY + (1. - q_lcl) * M_WET);

            rh_lcl = 100. * x_h20 * p_lcl / es;
            /*
             *  limit the RH that comes out and transfer
             */
            rh_lcl = (rh_lcl < 0) ? 0 : rh_lcl;
            rh_lcl = (rh_lcl > 100.) ? 100. : rh_lcl;

            *(rh + ival) = rh_lcl;
        } else {
            *(rh + ival) = BAD_FLT;
        }
    }
    return 0;
}

int met_cvt_rh_to_q( int nval, float *pres, int p_type, float *temp,
        int t_type, float *rh,  float *q, int q_type)
/*******************************************************************

   met_cvt_rh_to_q

   purpose: convert rh (relative humidity) into specific humidity 

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               nval             I      # of values in following 
                                                arrays to convert
      float *           pres             I      pressure array
      int               p_type           I      pressure type see met_convert.h
                                                for the types
      float *           temp             I      temperature array
      int               t_type           I      temperature type
      float *           rh               O      relative humidity in %
      float *           q                I      specific humidity array
      int               q_type           I      specific humidity type

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       22-Jan-2019     Original development

 *******************************************************************/ {
    int ival;
    double x_h20; /* volume mixing ratio */
    double es, es_wet, es_ice; /* saturation vapor pressure  */
    double p_lcl, t_lcl, q_lcl, rh_lcl, t_es, ph2o;
    /*
     *  loop through the inputs
     */
    for (ival = 0; ival < nval; ival++) {
        /*
         *  get local parameter values and convert them to 
         *  hPa (p), C (t), and kg kg^-1 (q)
         */
        p_lcl = *(pres + ival);
        t_lcl = *(temp + ival);
        rh_lcl = *(rh  + ival);
        if ((p_lcl != BAD_FLT) && (t_lcl != BAD_FLT) &&
                (rh_lcl != BAD_FLT)) {

            p_lcl = met_cvt_p_cvt(p_lcl, p_type, MET_UNITS__P_HPA);
            t_lcl = met_cvt_t_cvt(t_lcl, t_type, MET_UNITS__T_C);
            rh_lcl = ( rh_lcl > 100. ) ? 100. : rh_lcl;
            rh_lcl = ( rh_lcl < 0. ) ? 0. : rh_lcl;
            /*
             *  then, create the Q
             *  Note that for T < -50 C, the es is not good.  To have something
             *  compute es at -50 if t is lower
             */
            t_es = (t_lcl < -50.) ? -50. : t_lcl;
            es_wet = es_coef_wet[0] + t_es * (es_coef_wet[1] + t_es *
                    (es_coef_wet[2] + t_es *
                    (es_coef_wet[3] + t_es * (es_coef_wet[4] + t_es *
                    (es_coef_wet[5] + t_es * es_coef_wet[6])))));
            es_ice = es_coef_ice[0] + t_es * (es_coef_ice[1] + t_es *
                    (es_coef_ice[2] + t_es *
                    (es_coef_ice[3] + t_es * (es_coef_ice[4] + t_es *
                    (es_coef_ice[5] + t_es * es_coef_ice[6])))));
            es = (es_wet < es_ice) ? es_wet : es_ice;

            ph2o = es * rh_lcl / 100.;

            x_h20 = ph2o / p_lcl;

            q_lcl = x_h20 * M_WET / ( x_h20 * M_WET + ( 1. - x_h20 ) * M_DRY );
            q_lcl = met_cvt_q_cvt( q_lcl, MET_UNITS__Q_KG_KG, q_type );

            *(q + ival) = q_lcl;
        } else {
            *(q + ival) = BAD_FLT;
        }
    }
    return 0;
}

int met_cvt_ttd_to_rh(int nval, float *temp, int t_type, float *dwp_temp,
        int dwp_type, float *rh)
/*******************************************************************

   met_cvt_ttd_to_rh

   purpose: convert temperature and dewpoint temp into relative humidity

     Uses the Magnus formula for es from:
     Reference "The relationship between relative humidity and the dewpoint 
     temperature in moist air", Mark G. Lawrence, BAMS, Feb, 2005, 
     pp. 225 - 233, DOI:10.1175/BAMS-86-2-225

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               nval             I      # of values in following
                                                arrays to convert
      float *           temp             I      temperature array
      int               t_type           I      temperature type
      float *           dwp_temp         I      dewpoint temperature
      int               dwp_type         I      dewpoint type
      float *           rh               O      relative humidity in %

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       5 Aug 2013      Original development

 *******************************************************************/ {
    int ival;
    float rh_lcl;
    double td_lcl, t_lcl;
    for (ival = 0; ival < nval; ival++) {
        t_lcl = (double) *(temp + ival);
        td_lcl = (double) *(dwp_temp + ival);
        /*
         *  the equations use t in C, so make sure of this
         */
        t_lcl = met_cvt_t_cvt(t_lcl, t_type, MET_UNITS__T_C);
        td_lcl = met_cvt_t_cvt(td_lcl, dwp_type, MET_UNITS__T_C);
        rh_lcl = 100. * (exp((MAGNUS_A1 * td_lcl) / (MAGNUS_B1 + td_lcl))
                / exp((MAGNUS_A1 * t_lcl) / (MAGNUS_B1 + t_lcl)));

        rh_lcl = (rh_lcl < 0) ? 0 : rh_lcl;
        rh_lcl = (rh_lcl > 100.) ? 100. : rh_lcl;

        *(rh + ival) = rh_lcl;
    }
    return 0;
}

double met_cvt_p_cvt(double val, int in_type, int out_type)
/*******************************************************************

   met_cvt_p_cvt

   purpose: convert units for pressure 

   Returns double of new pressure

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      double            val              I      input pressure
      int               in_type          I      type of input pressure
      int               out_type         I      type to convert to

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       12-May-2009     Original development

 *******************************************************************/ {
    float mult;

   /* assume HPa is normal */
    if (in_type == out_type)
        return val;
    else {
        switch (in_type) {
        case MET_UNITS__P_PA:
            mult = .01;
            break;
        case MET_UNITS__P_HPA:
            mult = 1.;
            break;
        }
        switch (out_type) {
        case MET_UNITS__P_PA:
            mult *= 100.;
            break;
        case MET_UNITS__P_HPA:
            mult *= 1.;
            break;
        }
        return mult * val;
    }
}

double met_cvt_t_cvt(double val, int in_type, int out_type)
/*******************************************************************

   met_cvt_t_cvt

   purpose: convert units for temperature 

   Returns double of new temperature

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      double            val              I      input temperature 
      int               in_type          I      type of input temperature
      int               out_type         I      type to convert to

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       12-May-2009     Original development

 *******************************************************************/ {
    if (in_type == out_type)
        return val;
    else if (in_type == MET_UNITS__T_K)
        return val - C_IN_K;
    else
        return val + C_IN_K;
}

double met_cvt_q_cvt(double val, int in_type, int out_type)
/*******************************************************************

   met_cvt_q_cvt

   purpose: convert units for specific humidity 

   Returns double of new specific humidity

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      double            val              I      input specific humidity
      int               in_type          I      type of input specific humidity
      int               out_type         I      type to convert to

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       12-May-2009     Original development

 *******************************************************************/ {
    if (in_type == out_type)
        return val;
    else if (in_type == MET_UNITS__Q_KG_KG)
        return val / 1000.;
    else
        return val * 1000.;
}
