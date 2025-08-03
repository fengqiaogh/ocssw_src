#include "glint.h"
#include <genutils.h>

#define DSIGN(A, B) (B >= 0 ? fabs(A) : -fabs(A))

double acoss(double x) {
    if (x > 1.0)
        return 0;
    else if (x < -1.0)
        return OEL_PI;
    else
        return acos(x);
}

double asinn(double x) {
    if (x > 1.0)
        return OEL_PI / 2;
    else if (x < -1.0)
        return -OEL_PI / 2;
    else
        return asin(x);
}

/**
 * @brief C  inc_angle  Incident angle (radians)
C  Rs  Reflectivity s polarized
C  Rp  Reflectivity p polarized
C
C       n1 sin(inc_angle) = n2 sin(refract_angle)
C       Fresnel's equations
C                                   tan(inc_angle-refract_angle)**2
C       Reflectivity(P polarized) = -------------
C                                   tan(inc_angle+refract_angle)**2
C
C                                     sin(inc_angle-refract_angle)**2
C       Reflectivity(S polarized) =  -------------
C                                     sin(inc_angle+refract_angle)**2
C
C  Where:
C       inc_angle  Incident angle
C       n1  Index refraction of Air
C       x2  Refracted angle
C       n2  Index refraction of Water
 *
 * @param inc_angle Incident angle
 * @param Effective_Refl  Effective Reflectivity
 * @param BiRefl Bireflectivity
 */
void reflec_both(double inc_angle, double *Effective_Refl, double *BiRefl,double nw) {
    double ref = nw; // water ref index, no wavelength dependencies
    double refract_angle, Rs, Rp;
    if (inc_angle < 0.00001) {
        *Effective_Refl = .0204078;
        *BiRefl = 0.;
    } else {
        refract_angle = asin(sin(inc_angle) / ref);  // Snell's law
        Rs = (sin(inc_angle - refract_angle) / sin(inc_angle + refract_angle));
        Rs *= Rs;
        Rp = (tan(inc_angle - refract_angle) / tan(inc_angle + refract_angle));
        Rp *= Rp;
        *Effective_Refl = Rs + Rp;
        *Effective_Refl = *Effective_Refl / 2.;
        *BiRefl = -Rs + Rp;
        *BiRefl = *BiRefl / 2.;
    }
}

void reflec_(float *inc_angle, float *Effective_Refl) {
    double BiRefl, effective_Refl, inc_angle_;
    const double ref = 4.0e0 / 3.0e0;  // water ref index, no wavelength dependencies

    inc_angle_ = (double)(*inc_angle);
    reflec_both(inc_angle_, &effective_Refl, &BiRefl, ref);
    *Effective_Refl = (float)effective_Refl;
}

void getglint_iqu(float senz, float solz, float raz, float ws, float chi, float *glint_coef,
                  float *glint_coef_q, float *glint_coef_u,double nw) {
    const double deg_cuttof = 1e-7;
    // from Cox & Munk, 54
    const double ws_gl = .04964e0;
    float ws_ = fmax(ws, 0.001f);
    double senz_ = senz * OEL_DEGRAD;
    double solz_ = solz * OEL_DEGRAD;
    double raz_ = raz * OEL_DEGRAD;
    if (senz_ == 0.0)
        senz_ = deg_cuttof;
    if (solz_ == 0.0)
        solz_ = deg_cuttof;
    if (raz_ == 0.0)
        raz_ = deg_cuttof;
    double argument = cos(senz_) * cos(solz_) -
                      sin(senz_) * sin(solz_) *
                          cos(raz_);  // dot product of two vectors, the incident ray and the reflected ray
    // relaz = senaz - 180 - solaz
    double omega = acoss(argument) / 2.0e0;
    if (omega == 0.0)
        omega = deg_cuttof;
    argument = (cos(senz_) + cos(solz_)) / (2.0e0 * cos(omega));
    double beta = acoss(argument);
    if (beta == 0.0)
        beta = deg_cuttof;
    double sigc = ws_gl * sqrt(ws_);
    double expon = -tan(beta) * tan(beta) / 2. / sigc / sigc;
    if (expon < -30.)
        expon = -30.;  //   ! trap underflow
    if (expon > 30.)
        expon = +30.;  // ! trap overflow
    double prob = exp(expon) / (2. * OEL_PI * sigc * sigc);
    double BiRefl, effective_refl;
    reflec_both(omega, &effective_refl, &BiRefl, nw);
    double cs_beta2 = cos(beta) * cos(beta);
    double cs_beta4 = cs_beta2 * cs_beta2;
    *glint_coef = effective_refl * prob / (4.0e0 * cos(senz_) * cs_beta4);
    double rot_ang;
    if (omega > 0.0001) {
        double CR = (cos(solz_) - cos(2. * omega) * cos(senz_)) / (sin(2. * omega) * sin(solz_));
        double SR = sin(solz_) * sin(OEL_PI - raz_) / sin(2. * omega);
        rot_ang = DSIGN(1.0, CR) * asinn(SR);
    } else {
        rot_ang = OEL_PI / 2;
    }
    double c2r = cos(2. * rot_ang);
    double s2r = sin(2. * rot_ang);
    *glint_coef_q = c2r * BiRefl / effective_refl;
    *glint_coef_u = -s2r * BiRefl / effective_refl;
}

void getglint_(float *senz, float *solz, float *raz, float *ws, float *chi, float *glint_coef) {
    float glint_coef_q, glint_coef_u;
    const double ref = 4.0e0 / 3.0e0;  // water ref index, no wavelength dependencies

    getglint_iqu(*senz, *solz, *raz, *ws, *chi, glint_coef, &glint_coef_q, &glint_coef_u, ref);
}

void getglint_iqu_(float *senz, float *solz, float *raz, float *ws, float *chi, float *glint_coef,
                   float *glint_coef_q, float *glint_coef_u) {
    const double ref = 4.0e0 / 3.0e0;  // water ref index, no wavelength dependencies

    getglint_iqu(*senz, *solz, *raz, *ws, *chi, glint_coef, glint_coef_q, glint_coef_u, ref);
}