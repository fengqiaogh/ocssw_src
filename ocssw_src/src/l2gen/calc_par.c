/******************************************************************************
 * This routine computes daily Photosynthetically Active Radiation (PAR) just
 * above the surface.
 *
 * The atmosphere and surface are modeled as a 2-layer system. The first layer
 * contains molecules and aerosols and is located above the second layer, the
 * cloud/surface layer. The TOA radiances are transformed into reflectance
 * and corrected for gaseous absorption. The resulting reflectances are
 * further corrected for molecular and aerosol scattering and combined to
 * yield the instantaneous albedo of the cloud/surface layer in the PAR
 * interval, A. This albedo is then used to compute instantaneous PAR. Daily
 * PAR is finally obtained by integrating instantaneous PAR from sunrise to
 * sunset, taking into account statistically diurnal effects on A.
 *
 *      *** Important: This subroutine uses an external file containing
 *      aerosol optical properties.
 *
 *      Output:
 *      PAR: Photosynthetically Active Radiation (mW/cm^2/um)
 *
 *    Key changes in v 2.0:
 *      * wavelength integration weighting function added
 *
 *    Key changes in v 1.8:
 *      * AsatMax renamed to AcldMax to indicate "cloud" emphasis
 *      * Use (Asat-AsBar) rather than Asat in the test that determines if
 *        Asat is too large
 *      * Added 670 section to aerosols table since we need that to choose
 *        the right model in GetPhase subroutine, + modified code in that
 *        routine
 *      * Changed TauCons from 20 to 15
 *
 *    Key changes in v 1.7:
 *      * Corrected Tau term in spherical albedo expressions
 *      * Modified initial Asat definition with AsBar
 *      * Put 'cloudy' test after computation of new Asat
 *      * Replaced archaic intrinsic routines (COSD, SIND, ACOSD, JMOD, LONG)
 *
 *    Key changes in v 1.6:
 *      * Angstrom exponents are expected to be primarily positive,
 *        that is, in the range (-0.1, 2.0) ... initial test was removed.
 *      * Surface albedo tables were added, with interpolation subroutines
 *      * TauAer(500nm) calculation was added (for Surface albedo tables)
 *      * Replaced AsBar calculation with call to As interpolation subroutine
 *      * Added test to estimate if pixel is clear or cloudy
 *      * Replaced AsAvg calculation with call to As interpolation subroutine
 *      * Removed AsBar and AsAvg from calculation of Aavg_spectral, in integration
 *      * Added Abar_spectral vs. AsAvg test for "function" integrand term
 *
 *    Key changes in v 1.5:
 *      * addition of water vapor term in gaseous transmittance in integration
 *      * comparison test for Abar_spectral vs. AsAvg
 *
 *      Authors: Robert Frouin <RFrouin@ucsd.edu>, scientific algorithms
 *              John McPherson <JMcPherson@ucsd.edu>, program structure
 *
 *      Adapted from the original FORTRAN to C and better integrated with l2gen
 *      by Robert Lossing and Sean Bailey
 *
 * References:
 *
 * Frouin, R., and B. Chertock, 1992: A technique for global monitoring of net
 * solar irradiance at the ocean surface. Part I: Model. J. Appl. Meteo., 31,
 * 1056-1066.
 *
 * Frouin, R., D. W. Ligner, and C. Gautier, 1989: A Simple analytical formula
 * to compute clear sky total and photosynthetically available solar irradiance
 * at the ocean surface. J. Geophys. Res., 94, 9731-9742.
 *
 * Tanre, D., M. Herman, P.-Y. Deschamps, and A. De Leffe, 1979: Atmospheric
 * modeling for Space measurements of ground reflectances, including
 * bi-directional properties. Appl. Optics, 18, 21,3587-21,3597.
 *
 * Young, D. F., P. Minnis, D. R. Doelling, G. G. Gibson, and T. Wong:
 * Temporal interpolation methods for the Clouds and the Earth's Radiant
 * Energy System (CERES) Experiment. J. Appl. Meteo., 37, 572-590.
 *
 *****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <timeutils.h>
#include "genutils.h"
#include "l12_proto.h"
#include "l2_struc.h"
#include "par_utils.h"

static int ini = 0;
luts_par luts_data;
float *grid_ozone[2];
float *grid_watvap[2];
float *grid_tg[4], *grid_td[4], *grid_rho[6];
float *grid_scalar_par[3];
float *grid_scalar_inst_par[3];
static int16_t year, month, mday, doy;
static double sec;
static size_t t_step;
static float observed_time;
static float *t_range;
const float A = 8.435e-03;
const float B = -1.225e-04;
const float C = 1.40e-04;
float Ps0 = EARTH_SURF_PRESSURE;
const float delta = 0.0095;
const float Gamma0 = 1.f / 3.f;
const size_t N_5nm = 60;
const float F0_5nm[] = {163.33, 167.94, 174.72, 175.26, 170.75, 159.98, 163.31, 177.76, 186.23, 199.60,
                        204.19, 203.33, 202.92, 200.22, 200.56, 204.09, 197.04, 190.83, 195.62, 193.83,
                        191.92, 193.99, 188.39, 181.65, 183.47, 188.85, 190.71, 187.82, 186.81, 187.72,
                        187.50, 185.28, 183.95, 184.37, 185.20, 185.36, 184.33, 180.52, 178.12, 177.98,
                        176.71, 175.58, 171.55, 170.52, 169.82, 167.42, 166.57, 165.00, 163.21, 160.81,
                        155.44, 153.10, 155.10, 154.64, 152.72, 150.99, 148.61, 146.63, 145.06, 142.50};

const float wl_5nm[] = {402.5, 407.5, 412.5, 417.5, 422.5, 427.5, 432.5, 437.5, 442.5, 447.5, 452.5, 457.5,
                        462.5, 467.5, 472.5, 477.5, 482.5, 487.5, 492.5, 497.5, 502.5, 507.5, 512.5, 517.5,
                        522.5, 527.5, 532.5, 537.5, 542.5, 547.5, 552.5, 557.5, 562.5, 567.5, 572.5, 577.5,
                        582.5, 587.5, 592.5, 597.5, 602.5, 607.5, 612.5, 617.5, 622.5, 627.5, 632.5, 637.5,
                        642.5, 647.5, 652.5, 657.5, 662.5, 667.5, 672.5, 677.5, 682.5, 687.5, 692.5, 697.5};
static float weight_5nm[60];
static float Denom_5nm = 0;


void calc_scalar_inst_par(l2str *l2rec, int ip, float par_above_ins, float *par_scalar_ins) {
    l1str *l1rec = l2rec->l1rec;
    size_t dim_cot = luts_data.scalar_inst_luts.dim_cot;
    const float *cot_grid = luts_data.scalar_inst_luts.cot;
    float cf_planar_grid[dim_cot];
    float cf_scalar_grid[dim_cot];
    float second_derivative_cf_scalar_par[dim_cot];

    size_t dims[] = {
            luts_data.scalar_inst_luts.dim_wind_speed,
            luts_data.scalar_inst_luts.dim_solz,
            dim_cot,
    };
    for (size_t i = 0; i < dim_cot; i++) {
        float point[3] = {l1rec->ws[ip], l1rec->csolz[ip], cot_grid[i]};
        cf_planar_grid[i] = interp3d(dims, point, grid_scalar_inst_par, luts_data.scalar_inst_luts.cf_pard_p);
        cf_scalar_grid[i] = interp3d(dims, point, grid_scalar_inst_par, luts_data.scalar_inst_luts.cf_pard_m);
    }
    float point[2] = {l1rec->ws[ip], l1rec->csolz[ip]};

    spline(cf_planar_grid, cf_scalar_grid, dim_cot, 0, 0, second_derivative_cf_scalar_par);
    // need to define interp2d
    float inst_PARd_c = interpnd(2, dims, point, grid_scalar_inst_par, luts_data.scalar_inst_luts.pard_p_cs);
    float inst_PARo_c = interpnd(2, dims, point, grid_scalar_inst_par, luts_data.scalar_inst_luts.pard_m_cs);
    float inst_PARo_o = interpnd(2, dims, point, grid_scalar_inst_par, luts_data.scalar_inst_luts.pard_m_oc);
    float inst_cf_obs = par_above_ins / inst_PARd_c / EINSTEIN * 3600 * 24;
    float inst_cf_obs_spline;
    splint(cf_planar_grid, cf_scalar_grid, second_derivative_cf_scalar_par, dim_cot, inst_cf_obs, &inst_cf_obs_spline);
    // correct for 2d splie behavior when cloudy
    if(inst_cf_obs < cf_planar_grid[1] && inst_cf_obs_spline < inst_cf_obs){
        inst_cf_obs_spline = inst_cf_obs;
    }
    float inst_cf = inst_PARo_o/inst_PARo_c;
    float slope = (inst_PARo_c - inst_PARo_o)/(1-inst_cf);
    *par_scalar_ins = (slope * (inst_cf_obs_spline - inst_cf) + inst_PARo_o)/3600/24;
}

void calc_scalar_par_mean_cosine(l2str *l2rec, int ip, float par_above, float par_c, float *scalar_par,
                                 float *mean_cosine) {
    if (par_c == BAD_FLT) {
        *scalar_par = BAD_FLT;
        *mean_cosine = BAD_FLT;
        return;
    }
    
    l1str *l1rec = l2rec->l1rec;
    float point[] = {l1rec->ws[ip], (float) doy, l1rec->lat[ip]};
    size_t dims[] = {
            luts_data.scalar_luts.dim_wind_speed,
            luts_data.scalar_luts.dim_doy,
            luts_data.scalar_luts.dim_latitude,
    };
    float PARo_1 = interp3d(dims, point, grid_scalar_par, luts_data.scalar_luts.PARo);
    float PARo_2 = interp3d(dims, point, grid_scalar_par, luts_data.scalar_luts.PARo_c);
    float PARc_above = interp3d(dims, point, grid_scalar_par, luts_data.scalar_luts.PARc_above);
    float cf = interp3d(dims, point, grid_scalar_par, luts_data.scalar_luts.CF);
    float S1 = (PARo_2 - PARo_1) / (1. - cf);
    const float cf_obs = par_above / PARc_above;
    float PARo_est = S1 * (cf_obs - cf) + PARo_1;
    *scalar_par = PARo_est;

    float mu1 = interp3d(dims, point, grid_scalar_par, luts_data.scalar_luts.mu_cosine);
    float mu2 = interp3d(dims, point, grid_scalar_par, luts_data.scalar_luts.mu_cosine_c);

    S1 = (mu2 - mu1) / (1. - cf);
    float mu_est = S1 * (cf_obs - cf) + mu1;
    *mean_cosine = mu_est;
}

float calc_par(l2str *l2rec, int ip, int nbands, float *Lt, float taua, float angstrom, float *wl, float *fo,
               float *ko3, float *taumolbar) {
    int i, band;
    float r[nbands], csolz, transg;
    float weight[nbands];
    float rp[nbands], transa, sa, ps0, del;
    float asymfac, beta, tau[nbands], taumol[nbands], tauaer[nbands];
    float gamma, cos2x, pmol;
    float csenz, transav;

    float par0, integtrans, delt, saavg;
    float tmstep, sz, cossz, tg, td;
    float sola;
    float asavg, func[11], ra[nbands], rc[nbands];

    float sumSa, sumrc, sumtdir, sumtdif, denom, tatemp, ravgsat;
    float tdir, tdif;
    float q, q4, taucons, q4tau, f, asat, asbar, abar, acldmax;
    float fcsolz, fcsenz, fcossz;
    float tg_oz, tg_wv;
    float trise, tset;
    float par;

    int widx490, widx510;
    float ang500, ta500, slope;

    int cloudy = FALSE;

    static float kavg_oz = 5.2e-05;
    static float kavg_wv = 0.002;
    //	need some memory allocation
    float *aerphasefunc;
    float *omega;
    float *modelAngstrom;

    l1str *l1rec = l2rec->l1rec;

    int16_t year, month, mday, doy;
    double sec;
    unix2yds(l1rec->scantime, &year, &doy, &sec);
    unix2ymds(l1rec->scantime, &year, &month, &mday, &sec);

    aerphasefunc = (float *) allocateMemory(nbands * sizeof(float), "aerphasefunc");
    omega = (float *) allocateMemory(nbands * sizeof(float), "omega");
    modelAngstrom = (float *) allocateMemory(nbands * sizeof(float), "modelAngstrom");

    int asat_is_max = 0;

    /****************************************************************************
     * Perform data checks; assign defaults if necessary:
     ** If data invalid and uncorrectable, PAR set to BAD_FLT missing value.
     **
     ***************************************************************************/

    if (l1rec->solz[ip] < 0.0f || l1rec->solz[ip] > 90.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("SolZen (outside valid range 0-90) = %f\n", l1rec->solz[ip]);
        par = BAD_FLT;
        return par;
    }

    if (l1rec->senz[ip] < 0.0f || l1rec->senz[ip] > 90.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("ViewZen (outside valid range 0-90) = %f\n", l1rec->senz[ip]);
        par = BAD_FLT;
        return par;
    }

    if (l1rec->delphi[ip] < -180.0f || l1rec->delphi[ip] > 180.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("delphi (outside valid range -180:180) = %f\n", l1rec->delphi[ip]);
        par = BAD_FLT;
        return par;
    }

    float dobson = l1rec->oz[ip];

    csolz = l1rec->csolz[ip];
    csenz = l1rec->csenz[ip];
    cos2x = powf(cosf((l1rec->scattang[ip] * M_PI / 180.0)), 2.0);

    if (l1rec->scattang[ip] < 0.0 || l1rec->scattang[ip] > 180.0) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("scatang (outside valid range 0-180) = %f\n", l1rec->scattang[ip]);
        printf("scatang = -(cossz*zosvz)+(sinsz*sinvz*cosra)\n");
        par = BAD_FLT;
        return par;
    }

    if (taua <= 0.0f) {
        taua = 0.1f;
    }

    float surfpress = l1rec->pr[ip];
    if (surfpress <= 0.0f) {
        surfpress = EARTH_SURF_PRESSURE;
    }

    /* Need to use conventional Angstrom exponent.
     if it's missing GetAerPhase should work and return a default value */

    GetAerPhase(l2rec, ip, nbands, angstrom, aerphasefunc, omega, modelAngstrom);
    if (dobson < 0.0f) {
        dobson = EstimateDobson(year, month, mday, l1rec->lat[ip]);
    }

    float watvap = l1rec->wv[ip];
    if (watvap < 0.0f) {
        watvap = EstimateWatVap(year, month, mday, l1rec->lat[ip]);
    }

    /*****************************************************************************
     * Transform TOA radiances into reflectances: *
     ****************************************************************************/
    // TODO: consider passing in rhot instead of Lt
    float Airmass = (1.f / csolz) + (1.f / csenz);

    for (band = 0; band < nbands; band++) {
        r[band] = (M_PI * (Lt[band])) / (fo[band] * l1rec->fsol * csolz);
        // Compute gaseous transmittance:
        transg = expf((-ko3[band] * dobson * Airmass));
        rp[band] = r[band] / transg;
    }

    /****************************************************************************
     * Compute reflectances of the cloud-surface subsystem:
     **
     ***************************************************************************/

    /* Setup for transmittance equation:*/

    ps0 = EARTH_SURF_PRESSURE;
    del = 0.0095f;

    pmol = (2.0f * (1.0f - del) * 0.75f * (1.0f + cos2x) + 3.0f * del) / (2.0f + del);

    asymfac = 0.6666667f;
    beta = 0.5f * (1.0f + asymfac);
    gamma = 1.0f - asymfac;

    /******************************************************************************
     * Compute diffuse atmospheric transmittance:
     **
     ******************************************************************************/
    sumSa = 0.0f;
    sumrc = 0.0f;
    sumtdir = 0.0f;
    sumtdif = 0.0f;
    denom = 0.0f;

    weight[0] = (wl[0] - 400.0) + (wl[1] - wl[0]) / 2.0;
    for (band = 1; band < nbands - 1; band++) {
        weight[band] = (wl[band] - wl[band - 1]) / 2.0 + (wl[band + 1] - wl[band]) / 2.0;
    }
    weight[nbands - 1] = (700.0 - wl[nbands - 1]) + (wl[nbands - 1] - wl[nbands - 2]) / 2.0;

    for (band = 0; band < nbands; band++) {
        taumol[band] = taumolbar[band] * (surfpress / ps0);

        tauaer[band] = taua * powf(wl[band] / (float) input->aer_wave_long, modelAngstrom[band]);

        tau[band] = taumol[band] + gamma * tauaer[band];

        ra[band] = ((taumol[band] * pmol) + (omega[band] * tauaer[band] * aerphasefunc[band])) /
                   (4.0f * csolz * csenz);

        /* Compute diffuse atmospheric transmittance:*/

        transa = expf(-(taumol[band] + tauaer[band]) / csolz) *
                 expf((0.52f * taumol[band] + beta * tauaer[band]) / csolz);
        transav = expf(-(taumol[band] + tauaer[band]) / csenz) *
                  expf((0.52f * taumol[band] + beta * tauaer[band]) / csenz);

        // Get direct and diffuse transmittances
        tdir = expf(-tau[band] / csolz);

        tdif =
                (expf(-tau[band] / csolz)) * (expf((0.52f * taumol[band] + beta * tauaer[band]) / csolz) - 1.f);
        sumtdir += (tdir * fo[band]);
        sumtdif += (tdif * fo[band]);

        /* Compute spherical albedo of the atmosphere:*/
        sa = expf(-(taumol[band] + gamma * tauaer[band])) * (0.92f * taumol[band] + gamma * tauaer[band]);

        sumSa += sa * fo[band] * weight[band];
        denom += fo[band] * weight[band];

        rc[band] = (rp[band] - ra[band]) / ((transa * transav) + (sa * (rp[band] - ra[band])));
        sumrc += (rc[band] * fo[band] * weight[band]);
    }

    // Derive tauaer(500) from (469,555)
    //  indexes into modelAngstrom need to use bindex
    switch (l1rec->l1file->sensorID) {
        case MODIST:
        case MODISA:
            widx490 = windex(469., wl, nbands);
            widx510 = windex(555., wl, nbands);
            break;
        case VIIRSN:
        case VIIRSJ1:
        case VIIRSJ2:
            widx490 = windex(490., wl, nbands);
            widx510 = windex(555., wl, nbands);
            break;
        case OCTS:
            widx490 = windex(490., wl, nbands);
            widx510 = windex(520., wl, nbands);
            break;
        default:
            widx490 = windex(490., wl, nbands);
            widx510 = windex(510., wl, nbands);
            break;
    }

    slope = (modelAngstrom[widx510] - modelAngstrom[widx490]) / (wl[widx510] - wl[widx490]);
    ang500 = modelAngstrom[widx490] + slope * (wl[widx510] - wl[widx490]);
    ta500 = taua * powf((500.0 / (float) input->aer_wave_long), ang500);

    /******************************************************************************
     * Compute daily average PAR:
     **
     ******************************************************************************/

    /* Set up terms constant in the integration:*/
    // Sa_bar:
    saavg = sumSa / denom;
    // R_bar:
    ravgsat = sumrc / denom;

    /******************************************************************************
     *Get simple bi-directional reflectance factor:
     **
     ******************************************************************************/

    q = 1.0f / (3.0f * (1.0f - 0.853f));
    q4 = 4.0f * q;

    taucons = TAUCONS_LOW;
    q4tau = q4 / (taucons + q4);

    fcsenz = (3.f / 7.f) * (1.0f + 2.0f * csenz);
    fcsolz = (3.f / 7.f) * (1.0f + 2.0f * csolz);

    f = (1.0f - q4tau * fcsolz) /
        (0.49f * ((1.f + 4.f * csolz * csenz) / (csolz + csenz)) - q4tau * fcsolz * fcsenz);

    /******************************************************************************
     * Get As_bar:
     **
     ******************************************************************************/
    asbar = interp_as_taulow(csolz, ta500);

    /******************************************************************************
     * Get Albedo from R_bar:
     **
     ******************************************************************************/
    asat = (ravgsat - asbar) * f + asbar;

    /******************************************************************************
     * Test Asat to see if it's too large; adjust as necessary:
     **
     ******************************************************************************/

    acldmax = 1.0f - (q4 / (TAUCONS_HIGH + q4)) * fcsolz;

    if ((asat - asbar) >= acldmax) {
        asat_is_max = 1;
        taucons = TAUCONS_HIGH;
        q4tau = q4 / (taucons + q4);
    } else {
        asat_is_max = 0;
        taucons = q4 * ((fcsolz / (1.0f - (asat - asbar))) - 1.0f);

        /***************************************************************************
         * Compute new Asat if needed:
         **
         ***************************************************************************/

        if (taucons > TAUCONS_LOW) {
            q4tau = q4 / (taucons + q4);

            f = (1.0f - q4tau * fcsolz) /
                (0.49f * ((1.f + 4.f * csolz * csenz) / (csolz + csenz)) - q4tau * fcsolz * fcsenz);
            asat = (ravgsat - asbar) * f + asbar;

            /* New Test:*/

            if ((asat - asbar) >= acldmax) {
                asat_is_max = 1;
                taucons = TAUCONS_HIGH;
                q4tau = q4 / (taucons + q4);
            }
        } else {
            taucons = TAUCONS_LOW;
            q4tau = q4 / (taucons + q4);
        }
    }

    // Test to estimate if pixel is clear or cloudy:
    cloudy = FALSE;
    if (asat > asbar)
        cloudy = TRUE;

    /*****************************************************************************
     * Get parameters used in integration time-steps:
     **
     *****************************************************************************/

    /*  Get the rise and set time for this day (for integral):*/

    triseset(doy, l1rec->lon[ip], l1rec->lat[ip], &trise, &tset);

    /*  Set up for 10 intervals*/
    delt = (tset - trise) / 10.0f;

    /***************************************************************************
     *	Get parameters used in integration time-steps:
     **
     ***************************************************************************/

    for (i = 0; i < 11; i++) {
        tmstep = trise + i * delt;
        int intyear, intday;
        intyear = year;
        intday = doy;
        sunangs_(&intyear, &intday, &tmstep, l1rec->lon + ip, l1rec->lat + ip, &sz, &sola);
        cossz = cosf(sz * (float) M_PI / 180.0f);
        fcossz = (3.f / 7.f) * (1.0f + 2.0f * cossz);

        /***************************************************************************
         * Get gaseous transmittance term:
         **
         ***************************************************************************/
        if (sz >= MAX_SOLZEN) {
            tg_oz = 0.0f;
            tg_wv = 0.0f;
        } else {
            // kavg_oz (5.2E-5) is valid for dobson in true dobson units, not
            // the ozone units in l2gen
            tg_oz = expf(-kavg_oz * powf((dobson * 1000. / cossz), 0.99f));
            tg_wv = expf(-kavg_wv * powf((watvap / cossz), 0.87f));
        }

        tg = tg_oz * tg_wv;
        /***************************************************************************
         * Get average diffuse atmospheric transmittance and As_bar terms: *
         ***************************************************************************/
        td = 0.0f;

        if (sz < MAX_SOLZEN) {
            for (band = 0; band < nbands; band++) {
                tatemp = expf(-(taumol[band] + tauaer[band]) / cossz) *
                         expf((0.52f * taumol[band] + beta * tauaer[band]) / cossz);
                td += tatemp * fo[band] * weight[band];
            }
            td /= denom;
        }

        if (cossz < 0.05) {
            cossz = 0.05;
        }
        if (cloudy) {
            asavg = interp_as_tauhigh(cossz);
        } else {
            asavg = interp_as_taulow(cossz, ta500);
        }

        /*****************************************************************************
         *Get A_bar term:
         **
         *****************************************************************************/

        if (asat_is_max) {
            abar = 1.0f - q4tau * fcossz;
        } else {
            abar = asat * (1.0f - q4tau * fcossz) / (1.0f - q4tau * fcsolz);
        }

        /*new test:*/
        if (abar >= 1.0f) {
            taucons = TAUCONS_HIGH;
            q4tau = q4 / (taucons + q4);
            abar = 1.0f - q4tau * fcossz;
        }

        /****************************************************************************
         * Get "function" integrand term:
         **
         ****************************************************************************/

        if (abar <= asavg) {
            func[i] = cossz * tg * td / (1.0 - abar * saavg);
        } else {
            func[i] = cossz * tg * td * (1.0f - abar) / (1.0f - asavg) / (1.0f - abar * saavg);
        }
    }

    /*	Use trapezoidal rule for integration:*/

    integtrans = 0.0f;

    for (i = 0; i < 10; i++) {
        integtrans += delt * (func[i] + func[i + 1]) / 2.0f;
    }

    /* Finally, determine the PAR:*/
    // par0 determined as the mean value [400-700nm] for Thuillier 2003
    // from run/data/common/Thuillier_F0.dat
    par0 = PAR_0;
    par = par0 * l1rec->fsol * integtrans / 24.0f;

    free(aerphasefunc);
    free(omega);
    free(modelAngstrom);

    return par;
}

extern float *par_planar_a_inst, *par_planar_b_inst;

float calc_par_impl_of_2023(l2str *l2rec, int ip, int nbands, float *Lt, float taua, float angstrom,
                            float *wl, float *fo, float *ko3, float *taumolbar, float *parb, float *parc) {
    if (ini == 0) {
        get_luts_data(l2rec, &luts_data);
        ini = 1;
        unix2yds(l2rec->l1rec->scantime, &year, &doy, &sec);
        unix2ymds(l2rec->l1rec->scantime, &year, &month, &mday, &sec);
        grid_ozone[0] = luts_data.ozonedims.days;
        grid_ozone[1] = luts_data.ozonedims.latitude;
        grid_watvap[0] = luts_data.watvapdims.days;
        grid_watvap[1] = luts_data.watvapdims.latitude;

        grid_tg[0] = luts_data.tgdims.wavelength;
        grid_tg[1] = luts_data.tgdims.air_mass;
        grid_tg[2] = luts_data.tgdims.water_vapor_pressure;
        grid_tg[3] = luts_data.tgdims.ozone_concentration;

        grid_td[0] = luts_data.tddims.wavelength;
        grid_td[1] = luts_data.tddims.air_mass;
        grid_td[2] = luts_data.tddims.angstrom_coefficients;
        grid_td[3] = luts_data.tddims.optical_depth_at_550_nm;

        grid_rho[0] = luts_data.rhodims.solar_zenith_angle;
        grid_rho[1] = luts_data.rhodims.view_zenith_angle;
        grid_rho[2] = luts_data.rhodims.relative_azimuth_angle;
        grid_rho[3] = luts_data.rhodims.angstrom_coefficients;
        grid_rho[4] = luts_data.rhodims.optical_depth_at_550_nm;
        grid_rho[5] = luts_data.rhodims.wavelength;

        grid_scalar_par[0] = luts_data.scalar_luts.wind_speed;
        grid_scalar_par[1] = luts_data.scalar_luts.doy;
        grid_scalar_par[2] = luts_data.scalar_luts.latitude;

        grid_scalar_inst_par[0] = luts_data.scalar_inst_luts.wind_speed;
        grid_scalar_inst_par[1] = luts_data.scalar_inst_luts.cos_solz;
        grid_scalar_inst_par[2] = luts_data.scalar_inst_luts.cot;
        observed_time = sec / 3600;
        {
            size_t ntimes = l2rec->l1rec->cld_rad->ntimes;
            t_step = ntimes;
            t_range = l2rec->l1rec->cld_rad->timecldrange;
        }

        weight_5nm[0] = (wl_5nm[0] - 400.0) + (wl_5nm[1] - wl_5nm[0]) / 2.0;
        for (size_t j = 1; j < N_5nm - 1; j++)
            weight_5nm[j] = (wl_5nm[j] - wl_5nm[j - 1]) / 2.0 + (wl_5nm[j + 1] - wl_5nm[j]) / 2.0;

        weight_5nm[N_5nm - 1] = (700.0 - wl_5nm[N_5nm - 1]) + (wl_5nm[N_5nm - 1] - wl_5nm[N_5nm - 2]) / 2.0;
        for (size_t band = 0; band < N_5nm; band++)
            Denom_5nm += F0_5nm[band] * weight_5nm[band];
        // printf("Denom_5nm %f\n", Denom_5nm);
    }
    // calcs go here
    l1str *l1rec = l2rec->l1rec;
    float par = BAD_FLT;
    // int iscan = l2rec->l1rec->iscan;
    if (l1rec->solz[ip] < 0.0f || l1rec->solz[ip] > 90.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("SolZen (outside valid range 0-90) = %f\n", l1rec->solz[ip]);
        return par;
    }

    if (l1rec->senz[ip] < 0.0f || l1rec->senz[ip] > 90.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("ViewZen (outside valid range 0-90) = %f\n", l1rec->senz[ip]);
        return par;
    }

    if (l1rec->delphi[ip] < -180.0f || l1rec->delphi[ip] > 180.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("delphi (outside valid range -180:180) = %f\n", l1rec->delphi[ip]);
        return par;
    }

    float csolz = l1rec->csolz[ip];

    if (l1rec->scattang[ip] < 0.0 || l1rec->scattang[ip] > 180.0) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("scatang (outside valid range 0-180) = %f\n", l1rec->scattang[ip]);
        printf("scatang = -(cossz*zosvz)+(sinsz*sinvz*cosra)\n");
        return par;
    }
    float weight[nbands];
    weight[0] = (wl[0] - 400.0) + (wl[1] - wl[0]) / 2.0;
    for (size_t band = 1; band < nbands - 1; band++) {
        weight[band] = (wl[band] - wl[band - 1]) / 2.0 + (wl[band + 1] - wl[band]) / 2.0;
    }
    weight[nbands - 1] = (700.0 - wl[nbands - 1]) + (wl[nbands - 1] - wl[nbands - 2]) / 2.0;
    // replace optical depth, surface pressure and angstrom coefficients with
    // their default values

    if (taua <= 0.0f) {
        taua = TAU_550_LOW_BOUND;
        angstrom = ANGSTROM_DEFAULT_VALUE;
    }
    if (angstrom < -1.0f)
        angstrom = ANGSTROM_DEFAULT_VALUE;
    float surfpress = l1rec->pr[ip];
    if (surfpress <= 0.0f) {
        surfpress = EARTH_SURF_PRESSURE;
    }

    float dobson = l1rec->oz[ip] * 1000;
    float watvap = l1rec->wv[ip];
    if (dobson < 0) {
        float point[] = {doy, l1rec->lat[ip]};
        size_t dims[] = {luts_data.ozonedims.dimdays, luts_data.ozonedims.dimlatitude};
        dobson = interpnd(2, dims, point, grid_ozone, luts_data.lut_dobson);
    }
    if (watvap < 0) {
        {
            float point[] = {doy, l1rec->lat[ip]};
            size_t dims[] = {luts_data.watvapdims.dimdays, luts_data.watvapdims.dimlatitude};
            watvap = interpnd(2, dims, point, grid_watvap, luts_data.lut_watvap);
        }
    }

    float rho[nbands];
    float airmass = kasten_equation(l1rec->solz[ip]);
    float airmass_0 = kasten_equation(l1rec->senz[ip]);
    {
        size_t dims[] = {luts_data.tgdims.dimwavelength, luts_data.tgdims.dimair_mass,
                         luts_data.tgdims.dimwater_vapor_pressure, luts_data.tgdims.dimozone_concentration};
        for (size_t band = 0; band < nbands; band++) {
            rho[band] = (M_PI * (Lt[band])) / (fo[band] * l1rec->fsol * csolz);
            if (rho[band] < 0.0 || rho[band] > 1.0) {
                parb[ip] = BAD_FLT;
                parc[ip] = BAD_FLT;
                return BAD_FLT;
            }
            float point[] = {wl[band], airmass, watvap, dobson};
            float tg_solzen = interp4d(dims, point, grid_tg,
                                       luts_data.lut_tg);  // the same so it seems
            point[1] = airmass_0;
            float tg_viewzen = interp4d(dims, point, grid_tg, luts_data.lut_tg);
            rho[band] = rho[band] / tg_viewzen / tg_solzen;
        }
    }
    float Ra[nbands], TauMol[nbands], TauAer[nbands], TransA[nbands], TransAv[nbands], Sa[nbands];
    size_t dims[] = {luts_data.tddims.dimwavelength, luts_data.tddims.dimair_mass,
                     luts_data.tddims.dimangstrom_coefficients, luts_data.tddims.dimoptical_depth_at_550_nm};
    size_t dims_rho[] = {
            luts_data.rhodims.dimsolar_zenith_angle, luts_data.rhodims.dimview_zenith_angle,
            luts_data.rhodims.dimrelative_azimuth_angle, luts_data.rhodims.dimangstrom_coefficients,
            luts_data.rhodims.dimoptical_depth_at_550_nm, luts_data.rhodims.dimwavelength};

    float relaz = l1rec->sola[ip] - l1rec->sena[ip];
    if (relaz < 0)
        relaz += 360.0;
    if (relaz > 180.0)
        relaz = 360.0 - relaz;
    for (size_t band = 0; band < nbands; band++) {
        TauMol[band] = ((A / pow(wl[band] / 1000.0, 4)) + (B / pow(wl[band] / 1000.0, 5)) +
                        (C / pow(wl[band] / 1000.0, 6))) *
                       (surfpress / Ps0);
        TauAer[band] = taua * pow((550.0 / wl[band]), angstrom);

        float point[] = {wl[band], airmass, angstrom, taua};
        TransA[band] = interp4d(dims, point, grid_td, luts_data.lut_td);


        point[1] = airmass_0;
        TransAv[band] = interp4d(dims, point, grid_td, luts_data.lut_td);
        // Compute spherical albedo of the atmosphere:
        Sa[band] =
                exp(-(TauMol[band] + Gamma0 * TauAer[band])) * (0.92 * TauMol[band] + Gamma0 * TauAer[band]);

        float point_rho[] = {l1rec->solz[ip], l1rec->senz[ip], relaz, angstrom, taua, wl[band]};
        Ra[band] = interp6d(dims_rho, point_rho, grid_rho, luts_data.lut_rho);
    }
    float Sa_bar = 0.0f;
    float Denom = 0.0f;
    for (size_t band = 0; band < nbands; band++) {
        Sa_bar += Sa[band] * fo[band] * weight[band];
        Denom += fo[band] * weight[band];
    }
    Sa_bar /= Denom;

    float CF_obs, TauCld_obs;
    float albe_obs[nbands];

    getcldalbe(l1rec->cld_rad->taucld[ip], l1rec->cld_rad->cfcld[ip], csolz, observed_time, t_range, albe_obs,
               &TauCld_obs, &CF_obs, t_step, wl, nbands);
    // averaged observed albedo
    float AlbeCld_obs_bar = 0.0;
    for (size_t band = 0; band < nbands; band++)
        AlbeCld_obs_bar += albe_obs[band] / (float) nbands;
    // Get simple bi-directional reflectance factor, needs a refer

    const float q = 1.0 / (3.0 * (1.0 - 0.853));
    const float q4 = 4.0 * q;
    float TauCons[nbands];
    float q4tau[nbands];
    float F[nbands];
    const float cosvz = l1rec->csenz[ip];
    const float fCosVZ = (3. / 7.) * (1.0 + 2.0 * cosvz);
    const float fCosSolZen = (3. / 7.) * (1.0 + 2.0 * csolz);

    for (size_t band = 0; band < nbands; band++) {
        if (TauCld_obs >= 5.0) {
            TauCons[band] = TauCld_obs;
        } else {
            TauCons[band] = 5.0;
        }

        q4tau[band] = q4 / (TauCons[band] + q4);
        F[band] = (1.0 - q4tau[band] * fCosSolZen) /
                  (0.49 * ((1. + 4. * csolz * cosvz) / (csolz + cosvz)) - q4tau[band] * fCosSolZen * fCosVZ);
    }
    // # new parameterization of As

    float fr;
    float As[nbands], Asat[nbands];
    float As_bar = 0.0f, Asat_bar = 0.0f;
    const float chl = 0.1;

    const float airmass2 = kasten_equation(l1rec->solz[ip]) * kasten_equation(l1rec->senz[ip]);
    const float wind_speed = l1rec->ws[ip];
    const float AsatMax = 1.0 - (q4 / (300.0 + q4)) * fCosSolZen;
    // solving quadratic equation
    float a, b, c;

    for (size_t band = 0; band < nbands; band++) {
        fr = (1 - exp(-(TauMol[band] + TauAer[band]) / csolz) / TransA[band]) * (1. - CF_obs) + CF_obs;
        As[band] = getosa(wl[band], l1rec->solz[ip], wind_speed, chl, fr, &luts_data);
        a = -(1 / F[band] - 1) * exp(-(TauMol[band] + TauAer[band]) * airmass2) * Sa[band];
        b = (1 / F[band] - 1) * exp(-(TauMol[band] + TauAer[band]) * airmass2) * Sa[band] * As[band] +
            (1 / F[band] - 1) * exp(-(TauMol[band] + TauAer[band]) * airmass2) +
            TransA[band] * TransAv[band] + Sa[band] * (rho[band] - Ra[band]);
        c = -(1 / F[band] - 1) * exp(-(TauMol[band] + TauAer[band]) * airmass2) * As[band] -
            (rho[band] - Ra[band]);
        const float discr = b * b / 4. / a / a - c / a;
        Asat[band] = As[band];
        if (discr > 0.0 && isfinite(discr) && a != 0) {
            const float d = -sqrt(discr) - b / 2. / a;
            if ((d > 0) && (d < 1)) {
                Asat[band] = d;
            } else {
                Asat[band] = sqrt(discr) - b / 2. / a;
            }
        }
        if (Asat[band] < As[band])
            Asat[band] = As[band];
        As_bar += As[band] * fo[band] * weight[band] / Denom;
        Asat_bar += Asat[band] * fo[band] * weight[band] / Denom;
    }
    // instantaneous PAR
    float Func_inst = 0.0f, Funcb_inst  = 0.0f;
    {
        float Tg_bar_inst = 0.0;
        float Td_bar_inst = 0.0;
        float Td_inst[nbands];
        float Tg_inst[nbands];
        float CosSZ = l1rec->csolz[ip];
        float fCosSZ = (3. / 7.) * (1.0 + 2.0 * CosSZ);
        for (size_t band = 0; band < nbands; band++) {
            airmass = kasten_equation(l1rec->solz[ip]);
            float point_tg[] = {wl[band], airmass, watvap, dobson};
            float point_td[] = {wl[band], airmass, angstrom, taua};
            size_t dims_tg[] = {luts_data.tgdims.dimwavelength, luts_data.tgdims.dimair_mass,
                                luts_data.tgdims.dimwater_vapor_pressure,
                                luts_data.tgdims.dimozone_concentration};
            size_t dims_td[] = {luts_data.tddims.dimwavelength, luts_data.tddims.dimair_mass,
                                luts_data.tddims.dimangstrom_coefficients,
                                luts_data.tddims.dimoptical_depth_at_550_nm};
            Tg_inst[band] = interp4d(dims_tg, point_tg, grid_tg, luts_data.lut_tg);
            Td_inst[band] = interp4d(dims_td, point_td, grid_td, luts_data.lut_td);
            Tg_bar_inst += Tg_inst[band] * fo[band] * weight[band] / Denom;
            Td_bar_inst += Td_inst[band] * fo[band] * weight[band] / Denom;
        }
        float Abar = Asat_bar >= AsatMax ? AsatMax : Asat_bar;
        if (Abar >= 1.0)
            Abar = 1.0 - q4 / (300. + q4) * fCosSZ;
        if (Abar < 0)
            Abar = As_bar;
        if (Abar <= As_bar)
            Func_inst = CosSZ * Tg_bar_inst * Td_bar_inst / (1.0 - As_bar * Sa_bar);
        else
            Func_inst = CosSZ * Tg_bar_inst * Td_bar_inst * (1. - Abar) / (1.0 - As_bar) / (1.0 - Abar * Sa_bar);
        Funcb_inst= CosSZ * Tg_bar_inst * Td_bar_inst * (1. - Abar) / (1.0 - Abar * Sa_bar);  // # below surface
        
        const float fo_bar = PAR_0_2023;  // 176.41f;
        par_planar_a_inst[ip] = Func_inst * l1rec->fsol * fo_bar / 3600 / 24.0; ///*fo_bar / 3600;
        par_planar_b_inst[ip] = Funcb_inst * l1rec->fsol * fo_bar / 3600 / 24.0; ///*fo_bar / 3600;


    }


    // time integration starts here
    float Func_spectral[t_step][nbands], Funcb_spectral[t_step][nbands], Funcc_spectral[t_step][nbands];
    float Func[t_step], Funcb[t_step], Funcc[t_step];
    float trise, tset;
    triseset(doy, l1rec->lon[ip], l1rec->lat[ip], &trise, &tset);

    size_t st_int = t_step;
    size_t end_int = 0;
    for (size_t it = 0; it < t_step; it++) {  // need to check
        float sz = get_solz(doy, t_range[it], l1rec->lon[ip], l1rec->lat[ip]);
        if (sz >= 90.0 || t_range[it] < trise || t_range[it] > tset) {
            Func[it] = 0.0f;
            Funcb[it] = 0.0f;
            Funcc[it] = 0.0f;
            for (size_t band = 0; band < nbands; band++) {
                Func_spectral[it][band] = 0;
                Funcb_spectral[it][band] = 0;
                Funcc_spectral[it][band] = 0;
            }

            continue;
        }
        if (st_int > it)
            st_int = it;
        if (end_int < it)
            end_int = it;
        float As2_bar = 0.0f, Tg_bar = 0.0f, Td_bar = 0.0f, AlbeCld_bar = 0.0f;
        float Abar_spectral[nbands], Abar = 0.0f;
        float Aavg_spectral[nbands], Aavg = 0.0f;
        float Tg[nbands];
        float Td[nbands];
        // float temp1[nbands], temp2[nbands], temp3[nbands];
        float As2[nbands];
        float CosSZ = cos(sz * M_PI / 180.);
        float albe_cld[nbands], cf, tau;
        getcldalbe(l1rec->cld_rad->taucld[ip], l1rec->cld_rad->cfcld[ip], CosSZ, t_range[it], t_range,
                   albe_cld, &tau, &cf, t_step, wl,
                   nbands);  // something is wrong with getcldalbedo
        if (CosSZ < 0.05)
            CosSZ = 0.05;
        float fCosSZ = (3. / 7.) * (1.0 + 2.0 * CosSZ);
        size_t dims_tg[] = {luts_data.tgdims.dimwavelength, luts_data.tgdims.dimair_mass,
                            luts_data.tgdims.dimwater_vapor_pressure,
                            luts_data.tgdims.dimozone_concentration};
        size_t dims_td[] = {luts_data.tddims.dimwavelength, luts_data.tddims.dimair_mass,
                            luts_data.tddims.dimangstrom_coefficients,
                            luts_data.tddims.dimoptical_depth_at_550_nm};
        airmass = kasten_equation(sz);
        for (size_t band = 0; band < nbands; band++) {
            float point_tg[] = {wl[band], airmass, watvap, dobson};
            float point_td[] = {wl[band], airmass, angstrom, taua};
            Tg[band] = interp4d(dims_tg, point_tg, grid_tg, luts_data.lut_tg);
            Td[band] = interp4d(dims_td, point_td, grid_td, luts_data.lut_td);
            fr = (1. - exp(-(TauMol[band] + TauAer[band]) / CosSZ) / Td[band]) * (1. - cf) + cf;
            As2[band] = getosa(wl[band], sz, wind_speed, chl, fr, &luts_data);
            As2_bar += As2[band] * fo[band] * weight[band] / Denom;
            Tg_bar += Tg[band] * fo[band] * weight[band] / Denom;
            Td_bar += Td[band] * fo[band] * weight[band] / Denom;
            AlbeCld_bar += albe_cld[band] / (float) nbands;
        }
        // for (size_t band = 0; band < N_5nm; band++) {
        //     float point_tg[] = {wl_5nm[band], airmass, watvap, dobson};
        //     float point_td[] = {wl_5nm[band], airmass, angstrom, taua};
        //     Tg_bar +=
        //         interp4d(dims_tg, point_tg, grid_tg, luts_data.lut_tg) * F0_5nm[band] * weight_5nm[band];
        //     Td_bar +=
        //         interp4d(dims_td, point_td, grid_td, luts_data.lut_td) * F0_5nm[band] * weight_5nm[band];
        // }
        // Tg_bar /= Denom_5nm;
        // Td_bar /= Denom_5nm;
        for (size_t band = 0; band < nbands; band++) {
            Aavg_spectral[band] = Asat[band] >= AsatMax ? AsatMax : Asat[band];
            Abar_spectral[band] = Aavg_spectral[band] +
                                  (albe_cld[band] + As2[band] - 2 * albe_cld[band] * As2[band]) -
                                  (albe_obs[band] + As[band] - 2 * albe_obs[band] * As[band]);
            if (Abar_spectral[band] >= 1.0)
                Abar_spectral[band] = 1.0 - q4 / (300. + q4) * fCosSZ;
            if (Abar_spectral[band] < 0)
                Abar_spectral[band] = As2[band];

            if (Abar_spectral[band] <= As2[band])
                Func_spectral[it][band] = CosSZ * Tg[band] * Td[band] / (1.0 - As2[band] * Sa[band]);
            else
                Func_spectral[it][band] = CosSZ * Tg[band] * Td[band] * (1. - Abar_spectral[band]) /
                                          (1.0 - As2[band]) / (1.0 - Abar_spectral[band] * Sa[band]);
            Funcb_spectral[it][band] = CosSZ * Tg[band] * Td[band] * (1.0 - Abar_spectral[band]) /
                                       (1.0 - Abar_spectral[band] * Sa[band]);  // # below surface

            Funcc_spectral[it][band] = CosSZ * Tg[band] * Td[band] / (1.0 - As2[band] * Sa[band]);
        }
        Aavg = Asat_bar >= AsatMax ? AsatMax : Asat_bar;
        Abar = Aavg + (AlbeCld_bar + As2_bar - 2 * AlbeCld_bar * As2_bar) -
               (AlbeCld_obs_bar + As_bar - 2 * AlbeCld_obs_bar * As_bar);
        if (Abar >= 1.0)
            Abar = 1.0 - q4 / (300. + q4) * fCosSZ;
        if (Abar < 0)
            Abar = As2_bar;
        if (Abar <= As2_bar)
            Func[it] = CosSZ * Tg_bar * Td_bar / (1.0 - As2_bar * Sa_bar);
        else
            Func[it] = CosSZ * Tg_bar * Td_bar * (1. - Abar) / (1.0 - As2_bar) / (1.0 - Abar * Sa_bar);
        Funcb[it] = CosSZ * Tg_bar * Td_bar * (1. - Abar) / (1.0 - Abar * Sa_bar);  // # below surface
        Funcc[it] = CosSZ * Tg_bar * Td_bar / (1.0 - As2_bar * Sa_bar);
    }

    float PAR_spectral[nbands], PARb_spectral[nbands], PARc_spectral[nbands];
    float IntegTrans, IntegTransb, IntegTransc;
    float Par = 0, Parb = 0, Parc = 0;
    const float min_shift = t_range[st_int] - trise;
    const float max_shift = tset - t_range[end_int];
    size_t time_point_start = st_int;
    size_t time_point_end = end_int;

    // end debug print
    for (size_t band = 0; band < nbands; band++) {
        IntegTrans = 0;
        IntegTransb = 0;
        IntegTransc = 0;
        // trapezoid integration
        for (size_t it = st_int; it < end_int; it++) {
            IntegTrans += (Func_spectral[it][band] + Func_spectral[it + 1][band]) / 2.0;
            IntegTransb += (Funcb_spectral[it][band] + Funcb_spectral[it + 1][band]) / 2.0;
            IntegTransc += (Funcc_spectral[it][band] + Funcc_spectral[it + 1][band]) / 2.0;
        }
        IntegTrans += min_shift * Func_spectral[time_point_start][band] / 2.0 +
                      max_shift * Func_spectral[time_point_end][band] / 2.0;
        IntegTransb += min_shift * Funcb_spectral[time_point_start][band] / 2.0 +
                       max_shift * Funcb_spectral[time_point_end][band] / 2.0;
        IntegTransc += min_shift * Funcc_spectral[time_point_start][band] / 2.0 +
                       max_shift * Funcc_spectral[time_point_end][band] / 2.0;

        PAR_spectral[band] = IntegTrans / 24.0 * fo[band] * l1rec->fsol;
        PARb_spectral[band] = IntegTransb / 24.0 * fo[band] * l1rec->fsol;
        PARc_spectral[band] = IntegTransc / 24.0 * fo[band] * l1rec->fsol;
        Par += PAR_spectral[band] * weight[band] / 300.0;  // needs to modify to fit integration
        // step 300/5
        Parb += PARb_spectral[band] * weight[band] / 300.0;
        Parc += PARc_spectral[band] * weight[band] / 300.0;
    }
    IntegTrans = 0;
    IntegTransb = 0;
    IntegTransc = 0;
    for (size_t it = st_int; it < end_int; it++) {
        IntegTrans += (Func[it] + Func[it + 1]) / 2.0;
        IntegTransb += (Funcb[it] + Funcb[it + 1]) / 2.0;
        IntegTransc += (Funcc[it] + Funcc[it + 1]) / 2.0;
    }

    IntegTrans += min_shift * Func[time_point_start] / 2.0 + max_shift * Func[time_point_end] / 2.0;
    IntegTransb += min_shift * Funcb[time_point_start] / 2.0 + max_shift * Funcb[time_point_end] / 2.0;
    IntegTransc += min_shift * Funcc[time_point_start] / 2.0 + max_shift * Funcc[time_point_end] / 2.0;
    const float fo_bar = PAR_0_2023;  // 176.41f;
    Par = IntegTrans / 24.0 * fo_bar * l1rec->fsol;
    Parb = IntegTransb / 24.0 * fo_bar * l1rec->fsol;
    Parc = IntegTransc / 24.0 * fo_bar * l1rec->fsol;
    if (l1rec->glint_coef[ip] > MAXGLINT) {
        Par = BAD_FLT;
        Parb = BAD_FLT;
        Parc = BAD_FLT;
    }
    parb[ip] = Parb;
    parc[ip] = Parc;
    return Par;
}