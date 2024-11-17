#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#define NBIG 10U
#define NRAD 31U
#define NPHI 4U
#define NPHASE 16U
#define NWAVE 6U
#define NGAUS (2U * NRAD - 1U)
#define NNG 50U
#define NG (2U * NNG)
#define NUM 75U
#define NCHL 6U
#define NSUN 6U
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

static float PHSA[NPHASE][NWAVE][NBIG][NGAUS][NPHI], PHSR[NBIG][NGAUS][NPHI];
static float APHSRADA[NPHASE][NWAVE][NBIG][NGAUS][NPHI];
static float MU[NGAUS], PDIV[NG], PWT[NG], THETA[NGAUS];
static float PHSRADA[NGAUS][NGAUS][NPHI], PHSRADR[NGAUS][NGAUS][NPHI];
static float TAUA_RAT[NPHASE][NWAVE];
static float RAD1[NRAD][NPHI];
static float tst[2], tdf[2];
static float TSTARR[NWAVE][NRAD], TSTARA[NPHASE][NWAVE][NRAD];
namespace morel_vars {
    float BRDF[NWAVE][NSUN][NCHL][NRAD][NPHI];
    float SmBRDF[NWAVE][NSUN][NCHL][NBIG][NPHI];
    float Theta0[NSUN], Chl[NCHL], Wave[NWAVE], Thetav[NRAD];
    float LChl[NCHL];
}  // namespace morel_vars
bool initialized = false;

const float aindex = 1.334f;

float fresref(float muair, float index = aindex) {
    float theta1 = std::acos(muair);
    float stheta1 = std::sin(theta1);
    float stheta2 = stheta1 / index;
    float theta2 = std::asin(stheta2);
    float fresref = ((index - 1) / (index + 1)) * ((index - 1) / (index + 1));
    if (theta1 > 0.0) {
        float tan_ratio = std::tan(theta1 - theta2) / std::tan(theta1 + theta2);
        float sin_ratio = std::sin(theta1 - theta2) / std::sin(theta1 + theta2);
        tan_ratio *= tan_ratio;
        sin_ratio *= sin_ratio;
        fresref = 0.5f * (tan_ratio + sin_ratio);
    }
    return fresref;
}

void read_morel_data(const std::string &data_path) {
    using namespace morel_vars;
    std::string INFL_FOURIER31 =data_path + "/NEW_Morel_NRAD31-1-EDITED";
    std::string INFL_FOURIER10 = data_path + "/NEW_Morel_SMALL-1-EDITED";
    std::string line;

    {
        std::ifstream lut_file(INFL_FOURIER31.c_str(), std::ios::in);
        if (!lut_file.good()) {
            std::cerr << "-Error opening the Morel LUT txt file: " << INFL_FOURIER31 << std::endl;
            exit(EXIT_FAILURE);
        }
        for (size_t iwave = 0; iwave < NWAVE; iwave++) {
            for (size_t isun = 0; isun < NSUN; isun++) {
                for (size_t ichl = 0; ichl < NCHL; ichl++) {
                    std::getline(lut_file, line);
                    // std::cout << line << std::endl;
                    {
                        std::istringstream iss(line);
                        iss >> Wave[iwave] >> Theta0[isun] >> Chl[ichl];
                    }
                    for (size_t iphi = 0; iphi < NPHI; iphi++) {
                        std::getline(lut_file, line);
                        // std::cout << line << std::endl;
                        std::string number_line;
                        for (size_t il = 0; il < 7; il++) {
                            std::getline(lut_file, line);
                            number_line += line;
                        }
                        // std::cout << number_line << std::endl;
                        std::istringstream iss(number_line);
                        for (size_t irad = 0; irad < NRAD; irad++) {
                            iss >> BRDF[iwave][isun][ichl][irad][iphi];
                            // std::cout << BRDF[iwave][isun][ichl][irad][iphi]
                            //           << std::endl;
                        }
                    }
                }
            }
        }
        lut_file.close();
    }

    {
        std::ifstream lut_file(INFL_FOURIER10.c_str(), std::ios::in);
        if (!lut_file.good()) {
            std::cerr << "-Error opening the Morel LUT txt file: " << INFL_FOURIER10 << std::endl;
            exit(EXIT_FAILURE);
        }
        for (size_t iwave = 0; iwave < NWAVE; iwave++) {
            for (size_t isun = 0; isun < NSUN; isun++) {
                for (size_t ichl = 0; ichl < NCHL; ichl++) {
                    std::getline(lut_file, line);
                    // std::cout << line << std::endl;
                    {
                        std::istringstream iss(line);
                        iss >> Wave[iwave] >> Theta0[isun] >> Chl[ichl];
                    }
                    for (size_t iphi = 0; iphi < NPHI; iphi++) {
                        std::getline(lut_file, line);
                        // std::cout << line << std::endl;
                        std::string number_line;
                        for (size_t il = 0; il < 2; il++) {
                            std::getline(lut_file, line);
                            number_line += line;
                        }
                        // std::cout << number_line << std::endl;
                        std::istringstream iss(number_line);
                        for (size_t ibig = 0; ibig < NBIG; ibig++) {
                            iss >> SmBRDF[iwave][isun][ichl][ibig][iphi];
                            // std::cout <<
                            // SmBRDF[iwave][isun][ichl][ibig][iphi]
                            //           << std::endl;
                        }
                    }
                }
            }
        }
        lut_file.close();
    }
    for (size_t i = 0; i < NCHL; i++) {
        LChl[i] = std::log10(Chl[i]);
    }
}

void read_luts_txt() {
    std::string ocdataroot = std::getenv("OCDATAROOT");
    if (ocdataroot.empty()) {
        std::cerr << "-Error- : OCDATAROOT is not defined. Exiting ..." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string common_path = ocdataroot + "/eval" + "/common" + "/dtran_brdf";
    // file pathes
    std::string  infl_aer = common_path +"/Aerosols_Partial_Inegr.dat";
    std::string  infl_ray = common_path + "/Rayleigh_Partial_Inegr.dat";
    std::string  infl_rat = common_path + "/spec_var_EDITED.dat";
    //
    std::ifstream aer(infl_aer.c_str(), std::ios::in);
    if (!aer.good()) {
        std::cerr << "-Error-: opening the Aerosol LUT txt file: " << infl_aer << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    size_t line_counter = 0;
    // read the text Aerosol LUT file line by line.
    // let's convert it to netcdf
    for (auto &iphase: PHSA) {
        for (auto &iwave: iphase) {
            for (size_t m = 0; m < NPHI; m++) {
                for (size_t jup = 0; jup < NGAUS; jup++) {
                    line_counter++;
                    std::getline(aer, line);  // read the header
                    {
                        line_counter++;
                        std::getline(aer, line);
                        //                        std::cout << line << std::endl;
                        std::istringstream iss(line);
                        if (!(iss >> iwave[0][jup][m] >> iwave[1][jup][m] >> iwave[2][jup][m] >>
                                  iwave[3][jup][m] >> iwave[4][jup][m])) {
                            std::cerr << "-Error-: Couldn't parse the line # " << line_counter << " : "
                                      << line << std::endl;
                            exit(EXIT_FAILURE);
                        };
                    }
                    {
                        line_counter++;
                        std::getline(aer, line);
                        //                        std::cout << line << std::endl;
                        std::istringstream iss(line);
                        if (!(iss >> iwave[5][jup][m] >> iwave[6][jup][m] >> iwave[7][jup][m] >>
                                  iwave[8][jup][m] >> iwave[9][jup][m])) {
                            std::cerr << "-Error-: Couldn't parse the line # " << line_counter << " : "
                                      << line << std::endl;
                            exit(EXIT_FAILURE);
                        };
                    }
                }
            }
        }
    }
    aer.close();

    // read the text Rayleigh LUT file line by line.
    line_counter = 0;
    std::ifstream rayleigh(infl_ray.c_str(), std::ios::in);
    if (!rayleigh.good()) {
        std::cerr << "-Error-: opening the Rayleigh LUT txt file: " << infl_ray << std::endl;
        exit(EXIT_FAILURE);
    }

    for (size_t m = 0; m < NPHI; m++) {
        for (size_t jup = 0; jup < NGAUS; jup++) {
            line_counter++;
            std::getline(rayleigh, line);  // read the header
            {
                line_counter++;
                std::getline(rayleigh, line);
                //                std::cout << line << std::endl;
                std::istringstream iss(line);
                if (!(iss >> PHSR[0][jup][m] >> PHSR[1][jup][m] >> PHSR[2][jup][m] >> PHSR[3][jup][m] >>
                          PHSR[4][jup][m])) {
                    std::cerr << "-Error-: Couldn't parse the line # " << line_counter << " : " << line
                              << std::endl;
                    exit(EXIT_FAILURE);
                };
            }
            {
                line_counter++;
                std::getline(rayleigh, line);
                //                std::cout << line << std::endl;
                std::istringstream iss(line);
                if (!(iss >> PHSR[5][jup][m] >> PHSR[6][jup][m] >> PHSR[7][jup][m] >> PHSR[8][jup][m] >>
                          PHSR[9][jup][m])) {
                    std::cerr << "-Error-: Couldn't parse the line # " << line_counter << " : " << line
                              << std::endl;
                    exit(EXIT_FAILURE);
                };
            }
        }
    }
    rayleigh.close();
    // read the text Tau LUT file line by line.
    std::ifstream rat(infl_rat.c_str(), std::ios::in);
    if (!rat.good()) {
        std::cerr << "-Error-: opening the Tau LUT txt file: " << infl_rat << std::endl;
        exit(EXIT_FAILURE);
    }

    for (auto &iphase: TAUA_RAT) {
        std::getline(rat, line);  // read the header
        for (float &iwave: iphase) {
            std::getline(rat, line);
            //            std::cout << line << std::endl;
            std::istringstream iss(line);
            size_t wv;
            iss >> wv >> iwave;
            // std::cout << wv << " =  tau is " << TAUA_RAT[iphase][iwave]
            //           << std::endl;
        }
    }
    rat.close();

    // read transmittance LUTs

    std::string t_aer_path = common_path + "/tstar_aerosol.dat";
    std::string t_ray_path = common_path + "/tstar_rayleigh.dat";

    std::ifstream t_aer_file(t_aer_path.c_str(), std::ios::in);
    if (!t_aer_file.good()) {
        std::cerr << "-Error-: opening the transmittance LUT txt file: " << t_aer_path << std::endl;
        exit(EXIT_FAILURE);
    }

    for (auto &iphase: TSTARA) {
        for (auto &iwave: iphase) {
            std::getline(t_aer_file, line);
            // std::cout << line << std::endl;
            std::string number_line;
            for (size_t il = 0; il < 6; il++) {
                std::getline(t_aer_file, line);
                number_line += line;
            }
            // std::cout << number_line << std::endl;
            std::istringstream iss(number_line);
            for (size_t irad = 0; irad < NRAD - 3; irad++) {
                iss >> iwave[irad];
                std::cout << iwave[irad] << std::endl;
            }
        }
    }
    t_aer_file.close();

    std::ifstream t_ray_file(t_ray_path.c_str(), std::ios::in);
    if (!t_ray_file.good()) {
        std::cerr << "-Error-: opening the transmittance LUT txt file: " << t_ray_path << std::endl;
        exit(EXIT_FAILURE);
    }
    for (auto &iwave: TSTARR) {
        std::getline(t_ray_file, line);
        // std::cout << line << std::endl;
        std::string number_line;
        for (size_t il = 0; il < 6; il++) {
            std::getline(t_ray_file, line);
            number_line += line;
        }
        // std::cout << number_line << std::endl;
        std::istringstream iss(number_line);
        for (size_t irad = 0; irad < NRAD - 3; irad++) {
            iss >> iwave[irad];
            // std::cout << TSTARR[iwave][irad] << std::endl;
        }
    }
    t_ray_file.close();
    read_morel_data(common_path.c_str());
}

void Morel_BRDF(float Sun, float Chlor, float MBRDF[NWAVE][NRAD][NPHI], float SmMBRDF[NWAVE][NBIG][NPHI]) {
    using namespace morel_vars;
    int32_t isun = -1;
    int32_t ichl = -1;
    for (size_t i = 0; i < NSUN; i++) {
        if (Sun < Theta0[i]) {
            isun = i;
            break;
        }
    }
    for (size_t i = 0; i < NCHL; i++) {
        if (Chlor < Chl[i]) {
            ichl = i;
            break;
        }
    }
    float LChlor = std::log10(Chlor);
    //  C  Interpolate in the tables to get the water radiance

    float sun_interp = (Sun - Theta0[isun - 1]) / (Theta0[isun] - Theta0[isun - 1]);
    float chl_interp = (LChlor - LChl[ichl - 1]) / (LChl[ichl] - LChl[ichl - 1]);

    for (size_t iwave = 0; iwave < NWAVE; iwave++) {
        for (size_t i = 0; i < NRAD; i++) {
            for (size_t m = 0; m < NPHI; m++) {
                float interp_term =
                        (1. - sun_interp) * (1. - chl_interp) * BRDF[iwave][isun - 1][ichl - 1][i][m];
                interp_term =
                        interp_term + sun_interp * (1. - chl_interp) * BRDF[iwave][isun][ichl - 1][i][m];
                interp_term =
                        interp_term + (1. - sun_interp) * chl_interp * BRDF[iwave][isun - 1][ichl][i][m];
                interp_term = interp_term + sun_interp * chl_interp * BRDF[iwave][isun][ichl][i][m];
                MBRDF[iwave][i][m] = interp_term;
            }
        }
    }
    for (size_t iwave = 0; iwave < NWAVE; iwave++) {
        for (size_t i = 0; i < NBIG; i++) {
            for (size_t m = 0; m < NPHI; m++) {
                float interp_term =
                        (1. - sun_interp) * (1. - chl_interp) * SmBRDF[iwave][isun - 1][ichl - 1][i][m];
                interp_term =
                        interp_term + sun_interp * (1. - chl_interp) * SmBRDF[iwave][isun][ichl - 1][i][m];
                interp_term =
                        interp_term + (1. - sun_interp) * chl_interp * SmBRDF[iwave][isun - 1][ichl][i][m];
                interp_term = interp_term + sun_interp * chl_interp * SmBRDF[iwave][isun][ichl][i][m];
                SmMBRDF[iwave][i][m] = interp_term;
            }
        }
    }
}

extern "C" void diff_tran_corr_(int *iphase, float *solz, float *senz, float *phi, float *chl, float *taua,
                                float *correct) {
    std::cout
            << "\n\n\nDifftran BRDF exits; Uncomment to continue. Exiting ... Source file - dtranbrdf.cpp. Now we know the usage!\n\n\n\n"
            << std::endl;
    exit(EXIT_FAILURE);
    float MBRDF[NWAVE][NRAD][NPHI];
    float SmMBRDF[NWAVE][NBIG][NPHI];
    static float tr[] = {0.3132, 0.2336, 0.1547, 0.1330, 0.0957, 0.0446};
    int32_t iview = -1;
    if (!initialized) {
        read_luts_txt();
        for (size_t i = 0; i < NGAUS; i++) {
            THETA[i] = static_cast<float>(i) * 90.0f / static_cast<float>(NRAD - 1);
            MU[i] = std::cos(THETA[i] * M_PI / 180.0);
        }
        initialized = true;
    }
    Morel_BRDF(*solz, *chl, MBRDF, SmMBRDF);
    for (size_t i = 0; i < NRAD; i++) {
        if (*senz < THETA[i]) {
            iview = i;
            break;
        }
    }
    float Ta = *taua;
    size_t mgaus = NGAUS;
    size_t maxphi = NPHI;
    for (size_t iwave = 0; iwave < NWAVE; iwave++) {
        float tstar1, tstar2, ans1, ans2, tdiff1, tdiff2;
        for (size_t jup = iview - 1; jup <= iview; jup++) {
            //            float adelphi = (180.0 / M_PI) * (*phi);
            size_t jdn = mgaus + 1 - jup;
            float ya1 = 0.0f;
            float yr1 = 0.0f;
            for (size_t m = 0; m < maxphi; m++) {
                float order = static_cast<float>(m);
                float x1 = 0;
                float v1 = 0;
                for (size_t i = 0; i < NBIG; i++) {
                    float x2 = PHSA[*iphase][iwave][i][jup][m] * SmMBRDF[iwave][i][m];
                    float v2 = PHSR[i][jup][m] * SmMBRDF[iwave][i][m];
                    x1 += x2;
                    v1 += v2;
                }
                if (m > 0) {
                    x1 *= 2;
                    v1 *= 2;
                }
                ya1 += std::cos(order * *phi) * x1;
                yr1 += std::cos(order * *phi) * v1;
            }
            float za1 = 0.0f;
            float zr1 = 0.0f;
            for (size_t m = 0; m < maxphi; m++) {
                float order = static_cast<float>(m);
                float x1 = 0;
                float v1 = 0;
                for (size_t i = 0; i < NBIG; i++) {
                    float x2 = PHSA[*iphase][iwave][i][jdn][m] * SmMBRDF[iwave][i][m];
                    float v2 = PHSR[i][jdn][m] * SmMBRDF[iwave][i][m];
                    x1 += x2;
                    v1 += v2;
                }
                if (m > 0) {
                    x1 *= 2;
                    v1 *= 2;
                }
                za1 += std::cos(order * *phi) * x1;
                zr1 += std::cos(order * *phi) * v1;
            }
            float yl = 0.0f;
            for (size_t m = 0; m < maxphi; m++) {
                float order = static_cast<float>(m);
                float fac = m > 0 ? 2.0f : 1.0f;
                yl += MBRDF[iwave][jup][m] * fac * std::cos(order * *phi);
            }
            float rfres = fresref(MU[jup]);
            float aint_pl_a = (ya1 + rfres * za1) / yl / (1.0f - rfres);
            float aint_pl_r = (yr1 + rfres * zr1) / yl / (1.0f - rfres);
            float plterm_a = (1.0f - aint_pl_a / 2.0f) / MU[jup];
            float plterm_r = (1.0f - aint_pl_r / 2.0f) / MU[jup];
            float t_diff_a = std::exp(-plterm_a * Ta * TAUA_RAT[*iphase][iwave]);
            float t_diff_r = std::exp(-plterm_r * tr[iwave]);
            if (jup == iview - 1) {
                tstar1 =
                        TSTARR[iwave][jup] * std::pow(TSTARA[*iphase][iwave][jup], Ta * TAUA_RAT[*iphase][iwave]);
                tdiff1 = t_diff_a * t_diff_r;
                ans1 = (tdiff1 - tstar1) / (tstar1);
            } else {
                tstar2 =
                        TSTARR[iwave][jup] * std::pow(TSTARA[*iphase][iwave][jup], Ta * TAUA_RAT[*iphase][iwave]);
                tdiff2 = t_diff_a * t_diff_r;
                ans2 = (tdiff2 - tstar2) / (tstar2);
            }
        }
        float slope_ans = (ans2 - ans1) / (THETA[iview] - THETA[iview - 1]);
//        float slope_tst = (tstar2 - tstar1) / (THETA[iview] - THETA[iview - 1]);
//        float slope_tdf = (tdiff2 - tdiff1) / (THETA[iview] - THETA[iview - 1]);
        correct[iwave] = ans1 + slope_ans * (*senz - THETA[iview - 1]);
    }
};