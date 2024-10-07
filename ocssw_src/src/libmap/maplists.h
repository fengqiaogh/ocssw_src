#ifndef MAPLISTS_H /* avoid re-inclusion */
#define MAPLISTS_H

#define L3M_PARAMS      17

#define NLW_412         0
#define NLW_443         1
#define NLW_490         2
#define NLW_510         3
#define NLW_555         4
#define NLW_670         5
#define CZCS_PIGMENT    6
#define CHLOR_A         7
#define K_490           8
#define CHOR_A_K_490    9
#define ANGSTROM_510    10
#define EPSILON         11
#define TAU_865         12
#define PIXELS          13
#define SCENES          14
#define NDVI            15
#define BIOSPHERE       16

char *parmname_list[] = {
    "Normalized water-leaving radiance at 412 nm",
    "Normalized water-leaving radiance at 443 nm",
    "Normalized water-leaving radiance at 490 nm",
    "Normalized water-leaving radiance at 510 nm",
    "Normalized water-leaving radiance at 555 nm",
    "Normalized water-leaving radiance at 670 nm",
    "CZCS-like pigment concentration",
    "Chlorophyll a concentration",
    "Diffuse attenuation coefficient at 490 nm",
    "Ratio of chlor_a to K_490",
    "Epsilon of aerosol correction at 670 and 865 nm",
    "Angstrom coefficient at 510 nm",
    "Aerosol optical thickness at 865 nm",
    "Number of pixels per bin",
    "Number of scenes per bin",
    "Vegitation Index",
    "Chlorophyll a Concentration and Vegitation Index"
};

char *unit_list[] = {
    "mW cm^-2 um^-1 sr^-1", /* 412		*/
    "mW cm^-2 um^-1 sr^-1", /* 443		*/
    "mW cm^-2 um^-1 sr^-1", /* 490		*/
    "mW cm^-2 um^-1 sr^-1", /* 510		*/
    "mW cm^-2 um^-1 sr^-1", /* 555		*/
    "mW cm^-2 um^-1 sr^-1", /* 670		*/
    "mg m^-3", /* CZCS_pigment */
    "mg m^-3", /* chlor_a	*/
    "m^-1", /* Diffuse atte */
    "mg m^-2", /* chlor_a_K_490*/
    "", /* angstrom_510	*/
    "", /* epsilon	*/
    "", /* tau 865	*/
    "pixels", /* pixels	*/
    "scenes", /* scenes	*/
    "", /* ndvi  	*/
    "chl_a: mg m^-3, ndvi: dimensionless" /* biosphere  	*/
};

char *scaling_list[] = {
    "linear", /*  412		*/
    "linear", /*  443		*/
    "linear", /*  490		*/
    "linear", /*  510		*/
    "linear", /*  555		*/
    "linear", /*  670		*/
    "logarithmic", /*  CZCS_pigment*/
    "logarithmic", /*  chlor_a	*/
    "linear", /*  Diffuse att */
    "linear", /* chlor_a_K_490*/
    "linear", /*  angstrom 	*/
    "linear", /*  epsilon	*/
    "linear", /*  tau 865	*/
    "linear", /*  pixels	*/
    "linear", /*  scenes	*/
    "linear", /*  ndvi	*/
    "chl_a: logarithmic, ndvi: linear" /*  biosphere	*/
};

float32 slope_list[] = {
    0.02, /*  412  	*/
    0.02, /*  443		*/
    0.02, /*  490		*/
    0.02, /*  510		*/
    0.02, /*  555		*/
    0.004, /*  670		*/
    0.015, /*  CZCS-like	*/
    0.015, /*  chlor_a	*/
    0.004, /*  Diffuse att */
    0.1, /*  chlor_a_K_490*/
    0.02, /*  angstrom	*/
    0.01, /*  epsilon	*/
    0.005, /*  tau_865	*/
    1.0, /*  pixels	*/
    1.0, /*  scenes	*/
    0.005, /*  ndvi	*/
    0.0 /*  biosphere	*/
};

float32 intercept_list[] = {0, 0, 0, 0, 0, 0, -2.0, -2.0, 0, 0, -0.5, 0, 0, 0, 0, -0.35, 0.0};

float32 base = 10.0;

#endif /* MAPATTR_H */

