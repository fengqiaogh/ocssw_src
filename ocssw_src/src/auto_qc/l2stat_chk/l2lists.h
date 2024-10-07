#ifndef L2LISTS_H /* avoid re-inclusion */
#define L2LISTS_H

#define NLW_412         0
#define NLW_443         1
#define NLW_490         2
#define NLW_510         3
#define NLW_555         4
#define NLW_670         5
#define ANGSTROM        6
#define CHLOR_A         7
#define K_490           8
#define EPS_78          9
#define TAU_865         10

char *param_name[] = {
    "nLw_412", /* 412		*/
    "nLw_443", /* 443		*/
    "nLw_490", /* 490		*/
    "nLw_510", /* 510		*/
    "nLw_555", /* 555		*/
    "nLw_670", /* 670		*/
    "angstrom_510", /* Angstrom coeff at 510 - 865 nm */
    "chlor_a", /* chlor_a 	*/
    "K_490", /* Diffuse atte */
    "eps_78", /* epsilon	*/
    "tau_865", /* tau 865	*/
};

float32 slope[] = {
    0.001, /*  412  	*/
    0.001, /*  443		*/
    0.001, /*  490		*/
    0.001, /*  510		*/
    0.001, /*  555		*/
    0.0001, /*  670		*/
    0.0002, /*  angstrom_510 */
    0.001, /*  chlor_a	*/
    0.0002, /*  Diffuse att */
    0.01, /*  epsilon	*/
    0.0001, /*  tau_865	*/
};

float32 intercept[] = {0, 0, 0, 0, 0, 0, 0, 32.0, 0, 0, 0};

char *flag_names[] = {
    "ATMFAIL", /* 0 */
    "LAND", /* 1 */
    "BADANC", /* 2 */
    "HIGLINT", /* 3 */
    "HILT", /* 4 */
    "HISATZEN", /* 5 */
    "COASTZ", /* 6 */
    "NEGLW", /* 7 */
    "STRAYLIGHT", /* 8 */
    "CLDICE", /* 9 */
    "COCCOLITH", /* 10 */
    "TURBIDW", /* 11 */
    "HISOLZEN", /* 12 */
    "HITAU", /* 13 */
    "LOWLW", /* 14 */
    "CHLFAIL", /* 15 */
    "NAVWARN", /* 16 */
    "ABSAER", /* 17 */
    "TRICHO", /* 18 */
    "MAXAERITER", /* 19 */
    "MODGLINT", /* 20 */
    "CHLWARN", /* 21 */
    "ATMWARN", /* 22 */
    "DARKPIXEL", /* 23 */
    "SPARE", /* 24 */
    "SPARE", /* 25 */
    "SPARE", /* 26 */
    "SPARE", /* 27 */
    "SPARE", /* 28 */
    "SPARE", /* 29 */
    "SPARE", /* 30 */
    "SPARE", /* 31 */
};

#endif /* L2LISTS_H */

