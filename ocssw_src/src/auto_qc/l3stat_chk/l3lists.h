/*
 W. Robinson, SAIC, 21 Apr 2005  eliminate chlor_a_K_490
 */
#ifndef L2LISTS_H /* avoid re-inclusion */
#define L2LISTS_H

#define L3M_PARAMS      14

#define NLW_412         0
#define NLW_443         1
#define NLW_490         2
#define NLW_510         3
#define NLW_555         4
#define NLW_670         5
#define ANGSTROM        6
#define CHLOR_A         7
#define K_490           8
#define CHLOR_A_K490    9
#define EPS_78          10
#define TAU_865         11

char *param_name[] = {
    "nLw_412", /* 412		*/
    "nLw_443", /* 443		*/
    "nLw_490", /* 490		*/
    "nLw_510", /* 510		*/
    "nLw_555", /* 555		*/
    "nLw_670", /* 670		*/
    "angstrom_510", /* Angstrom at 510 */
    "chlor_a", /* chlor_a 	*/
    "K_490", /* Diffuse atte */
    "eps_78", /* epsilon	*/
    "tau_865" /* tau 865	*/
};

char *param_flds[] = {
    "nLw_412_sum,nLw_412_sum_sq", /* 412          */
    "nLw_443_sum,nLw_443_sum_sq", /* 443          */
    "nLw_490_sum,nLw_490_sum_sq", /* 490          */
    "nLw_510_sum,nLw_510_sum_sq", /* 510          */
    "nLw_555_sum,nLw_555_sum_sq", /* 555          */
    "nLw_670_sum,nLw_670_sum_sq", /* 670          */
    "angstrom_510_sum,angstrom_510_sum_sq", /* Angstrom at 510 */
    "chlor_a_sum,chlor_a_sum_sq", /* chlor_a      */
    "K_490_sum,K_490_sum_sq", /* Diffuse atte */
    "eps_78_sum,eps_78_sum_sq", /* epsilon      */
    "tau_865_sum,tau_865_sum_sq" /* tau 865      */
};

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

