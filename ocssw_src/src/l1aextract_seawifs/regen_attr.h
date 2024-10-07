#ifndef REGEN_ATTR_H
#define REGEN_ATTR_H

#define MAXVAL  255
#define NBANDS  8
#define VECTORS  3
#define NAV_COEFS 6
#define NAV_FLAGS 8
#define ENG_FLAGS 4

#define NADIR  0
#define FORWARD  1
#define AFT  2
#define CHANGE  3

#define PROD_NAME "Product Name"
#define LONGNAME "long_name"
#define SLOPE  "slope"
#define INTERCEPT "intercept"
#define RANGE  "valid_range"
#define UNITS  "units"

/*** Scan-Line Attribute data set names */
#define MSEC  "msec"
#define ENGQUAL  "eng_qual"
#define SFLAGS  "s_flags"
#define SLAT  "slat"
#define SLON  "slon"
#define CLAT  "clat"
#define CLON  "clon"
#define ELAT  "elat"
#define ELON  "elon"
#define CSOLZ  "csol_z"
#define TILT  "tilt"

#define SSATP  "s_satp"
#define SZEROP  "s_zerop"

/*** Geophysical Data set names */
#define NLW412  "nLw_412"
#define NLW443  "nLw_443"
#define NLW490  "nLw_490"
#define NLW510  "nLw_510"
#define NLW555  "nLw_555"
#define LA670  "La_670"
#define LA865  "La_865"
#define CZCSPIGMENT "CZCS_pigment"
#define CHLORA  "chlor_a"
#define K490  "K_490"
#define EPSILON  "eps_68"
#define TAU865  "tau_865"
#define L2FLAGS  "l2_flags"

/*** Raw SeaStar Data */
#define GAIN  "gain"
#define TDI  "tdi"

#define SCID  "sc_id"
#define SCTTAG  "sc_ttag"
#define SCSOH  "sc_soh"
#define INSTTLM  "inst_tlm"
#define L1ADATA  "l1a_data"
#define STARTSYN "start_syn"
#define STOPSYN  "stop_syn"
#define DARKREST "dark_rest"

/*** Converted Telemetry dataset names */
#define INSTANA  "inst_ana"
#define INSTDIS  "inst_dis"
#define SCANA  "sc_ana"
#define SCDIS  "sc_dis"
#define SCANTEMP "scan_temp"
#define SIDE  "side"

/*** Navigation dataset names */
#define ORBVEC  "orb_vec"
#define LVERT  "l_vert"
#define SUNREF  "sun_ref"
#define ATTANG  "att_ang"
#define SENMAT  "sen_mat"
#define SCANELL  "scan_ell"
#define NFLAG  "nflag"

/*** Sensor Tilt dataset names */
#define NTILTS   "ntilts"
#define TILT_FLAGS "tilt_flags"
#define TILT_RANGES "tilt_ranges"
#define TILT_LATS "tilt_lats"
#define TILT_LONS "tilt_lons"

/*** Sensor Calibration dataset names */
#define ENTRYYEAR "entry_year"
#define ENTRYDAY "entry_day"
#define MIRROR  "mirror"
#define TIMEFACTOR "time_factor"
#define COUNTS  "counts"
#define RADS  "rads"

#endif /* REGEN_ATTR_H */
