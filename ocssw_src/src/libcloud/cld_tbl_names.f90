module cld_tbl_names
  !  cld_tbl_names is devoted now to defining easy names for the index
  !  for a table name in the cloud_table_names.txt file in the share/
  !   shared/cloud/<inst> area
  !
  !  W. Robinson, SAIC, Sep 2021
  !
  !   common area (../common) files
  integer, parameter :: FAST_TRANS_COEFF_DRY = 0
  integer, parameter :: FAST_TRANS_COEFF_OZO = 1
  integer, parameter :: FAST_TRANS_COEFF_WCO = 2
  integer, parameter :: FAST_TRANS_COEFF_WTL = 3
  integer, parameter :: FAST_TRANS_COEFF_WTS = 4
  integer, parameter :: SNOW_ALBEDO = 5
  integer, parameter :: ECOLOGY_MAP = 6
  integer, parameter :: GLOBAL_EMISS = 7
  integer, parameter :: ALBEDO_CLIM = 8
  !
  !  sensor-specific area
  integer, parameter :: TRANSMITTANCE_COEFFS = 9
  integer, parameter :: SINGLE_SCAT = 10
  integer, parameter :: SCAT_ICE = 11
  integer, parameter :: SCAT_ICE_03_M_S = 12
  integer, parameter :: SCAT_ICE_07_M_S = 13
  integer, parameter :: SCAT_ICE_15_M_S = 14
  integer, parameter :: SCAT_ICE_SD_03_M_S = 15 
  integer, parameter :: SCAT_ICE_SD_07_M_S = 16
  integer, parameter :: SCAT_ICE_SD_15_M_S = 17
  integer, parameter :: SCAT_WATER = 18
  integer, parameter :: SCAT_WATER_03_M_S = 19
  integer, parameter :: SCAT_WATER_07_M_S = 20
  integer, parameter :: SCAT_WATER_15_M_S = 21
  integer, parameter :: SCAT_WATER_SD_03_M_S = 22
  integer, parameter :: SCAT_WATER_SD_07_M_S = 23
  integer, parameter :: SCAT_WATER_SD_15_M_S = 24

end module cld_tbl_names
