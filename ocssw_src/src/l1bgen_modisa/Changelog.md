
# l1bgen_modisa Changelog

## 6.2.1 - 2019-06-18
 - moved large static arrays off the stack

### Source Changed
 * CMakeLists.txt
 * Preprocess.c
 * L1B_Setup.c

## 6.2 - 2018-04-09
 - isolated MODIS extract metadata code
 - removed extra libs from modis CMakeLists.txt
 - Collection "6.1" of MODIS L1 processing code: l1agen_modis v6.0.6, geogen_modis v6.1.0, l1bgen_modisa v6.2.1, l1bgen_modist v6.2.2.
 - Cleaned up our changes to MODIS L1 code
 - Updated MODIS L1 code to be compatible with opt (aka lib3) r14227, and simplified changes wrt SDST delivery so as to make future integration easier.
 - moved MODIS L1 files out of subdirectory
  
### Source Changed

  * Changelog.md
  * CMakeLists.txt
  * Emissive_Cal.c
  * error_messages.txt
  * from_SDST/HISTORY.txt
  * from_SDST/MOD_PR02AQUA.mk
  * from_SDST/MOD_PR02AQUA_pr.txt
  * from_SDST/MYD021KM.fs
  * from_SDST/MYD021KM.mcf
  * from_SDST/MYD02HKM.fs
  * from_SDST/MYD02HKM.mcf
  * from_SDST/MYD02OBC.fs
  * from_SDST/MYD02OBC.mcf
  * from_SDST/MYD02QKM.fs
  * from_SDST/MYD02QKM.mcf
  * from_SDST/PACKINGLIST.txt
  * from_SDST/PGE02_HISTORY.txt
  * from_SDST/README.txt
  * Granule.c
  * Granule.h
  * HDF_Lib.c
  * HISTORY.txt
  * L1B.c
  * L1B_Setup.c
  * L1B_Tables.c
  * L1B_Tables.h
  * Metadata.c
  * Metadata.h
  * MetadataP.h
  * MOD_PR02AQUA.mk
  * MOD_PR02AQUA.pcf
  * MOD_PR02AQUA_pr.txt
  * MYD021KM.fs
  * MYD021KM.mcf
  * MYD02HKM.fs
  * MYD02HKM.mcf
  * MYD02OBC.fs
  * MYD02OBC.mcf
  * MYD02QKM.fs
  * MYD02QKM.mcf
  * MYD_PR02.pff.doc
  * PACKINGLIST.txt
  * PGE02_HISTORY.txt
  * PGS_MODIS_36110.h
  * Preprocess.c
  * PreprocessP.h
  * README.txt
  * Reflective_Cal.c

## <VERSION> - 2015-03-12
### Added
  * Changelog.md
