
# l1bgen_modist Changelog

## 6.2.1 - 2019-06-18
 - moved large static arrays off the stack

### Source Changed
 * CMakeLists.txt
 * Preprocess.c
 * L1B_Setup.c

## 6.2 - 2018-04-09

 - isolated MODIS extract metadata code
 - removed extra libs from modis CMakeLists.txt
 - Cleaned up our changes to MODIS L1 code
 - Testing MODIS Terra L1B TEB crosstalk correction (C6.1)
 - Updated MODIS L1 code to be compatible with opt (aka lib3) r14227, and simplified changes wrt SDST delivery so as to make future integration easier.
 - moved MODIS L1 files out of subdirectory

### Source Changed:
  * Changelog.md
  * CMakeLists.txt
  * Emissive_Cal.c
  * error_messages.txt
  * from_SDST/HISTORY.txt
  * from_SDST/MOD021KM.fs
  * from_SDST/MOD021KM.mcf
  * from_SDST/MOD02HKM.fs
  * from_SDST/MOD02HKM.mcf
  * from_SDST/MOD02OBC.fs
  * from_SDST/MOD02OBC.mcf
  * from_SDST/MOD02QKM.fs
  * from_SDST/MOD02QKM.mcf
  * from_SDST/MOD_PR02.mk
  * from_SDST/MOD_PR02_pr.txt
  * from_SDST/PACKINGLIST.txt
  * from_SDST/PGE02_HISTORY.txt
  * from_SDST/PGE02.README.txt
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
  * MOD021KM.fs
  * MOD021KM.mcf
  * MOD02HKM.fs
  * MOD02HKM.mcf
  * MOD02OBC.fs
  * MOD02OBC.mcf
  * MOD02QKM.fs
  * MOD02QKM.mcf
  * MOD_PR02.mk
  * MOD_PR02.pcf
  * MOD_PR02.pff.doc
  * MOD_PR02_pr.txt
  * PACKINGLIST.txt
  * PGE02_HISTORY.txt
  * PGE02.README.txt
  * Preprocess.c
  * Preprocess.h
  * PreprocessP.h
  * README.txt
  * Reflective_Cal.c

## <VERSION> - 2015-03-12
### Added
  * Changelog.md
