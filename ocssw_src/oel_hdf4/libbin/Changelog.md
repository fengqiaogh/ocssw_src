
# libbin Changelog

## <VERSION> - 2019-02-01

 - moved files from bin++ into this library and deleted bin++


## <VERSION> - 2018-04-10

 - first reorg of subsensor defaults files
 - use new sensor info functions
 - implemented OBPG standard code formatting
 - Allowed the eclipse IDE to reformat code according to K&R rules, and made a few other format changes to make the 32-bit and 64-bit versions as parallel as possible.
 - Changed HDF4 types to standard types:     int16, int32, uint8, uint16 -> same, with _t suffix     float32 -> float; float64 -> double
 - get ready for shared libs
 - made units bigger and check for buffer overrun
 - modified deflate and chunking code
 - got rid of compiler warnings

### Source Changed
  * bin_csub64.c
  * bin_csub.c
  * Changelog.md
  * CMakeLists.txt
  * get_beg_ext64.c
  * get_beg_ext.c
  * getl3b64.h
  * getl3b.h
  * getl3b_misc64.c
  * getl3b_misc.c
  * l3b_misc64.c
  * l3b_misc.c
  * meta_l3b64.h
  * meta_l3b.h
  * ncdfbin_utils.c
  * ncdfbin_utils.h
  * put_l3b64.c
  * put_l3b.c
  * read_l3b_meta64.c
  * read_l3b_meta.c
  * seabin64.h
  * seabin.h
  * seaproto64.h
  * seaproto.h
  * seaprotoi64.h
  * seaprotoi.h
  * write_l3b_meta64.c
  * write_l3b_meta.c
  * wr_vdata.c

## 2.0 - 2018-02-14
### only kept the 64bit version of the code. deleted 32 bit version.

## <VERSION> - 2015-03-12
### Added
  * Changelog.md
