
# libgenutils Changelog

## <VERSION> - 2019
  * changed alloc_2d to allocate2d so the order of the arguments
    matches the order of accessing the array.

## argpar 2.1.1 - 2018-07-27
### Fixed
  * More thoroughly stripped leading/trailing quotes/spaces from key and value from par files

## <VERSION> - 2018-04-10

 - Created script to generate alloc_2d and ran it.
 - Added tests, reorganized alloc2d
 - argpar features added, small changes elsewhere
 - Moved olog into its proper home.
 - added a switch to turn off the par file option
 - implemented OBPG standard code formatting.
 - get ready for shared libs
 - removed roxml
 - Added small libraries formerly part of val-extract
 - made clo ignore dashes on the begining of an option key
 - Removed dashes from common options in order to be compatible with gnu style Dashes are stripped on the command line.
 - Fixed addAlias to check if NOT option exists.
 - Added par option to clo library
 - removed dependancy on timeutils

### Source Changed

  * alloc_2d.c
  * alloc_2d.h
  * alloc_2d.pl
  * allocateMemory.c
  * argpar.c
  * argpar.h
  * argpar-help.c
  * argpar-json.c
  * Changelog.md
  * clo.c
  * clo.h
  * CMakeLists.txt
  * endianess.c
  * fileFormatUtils.c
  * filesize.c
  * fread_swap.c
  * genutils_globals.c
  * genutils.h
  * getlut_file.c
  * isValidInt.c
  * lenstr.f
  * lowcase.c
  * lspline.c
  * nr_spline.c
  * olog/buffer.c
  * olog/buffer.h
  * olog.c
  * olog/file.c
  * olog/file.h
  * olog.h
  * olog_loader.c
  * olog/loader.h
  * olog-main.c
  * olog/stream.c
  * olog/streamf.c
  * olog/streamf.h
  * olog/stream.h
  * parse_file_name.c
  * passthebuck.h
  * phash.c
  * phash.h
  * pqueue.h
  * pqueue-ll.c
  * shash.c
  * shash.h
  * swapc_bytes.c
  * table_io.cpp
  * table_io.h
  * table_io_wrapper.h
  * trimBlanks.c
  * upcase.c
  * vincenty.c
  * vincenty.h
  * xmlUtils.c
  * xmlUtils.h

## <VERSION> - 2015-03-12
### Added
  * Changelog.md
