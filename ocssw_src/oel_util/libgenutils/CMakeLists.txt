cmake_minimum_required(VERSION 3.12)

find_package(Jansson REQUIRED)
include_directories(SYSTEM ${JANSSON_INCLUDE_DIRS})

set(LIBS pthread m ${JANSSON_LIBRARIES})

ADD_LIBRARY(genutils
  allocate2d.c
  allocate3d.c
  allocate4d.c
  allocate5d.c
  endianess.c
  parse_file_name.c
  swapc_bytes.c
  fread_swap.c
  genutils_globals.c
  clo.c
  table_io.cpp
  nr_spline.c
  lspline.c
  lenstr.f
  filesize.c
  getlut_file.c
  lowcase.c
  upcase.c
  allocateMemory.c
  trimBlanks.c
  isValidInt.c
  fileFormatUtils.c

  argpar.c
  argpar-help.c
  argpar-json.c
  phash.c
  shash.c
  pqueue-ll.c
  vincenty.c

  olog.c
  olog_loader.c
  olog/buffer.c
  olog/file.c
  olog/stream.c
  olog/streamf.c

  replace_ocroots.cpp
  resolution_utils.cpp
  readFileList.cpp
)

target_link_libraries(genutils ${LIBS})
target_include_directories(genutils PUBLIC ".")

install (TARGETS genutils DESTINATION lib)
