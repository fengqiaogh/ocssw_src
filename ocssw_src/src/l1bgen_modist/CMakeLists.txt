cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

find_package(PGSTK REQUIRED)
include_directories(${PGSTK_INCLUDE_DIRS})
list(APPEND LIBS ${PGSTK_LIBRARIES})

list(APPEND LIBS h4metadata)

add_executable(l1bgen_modist
  L1B.c
  Preprocess.c
  L1B_Tables.c
  Granule.c
  HDF_Lib.c
  Metadata.c
  Reflective_Cal.c
  Emissive_Cal.c
  L1B_Setup.c
)

add_definitions(-DALLOW_MISSING_GRANULES -DNOCHECKLUT -DSEADAS -Wno-unused-local-typedefs)

target_link_libraries(l1bgen_modist
  ${LIBS}
)

install (TARGETS l1bgen_modist DESTINATION bin)

