cmake_minimum_required(VERSION 3.12)

find_package(PGSTK REQUIRED)
include_directories(${PGSTK_INCLUDE_DIRS})
list(APPEND LIBS ${PGSTK_LIBRARIES})

add_executable(l0cnst_write_modis
  l0cnst_write_modis.c
)

target_link_libraries(l0cnst_write_modis
  ${LIBS}
)

install (TARGETS l0cnst_write_modis DESTINATION bin)

