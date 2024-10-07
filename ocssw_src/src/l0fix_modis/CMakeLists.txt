cmake_minimum_required(VERSION 3.12)

find_package(PGSTK REQUIRED)
include_directories(${PGSTK_INCLUDE_DIRS})
set(LIBS ${PGSTK_LIBRARIES})

add_executable(l0fix_modis
  l0fix_modis.c
)

target_link_libraries(l0fix_modis
  ${LIBS}
)

install (TARGETS l0fix_modis DESTINATION bin)

