cmake_minimum_required(VERSION 3.12)

find_package(PGSTK REQUIRED)
include_directories(${PGSTK_INCLUDE_DIRS})
set(LIBS ${PGSTK_LIBRARIES})

add_executable(l0chunk_modis
  l0chunk_modis.c
)

target_link_libraries(l0chunk_modis
  ${LIBS}
)

install (TARGETS l0chunk_modis DESTINATION bin)

