cmake_minimum_required(VERSION 3.12)

find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
  
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
add_definitions(-D_DEFAULT_SOURCE -D_BSD_SOURCE)
set(LIBS ${LibGenutils_LIBRARIES})

add_library(valextract val_extract.c)
target_include_directories(valextract PUBLIC ".")
target_link_libraries(valextract m genutils ${NETCDF_LIBRARIES})
install (TARGETS valextract DESTINATION lib)

add_executable(val_extract val_extract-main.c)
target_link_libraries(val_extract valextract ${LibGenutils_LIBRARIES})

add_executable(sstval_extract sstval_extract-main.c)
target_link_libraries(sstval_extract valextract ${LibGenutils_LIBRARIES})

install(TARGETS val_extract sstval_extract DESTINATION bin)

