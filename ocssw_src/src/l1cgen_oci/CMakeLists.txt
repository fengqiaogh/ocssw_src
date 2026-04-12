cmake_minimum_required(VERSION 3.12)

set(LIBS filetype l1)
set(CMAKE_CXX_STANDARD 17)
find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})
find_package(Boost REQUIRED)
find_package(CGAL REQUIRED)
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_compile_definitions(GIT_COMMIT_HASH="${GIT_HASH}")
add_executable(l1cgen_oci search_l1c.cpp l1c_processor.cpp main_l1cgen.cpp l1c_input.cpp bin_l1b.cpp utils.cpp surface_intersector.cpp area_weighting.cpp surface_correction.cpp l1b_file.cpp l1c_file.cpp)

target_link_libraries(l1cgen_oci PRIVATE
        ${LIBS} Boost::headers CGAL::CGAL)  

install(TARGETS l1cgen_oci DESTINATION bin)
