cmake_minimum_required(VERSION 3.12)

set(LIBS l2 bin dataday)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

set(CMAKE_CXX_STANDARD 17)

# should be removed after boost update
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
        set(CMAKE_CXX_FLAGS "-DMACINTOSH ${CMAKE_CXX_FLAGS} -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION -DMAC_MEM")
endif ()

# for profiling with debug symbols
if(BUILD_L2BIN_PROFILE)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g")
endif ()
# use address sanitizer
if(BUILD_L2BIN_ASAN)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -fsanitize=address  -g")
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}  -fsanitize=address  -g")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif ()

add_executable(l2bin
        l2bin.cpp
        l2bin_input.cpp area_weighting.cpp l2bin_utils.cpp
)


target_link_libraries(l2bin
        ${LIBS}
        ${LibBin_LIBRARIES}
)


install(TARGETS l2bin  DESTINATION bin)
