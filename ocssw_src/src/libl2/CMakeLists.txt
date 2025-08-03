cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBS l1 dfutils)

set(CMAKE_CXX_STANDARD 17)
# print more debug info
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Debug set")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG_LOG=1")
endif ()

# for profiling with debug symbols
if(BUILD_L2BIN_PROFILE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g")
endif ()
# use address sanitizer
if(BUILD_L2BIN_ASAN)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -fsanitize=address -g")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}  -fsanitize=address  -g")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif ()

# should be removed after boost update
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
        set(CMAKE_CXX_FLAGS "-DMACINTOSH ${CMAKE_CXX_FLAGS} -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION")
endif ()
ADD_LIBRARY(l2
        setupflags.cpp
        expand3D.cpp
        get_geospatial.cpp
        find_variable.c
        l2_reader.cpp
        l2_metadata.cpp
        l2_variable.cpp
        l2_utils.cpp
        l2_wrapper.cpp
)

add_executable(test_l2 test_code.cpp)

target_link_libraries(l2 ${LIBS})
target_link_libraries(test_l2 l2)
target_include_directories(l2 PUBLIC ".")

install(TARGETS l2 DESTINATION lib)
