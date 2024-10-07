cmake_minimum_required(VERSION 3.12)

# set these to make CGAL run in a reaesonable amount of time even when
# the rest of the project is built in Debug mode
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")

set(LIBS "")

find_package(CGAL)
include_directories(${CGAL_INCLUDE_DIRS})
list(APPEND LIBS ${CGAL_LIBRARIES})

ADD_LIBRARY(cgal_interp
  cgal_interp.cpp
)

target_link_libraries(cgal_interp ${LIBS})
target_include_directories(cgal_interp PUBLIC ".")

install (TARGETS cgal_interp DESTINATION lib)
