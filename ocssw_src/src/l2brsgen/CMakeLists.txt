cmake_minimum_required(VERSION 3.12)

set(LIBS genutils l2 timeutils)

if(BUILD_HISTORICAL)
    set(CMAKE_C_FLAGS "-D BUILD_HISTORICAL")
    set(CMAKE_CXX_FLAGS "-D BUILD_HISTORICAL")
endif(BUILD_HISTORICAL)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})
list(APPEND LIBS ${PNG_LIBRARIES})

add_executable(l2brsgen
  main_l2brsgen.c
  l2brsgen_input.c
  put_l2brs.c
)

target_link_libraries(l2brsgen ${LIBS})

install (TARGETS l2brsgen DESTINATION bin)
