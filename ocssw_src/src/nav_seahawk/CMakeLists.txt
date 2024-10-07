cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

ADD_LIBRARY(geomag70
  geomag70_py.c
)
#target_link_libraries(l2 ${LIBS})
target_include_directories(geomag70 PUBLIC ".")
set_target_properties(geomag70 PROPERTIES PREFIX "")
set_property(TARGET geomag70 PROPERTY POSITION_INDEPENDENT_CODE ON)


install (TARGETS geomag70 DESTINATION opt/lib)
