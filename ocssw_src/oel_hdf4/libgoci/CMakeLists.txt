cmake_minimum_required(VERSION 3.12)

#set_ocssw_policy()

find_package(HDF5 REQUIRED COMPONENTS C HL)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})

find_package(PROJ REQUIRED CONFIG)
list(APPEND LIBS PROJ::proj)

ADD_LIBRARY(goci
  goci.c
  goci_slot.c
)

target_link_libraries(goci ${LIBS})
target_include_directories(goci PUBLIC ".")

install (TARGETS goci DESTINATION lib)
