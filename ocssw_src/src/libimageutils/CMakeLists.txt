cmake_minimum_required(VERSION 3.12)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})
set(LIBS ${PNG_LIBRARIES})

add_library (imageutils
  color_quant.c
)

target_link_libraries(imageutils ${LIBS} )
target_include_directories(imageutils PUBLIC ".")

add_executable(imageExample
    imageExample.c
)

target_link_libraries(imageExample
    imageutils
)

install (TARGETS imageutils DESTINATION lib)
#install (TARGETS imageExample DESTINATION bin)
