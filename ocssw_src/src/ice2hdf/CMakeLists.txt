cmake_minimum_required(VERSION 3.12)

add_executable(ice2hdf
  ice2hdf.c
)

target_link_libraries(ice2hdf 
    hdt4utils
)

install (TARGETS ice2hdf DESTINATION bin)
