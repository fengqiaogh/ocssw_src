cmake_minimum_required (VERSION 3.12)

project (OEL_HDF4)

file(GLOB files "lib*")

foreach(dir ${files})
if (BUILD_HISTORICAL OR NOT((${dir} MATCHES "seawifs") OR (${dir} MATCHES "czcs")))
  add_subdirectory(${dir})
   endif()
endforeach(dir)

add_subdirectory(filetype)
add_subdirectory(h5io)
