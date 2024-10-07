cmake_minimum_required (VERSION 3.12)

project (OEL_HDF4)

file(GLOB files "lib*")

foreach(dir ${files})
  add_subdirectory(${dir})
endforeach(dir)

add_subdirectory(filetype)
add_subdirectory(h5io)
