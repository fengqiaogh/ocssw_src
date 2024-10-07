# LibNovas_FOUND - true if library and headers were found
# LibNovas_INCLUDE_DIRS - include directories
# LibNovas_LIBRARIES - library directories

set(LibNovas_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libnovas${CMAKE_SHARED_LIBRARY_SUFFIX})
set(LibNovas_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/novas)

set(LibNovas_FOUND TRUE)
set(LibNovas_LIBRARIES ${LibNovas_LIBRARY} m)
set(LibNovas_INCLUDE_DIRS ${LibNovas_INCLUDE_DIR})

mark_as_advanced(FORCE LibNovas_INCLUDE_DIR LibNovas_LIBRARY)
