# ViirsCal_FOUND - true if library and headers were found
# LibViirsCal_INCLUDE_DIRS - include directories
# LibViirsCal_LIBRARIES - library directories
set(LibViirsCal_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libViirsCal${CMAKE_SHARED_LIBRARY_SUFFIX})
set(LibViirsCal_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/ViirsCal)
set(LibViirsCal_LIBRARIES ${LibViirsCal_LIBRARY})
set(LibViirsCal_INCLUDE_DIRS ${LibViirsCal_INCLUDE_DIR})

find_package(LibViirsCmn REQUIRED)
list(APPEND LibViirsCal_LIBRARIES ${LibViirsCmn_LIBRARIES})
list(APPEND LibViirsCal_INCLUDE_DIRS ${LibViirsCmn_INCLUDE_DIRS})

set(LibViirsCal_FOUND TRUE)

mark_as_advanced(FORCE LibViirsCal_INCLUDE_DIR LibViirsCal_LIBRARY)
