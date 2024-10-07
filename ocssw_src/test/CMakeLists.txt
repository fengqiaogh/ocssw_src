
set(FUZZ "0%")
if(DEFINED ENV{OCTEST_TOLERANCE})
  set(FUZZ "6%")
endif(DEFINED ENV{OCTEST_TOLERANCE})

enable_testing()

include(CData)
include(GetSubdirs)

set(CDATA_DEFAULT_URL "https://oceandata.sci.gsfc.nasa.gov/ocssw_test/")
set(CDATA_DEFAULT_STRIP_PARTS 1)
set(CDATA_DEFAULT_BASE "${CMAKE_SOURCE_DIR}/testdata")


# add all of the sub directories in the test directory
get_subdirs(SUBDIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(subdir ${SUBDIRS})
    add_subdirectory(${subdir})
endforeach(subdir)

