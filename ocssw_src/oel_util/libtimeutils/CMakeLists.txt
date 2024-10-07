cmake_minimum_required(VERSION 3.12)

find_package(Boost REQUIRED)
set(LIBS genutils m)

ADD_LIBRARY(timeutils
  addmsec.c
  date2ydmsec.c
  gmt_offset.c
  now.c
  unix2isodate.c
  unix2yds.c
  unix2ydhmsf.c
  unix2ymds.c
  unix2ymdhms.c
  yd2md.c
  ydhms2ydmsec.c
  ydhmsf.c
  yds2unix.c
  ymds2unix.c
  zulu2unix.c
  isleap.c
  leapsecond.c
  get_time.c
  jdate.c
  ccsds_to_yds.c
  yds2tai.c
  jday.c
  udunits2unix.cpp  
)

target_link_libraries(timeutils ${LIBS})
target_include_directories(timeutils PUBLIC "." PUBLIC ${Boost_INCLUDE_DIRS})

install (TARGETS timeutils DESTINATION lib)
