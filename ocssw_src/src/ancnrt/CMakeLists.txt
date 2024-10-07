cmake_minimum_required(VERSION 3.12)

#set_ocssw_policy()

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})

find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIRS})

add_executable(ancnrt
  ancnrt.c
  readgrib.c
  fillenv.c
  ANCroutines.c
  rdattr.c
  countann.c
  rdsdsid.c
  pexit.c
  resize_oz.c
)

add_executable(ancnrt_2p5
  ancnrt_2p5.c
  readgrib.c
  fillenv.c
  ANCroutines.c
  rdattr.c
  countann.c
  rdfiles.f
  rdsdsid.c
  pexit.c
  julian.f
  lnstrg.f
  leap.f
  resize_oz.c
)

add_executable(o3nrt_omi
  o3nrt_omi.c
  ANCroutines.c
  rdgrid.f
  countann.c
  fillenv.c
  pexit.c
  gregor.f
  leap.f
  julian.f
  getfn.f
  lnstrg.f
  epochbk.c
  julday.f
  rdattr.c
  rdfiles.f
  rdsdsid.c
  world_avg.c
)

add_executable(o3_toms
  o3_toms.c
  ANCroutines.c
  rd_toms_ascii.c
  countann.c
  fillenv.c
  pexit.c
  epochbk.c
  rdattr.c
  rdsdsid.c
  fill_smooth.c
)

if(BUILD_AQUARIUS)
  add_executable(anc_seaice
    anc_seaice.c
    readgrib.c
    fillenv.c
    ANCroutines.c
    rdattr.c
    countann.c
    rdfiles.f
    rdsdsid.c
    pexit.c
    julian.f
    lnstrg.f
    leap.f
    resize_oz.c
    mk_smooth_ice_map.f
  )

  add_executable(ancnrt_aquarius
    ancnrt_aquarius.c
    readgrib.c
    fillenv.c
    ANCroutines.c
    rdattr.c
    countann.c
    rdfiles.f
    rdsdsid.c
    pexit.c
    julian.f
    lnstrg.f
    leap.f
    resize_oz.c
  )

  add_executable(solar_xray
   solar_xray.c
  )

  target_link_libraries(ancnrt_aquarius
    timeutils
    genutils
    ${HDF4_LIBRARIES}
  )

  target_link_libraries(anc_seaice
    aquarius
    genutils
    ${HDF4_LIBRARIES}
  )

  target_link_libraries(solar_xray
    genutils
    ${HDF5_LIBRARIES}
  )

  if(BUILD_JPL) 

    add_executable(interp_hycom
     interp_hycom.f
    )

    add_executable(interp_hycom_ascii
     interp_hycom_ascii.f
    )

    add_executable(tec
      tec.c
    )

    add_executable(swh
      swh.c
    )

    target_link_libraries(interp_hycom
      netcdff
      dfutils
      ${HDF5_LIBRARIES}
    )

    target_link_libraries(interp_hycom_ascii
      dfutils
      netcdff
      z
    )

    target_link_libraries(tec
      ${HDF5_LIBRARIES}
    )

    target_link_libraries(swh
      ${HDF5_LIBRARIES}
    )

  endif(BUILD_JPL)

endif(BUILD_AQUARIUS)

target_link_libraries(ancnrt
  timeutils
  genutils
  ${HDF4_LIBRARIES}
)

target_link_libraries(ancnrt_2p5
  timeutils
  genutils
  ${HDF4_LIBRARIES}
)

target_link_libraries(o3nrt_omi
  genutils
  ${HDF4_LIBRARIES}
)

target_link_libraries(o3_toms
  genutils
  h5io
  ${HDF4_LIBRARIES}
  ${HDF5_LIBRARIES}
)

if(BUILD_AQUARIUS) 
  install (TARGETS ancnrt_aquarius anc_seaice solar_xray DESTINATION bin)
  if (BUILD_JPL)
    install (TARGETS tec interp_hycom interp_hycom_ascii swh DESTINATION bin)
  endif(BUILD_JPL)
endif(BUILD_AQUARIUS)

install (TARGETS ancnrt ancnrt_2p5 o3nrt_omi o3_toms DESTINATION bin)
