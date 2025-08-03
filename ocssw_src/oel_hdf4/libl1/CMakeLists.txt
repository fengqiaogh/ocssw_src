cmake_minimum_required(VERSION 3.12)

#set_ocssw_policy()

set(LIBS "")
if (BUILD_HISTORICAL)
    set(historical_src "l1_czcs.c" "l1_seawifs.c" "l1_aci.c" "l1_seawifs_netcdf.cpp" "l1_czcs_netcdf.cpp"
    avhrrsub5h.f
    airmass_avhrr.f
    hmf8.f
    etbpsub.f
    intpos.f
    raygetpol.f
    simpsn.f)
    set(historical_libs czcs seawifs)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D BUILD_HISTORICAL")
endif(BUILD_HISTORICAL)

list(APPEND LIBS  piutils dfutils genutils ${historical_libs}  goci  meris oli )

find_package(GeographicLib REQUIRED)
include_directories(${GeographicLib_INCLUDE_DIRS})
list(APPEND LIBS ${GeographicLib_LIBRARIES})

find_package(GeoTIFF REQUIRED)
include_directories(${GeoTIFF_INCLUDE_DIRS})
list(APPEND LIBS ${GeoTIFF_LIBRARIES})

find_package(OpenJP2 REQUIRED)
include_directories(${OpenJP2_INCLUDE_DIRS})
list(APPEND LIBS ${OpenJP2_LIBRARIES})

find_package(NetCDF REQUIRED COMPONENTS CXX F77)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

list(APPEND LIBS timeutils)
list(APPEND LIBS h5io)
list(APPEND LIBS filetype)
list(APPEND LIBS bilio)

add_library(l1
    l1_options.c
    filehandle_init.c
    l1_io.c
    alloc_l1.c
    init_l1.c
    windex.c
    l1subpix.c
    l1_generic_write.c
#    l1_hdf_generic_write.c
    l1_hdf_generic_read.c
    l1_nc_generic_read.c
    l1_xcal.c
    windex.c
    brightness.c
    geom_per_band.c
    xcal.c
    get_f0.c
    setflags_l1.c
    cloud_flag.c
    scene_meta.c
    ll2vec.c
    ftrim.c
    l1_octs.c
    l1_octs_netcdf.cpp
    getglint.c
    l1_meris_N1.c
    l1_safe.cpp
    smile.c
    l1_meris_CC.c
    l1_aviris.c
    l1_aviris_nc.c
    l1_prism.c
    l1_viirs_h5.c
    l1_hico.c
    l1_goci.c
    l1_sgli.c
    l1_modis.c
    l1_viirs_l1b.c
    ${historical_src}
    l1_ocm.c
    l1_ocm2.c
    l1_ocmdb.c
    l1_oli.c
    l1_l5tm.c
    l1_l7etm.c
    l1_msi.cpp
    l1_hawkeye.c
    l1_misr.c
    l1_oci.c
    l1_seabass.cpp
    l1_l1c.c
    l1_l1c_anc.c
    l1_spexone.c
    uncertainty.c
)

target_link_libraries(l1 ${LIBS})
target_include_directories(l1 PUBLIC ".")

install (TARGETS l1 DESTINATION lib)
