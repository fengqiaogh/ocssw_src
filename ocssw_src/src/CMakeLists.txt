cmake_minimum_required(VERSION 3.12)

# NOTE, the order matters... when adding new code, put elements that depend on
# others AFTER the dependency! (e.g. libbin needs libhdfutils, so libhdfutils
# goes first...)


# libraries
add_subdirectory(libanc)

add_subdirectory(libl2)
add_subdirectory(libdataday)
add_subdirectory(libimageutils)
add_subdirectory(libterrain)

# Common progs
add_subdirectory(l2gen)
add_subdirectory(l2bin)
add_subdirectory(l3bin)
add_subdirectory(l1aextract_modis)
add_subdirectory(l1bextract_meris)
add_subdirectory(l2extract)
add_subdirectory(l2mapgen)
add_subdirectory(l3bindump)
add_subdirectory(get_product_info)
add_subdirectory(l3mapgen)
add_subdirectory(nccmp)
add_subdirectory(val_extract)
add_subdirectory(matchup_tools)
add_subdirectory(l1agen_hawkeye)
add_subdirectory(l1agen_oci)
add_subdirectory(l1bgen_oci)
add_subdirectory(manifest)
add_subdirectory(l1cgen)
add_subdirectory(ancgen)
add_subdirectory(geo_eval)
# MODIS-specific
add_subdirectory(l1agen_modis)
add_subdirectory(l1bgen_modist)
add_subdirectory(l1bgen_modisa)
add_subdirectory(geogen_modis)
add_subdirectory(l0fix_modis)
add_subdirectory(l0cnst_write_modis)
add_subdirectory(l0chunk_modis)
add_subdirectory(l1asubset_modis)
add_subdirectory(atteph_info_modis)
add_subdirectory(geocheck_modis)
add_subdirectory(pdsmerge)
if (BUILD_MISC)
    add_subdirectory(l0cnst_read_modis)
endif ()


add_subdirectory(scripts)
add_subdirectory(l1aqc_oci)

add_subdirectory(aviris2oci)
add_subdirectory(aviris2ncdf)
add_subdirectory(prism2oci)
add_subdirectory(prismbil2oci)
add_subdirectory(l1bgen_hico)
add_subdirectory(metadata_gibs)

# VIIRS QC programs
add_subdirectory(l0qc_viirs)
add_subdirectory(l2qc_viirs)

if (BUILD_MISC)
  # Ancillary-specific
  add_subdirectory(met_reverse)
  add_subdirectory(auto_qc)
  add_subdirectory(oisst_info)
  add_subdirectory(ancnrt)
  add_subdirectory(gmao)
  #add_subdirectory(ice2hdf)
  #  Miscellaneous
  add_subdirectory(landbin)
  add_subdirectory(bindepths)
  #  add_subdirectory(biosmap)
  #  add_subdirectory(l3despeckle)
  add_subdirectory(get_utm_zone)
  #add_subdirectory(darktargetc)
  add_subdirectory(anctools)
  add_subdirectory(l3bstats)
  add_subdirectory(get_dataday)
endif(BUILD_MISC)

if (BUILD_HISTORICAL)
  add_subdirectory(l1agen_czcs)
  add_subdirectory(l1agen_seawifs)
  add_subdirectory(orbnum_seawifs)
  add_subdirectory(gpsextract_seawifs)
  add_subdirectory(neb2raw_seawifs)
  add_subdirectory(raw2l0_seawifs)
  add_subdirectory(ql02raw_seawifs)
  add_subdirectory(tlmgen_seawifs)
  add_subdirectory(l0regen_seawifs)
  add_subdirectory(viirs_sim_sdr)
  add_subdirectory(modis_timestamp)
  add_subdirectory(l1aextract_seawifs)
  add_subdirectory(libswfnav)
  # These are commented out since they are the only programs that use
  # the HDF4 Fortran API.  To build them you need to enable Fortran
  # when building HDF4
  #
  #add_subdirectory(interpnav_seawifs)
  #add_subdirectory(navassess_seawifs)

endif(BUILD_HISTORICAL)

if (BUILD_CLOUD)
  add_subdirectory(libcloud)
endif(BUILD_CLOUD)
