enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l3mapgen-data"
  VERBOSE
  TEST l3mapgen-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l3mapgen-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

# hdf.input.nc.par - HDF4 input, netCDF output, one product
# multiple.nc.par - multiple output products, area interp, smi projection
# out.hdf.par - netCDF input, HDF4 output
# out.nc.par - netCDF output
# out.png.par - PNG output
# out.ppm.par - PPM output
# out.tiff.par - GEO TIFF output
# projection.nc.par - mollweide projection
# resolution.nc.par - 90km resolution pixel
# rgb.png.par - make a RGB true color image
#


# get a list of the netCDF par files in this directory
file(GLOB files "${dataDir}/*.nc.par")

# loop through all netCDF par files found and make tests
foreach(filename ${files})

  GET_FILENAME_COMPONENT(parfile ${filename} NAME)
  STRING(REGEX REPLACE ".par" "" outfile ${parfile})

  add_test(NAME "l3mapgen-${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
l3mapgen par=${parfile} && \
nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} ${outfile} output/${outfile}")

endforeach(filename)


# test for HDF4 SMI output
add_test(NAME "l3mapgen-hdf4"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l3mapgen par=out.hdf.par && \
compare_hdf4 -G \"Software Version,Processing Time,Processing Control,_lastModified\" out.hdf output/out.hdf")


# test PNG output
add_test(NAME "l3mapgen-png"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l3mapgen par=out.png.par && \
compare -metric AE -fuzz ${FUZZ} out.png output/out.png /dev/null")


# test PNG transparent missing data output
add_test(NAME "l3mapgen-transparent-png"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l3mapgen par=out.transparent.png.par && \
compare -metric AE -fuzz ${FUZZ} out.transparent.png output/out.transparent.png /dev/null")


# test PPM output
add_test(NAME "l3mapgen-ppm"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l3mapgen par=out.ppm.par && \
compare -metric AE -fuzz ${FUZZ} out.ppm output/out.ppm /dev/null")


# need to add a test for the geotiff metadata tags
#
add_test(NAME "l3mapgen-tiff"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l3mapgen par=out.tiff.par && \
tiffcmp out.tiff output/out.tiff")

# test geotiff transparent missing data
add_test(NAME "l3mapgen-transparent-tiff"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l3mapgen par=out.transparent.tiff.par && \
tiffcmp out.transparent.tiff output/out.transparent.tiff")

# test RGB PNG output
add_test(NAME "l3mapgen-rgb"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l3mapgen par=rgb.png.par && \
compare -metric AE -fuzz ${FUZZ} rgb.png output/rgb.png /dev/null")

# add a test for multiple tiff files
add_test(NAME "l3mapgen-multiple-tiff"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l3mapgen par=multiple_PRODUCT.tiff.par && \
tiffcmp multiple_Rrs_443.tiff output/multiple_Rrs_443.tiff && \
tiffcmp multiple_Rrs_555.tiff output/multiple_Rrs_555.tiff")
