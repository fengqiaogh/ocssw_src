enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l2brsgen-data"
  VERBOSE
  TEST l2brsgen-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l2brsgen-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

# get a list of the HDF par files in this directory
file(GLOB files "${dataDir}/*.hdf.par")

# loop through all HDF par files found and make tests
foreach(filename ${files})

  GET_FILENAME_COMPONENT(parfile ${filename} NAME)
  STRING(REGEX REPLACE ".par" "" outfile ${parfile})

  add_test(NAME "l2brsgen-${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
l2brsgen par=${parfile} && \
compare_hdf4 -G \"Processing Time,Processing Control\" ${outfile} output/${outfile}")

endforeach(filename)


# get a list of the png and ppm par files in this directory
file(GLOB files "${dataDir}/*.ppm.par")
file(GLOB files2 "${dataDir}/*.png.par")
list(APPEND files ${files2})

# loop through all HDF par files found and make tests
foreach(filename ${files})

  GET_FILENAME_COMPONENT(parfile ${filename} NAME)
  STRING(REGEX REPLACE ".par" "" outfile ${parfile})

  add_test(NAME "l2brsgen-${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
l2brsgen par=${parfile} && \
compare -metric AE -fuzz ${FUZZ} ${outfile} output/${outfile} /dev/null")

endforeach(filename)
