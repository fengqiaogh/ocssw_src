enable_testing()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l1agen_seawifs-data"
  VERBOSE
  TEST l1agen_seawifs-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l1agen_seawifs-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# loop through L0 files
################################################################################
file(GLOB files "${dataDir}/*.L0_*")
foreach(filename ${files})

  GET_FILENAME_COMPONENT(l0file ${filename} NAME)
  STRING(REGEX REPLACE "L0_" "L1A_" l1file ${l0file})

  add_test(NAME "l1agen_seawifs-${l1file}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
l1agen_seawifs -n sdps ${l0file} output/${l1file} && \
compare_hdf4 -G \"Processing Time,Processing Control\" ${l1file} output/${l1file}")

endforeach(filename)
