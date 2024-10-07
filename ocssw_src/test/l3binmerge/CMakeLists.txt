enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l3binmerge-data"
  VERBOSE
  TEST l3binmerge-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l3binmerge-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

set(filename "l3binmerge_Rrs_443_555")
set(outfile "${filename}.nc")
set(parfile "${filename}.par")
add_test(NAME "l3binmerge-${filename}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
l3binmerge par=${parfile} && \
nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} ${outfile} output/${outfile}")
