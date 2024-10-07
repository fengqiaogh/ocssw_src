enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/ancgen-data"
  VERBOSE
  TEST ancgen-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "ancgen-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# PACE HKT tests
################################################################################

# Crossing the dateline
set(base_file "PACE.20240515T003705.L1C.ANC")
set(anc_file "${base_file}.nc")
set(par_file "${base_file}.par")
add_test(
  NAME "ancgen-${anc_file}"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
ancgen par=${par_file} && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${anc_file} output/${anc_file}"
)
