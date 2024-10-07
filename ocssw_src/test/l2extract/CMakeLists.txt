enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l2extract-data"
  VERBOSE
  TEST l2extract-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l2extract-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# PACE HKT tests
################################################################################

# normal extraction
set(base_file "PACE_OCI.20240324T011051.L2_OC")
set(l2_file "${base_file}_extract.nc")
set(par_file "${base_file}.par")
add_test(
  NAME "l2extract-${base_file}"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l2extract par=${par_file} && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${l2_file} output/${l2_file}"
)

# 3D extraction
set(base_file "PACE_OCI.20240324T011051.L2_OC_3D")
set(l2_file "${base_file}_extract.nc")
set(par_file "${base_file}.par")
add_test(
  NAME "l2extract-${base_file}"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l2extract par=${par_file} && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${l2_file} output/${l2_file}"
)
