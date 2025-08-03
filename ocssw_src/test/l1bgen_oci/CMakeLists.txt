enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l1bgen_oci-data"
  VERBOSE
  TEST l1bgen_oci-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l1bgen_oci-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# PACE l1bgen tests
################################################################################

# first test case
set(base_file "PACE_OCI.20240404T141256")
set(l1b_file "${base_file}.L1B.nc")
set(par_file "${base_file}.L1B.par")
add_test(
  NAME "l1bgen_oci-${l1b_file}"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l1bgen_oci par=${par_file} && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${l1b_file} output/${l1b_file}"
)

# same test case using definftive ephem
set(base_file "PACE_OCI.20240404T141256")
set(l1b_file "${base_file}.L1B.eph.nc")
set(par_file "${base_file}.L1B.eph.par")
add_test(
  NAME "l1bgen_oci-${l1b_file}"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l1bgen_oci par=${par_file} && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${l1b_file} output/${l1b_file}"
)

# test case for solar cal data
set(base_file "PACE_OCI_SOL.20250206T065825")
set(l1b_file "${base_file}.L1B.nc")
set(par_file "${base_file}.L1B.par")
add_test(
  NAME "l1bgen_oci-${l1b_file}"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l1bgen_oci par=${par_file} && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${l1b_file} output/${l1b_file}"
)
