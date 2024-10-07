enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l1cgen-data"
  VERBOSE
  TEST l1cgen-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l1cgen-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# OCI L1C tests
################################################################################

# Make the L1C grid
set(l1cfile "PACE.20240320T103854.L1C.nc")
set(listfile "l1cgen-grid.list")
add_test(
  NAME "l1cgen-grid"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1cgen par=../l1cgen-grid.par && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${l1cfile} output/${l1cfile} && \
diff ${listfile} output/${listfile}"
)

# make a full OCI L1C file
set(l1cfile "PACE_OCI.20240320T103854.L1C.nc")
add_test(
  NAME "l1cgen-full"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l1cgen par=l1cgen.par && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${l1cfile} output/${l1cfile}"
)

