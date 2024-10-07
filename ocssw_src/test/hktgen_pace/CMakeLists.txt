enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/hktgen_pace-data"
  VERBOSE
  TEST hktgen_pace-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "hktgen_pace-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# PACE HKT tests
################################################################################

# Crossing the dateline
set(hskfile "HSK000011063_20240320T112501.HSK")
set(hktfile "PACE.20240320T095839.HKT.nc")
add_test(
  NAME "hktgen_pace-${hktfile}"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
hktgen_pace --pversion my_version --doi my_doi -o output/${hktfile} ${hskfile} && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${hktfile} output/${hktfile}"
)
