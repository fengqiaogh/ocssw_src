enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l1agen_oci-data"
  VERBOSE
  TEST l1agen_oci-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l1agen_oci-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# OCI L1A tests
################################################################################

# Crossing the dateline
set(datetime "20240430T012933")
set(l1file "PACE_OCI.${datetime}.L1A.nc")
set(ifile "OCI000037515.oci")
add_test(
  NAME "l1agen_oci-cross_dateline"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci ifile=../cross/${ifile} granule_len=1 -t=${datetime} -k=../cross/hktlist.txt -o=l1agen_oci-cross_out.txt && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} cross/${l1file} output/${l1file}"
)

# First granule in an orbit
# also test pversion
set(datetime "20240430T074140")
set(l1file "PACE_OCI.${datetime}.L1A.nc")
set(ifile "OCI000037615.oci")
add_test(
  NAME "l1agen_oci-first_granule"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci ifile=../firstGranule/${ifile} granule_len=1 -t=${datetime} -k=../firstGranule/hktlist.txt -o=l1agen_oci-first_out.txt -pversion=CoolVersion && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} firstGranule/${l1file} output/${l1file}"
)

# Last granule in an orbit
set(datetime "20240430T082957")
set(l1file "PACE_OCI.${datetime}.L1A.nc")
set(ifile "OCI000037642.oci")
add_test(
  NAME "l1agen_oci-last_granule"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci ifile=../lastGranule/${ifile} granule_len=1 -t=${datetime} -k=../lastGranule/hktlist.txt -o=l1agen_oci-last_out.txt && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} lastGranule/${l1file} output/${l1file}"
)

# Lunar calibration
set(datetime "20240424T132451")
set(l1file "PACE_OCI_LUN.${datetime}.L1A.nc")
set(ifile "OCI000035266.oci")
add_test(
  NAME "l1agen_oci-lun"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci ifile=../lunar/${ifile} granule_len=1 -t=${datetime} -k=../lunar/hktlist.txt -o=l1agen_oci-lun_out.txt && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} lunar/${l1file} output/${l1file}"
)

# Daily solar calibration
set(datetime "20240501T072601")
set(l1file "PACE_OCI_SOL.${datetime}.L1A.nc")
set(ifile "OCI000038032.oci")
add_test(
  NAME "l1agen_oci-sol"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci ifile=../solar/${ifile} granule_len=1 -t=${datetime} -k=../solar/hktlist.txt -o=l1agen_oci-sol_out.txt && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} solar/${l1file} output/${l1file}"
)

# Earth spectral
set(datetime "20240329T005320")
set(l1file "PACE_OCI_SPEC.${datetime}.L1A.nc")
set(ifile "OCI000024932.oci")
add_test(
  NAME "l1agen_oci-spec"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci ifile=../spectral/${ifile} granule_len=1 -t=${datetime} -k=../spectral/hktlist.txt -o=l1agen_oci-spec_out.txt && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} spectral/${l1file} output/${l1file}"
)

# Response curve (LIN)
set(datetime "20240420T041922")
set(l1file "PACE_OCI_LIN.${datetime}.L1A.nc")
set(ifile "ocilist.txt")
add_test(
  NAME "l1agen_oci-lin"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci ifile=../responseCurve/${ifile} granule_len=1 -t=${datetime} -k=../responseCurve/hktlist.txt -o=l1agen_oci-lin_out.txt && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} responseCurve/${l1file} output/${l1file}"
)

# Monthly solar calibration (SPCA)
set(datetime "20240509T035619")
set(l1file "PACE_OCI_SPCA.${datetime}.L1A.nc")
set(ifile "ocilist.txt")
add_test(
  NAME "l1agen_oci-spca"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci ifile=../monthlySolarCal/${ifile} granule_len=1 -t=${datetime} -k=../monthlySolarCal/hktlist.txt -o=l1agen_oci-lin_out.txt && \
cd .. && \
nccmp -mgdf -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} monthlySolarCal/${l1file} output/${l1file}"
)
