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

# global attributes to ignore because it will change based on how the files are set up or ran
# ignore these to avoid failure when comparing data
set(ignoreGlobalAtts "date_created,software_version")


# Full Granule
# One full run without any file splitting or other data types in between.
set(l1aFile "PACE_OCI.20250909T080250.L1A.nc")
add_test(
  NAME "l1agen_oci-full-granule"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci par=../fullGranule/par.par && cd .. && \
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} fullGranule/${l1aFile} output/${l1aFile}"
)


# HKT Day Rollover
# HKT time < Requested Time and it crosses to a new day
set(l1aFile "PACE_OCI.20250809T000034.L1A.nc")
add_test(
  NAME "l1agen_oci-hkt-day-rollover"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci par=../hktDayRollover/par.par && cd .. && \
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} hktDayRollover/${l1aFile} output/${l1aFile}"
)

# Negative Navigation
# Time is < 10 seconds into a new day. The code has a 10 second padding between
# the start and end times, so time will contain ~6 seconds from the previous day
# which will make the navigation time negative
set(l1aFile "PACE_OCI.20250712T000004.L1A.nc")
add_test(
  NAME "l1agen_oci-negative-navigation-time"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci par=../negativeNavigation/par.par && cd .. && \
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} negativeNavigation/${l1aFile} output/${l1aFile}"
)

# File Split
# _SPEC file is in-between science data. Second science data is merged into the first one.
# Tests **BOTH** science file and _SPEC file
set(sciFile "PACE_OCI.20250904T201240.L1A.nc")          # expected science file name
set(specFile "PACE_OCI_SPEC.20250904T201622.L1A.nc")    # expected spec file name
add_test(
  NAME "l1agen_oci-file-split"
  WORKING_DIRECTORY ${dataDir}/output
  COMMAND bash -c "set -o xtrace; \
l1agen_oci par=../fileSplit/par.par && cd .. && \
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} fileSplit/${sciFile} output/${sciFile} && \
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} fileSplit/${specFile} output/${specFile}"
)


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
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} cross/${l1file} output/${l1file}"
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
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} lunar/${l1file} output/${l1file}"
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
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} solar/${l1file} output/${l1file}"
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
nccmp -mgdf -C 10 -G ${ignoreGlobalAtts} $ENV{OCTEST_TOLERANCE} responseCurve/${l1file} output/${l1file}"
)
