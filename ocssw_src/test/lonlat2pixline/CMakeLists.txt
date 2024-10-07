enable_testing()

#include(PrintAllVariables)
#print_all_variables()

# download lonlat2pixline's data
cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/lonlat2pixline-data"
  VERBOSE
  TEST lonlat2pixline-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "lonlat2pixline-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# SeaWiFS find point
################################################################################
set(infile S2002079134718.L1A_GAC.middle)
set(outfile ${infile}_point.txt)
add_test(NAME "lonlat2pixline_${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -o output/${outfile} ${infile} -34 -8 && \
diff -u ${outfile} output/${outfile}")

################################################################################
# SeaWiFS find box
################################################################################
set(outfile ${infile}_box.txt)
add_test(NAME "lonlat2pixline_${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -o output/${outfile} ${infile} -34 -13 -30 -8 && \
diff -u ${outfile} output/${outfile}")

################################################################################
# SeaWiFS find point failure, exit=100
################################################################################
add_test(NAME "lonlat2pixline_${infile}_point-fail"
    WORKING_DIRECTORY ${dataDir}
    COMMAND test_exit_code 100 lonlat2pixline ${infile} 0 0 )

################################################################################
# SeaWiFS find box failure, exit=1
################################################################################
add_test(NAME "lonlat2pixline_${infile}_box-fail"
    WORKING_DIRECTORY ${dataDir}
    COMMAND test_exit_code 100 lonlat2pixline ${infile} 0 0 2 2 )

################################################################################
# SeaWiFS partial box extracted, exit=0
################################################################################
set(outfile ${infile}_partial.txt)
add_test(NAME "lonlat2pixline_${outfile}"
  WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -o output/${outfile} ${infile} -34 -8 50 50 && \
diff -u ${outfile} output/${outfile}")

################################################################################
# SeaWiFS full box not extracted, exit=110
################################################################################
set(outfile ${infile}_full.txt)
add_test(NAME "lonlat2pixline_${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -F -o output/${outfile} ${infile} -34 -8 50 50 ; \
(test $? -eq 110 || exit -1) && \
diff -u ${outfile} output/${outfile}")

################################################################################
# SeaWiFS whole file is inside box, exit=120
################################################################################
set(outfile ${infile}_whole.txt)
add_test(NAME "lonlat2pixline_${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -F -o output/${outfile} ${infile} -50 -50 50 50 ; \
(test $? -eq 120 || exit -1) && \
diff -u ${outfile} output/${outfile}")

################################################################################
# Aqua find point
################################################################################
set(infile A2008080195500.L1B_LAC.subline)
set(geofile A2008080195500.GEO.subline)
set(outfile ${infile}_point.txt)
add_test(NAME "lonlat2pixline_${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -o output/${outfile} ${infile} ${geofile} -92 -47 && \
diff -u ${outfile} output/${outfile}")

################################################################################
# Aqua find box
################################################################################
set(outfile ${infile}_box.txt)
add_test(NAME "lonlat2pixline_${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -o output/${outfile} ${infile} ${geofile} -92 -47 -91 -46 && \
diff -u ${outfile} output/${outfile}")

################################################################################
# Meris find point
################################################################################
set(infile MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix)
set(outfile ${infile}_point.txt)
add_test(NAME "lonlat2pixline_${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -o output/${outfile} ${infile} -40 42 && \
diff -u ${outfile} output/${outfile}")

################################################################################
# Meris find box
################################################################################
set(outfile ${infile}_box.txt)
add_test(NAME "lonlat2pixline_${outfile}"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
lonlat2pixline -o output/${outfile} ${infile} -40 41.2 -36.5 43 && \
diff -u ${outfile} output/${outfile}")
