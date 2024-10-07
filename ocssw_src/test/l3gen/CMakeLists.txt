enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l3gen-data"
  VERBOSE
  TEST l3gen-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l3gen-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

################################################################################
# chl creation - Aqua
################################################################################
add_test(NAME "l3gen_A2003060.L3b_DAY_CHL"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
l3gen par=A2003060.L3b_DAY_CHL.par && \
nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} \
A2003060.L3b_DAY_CHL.nc output/A2003060.L3b_DAY_CHL.nc")

################################################################################
# Hu chl creation - SeaWiFS
################################################################################
add_test(NAME "l3gen_S2002002.L3b_DAY_RRS_CHLHU"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
l3gen par=S2002002.L3b_DAY_RRS_CHLHU.par && \
nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} \
S2002002.L3b_DAY_RRS_CHLHU.nc output/S2002002.L3b_DAY_RRS_CHLHU.nc")

################################################################################
# (G)IOP creation - SeaWiFS
################################################################################
add_test(NAME "l3gen_S2002002.L3b_DAY_RRS_GIOP"
    WORKING_DIRECTORY ${dataDir}
    COMMAND bash -c "set -o xtrace; \
l3gen par=S2002002.L3b_DAY_RRS_GIOP.par && \
nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} \
S2002002.L3b_DAY_RRS_GIOP.nc output/S2002002.L3b_DAY_RRS_GIOP.nc")

################################################################################
# test out the Virtual Constellation stuff - netCDF4
################################################################################

## Generate the Rrs_vc product
#add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc.nc"
#    WORKING_DIRECTORY ${dataDir}
#    COMMAND l3gen par=A2003060.L3b_DAY_RRS_small_vc.nc.par)
#
## convert to OCRVC sensor
#add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc.nc-convert"
#    WORKING_DIRECTORY ${dataDir}
#    COMMAND convert2ocrvc.py output/A2003060.L3b_DAY_RRS_small_vc.L3b.nc)
#
#add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc.nc-check"
#    WORKING_DIRECTORY ${dataDir}
#    COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} A2003060.L3b_DAY_RRS_small_vc.L3b.nc output/A2003060.L3b_DAY_RRS_small_vc.L3b.nc )
#
#add_test(NAME "l3gen-A2003060.L3b_DAY_RRS_small_vc.L3b.nc-check2"
#    WORKING_DIRECTORY ${dataDir}
#    COMMAND compare_l3b.sh A2003060.L3b_DAY_RRS_small_vc.L3b.nc output/A2003060.L3b_DAY_RRS_small_vc.L3b.nc)
#
#
## Generate the VC product
#add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc2.nc"
#    WORKING_DIRECTORY ${dataDir}
#    COMMAND l3gen par=A2003060.L3b_DAY_RRS_small_vc2.nc.par)
#
#add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc2.nc-check"
#    WORKING_DIRECTORY ${dataDir}
#    COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} A2003060.L3b_DAY_RRS_small_vc2.L3b.nc output/A2003060.L3b_DAY_RRS_small_vc2.L3b.nc )
#
#add_test(NAME "l3gen-A2003060.L3b_DAY_RRS_small_vc2.nc-check2"
#    WORKING_DIRECTORY ${dataDir}
#    COMMAND compare_l3b.sh A2003060.L3b_DAY_RRS_small_vc2.L3b.nc output/A2003060.L3b_DAY_RRS_small_vc2.L3b.nc)

