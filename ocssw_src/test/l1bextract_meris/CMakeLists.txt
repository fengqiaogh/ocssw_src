enable_testing()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l1bextract_meris-data"
  VERBOSE
  TEST l1bextract_meris-data
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l1bextract_meris-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

# test the SAFE extractor
add_test(NAME "l1bextract_meris_safe"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "set -o xtrace; \
l1bextract_safe_nc --spixl 300 --epixl 900 --sline 6000 --eline 6200 \
ENV_ME_1_RRG____20061115T162732_20061115T171127_________________2636_053_026______DSI_R_NT____.SEN3 \
output/safe.extract && \
cd ${dataDir}/output/safe.extract && \
l2gen atmocor=off l2prod=Lt_490 ifile=xfdumanifest.xml ofile=meris.safe.extract.L2.nc && \
cd ${dataDir} && \
nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} \
meris.safe.extract.L2.nc output/safe.extract/meris.safe.extract.L2.nc")
