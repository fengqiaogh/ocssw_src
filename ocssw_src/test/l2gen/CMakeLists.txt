enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l2gen-data"
  VERBOSE
  TEST l2gen-data
  TARGET l2gen-data-t
  SAVE_DIR dataDir)

# refresh output directory
add_test(NAME "l2gen-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND bash -c "rm -rf output; mkdir -p output")

set(MISSIONS
     czcs
     octs
     seawifs
     aqua
     terra
     viirs
     hico
     olci
     oli
     ocm
     ocm2
)


################################################################################
# loop through l1info tests
################################################################################
foreach(mission ${MISSIONS})
  file(GLOB files "${dataDir}/${mission}/*.l1info")

  foreach(filename ${files})
    GET_FILENAME_COMPONENT(infofile ${filename} NAME)
    STRING(REGEX REPLACE ".l1info" "" l1file ${infofile})

    set(geofile "")
    if(${mission} MATCHES aqua|terra)
      STRING(REGEX REPLACE "L1B_LAC" "GEO" geofile ${l1file})
    endif(${mission} MATCHES aqua|terra)

    if(${mission} MATCHES viirs)
      STRING(REGEX REPLACE "SVM01" "GMTCO" geofile ${l1file})
    endif(${mission} MATCHES viirs)

    add_test(NAME "l1info-${l1file}"
      WORKING_DIRECTORY ${dataDir}/${mission}
      COMMAND bash -c "set -o xtrace; \
l1info -i 50 -o ../output/${infofile} ${l1file} ${geofile} && \
cd ${dataDir} && \
diff ${mission}/${infofile} output/${infofile}")

  endforeach(filename)
endforeach(mission)


################################################################################
# loop through l2gen tests
################################################################################
foreach(mission ${MISSIONS})
    file(GLOB files "${dataDir}/${mission}/*L2*.par")
    foreach(filename ${files})

      GET_FILENAME_COMPONENT(L2PAR ${filename} NAME)
      STRING(REGEX REPLACE ".par" ".nc" L2FILE ${L2PAR})

      add_test(NAME "l2gen-${L2FILE}"
             WORKING_DIRECTORY ${dataDir}/${mission}
             COMMAND bash -c "set -o xtrace; \
l2gen par=${L2PAR} && \
cd ${dataDir} && \
nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${mission}/${L2FILE} output/${L2FILE}")

    endforeach(filename)
endforeach(mission)
