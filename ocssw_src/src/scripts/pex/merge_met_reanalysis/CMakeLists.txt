cmake_minimum_required(VERSION 3.12)
find_program(PYTHON "python")

set (CMAKE_INSTALL_PREFIX $ENV{OCSSWROOT})
set (PACKAGE_VERSION "1.1")

#set (CMAKE_PREFIX_PATH $ENV{LIB3_DIR})
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/oel_hdf4/cmake/Modules)
#message("MODULE PATH: ${CMAKE_MODULE_PATH}")

#find_package(HDF4 REQUIRED)
#include_directories(${HDF4_INCLUDE_DIRS})
#list(APPEND LIBS ${HDF4_LIBRARIES})

if (PYTHON)
    set(TGT         "merge_met_reanalysis")
    set(OUTPUT      "${CMAKE_INSTALL_PREFIX}/bin/${TGT}")
    set(COPY_TO_PATH "${MANIFESTDISTRIB}/${TGT}")
    set(PYREQUIRE   "${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt")
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(PYVER       "--python $ENV{VIRTUAL_ENV}/bin/python")
    set(SHEBANG     "--python-shebang=/usr/bin/env\\ python3.8")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    message("${PYREQUIRE}")
#add_custom_command(OUTPUT ALL
#    COMMAND ${CMAKE_COMMAND} -E env "INCLUDE_DIRS=${HDF4_INCLUDE_DIRS}" 
#    COMMAND ${CMAKE_COMMAND} -E env "LIBRARY_DIRS=${HDF4_LIBRARIES}"
#    COMMAND  pex ${CMAKE_CURRENT_BINARY_DIR} -v -r ${PYREQUIRE}  -o ${OUTPUT} -c merge_met_reanalysis.py --no-compile --disable-cache --inherit-path
#)

    install(CODE "execute_process(COMMAND pex ${CMAKE_CURRENT_BINARY_DIR} --no-index  -f ${CMAKE_INSTALL_PREFIX}/src/scripts/pex/wheel -v -r ${PYREQUIRE}  -o ${OUTPUT} -c merge_met_reanalysis.py --no-compile --disable-cache --inherit-path ${PYVER} ${SHEBANG})")

    install(CODE "MESSAGE(\"\nIf this failed, you need to:\nexport LIBRARY_DIRS=${CMAKE_INSTALL_PREFIX}/opt/lib\nexport INCLUDE_DIRS=${CMAKE_INSTALL_PREFIX}/opt/include\nand rerun: make install\")")

    install(CODE "MESSAGE(\"Copying ${TGT} binary to ${COPY_TO_PATH}\")")
    install(CODE "execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${OUTPUT}
${COPY_TO_PATH})")

endif()
