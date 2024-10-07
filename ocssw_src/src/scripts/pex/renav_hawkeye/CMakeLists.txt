cmake_minimum_required(VERSION 3.12)
project(l1bgen_hico)

find_program(PYTHON "python")

set (CMAKE_INSTALL_PREFIX $ENV{OCSSWROOT})
set (PACKAGE_VERSION "1.0.0")

if (PYTHON)
    set(TGT         "renav_hawkeye")
    set(OUTPUT      "${CMAKE_INSTALL_PREFIX}/bin/${TGT}")
    set(COPY_TO_PATH "${MANIFESTDISTRIB}/${TGT}")
    set(HAWKNAVDIR  "${CMAKE_INSTALL_PREFIX}/src/nav_seahawk/hawknav")
    set(IGRFDIR     "${CMAKE_INSTALL_PREFIX}/src/nav_seahawk/pyIGRF")
    set(PYREQUIRE   "${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt")
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(HAWKNAV     "${CMAKE_CURRENT_BINARY_DIR}/hawknav")
    set(IGRF        "${CMAKE_CURRENT_BINARY_DIR}/pyIGRF")
    set(PYVER       "--python $ENV{VIRTUAL_ENV}/bin/python")
    set(SHEBANG     "--python-shebang=/usr/bin/env\\ python3")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    file(MAKE_DIRECTORY ${HAWKNAV})
    #file(MAKE_DIRECTORY ${IGRF})
    file(COPY "${HAWKNAVDIR}" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    #file(COPY "${IGRFDIR}" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    install(CODE "execute_process(COMMAND pex ${CMAKE_CURRENT_BINARY_DIR} --no-index  -f ${CMAKE_INSTALL_PREFIX}/src/scripts/pex/wheel -v -r ${PYREQUIRE} -c renav_hawkeye.py -o ${OUTPUT} -D ${CMAKE_CURRENT_BINARY_DIR} --no-compile --disable-cache ${PYVER} ${SHEBANG})")

    install(CODE "MESSAGE(\"Copying ${TGT} binary to ${COPY_TO_PATH}\")")
    install(CODE "execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${OUTPUT}
${COPY_TO_PATH})")

endif()
