cmake_minimum_required(VERSION 3.12)
project(l1cextract)

find_program(PYTHON "python")

set (CMAKE_INSTALL_PREFIX $ENV{OCSSWROOT})
set (PACKAGE_VERSION "1.1.0")

if (PYTHON)
    set(TGT         "l1bextract_oci")
    set(OUTPUT      "${CMAKE_INSTALL_PREFIX}/bin/${TGT}")
    set(COPY_TO_PATH "${MANIFESTDISTRIB}/${TGT}")
    set(SCRIPTSDIR "${CMAKE_INSTALL_PREFIX}/src/scripts/seadasutils")
    set(PYREQUIRE   "${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt")
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(MODULES     "${CMAKE_CURRENT_BINARY_DIR}/seadasutils")
    set(PYVER       "--python $ENV{VIRTUAL_ENV}/bin/python")
    set(SHEBANG     "--python-shebang=/usr/bin/env\\ python3")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    file(MAKE_DIRECTORY ${MODULES})
    configure_file("${SCRIPTSDIR}/setupenv.py" ${MODULES} COPYONLY)
    configure_file("${SCRIPTSDIR}/netcdf_utils.py" ${MODULES} COPYONLY)
    configure_file("${SCRIPTSDIR}/pixlin_utils.py" ${MODULES} COPYONLY)
#    configure_file("${SCRIPTSDIR}/MetaUtils.py" ${MODULES} COPYONLY)
#    configure_file("${SCRIPTSDIR}/ProcUtils.py" ${MODULES} COPYONLY)
#    configure_file("${SCRIPTSDIR}/DictUtils.py" ${MODULES} COPYONLY)
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E touch "${MODULES}/__init.py__")

    install(CODE "execute_process(COMMAND pex ${CMAKE_CURRENT_BINARY_DIR} --no-index  -f ${CMAKE_INSTALL_PREFIX}/src/scripts/pex/wheel -v -r ${PYREQUIRE} -c l1bextract_oci.py -o ${OUTPUT} -D ${CMAKE_CURRENT_BINARY_DIR} --no-compile --disable-cache ${PYVER} ${SHEBANG})")

    install(CODE "MESSAGE(\"Copying ${TGT} binary to ${COPY_TO_PATH}\")")
    install(CODE "execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${OUTPUT}
${COPY_TO_PATH})")

endif()
