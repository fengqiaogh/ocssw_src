cmake_minimum_required(VERSION 3.12)
find_program(PYTHON "python")

set (CMAKE_INSTALL_PREFIX $ENV{OCSSWROOT})
set (PACKAGE_VERSION "1.1")

if (PYTHON)
    set(TGT         "gmao_aot_subset")
    set(OUTPUT      "${CMAKE_INSTALL_PREFIX}/bin/${TGT}")
    set(COPY_TO_PATH "${MANIFESTDISTRIB}/${TGT}")
    set(PYREQUIRE   "${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt")
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(PYVER       "--python $ENV{VIRTUAL_ENV}/bin/python")
    set(SHEBANG     "--python-shebang=/usr/bin/env\\ python3.8")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    install(CODE "execute_process(COMMAND pex ${CMAKE_CURRENT_BINARY_DIR} --no-index  -f ${CMAKE_INSTALL_PREFIX}/src/scripts/pex/wheel -v -r ${PYREQUIRE}  -o ${OUTPUT} -c gmao_aot_subset --no-compile --disable-cache ${PYVER} ${SHEBANG})")

    install(CODE "MESSAGE(\"Copying ${TGT} binary to ${COPY_TO_PATH}\")")
    install(CODE "execute_process( COMMAND ${CMAKE_COMMAND} -E copy ${OUTPUT}
${COPY_TO_PATH})")

endif()
