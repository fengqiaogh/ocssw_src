
cmake_minimum_required(VERSION 3.12)

set(LIBS genutils piutils)

add_library (gring
  gring.cpp
)

add_executable ( test_gring
  test.cpp
)

target_link_libraries(test_gring gring ${LIBS})

target_link_libraries(gring ${LIBS})
target_include_directories(gring PUBLIC ".")

install (TARGETS gring DESTINATION lib) 