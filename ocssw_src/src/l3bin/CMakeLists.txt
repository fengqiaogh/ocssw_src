cmake_minimum_required(VERSION 3.12)

# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -fsanitize=leak")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -fsanitize=leak")
add_executable(l3bin
  l3bin.cpp
  l3bin_input.cpp
)

add_executable(l3binmerge
  l3bin_input.cpp
  l3binmerge.cpp
)

target_link_libraries(l3bin
  bin
)

target_link_libraries(l3binmerge
  bin
)

install (TARGETS l3bin l3binmerge DESTINATION bin)
