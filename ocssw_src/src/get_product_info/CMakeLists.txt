cmake_minimum_required(VERSION 3.12)

add_executable(get_product_info
  get_product_info.c
)

target_link_libraries(get_product_info
  l2
)

install (TARGETS get_product_info DESTINATION bin)
