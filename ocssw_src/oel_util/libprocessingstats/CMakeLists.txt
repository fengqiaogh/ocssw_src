
cmake_minimum_required(VERSION 3.12)

add_library (processing_tracker
  processing_tracker.cpp
)

target_include_directories(processing_tracker PUBLIC ".")

install (TARGETS processing_tracker DESTINATION lib) 