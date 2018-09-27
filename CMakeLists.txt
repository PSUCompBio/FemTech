cmake_minimum_required(VERSION 3.0)

PROJECT (digitalbrain)
include_directories(include)
ADD_SUBDIRECTORY (src)

#===============================================================================
# Examples and tests
#===============================================================================
option(EXAMPLES "Activate examples" OFF)
if(EXAMPLES)
  add_subdirectory(examples)
endif()
