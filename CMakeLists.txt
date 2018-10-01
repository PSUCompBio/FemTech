cmake_minimum_required(VERSION 3.0)

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()


# set install rules (windows only)
if(WIN32)
		add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

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
