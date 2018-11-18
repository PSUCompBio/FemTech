cmake_minimum_required(VERSION 3.9)

PROJECT (digitalbrain)

# set install rules (windows only)
if(WIN32)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

include_directories(include)

option(ENABLE_MPI "Enable MPI parallelization" OFF)

SET(PARALLEL 0 CACHE INTEGER "Defines a Parallel Build")

if(ENABLE_MPI)
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})
	set(CMAKE_C_FLAGS "${CMAKE_FLAGS} ${MPI_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")
	add_subdirectory(third-party/parmetis-4.0.3)
	include_directories(third-party/parmetis-4.0.3/include)
	include_directories(third-party/parmetis-4.0.3/metis/include)
	unset(PARALLEL CACHE)
	SET(PARALLEL 1 CACHE INTEGER "Defines a Parallel Build")
	add_definitions( -DPARALLEL=${PARALLEL} )
endif(ENABLE_MPI)

ADD_SUBDIRECTORY (src)

#===============================================================================
# Examples and tests
#===============================================================================
option(EXAMPLES "Activate examples" OFF)
if(EXAMPLES)
  add_subdirectory(examples)
endif()
