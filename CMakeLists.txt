#To compile in debug mode                    : -DCMAKE_BUILD_TYPE=Debug
#Default behaviour is to use Release mode.

cmake_minimum_required(VERSION 3.3)

PROJECT (FemTech)

#============================================================================
# BUILD TYPE SPECIFIC PARTS
#==============================================================================
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_NO_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Build type (default Release)" FORCE)
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL DEBUG)
  message("Adding Debug flag...")
  add_definitions(-DDEBUG)
  # Add compile definition requires cmake 3.12+
  # add_compile_definitions(DEBUG)
endif()

#Setting c++11 parameters
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

message(STATUS "Configuring application: ${PROJECT_NAME}")
message(STATUS "Install Directory : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Current built type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Current build type flags: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
# Set additional compiler flags
STRING(TOLOWER "${PROC_ARCH}" PROC_ARCH)
if (NOT PROC_ARCH)
  message(STATUS "No processor instruction set selected, default to Native")
  set(PROC_ARCH "native" CACHE STRING "Instruction set (default native)" FORCE)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=${PROC_ARCH} -mtune=${PROC_ARCH}")

# set install rules (windows only)
if(WIN32)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

include_directories(include)

option(ENABLE_MPI "Enable MPI parallelization" ON)

#SET(PARALLEL 0 CACHE INTEGER "Defines a Parallel Build")

# Define a set of values to identify each MPI implementation supported (when enabled)
#set(MS_MPI "Microsoft")
#set(INTEL_MPI "Intel")
#set(DEFAULT_MPI "Default")  # Let's cmake automatically pick one

# Define the MPI to be used at each environment (when enabled ENABLE_MPI)
#set(SELECTED_MPI ${DEFAULT_MPI})
# At windows use MS one by default
#if(WIN32)
#	set(SELECTED_MPI ${MS_MPI})
#        message("Set MS_MPI") # debugging purposes
#endif(WIN32)

if(ENABLE_MPI)
#	set(MPI_IMPLEMENTATION ${SELECTED_MPI} CACHE STRING "MPI Implementation used by the software")
#	set_property(CACHE MPI_IMPLEMENTATION PROPERTY STRINGS ${MS_MPI} ${INTEL_MPI} ${DEFAULT_MPI})

#	message("MPI_IMPLEMENTATION: ${MPI_IMPLEMENTATION}") # debugging purposes

	find_package(MPI REQUIRED)

#	if("${MPI_IMPLEMENTATION}" STREQUAL "${MS_MPI}")
#		# Modify the compiler obtained by find_package with the user defined settings
#		set(MS_MPI_EXECUTABLE "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe")
#		set(MPIEXEC_EXECUTABLE ${MS_MPI_EXECUTABLE})
		#message("MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}") # debugging purposes
#	endif()

	include_directories(${MPI_INCLUDE_PATH})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
	add_subdirectory(third-party/parmetis-4.0.3)
	include_directories(third-party/parmetis-4.0.3/include)
	include_directories(third-party/parmetis-4.0.3/metis/include)
endif(ENABLE_MPI)

if(NOT WIN32)
	#Compiler and Linker parameters common to all compilers
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
endif()

#Get git related parameters
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating gitbranch.h")
configure_file(
  ${CMAKE_SOURCE_DIR}/include/gitbranch.h.in
  ${CMAKE_BINARY_DIR}/generated/gitbranch.h
)
include_directories(${CMAKE_BINARY_DIR}/generated)

#Add external packages and headers
add_subdirectory(third-party)

#Add JsonCPP Include and Link Directories
ExternalProject_Get_Property(JsonCpp source_dir)
include_directories(${source_dir}/include)

ExternalProject_Get_Property(JsonCpp binary_dir)
link_directories(${binary_dir}/src/lib_json)

#Add Boost Include Directories
ExternalProject_Get_Property(Boost source_dir)
include_directories(${source_dir})



#============================================================================
# BLAS/Lapack
#==============================================================================
if(NOT WIN32)
	find_package(BLAS REQUIRED)
	find_package(LAPACK REQUIRED)
	if(LAPACK_FOUND)
	  message("Found Lapack")
	  message("Lapack link library is ${LAPACK_LIBRARIES}")
	endif()
	if(BLAS_FOUND)
	  message("Found BLAS")
	  message("BLAS link library is ${BLAS_LIBRARIES}")
	endif()
	if(LAPACK_FOUND AND BLAS_FOUND)
	   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
	endif()
elseif(WIN32)
	add_library(Blas SHARED IMPORTED)
	set_target_properties(Blas PROPERTIES
		IMPORTED_LOCATION "C:/Users/rhk12/Documents/lapack-3.10.1/build/bin/libblas.dll"
		IMPORTED_IMPLIB "C:/Users/rhk12/Documents/lapack-3.10.1/build/lib/libblas.lib"
		INTERFACE_INCLUDE_DIRECTORIES "C:/Users/rhk12/Documents/lapack-3.10.1/build/include"
	)
	add_library(Lapack SHARED IMPORTED)
	set_target_properties(Lapack PROPERTIES
		IMPORTED_LOCATION "C:/Users/rhk12/Documents/lapack-3.10.1/build/bin/liblapack.dll"
		IMPORTED_IMPLIB "C:/Users/rhk12/Documents/lapack-3.10.1/build/lib/liblapack.lib"
		INTERFACE_INCLUDE_DIRECTORIES "C:/Users/rhk12/Documents/lapack-3.10.1/build/include"
	)
endif()


ADD_SUBDIRECTORY (src)

#===============================================================================
# Examples and tests
#===============================================================================
option(EXAMPLES "Activate examples" OFF)
if(EXAMPLES)
  add_subdirectory(examples)
endif()
