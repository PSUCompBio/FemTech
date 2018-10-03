cmake_minimum_required(VERSION 3.0)

PROJECT (digitalbrain)

# set install rules (windows only)
if(WIN32)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

include_directories(include)

ADD_SUBDIRECTORY (src)

option(MPI "Parallel MPI" OFF)
if(MPI)
	if(WIN32)
			find_package(MPI REQUIRED)
			if(MPI_FOUND)
	      message("Using MPI")
				
				include_directories(${MPI_INCLUDE_PATH})

				add_definitions(-DMPI_C_COMPILER=C:/Program Files/Microsoft MPI/Bin/mpicc -DMPI_CXX_COMPILER=C:/Program Files/Microsoft MPI/Bin/mpicxx)

				if(MPI_COMPILE_FLAGS)
		  			set_target_properties(digitalbrain PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
				endif(MPI_COMPILE_FLAGS)

				if(MPI_LINK_FLAGS)
		  			set_target_properties(digitalbrain PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
				endif(MPI_LINK_FLAGS)

				#	FIND_PATH(MPI_M_C mpicc $ENV{PATH})
				#	FIND_PATH(MPI_INTEL_CXX mpicxx $ENV{PATH})
				#	MESSAGE(STATUS "MPI_INTEL_C: ${MPI_BIN}")
				#	MESSAGE(STATUS "MPI_INTEL_CXX: ${MPI_INTEL_CXX}")

			endif(MPI_FOUND)
   	endif(WIN32)
endif(MPI)

#===============================================================================
# Examples and tests
#===============================================================================
option(EXAMPLES "Activate examples" OFF)
if(EXAMPLES)
  add_subdirectory(examples)
endif()
