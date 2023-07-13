SET(ex6_src ex6.cpp)

ADD_EXECUTABLE(ex6 ${ex6_src})

if(ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex6 FemTech parmetis ${lapackblas_libraries})
elseif(NOT ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex6 FemTech)
elseif(ENABLE_MPI AND WIN32)
	LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
	target_link_libraries(ex6 FemTech parmetis Blas Lapack ${MPI_LIBRARIES})
	# copy dlls to build directory
	file(COPY ${LAPACK}/bin/libcblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex6)
	file(COPY ${LAPACK}/bin/libblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex6)
	file(COPY ${LAPACK}/bin/liblapack.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex6)
	file(COPY ${LAPACK}/bin/liblapacke.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex6)
endif()

file(COPY 10e.k DESTINATION .)
file(COPY materials.dat DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 2 \"$(TargetPath)\" 10e.k" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex6/ex6.vcxproj.user @ONLY )
endif()
