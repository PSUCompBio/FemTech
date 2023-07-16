SET(ex15_src ex15.cpp)

ADD_EXECUTABLE(ex15 ${ex15_src})

if(ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex15 FemTech parmetis ${lapackblas_libraries})
elseif(NOT ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex15 FemTech)
elseif(ENABLE_MPI AND WIN32)
	LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
	target_link_libraries(ex15 FemTech parmetis Blas Lapack ${MPI_LIBRARIES})
	# copy dlls to build directory
	file(COPY ${LAPACK}/bin/libcblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex15)
	file(COPY ${LAPACK}/bin/libblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex15)
	file(COPY ${LAPACK}/bin/liblapack.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex15)
	file(COPY ${LAPACK}/bin/liblapacke.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex15)
endif()

file(COPY brain.k DESTINATION .)
file(COPY materials.dat DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 1 \"$(TargetPath)\" brain.k" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex15/ex15.vcxproj.user @ONLY )
endif()
