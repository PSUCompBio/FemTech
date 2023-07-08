SET(ex3_src ex3.cpp)

ADD_EXECUTABLE(ex3 ${ex3_src})

if(ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex3 FemTech parmetis ${lapackblas_libraries})
elseif(NOT ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex3 FemTech)
elseif(ENABLE_MPI AND WIN32)
	LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
	target_link_libraries(ex3 FemTech parmetis Blas Lapack ${MPI_LIBRARIES})
	# copy dlls to build directory
	file(COPY ${CMAKE_SOURCE_DIR}/lib/libcblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex3)
	file(COPY ${CMAKE_SOURCE_DIR}/lib/libblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex3)
	file(COPY ${CMAKE_SOURCE_DIR}/lib/liblapack.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex3)
	file(COPY ${CMAKE_SOURCE_DIR}/lib/liblapacke.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex3)
endif()

file(COPY multipart-simple-tet.inp DESTINATION .)
file(COPY materials.dat DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 2 \"$(TargetPath)\" multipart-simple-tet.inp" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex3/ex3.vcxproj.user @ONLY )
endif()
