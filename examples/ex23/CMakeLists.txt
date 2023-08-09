SET(ex23_src ex23.cpp)

ADD_EXECUTABLE(ex23 ${ex23_src})
ADD_DEPENDENCIES(ex23 JsonCpp)

if(ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex23 FemTech parmetis ${lapackblas_libraries})
elseif(NOT ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex23 FemTech)
elseif(ENABLE_MPI AND WIN32)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
	target_link_libraries(ex23 FemTech parmetis Blas Lapack ${MPI_LIBRARIES})
	# copy dlls to build directory
	file(COPY ${LAPACK}/bin/libcblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex23)
	file(COPY ${LAPACK}/bin/libblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex23)
	file(COPY ${LAPACK}/bin/liblapack.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex23)
	file(COPY ${LAPACK}/bin/liblapacke.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex23)
endif()

file(INSTALL brainfiber.inp DESTINATION .)
file(INSTALL brainfibermap.txt DESTINATION .)
file(INSTALL fiber_test.json DESTINATION .)
file(INSTALL materials.dat DESTINATION .)
file(INSTALL updateOutputJson.py DESTINATION .)
file(INSTALL gnuplot.script DESTINATION .)
file(INSTALL gnuplotEnergy.script DESTINATION .)
file(INSTALL README.md DESTINATION .)
file(INSTALL simulationMovie.py DESTINATION .)
file(INSTALL mps95Movie.py DESTINATION .)
file(INSTALL addGraph.py DESTINATION .)
file(INSTALL test_ex23.sh DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 2 \"$(TargetPath)\" brain.inp" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex23/ex23.vcxproj.user @ONLY )
endif()
