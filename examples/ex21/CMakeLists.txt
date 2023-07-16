SET(ex21_src ex21.cpp)

ADD_EXECUTABLE(ex21 ${ex21_src})
ADD_DEPENDENCIES(ex21 JsonCpp)

if(ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex21 FemTech parmetis  ${lapackblas_libraries} jsoncpp)
elseif(NOT ENABLE_MPI AND NOT WIN32)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex21 FemTech)
elseif(ENABLE_MPI AND WIN32)
	LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
	target_link_libraries(ex21 FemTech parmetis Blas Lapack ${MPI_LIBRARIES})
	# copy dlls to build directory
	file(COPY ${LAPACK}/bin/libcblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex21)
	file(COPY ${LAPACK}/bin/libblas.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex21)
	file(COPY ${LAPACK}/bin/liblapack.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex21)
	file(COPY ${LAPACK}/bin/liblapacke.dll DESTINATION ${CMAKE_BINARY_DIR}/examples/ex21)
endif()

file(INSTALL fine_brain.inp DESTINATION .)
file(INSTALL coarse_brain.inp DESTINATION .)
file(INSTALL materials.dat DESTINATION .)
file(INSTALL pressureInput.json DESTINATION .)
file(INSTALL README.md DESTINATION .)
file(INSTALL pressureBC.json DESTINATION .)
file(INSTALL mapPressureBCtoScaledMesh.py DESTINATION .)
file(INSTALL pressureToNodeIDMap.json DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 2 \"$(TargetPath)\" fine_brain.inp" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex5/ex5.vcxproj.user @ONLY )
endif()
