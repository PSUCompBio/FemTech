SET(ex1_src ex1.cpp)

ADD_EXECUTABLE(ex1 ${ex1_src})

if(ENABLE_MPI)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex1 FemTech parmetis ${lapackblas_libraries})
else()
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex1 FemTech)
endif(ENABLE_MPI)

# NOTE: once Cmake release version 3.13 the following
# option can be used to set up properties of ex1 in
# Visual Studio. Only the 1st one (VS_DEBUGGER_COMMAND)
# works with Cmake version 3.12. version 3.13 supports the
# second option. for more information see:
# https://stackoverflow.com/questions/52649753/how-to-send-visual-studio-configuration-options-from-cmake
#
#if(WIN32)
	#set_target_properties(ex1 PROPERTIES
	#	VS_DEBUGGER_COMMAND "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe"
	#	VS_DEBUGGER_COMMAND_ARGUMENTS "-n 2 \"$(TargetPath)\" 1-elt-cube.k")
#endif(WIN32)

file(COPY mixed-hex-tet.k DESTINATION .)
file(COPY materials.dat DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 2 \"$(TargetPath)\" mixed-hex-tet.k" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex1/ex1.vcxproj.user @ONLY )
endif()
