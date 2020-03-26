SET(ex18_src ex18.cpp)

ADD_EXECUTABLE(ex18 ${ex18_src})
ADD_DEPENDENCIES(ex18 JsonCpp)

if(ENABLE_MPI)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex18 FemTech parmetis  ${lapackblas_libraries} jsoncpp)
else()
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex18 FemTech)
endif(ENABLE_MPI)

file(INSTALL coarse_brain.inp DESTINATION .)
file(INSTALL materials.dat DESTINATION .)
file(INSTALL input_case9218Jietal2015.json DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 2 \"$(TargetPath)\" brain.inp" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex18/ex18.vcxproj.user @ONLY )
endif()
