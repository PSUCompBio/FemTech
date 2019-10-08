SET(ex14_src ex14.cpp)

ADD_EXECUTABLE(ex14 ${ex14_src})
ADD_DEPENDENCIES(ex14 JsonCpp)

if(ENABLE_MPI)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex14 FemTech parmetis ${lapackblas_libraries} jsoncpp)
else()
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex14 FemTech)
endif(ENABLE_MPI)

file(INSTALL sphereinsphere.inp DESTINATION .)
file(INSTALL sphere-in-sphere-hex.inp DESTINATION .)
file(INSTALL spherefromansys.inp DESTINATION .)
file(INSTALL materials.dat DESTINATION .)
file(INSTALL input.json DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 1 \"$(TargetPath)\" sphereinsphere.inp" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex14/ex14.vcxproj.user @ONLY )
endif()
