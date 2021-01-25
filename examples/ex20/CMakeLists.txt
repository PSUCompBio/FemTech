SET(ex20_src ex20.cpp)

ADD_EXECUTABLE(ex20 ${ex20_src})

if(ENABLE_MPI)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(ex20 FemTech parmetis ${lapackblas_libraries})
else()
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(ex20 FemTech)
endif(ENABLE_MPI)

file(INSTALL 1-elt-cube.k DESTINATION .)
file(INSTALL 1shell.inp DESTINATION .)
file(INSTALL 2shells.inp DESTINATION .)
file(INSTALL materials.dat DESTINATION .)
file(INSTALL compareResults.py DESTINATION .)
file(INSTALL gnuplot.script DESTINATION .)
file(INSTALL gnuplotEnergy.script DESTINATION .)
file(INSTALL README.md DESTINATION .)
file(INSTALL LSDYNA DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 1 \"$(TargetPath)\" 1-elt-cube.k" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/ex20/ex20.vcxproj.user @ONLY )
endif()
