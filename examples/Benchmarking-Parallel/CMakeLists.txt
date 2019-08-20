SET(Benchmarking-Parallel_src Benchmarking-Parallel.cpp)

ADD_EXECUTABLE(Benchmarking-Parallel ${Benchmarking-Parallel_src})

if(ENABLE_MPI)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib ${CMAKE_BINARY_DIR}/third-party/parmetis-4.0.3/libparmetis)
  TARGET_LINK_LIBRARIES(Benchmarking-Parallel FemTech parmetis ${lapackblas_libraries})
else()
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
  TARGET_LINK_LIBRARIES(Benchmarking-Parallel FemTech)
endif(ENABLE_MPI)

file(INSTALL 10elements.inp DESTINATION .)
file(INSTALL 20elements.inp DESTINATION .)
file(INSTALL 27elements.inp DESTINATION .)
file(INSTALL materials.dat DESTINATION .)
file(INSTALL README.md DESTINATION .)
file(INSTALL abaqus DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)

# Define debugger options when working with windows & visual studio generator
# We set them here but could be passed as options if required
# We use a template for the same, replacing the options we need
if(WIN32 AND MSVC)
	set( DEBUG_COMMAND "C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 1 \"$(TargetPath)\" 10elements.inp" )
  set( DEBUG_COMMAND_ARGUMENTS "-n 1 \"$(TargetPath)\" 20elements.inp" )
	set( DEBUG_COMMAND_ARGUMENTS "-n 1 \"$(TargetPath)\" 27elements.inp" )
	configure_file( template.vcxproj.user.in ${CMAKE_BINARY_DIR}/examples/Benchmarking-Parallel/Benchmarking-Parallel.vcxproj.user @ONLY )
endif()
