option (EXAMPLE1  "Mixed Hex Tet Sim" OFF)
option (EXAMPLE2  "Sphere Mesh Sim" OFF)
option (EXAMPLE3  "Multipart tet sim using Abaqus INP" OFF)
option (EXAMPLE4  "3-cube multipart tet sim using Abaqus INP" OFF)
option (EXAMPLE5  "Brain Sim" OFF)
option (EXAMPLE6  "Larger 10e Dyna Input" OFF)
option (EXAMPLE7  "2 element hex model" OFF)

if (EXAMPLE1)
	add_subdirectory(ex1)
	# Define the startup project based on the example project that is enabled
	set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT "ex1")
endif (EXAMPLE1)

if (EXAMPLE2)
	add_subdirectory(ex2)
endif (EXAMPLE2)

if (EXAMPLE3)
	add_subdirectory(ex3)
endif (EXAMPLE3)

if (EXAMPLE4)
	add_subdirectory(ex4)
endif (EXAMPLE4)

if (EXAMPLE5)
	add_subdirectory(ex5)
endif (EXAMPLE5)

if (EXAMPLE6)
	add_subdirectory(ex6)
endif (EXAMPLE6)

if (EXAMPLE7)
	add_subdirectory(ex7)
endif (EXAMPLE7)
