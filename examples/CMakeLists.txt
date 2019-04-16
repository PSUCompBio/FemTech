option (EXAMPLE1  "3D Mixed Hex Tet Sim" OFF)
option (EXAMPLE2  "3D Sphere Mesh Sim" OFF)
option (EXAMPLE3  "3D Multipart tet sim using Abaqus INP" OFF)
option (EXAMPLE4  "3D 3-cube multipart tet sim using Abaqus INP" OFF)
option (EXAMPLE5  "3D Brain Sim" OFF)
option (EXAMPLE6  "3D Larger 10e Dyna Input" OFF)
option (EXAMPLE7  "3D 2-element hex model" OFF)
option (EXAMPLE8  "2D 2-element quad model" OFF)
option (EXAMPLE9  "3D 1-element hex model" OFF)
option (EXAMPLE10  "3D 1-element tet model" OFF)
option (EXAMPLE11  "3D 1-element hex model" OFF)
option (EXAMPLE12  "3D 1-element truss model - lsdyna k format" OFF)
option (EXAMPLE13  "3D brain tractography truss model - abaqus inp format" OFF)

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

if (EXAMPLE8)
	add_subdirectory(ex8)
endif (EXAMPLE8)

if (EXAMPLE9)
	add_subdirectory(ex9)
endif (EXAMPLE9)

if (EXAMPLE10)
	add_subdirectory(ex10)
endif (EXAMPLE10)

if (EXAMPLE11)
	add_subdirectory(ex11)
endif (EXAMPLE11)

if (EXAMPLE12)
	add_subdirectory(ex12)
endif (EXAMPLE12)

if (EXAMPLE13)
	add_subdirectory(ex13)
endif (EXAMPLE13)
