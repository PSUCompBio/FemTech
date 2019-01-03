option (EXAMPLE1  "Example of a Mixed Hex Tet Sim" OFF)
option (EXAMPLE2  "Example of a Sphere Mesh Sim" OFF)

if (EXAMPLE1)
   add_subdirectory(ex1)
endif (EXAMPLE1)

if (EXAMPLE2)
	 add_subdirectory(ex2)
endif (EXAMPLE2)
