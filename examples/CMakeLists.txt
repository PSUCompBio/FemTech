option (EXAMPLE1  "Example of a Mixed Hex Tet Sim" OFF)
option (EXAMPLE2  "Example of a Sphere Mesh Sim" OFF)

if (EXAMPLE1)
	add_subdirectory(ex1)
	# Now define the startup project based on the example project that is enabled
	set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT "ex1")
endif (EXAMPLE1)

if (EXAMPLE2)
	 add_subdirectory(ex2)
endif (EXAMPLE2)
