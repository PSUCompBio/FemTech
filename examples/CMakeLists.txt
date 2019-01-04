option (EXAMPLE1  "Build Example 1" OFF)
if (EXAMPLE1)
	add_subdirectory(ex1)
	# Now define the startup project based on the example project that is enabled
	set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT "ex1")
endif (EXAMPLE1)
