Set(timestep_src
    CalculateTimeStep.cpp
		StableTimeStep.cpp
    )

add_library(timestep OBJECT ${timestep_src})
add_dependencies(timestep JsonCpp)
