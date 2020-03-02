Set(math_src
		DeterminateF.cpp
		InverseF.cpp
		math.cpp
    Geometry.cpp
    )

add_library(math OBJECT ${math_src})
add_dependencies(math JsonCpp)
