Set(fem_src
    fem.cpp
	ShapeFunctions3d.cpp
    )

add_library(fem OBJECT ${fem_src})
