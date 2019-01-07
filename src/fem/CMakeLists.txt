Set(fem_src
    fem.cpp
		ShapeFunctions.cpp
		ShapeFunction_C3D8.cpp
		GaussQuadrature3D.cpp
    )

add_library(fem OBJECT ${fem_src})
