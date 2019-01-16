Set(fem_src
    fem.cpp
		ShapeFunctions.cpp
		GaussQuadrature3D.cpp
		ShapeFunction_C3D8.cpp
		ShapeFunction_C3D4.cpp
    )

add_library(fem OBJECT ${fem_src})
