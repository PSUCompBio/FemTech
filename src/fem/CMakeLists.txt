Set(fem_src
		Assembly.cpp
    BoundaryConditions.cpp
		AllocateArrays.cpp
		ShapeFunctions/ShapeFunctions.cpp
		ShapeFunctions/GaussQuadrature3D.cpp
		ShapeFunctions/ShapeFunction_C3D8.cpp
		ShapeFunctions/ShapeFunction_C3D4.cpp
		Mass/Mass3D.cpp
		Stiffness/Stiffness3D.cpp
    Solver/SteadyImplicit.cpp
    )

add_library(fem OBJECT ${fem_src})
