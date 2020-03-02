Set(fem_src
		Assembly.cpp
    BoundaryConditions.cpp
		AllocateArrays.cpp
		ShapeFunctions/ShapeFunctions.cpp
		ShapeFunctions/GaussQuadrature3D.cpp
		ShapeFunctions/ShapeFunction_C3D8.cpp
		ShapeFunctions/ShapeFunction_C3D4.cpp
		ShapeFunctions/ShapeFunction_T3D2.cpp
		Mass/Mass3D.cpp
		Stiffness/Stiffness3D.cpp
    Solver/SteadyImplicit.cpp
    Solver/UnsteadyNewmarkImplicit.cpp
		Solver/ExplicitDynamics.cpp
		SolidMechanics/GetForce.cpp
		SolidMechanics/GetForce_3D.cpp
		SolidMechanics/CalculateAcclerations.cpp
		SolidMechanics/CalculateFR.cpp
		SolidMechanics/StressUpdate.cpp
		SolidMechanics/InternalForceUpdate.cpp
		SolidMechanics/TrussStressForceUpdate.cpp
    SolidMechanics/CheckEnergy.cpp
    )

add_library(fem OBJECT ${fem_src})
add_dependencies(fem JsonCpp)
