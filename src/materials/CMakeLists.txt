Set(materials_src
    materials.cpp
		CompressibleNeoHookean.cpp
    NeoHookeanAbaqus.cpp
		StVenantKirchhoff.cpp
    LinearElastic.cpp
    HGOIsotropic.cpp
    HGOIsotropicViscoelastic.cpp
    Viscoelastic.cpp
    Hourglass.cpp
    )

add_library(materials OBJECT ${materials_src})
add_dependencies(materials JsonCpp)
