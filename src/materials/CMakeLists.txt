Set(materials_src
    materials.cpp
		CompressibleNeoHookean.cpp
		StVenantKirchhoff.cpp
    LinearElastic.cpp
    HGOIsotropic.cpp
    HGOIsotropicViscoelastic.cpp
    )

add_library(materials OBJECT ${materials_src})
