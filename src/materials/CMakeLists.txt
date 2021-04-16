Set(materials_src
    materials.cpp
		CompressibleNeoHookean.cpp
    CompressibleNeoHookeanAbaqus.cpp
		StVenantKirchhoff.cpp
    LinearElastic.cpp
    HGOIsotropic.cpp
    HGOIsotropicViscoelastic.cpp
    Ogden.cpp
    OgdenViscoelastic.cpp
    lsDynaKMEquivalent.cpp
    lsDynaKMEquivalentNew.cpp
    )

add_library(materials OBJECT ${materials_src})
add_dependencies(materials JsonCpp)
