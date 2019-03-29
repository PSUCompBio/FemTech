Set(materials_src
    materials.cpp
		CompressibleNeoHookean.cpp
		StVenantKirchhoff.cpp
    )

add_library(materials OBJECT ${materials_src})
