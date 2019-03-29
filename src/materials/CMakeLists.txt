Set(materials_src
    materials.cpp
		CompressibleNeoHookean.cpp
    )

add_library(materials OBJECT ${materials_src})
