Set(io_src
    input/ReadInputFile.cpp
    output/WriteVTU.cpp
		output/WritePVD.cpp
    output/FreeArrays.cpp
    PartitionMesh.cpp
    input/ReadLsDyna.cpp
    input/ReadAbaqus.cpp
		input/ReadMaterials.cpp
    input/ReadJSON.cpp
    InitFinalizeFemTech.cpp
    parallel_log.cpp
    input/ReadMapping.cpp
    )

add_library(io OBJECT ${io_src})
add_dependencies(io JsonCpp)
