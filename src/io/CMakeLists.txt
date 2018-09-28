Set(io_src
    input/ReadInputFile.cpp
    output/WriteVTU.cpp
    )

add_library(io OBJECT ${io_src})
