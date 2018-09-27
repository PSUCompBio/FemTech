SET(ex1_src ex1.cpp)

ADD_EXECUTABLE(db ${ex1_src})

TARGET_LINK_LIBRARIES(db digitalbrain)

#file(COPY input.inp DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)
