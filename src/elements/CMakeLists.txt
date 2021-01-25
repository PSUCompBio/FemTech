Set(elements_src
    Elements.cpp
    CorotationalSystem.cpp
    CorotationalSystemTransformation.cpp
    ShellHourglassStrain.cpp
		ElementCalculations/CalculateDeformationGradient.cpp
		ElementCalculations/StrainDisplacementMatrix.cpp
    ElementCalculations/StrainDisplacementMatrixShell.cpp
    ElementCalculations/ShellArea.cpp
    CharacteristicLength/CalculateCharacteristicLength.cpp
    CharacteristicLength/CalculateCharacteristicLength_C3D4.cpp
    CharacteristicLength/CalculateCharacteristicLength_C3D8.cpp
    CharacteristicLength/CalculateCharacteristicLength_S4.cpp
    ElementCalculations/CalculateStrain.cpp
    CenterOfMass/CalculateCentroidAndVolume.cpp
    CenterOfMass/CalculateCentroid_C3D4.cpp
    CenterOfMass/CalculateCentroidAndVolume_C3D8.cpp
    )

add_library(elements OBJECT ${elements_src})
add_dependencies(elements JsonCpp)
