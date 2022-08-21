Set(elements_src
    Elements.cpp
    ElementCalculations/CalculateDeformationGradient.cpp
    ElementCalculations/CalculateF_XiAndInverse.cpp
    ElementCalculations/ComputeH.cpp
    ElementCalculations/ComputeHG.cpp
    ElementCalculations/StrainDisplacementMatrix.cpp
    CharacteristicLength/CalculateCharacteristicLength.cpp
    CharacteristicLength/CalculateCharacteristicLength_C3D4.cpp
    CharacteristicLength/CalculateCharacteristicLength_C3D8.cpp
    ElementCalculations/CalculateStrain.cpp
    ElementCalculations/CalculateElementStress.cpp
    ElementCalculations/CalculateElementPressure.cpp
    CenterOfMass/CalculateCentroidAndVolume.cpp
    CenterOfMass/CalculateCentroid_C3D4.cpp
    CenterOfMass/CalculateCentroidAndVolume_C3D8.cpp
    )

add_library(elements OBJECT ${elements_src})
add_dependencies(elements JsonCpp)
