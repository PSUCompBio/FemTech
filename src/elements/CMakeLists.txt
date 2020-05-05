Set(elements_src
    Elements.cpp
		ElementCalculations/CalculateDeformationGradient.cpp
		ElementCalculations/StrainDisplacementMatrix.cpp
    CharacteristicLength/CalculateCharacteristicLength.cpp
    CharacteristicLength/CalculateCharacteristicLength_C3D4.cpp
    CharacteristicLength/CalculateCharacteristicLength_C3D8.cpp
    ElementCalculations/CalculateStrain.cpp
    ElementCalculations/Calculate_lnV.cpp
    CenterOfMass/CalculateCentroidAndVolume.cpp
    CenterOfMass/CalculateCentroid_C3D4.cpp
    CenterOfMass/CalculateCentroidAndVolume_C3D8.cpp

    )

add_library(elements OBJECT ${elements_src})
