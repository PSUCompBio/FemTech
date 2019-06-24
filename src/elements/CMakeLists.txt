Set(elements_src
    elements.cpp
		ElementCalculations/CalculateDeformationGradient.cpp
		ElementCalculations/StrainDisplacementMatrix.cpp
    CharacteristicLength/CalculateCharacteristicLength.cpp
    CharacteristicLength/CalculateCharacteristicLength_C3D4.cpp
    CharacteristicLength/CalculateCharacteristicLength_C3D8.cpp
    ElementCalculations/CalculateStrain.cpp
    ElementCalculations/SumOfDeformationGradient.cpp
    )

add_library(elements OBJECT ${elements_src})
