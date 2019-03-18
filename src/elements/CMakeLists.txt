Set(elements_src
    elements.cpp
		ElementCalculations/CalculateDeformationGradient.cpp
    )

add_library(elements OBJECT ${elements_src})
