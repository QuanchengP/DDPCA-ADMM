#ifndef TriangularQuadrature_hpp
#define TriangularQuadrature_hpp

#include "../General/General.hpp"

#include <vector>
#include <cmath>

namespace Ddpca {

/****************************************************************************************************/
// Implementation of numerical Gaussian quadrature on triangular elements
// integral boundary transformation: from [-1, 1] * [-1, 1] to [0, 1] * [0, 1-x]
// H.T. Rathod, K.V. Nagaraja, B. Venkatesudu, N.L. Ramesh.
// Gauss Legendre quadrature over a triangle. IIS, 2004, 84: 183-188.
// Template class definition with default 2-point Gaussian quadrature
template <I64 NumbPoints = 2>
class TriangularQuadrature {
public:
    // Member variable declarations
    I64 numbGaussPoints;                      // Number of Gaussian integration points
    std::vector<std::array<Real, 3>> gaussPoints; // Natural coordinates of Gaussian integration points
    std::vector<Real> weights;                   // Gaussian integration weight factors
    
    // Validate the validity of integration point count
    static_assert(NumbPoints == 2 || NumbPoints == 3,
                 "TriangularQuadrature only supports 2 or 3 Gauss points");
    
public:
    // Constructor - Initialize Gaussian integration points for triangular element
    TriangularQuadrature() {
        // Get Gaussian quadrature points and weights based on template parameter
        constexpr auto gaussPoints1D = GaussData<NumbPoints>::GetPoints();
        constexpr auto gaussWeights1D = GaussData<NumbPoints>::GetWeights();
        
        // Calculate number of 2D Gaussian integration points and initialize data structures
        I64 gapoSize = gaussPoints1D.size();
        numbGaussPoints = gapoSize * gapoSize;
        gaussPoints.resize(numbGaussPoints);
        weights.resize(numbGaussPoints);
        
        // Generate Gaussian integration points and weights for triangular element
        for(I64 ti = 0; ti < gapoSize; ++ti) {
            for(I64 tj = 0; tj < gapoSize; ++tj) {
                const I64 index = ti * gapoSize + tj;
                gaussPoints[index][0] = (1.0 + gaussPoints1D[ti]) / 2.0;
                gaussPoints[index][1] = 
                    (1.0 - gaussPoints1D[ti]) * (1.0 + gaussPoints1D[tj]) / 4.0;
                gaussPoints[index][2] = 
                    1.0 - gaussPoints[index][0] - gaussPoints[index][1];
                weights[index] = 
                    (1.0 - gaussPoints1D[ti]) / 8.0 * gaussWeights1D[ti] * gaussWeights1D[tj];
            }
        }
    }
    
    // Destructor
    ~TriangularQuadrature() = default;
};

// Global instance - Default with 2-point Gaussian quadrature
inline TriangularQuadrature<2>& GetTriangularQuadrature() {
    static TriangularQuadrature<2> instance;
    return instance;
}

// Global instance - Customizable integration point count
template <I64 NumbPoints>
inline TriangularQuadrature<NumbPoints>& GetTriangularQuadrature() {
    static TriangularQuadrature<NumbPoints> instance;
    return instance;
}

} // namespace Ddpca

#endif // TriangularQuadrature_hpp