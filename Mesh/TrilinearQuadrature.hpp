#ifndef TrilinearQuadrature_hpp
#define TrilinearQuadrature_hpp

#include "../General/General.hpp"
#include "../General/DenseMatrix.hpp"

#include <vector>

namespace Ddpca{

/****************************************************************************************************/
// Implementation of numerical Gaussian quadrature on trilinear elements
// Template class definition with default 2-point Gaussian quadrature
template <I64 NumbPoints = 2> 
class TrilinearQuadrature {
public:
    // Member variable declarations
    I64 numbGaussPoints;                    // Number of Gaussian integration points
    std::vector<std::array<Real, 3>> gaussPoints; // Natural coordinates of Gaussian integration points
    std::vector<Real> weights;                   // Gaussian integration weight factors
    std::vector<std::array<Real, 3>> cornerNodes; // Natural coordinates of eight corner nodes
    
    // Shape function values at integration points (numbGaussPoints × 8 matrix)
    DenseMatrix shapeFunctions; 
    // Partial derivatives of shape functions with respect to ξ/η/ζ (3 × 8 matrix for each point)
    std::vector<DenseMatrix> shapeDerivatives; 
    
    // Validate the validity of integration point count
    static_assert(NumbPoints == 2 || NumbPoints == 3, 
                 "TrilinearQuadrature only supports 2 or 3 Gauss points");
    
    // Private methods
    // Initialize natural coordinates of eight corner nodes
    void InitializeCornerNodes();
    
    // Generate 3D Gaussian integration points and weights
    void GenerateGaussPointsAndWeights(const std::array<Real, NumbPoints>& points1D, 
                                     const std::array<Real, NumbPoints>& weights1D);
    
    // Compute shape functions and their derivatives
    void ComputeShapeFunctions();

public:
    // Constructor - Initialize Gaussian integration points for trilinear element
    TrilinearQuadrature();
    
    // Destructor
    ~TrilinearQuadrature() = default;
};

// Constructor implementation
template <I64 NumbPoints>
TrilinearQuadrature<NumbPoints>::TrilinearQuadrature() {
    // Get Gaussian quadrature points and weights based on template parameter
    constexpr auto gaussPoints1D = GaussData<NumbPoints>::GetPoints();
    constexpr auto gaussWeights1D = GaussData<NumbPoints>::GetWeights();
    
    // Initialize natural coordinates of eight corner nodes
    InitializeCornerNodes();
    
    // Calculate number of 3D Gaussian integration points and initialize data structures
    numbGaussPoints = gaussPoints1D.size() * gaussPoints1D.size() * gaussPoints1D.size();
    gaussPoints.resize(numbGaussPoints);
    weights.resize(numbGaussPoints);
    
    // Now using the updated DenseMatrix which supports arbitrary dimensions
    // Create a matrix of size numbGaussPoints x 8 for shape functions
    shapeFunctions.Resize(numbGaussPoints, 8);
    shapeDerivatives.resize(numbGaussPoints);
    
    // Create a 3x8 matrix for each integration point's derivatives
    for (I64 i = 0; i < numbGaussPoints; ++i) {
        shapeDerivatives[i].Resize(3, 8);
    }
    
    // Generate 3D Gaussian integration points and weights
    GenerateGaussPointsAndWeights(gaussPoints1D, gaussWeights1D);
    
    // Compute shape functions and their derivatives
    ComputeShapeFunctions();
}

// Initialize natural coordinates of eight corner nodes
template <I64 NumbPoints>
void TrilinearQuadrature<NumbPoints>::InitializeCornerNodes() {
    cornerNodes = {
        {-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0},
        {1.0, 1.0, -1.0}, {-1.0, 1.0, -1.0},
        {-1.0, -1.0, 1.0}, {1.0, -1.0, 1.0},
        {1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}
    };
}

// Generate 3D Gaussian integration points and weights
template <I64 NumbPoints>
void TrilinearQuadrature<NumbPoints>::GenerateGaussPointsAndWeights(
    const std::array<Real, NumbPoints>& points1D, 
    const std::array<Real, NumbPoints>& weights1D) {
    const I64 dim = points1D.size();
    
    for (I64 i = 0; i < dim; ++i) {
        for (I64 j = 0; j < dim; ++j) {
            for (I64 k = 0; k < dim; ++k) {
                // Avoid repeated calculations when computing the index
                const I64 index = i * dim * dim + j * dim + k;
                gaussPoints[index][0] = points1D[i];
                gaussPoints[index][1] = points1D[j];
                gaussPoints[index][2] = points1D[k];
                weights[index] = weights1D[i] * weights1D[j] * weights1D[k];
            }
        }
    }
}

// Compute shape functions and their derivatives
template <I64 NumbPoints>
void TrilinearQuadrature<NumbPoints>::ComputeShapeFunctions() {
    // 预计算常用值以减少重复计算
    const Real inv8 = 1.0 / 8.0;
    
    for (I64 i = 0; i < numbGaussPoints; ++i) {
        const auto& gp = gaussPoints[i];
        auto& derivs = shapeDerivatives[i];
        
        for (I64 j = 0; j < 8; ++j) {
            const auto& cn = cornerNodes[j];
            // 预计算公共因子
            const Real xi_term = 1.0 + cn[0] * gp[0];
            const Real eta_term = 1.0 + cn[1] * gp[1];
            const Real zeta_term = 1.0 + cn[2] * gp[2];
            
            // 使用预计算的值
            derivs(0, j) = cn[0] * eta_term * zeta_term * inv8;
            derivs(1, j) = xi_term * cn[1] * zeta_term * inv8;
            derivs(2, j) = xi_term * eta_term * cn[2] * inv8;
            
            shapeFunctions(i, j) = xi_term * eta_term * zeta_term * inv8;
        }
    }
}

// Global instance - Default with 2-point Gaussian quadrature
inline TrilinearQuadrature<2>& GetTrilinearQuadrature() {
    static TrilinearQuadrature<2> instance;
    return instance;
}

// Global instance - Customizable integration point count
template <I64 NumbPoints>
inline TrilinearQuadrature<NumbPoints>& GetTrilinearQuadrature() {
    static TrilinearQuadrature<NumbPoints> instance;
    return instance;
}

} // namespace Ddpca

#endif // TrilinearQuadrature_hpp