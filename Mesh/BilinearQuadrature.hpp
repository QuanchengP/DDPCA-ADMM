#ifndef BilinearQuadrature_hpp
#define BilinearQuadrature_hpp

#include "../General/General.hpp"
#include "../General/DenseMatrix.hpp"
#include "../Mesh/Coordinate.hpp"

#include <array>
#include <vector>
#include <cmath>

namespace Ddpca{

/****************************************************************************************************/
// Implementation of numerical Gaussian quadrature on bilinear elements
// Template class definition with default 2-point Gaussian quadrature
template <I64 NumbPoints = 2> 
class BilinearQuadrature {
public:
    // Member variable declarations
    I64 numbGaussPoints;                    // Number of Gaussian integration points
    std::vector<std::array<Real, 2>> gaussPoints; // Natural coordinates of Gaussian integration points
    std::vector<Real> weights;                   // Gaussian integration weight factors
    std::vector<std::array<Real, 2>> cornerNodes; // Natural coordinates of four corner nodes
    
    // Shape function values at integration points (numbGaussPoints × 4 matrix)
    DenseMatrix shapeFunctions;
    // Partial derivatives of shape functions with respect to ξ/η (2 × 4 matrix for each point)
    std::vector<DenseMatrix> shapeDerivatives; 

public:
    // Constructor - Initialize Gaussian integration points for bilinear element
    BilinearQuadrature();
    
    // Destructor
    ~BilinearQuadrature() = default;
    
    // Validate the validity of integration point count
    static_assert(NumbPoints == 2 || NumbPoints == 3, 
                 "BilinearQuadrature only supports 2 or 3 Gauss points");
    
    // Private methods
    // Initialize natural coordinates of four corner nodes
    void InitializeCornerNodes();
    
    // Generate 2D Gaussian integration points and weights
    void GenerateGaussPointsAndWeights(
        const std::array<Real, NumbPoints>& points1D, 
        const std::array<Real, NumbPoints>& weights1D);
    
    // Compute shape functions and their derivatives
    void ComputeShapeFunctions();
};

// Constructor implementation
template <I64 NumbPoints>
BilinearQuadrature<NumbPoints>::BilinearQuadrature() {
    // Get Gaussian quadrature points and weights based on template parameter
    constexpr auto gaussPoints1D = GaussData<NumbPoints>::GetPoints();
    constexpr auto gaussWeights1D = GaussData<NumbPoints>::GetWeights();
    
    // Initialize natural coordinates of four corner nodes
    InitializeCornerNodes();
    
    // Calculate number of 2D Gaussian integration points and initialize data structures
    numbGaussPoints = gaussPoints1D.size() * gaussPoints1D.size();
    gaussPoints.resize(numbGaussPoints);
    weights.resize(numbGaussPoints);
    // Create a matrix of size numbGaussPoints x 4 for shape functions
    shapeFunctions.Resize(numbGaussPoints, 4);
    shapeDerivatives.resize(numbGaussPoints);
    // Create a 2x4 matrix for each integration point's derivatives
    for (I64 i = 0; i < numbGaussPoints; ++i) {
        shapeDerivatives[i].Resize(2, 4);
    }
    
    // Generate 2D Gaussian integration points and weights
    GenerateGaussPointsAndWeights(gaussPoints1D, gaussWeights1D);
    
    // Compute shape functions and their derivatives
    ComputeShapeFunctions();
}

// Initialize natural coordinates of four corner nodes
template <I64 NumbPoints>
void BilinearQuadrature<NumbPoints>::InitializeCornerNodes() {
    cornerNodes = {
        {-1.0, -1.0}, {1.0, -1.0},
        {1.0, 1.0}, {-1.0, 1.0}
    };
}

// Generate 2D Gaussian integration points and weights
template <I64 NumbPoints>
void BilinearQuadrature<NumbPoints>::GenerateGaussPointsAndWeights(
    const std::array<Real, NumbPoints>& points1D, 
    const std::array<Real, NumbPoints>& weights1D) {
    //
    const I64 dim = points1D.size();
    
    for (I64 i = 0; i < dim; ++i) {
        for (I64 j = 0; j < dim; ++j) {
            // Avoid repeated calculations when computing the index
            const I64 index = i * dim + j;
            gaussPoints[index][0] = points1D[i];
            gaussPoints[index][1] = points1D[j];
            weights[index] = weights1D[i] * weights1D[j];
        }
    }
}

// Compute shape functions and their derivatives
template <I64 NumbPoints>
void BilinearQuadrature<NumbPoints>::ComputeShapeFunctions() {
    //
    const Real inv4 = 1.0 / 4.0;
    
    for (I64 i = 0; i < numbGaussPoints; ++i) {
        const auto& gp = gaussPoints[i];
        auto& derivs = shapeDerivatives[i];
        
        for (I64 j = 0; j < 4; ++j) {
            const auto& cn = cornerNodes[j];
            // Precompute common factors
            const Real xi_term = 1.0 + cn[0] * gp[0];
            const Real eta_term = 1.0 + cn[1] * gp[1];
            
            // Calculate partial derivatives
            derivs(0, j) = cn[0] * eta_term * inv4;  // Partial derivative with respect to xi
            derivs(1, j) = xi_term * cn[1] * inv4;  // Partial derivative with respect to eta
            
            // Calculate shape function values
            shapeFunctions(i, j) = xi_term * eta_term * inv4;
        }
    }
}

// Global instance - Default with 2-point Gaussian quadrature
inline BilinearQuadrature<2>& GetBilinearQuadrature() {
    static BilinearQuadrature<2> instance;
    return instance;
}

// Global instance - Customizable integration point count
template <I64 NumbPoints>
inline BilinearQuadrature<NumbPoints>& GetBilinearQuadrature() {
    static BilinearQuadrature<NumbPoints> instance;
    return instance;
}

// Function to calculate bilinear element integration weight
template <I64 NumbPoints>
Real BilinearQuadratureJacobian(
    const std::array<Real, 2> &xiEtaParams, 
    const std::vector<Coordinate> &elementCoordinates) {      
    // Create two 3D vectors representing partial derivatives with respect to xi and eta
    Coordinate dXdxi;
    Coordinate dXdeta;
    
    const Real inv4 = 1.0 / 4.0;
    
    // Use global instance with 2-point Gaussian quadrature
    const auto& biliQuad = GetBilinearQuadrature();
    
    for(I64 ti = 0; ti < 3; ++ti) {
        for(I64 tj = 0; tj < 4; ++tj) {
            const auto& cn = biliQuad.cornerNodes[tj];
            // Fill partial derivative vector with respect to xi
            dXdxi[ti] += elementCoordinates[tj][ti] * 
                (cn[0] * inv4 + cn[0] * cn[1] * xiEtaParams[1] * inv4); 
            // Fill partial derivative vector with respect to eta
            dXdeta[ti] += elementCoordinates[tj][ti] * 
                (cn[1] * inv4 + cn[0] * cn[1] * xiEtaParams[0] * inv4); 
        } 
    } 
    
    // Calculate cross product of the two partial derivative vectors
    Coordinate crossProduct = dXdxi.Cross(dXdeta);
    
    // Calculate L2 norm of cross product as weight factor
    return NRM2(crossProduct.data);
}

} // namespace Ddpca

#endif // BilinearQuadrature_hpp