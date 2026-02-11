#ifndef _FemStiffnessMatrix_hpp
#define _FemStiffnessMatrix_hpp

#include "Mesh.hpp"
#include "../General/DenseMatrix.hpp"
#include "../General/SparseMatrix.hpp"
#include "TrilinearQuadrature.hpp"

#include <cassert>

namespace Ddpca {

void CalculateElasticity(
    const Real materialElasticity, 
    const Real materialPoisson, 
    DenseMatrix& outputElasticity);

void FemStiffnessMatrix(
    Mesh& inputMesh, 
    const Real materialElasticity, 
    const Real materialPoisson, 
    SparseMatrix& outputMatrix, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel);

}// namespace Ddpca

#endif //_FemStiffnessMatrix_hpp