#ifndef _RJacobian_hpp
#define _RJacobian_hpp

#include "../General/SparseMatrix.hpp"

namespace Ddpca {

//relaxation: jacobian
class RJacobian{

public:

    SparseMatrix hierarchyDiagonal;//1 ~ maxLevel
    AlignedVectorRx residual;

    void Establish(const SparseMatrix& hierarchyStiffness);

    void Apply(
        const SparseMatrix& hierarchyStiffness, 
        const AlignedVectorRx& rhs, 
        AlignedVectorRx& solution, 
        const I64 times, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

}; // class RJacobian

} // namespace Ddpca

#endif // _RJacobian_hpp