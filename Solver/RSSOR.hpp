#ifndef _RSSOR_hpp
#define _RSSOR_hpp

#include "../General/SparseMatrix.hpp"

namespace Ddpca{

//relaxation: symmetric successive over-relaxation
class RSSOR{
public:

    Real omega = 0.9;
    
    //Cache-friendly: help much
    SparseMatrix hierStriLowe;//1 ~ maxLevel
    SparseMatrix hierStriUppe;//1 ~ maxLevel
    SparseMatrix hierarchyDiagonal;//1 ~ maxLevel

    AlignedVectorRx p_0;
    AlignedVectorRx bp_0;
    AlignedVectorRx p_1;
    AlignedVectorRx x_1;

public:

    void Establish(const SparseMatrix& hierarchyStiffness);

    void Apply(
        const AlignedVectorRx& rhs, 
        AlignedVectorRx& solution, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

}; // class RSSOR

} // namespace Ddpca

#endif //_RSSOR_hpp