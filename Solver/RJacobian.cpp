#include "RJacobian.hpp"
#include "../Mesh/Coordinate.hpp"

namespace Ddpca {

void RJacobian::Establish(const SparseMatrix& hierarchyStiffness){
    //
    hierarchyDiagonal.M = hierarchyStiffness.M;
    hierarchyDiagonal.N = hierarchyStiffness.N;
    hierarchyDiagonal.nnz = hierarchyStiffness.M;
    hierarchyDiagonal.row_ptr = AlignedAllocate<I64>(hierarchyDiagonal.M + 1);
    hierarchyDiagonal.row_ptr[0] = 0;
    hierarchyDiagonal.col_ind = AlignedAllocate<I64>(hierarchyDiagonal.nnz);
    hierarchyDiagonal.val = AlignedAllocate<Real>(hierarchyDiagonal.nnz);
    //
    for(I64 ti = 0; ti < hierarchyStiffness.M; ++ ti){
        I64 tj_start = hierarchyStiffness.row_ptr[ti];
        I64 tj_end = hierarchyStiffness.row_ptr[ti + 1];
        for(I64 tj = tj_start; tj < tj_end; ++ tj){
            I64 col_tj = hierarchyStiffness.col_ind[tj];
            Real val_tj = hierarchyStiffness.val[tj];
            if(ti == col_tj){
                hierarchyDiagonal.row_ptr[ti + 1] = ti + 1;
                hierarchyDiagonal.col_ind[ti] = col_tj;
                hierarchyDiagonal.val[ti] = val_tj;
                break;
            }
        }
    }
    residual.resize(hierarchyStiffness.M);
}

void RJacobian::Apply(
    const SparseMatrix& hierarchyStiffness, 
    const AlignedVectorRx& rhs, 
    AlignedVectorRx& solution, 
    const I64 times, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    for(I64 ti = 0; ti < times; ++ ti){
        MV(-1.0, hierarchyStiffness, solution, 1.0, rhs, residual, threadTask, nestLevel);
        std::transform(std::execution::unseq, 
            residual.begin(), residual.end(), 
            hierarchyDiagonal.val, 
            residual.begin(), 
            [](const Real x, const Real y){return x / y;});
        XPEY(solution, residual);
    }
}

} // namespace Ddpca