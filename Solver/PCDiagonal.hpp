#ifndef _PCDiagonal_hpp
#define _PCDiagonal_hpp

#include "../General/SparseMatrix.hpp"

#include <cassert>

namespace Ddpca {

//diagonal preconditioner
class PCDiagonal{
public:

    SparseMatrix inverseDiagonal;

    void Establish(const SparseMatrix& tempStiff, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel){
        //
        inverseDiagonal.M = tempStiff.M;
        inverseDiagonal.N = tempStiff.N;
        inverseDiagonal.nnz = tempStiff.M;
        inverseDiagonal.row_ptr = AlignedAllocate<I64>(inverseDiagonal.M + 1);
        inverseDiagonal.row_ptr[0] = 0;
        inverseDiagonal.col_ind = AlignedAllocate<I64>(inverseDiagonal.nnz);
        inverseDiagonal.val = AlignedAllocate<Real>(inverseDiagonal.nnz);
        //
        I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
        std::vector<I64> startIndex(numbPart), endIndex(numbPart);
        EvenlyDistribute(tempStiff.M, numbPart, startIndex, endIndex);
        std::function<void(I64,I64)> taskFunction = 
            [&](I64 tp, I64){
                I64 start_tp = startIndex[tp];
                I64 end_tp = endIndex[tp];
                for(I64 ti = start_tp; ti < end_tp; ++ ti){
                    I64 tj_start = tempStiff.row_ptr[ti];
                    I64 tj_end = tempStiff.row_ptr[ti + 1];
                    for(I64 tj = tj_start; tj < tj_end; ++ tj){
                        I64 col_tj = tempStiff.col_ind[tj];
                        if(ti == col_tj){
                            inverseDiagonal.row_ptr[ti + 1] = ti + 1;
                            inverseDiagonal.col_ind[ti] = col_tj;
                            inverseDiagonal.val[ti] = 1.0 / tempStiff.val[tj];
                            break;
                        }
                    }
                }
            };
    }

    void Apply(const AlignedVectorRx& rhs, 
        AlignedVectorRx& solution, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel){
        //
        MV(inverseDiagonal, rhs, solution, threadTask, nestLevel);
    }

}; // class PCDiagonal

} // namespace Ddpca

#endif // _PCDiagonal_hpp