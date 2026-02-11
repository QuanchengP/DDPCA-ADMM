#include "RSSOR.hpp"
#include "../Mesh/Coordinate.hpp"

namespace Ddpca{

void RSSOR::Establish(const SparseMatrix& hierarchyStiffness){
    //
    hierStriLowe.M = hierarchyStiffness.M;
    hierStriLowe.N = hierarchyStiffness.N;
    hierStriLowe.nnz = (hierarchyStiffness.nnz - hierStriLowe.M) / 2;
    hierStriLowe.row_ptr = AlignedAllocate<I64>(hierStriLowe.M + 1);
    hierStriLowe.col_ind = AlignedAllocate<I64>(hierStriLowe.nnz);
    hierStriLowe.val = AlignedAllocate<Real>(hierStriLowe.nnz);
    hierStriLowe.row_ptr[0] = 0;
    //
    hierStriUppe.M = hierarchyStiffness.M;
    hierStriUppe.N = hierarchyStiffness.N;
    hierStriUppe.nnz = hierStriLowe.nnz;
    hierStriUppe.row_ptr = AlignedAllocate<I64>(hierStriUppe.M + 1);
    hierStriUppe.col_ind = AlignedAllocate<I64>(hierStriUppe.nnz);
    hierStriUppe.val = AlignedAllocate<Real>(hierStriUppe.nnz);
    hierStriUppe.row_ptr[0] = 0;
    //
    hierarchyDiagonal.M = hierarchyStiffness.M;
    hierarchyDiagonal.N = hierarchyStiffness.N;
    hierarchyDiagonal.nnz = hierarchyStiffness.M;
    hierarchyDiagonal.row_ptr = AlignedAllocate<I64>(hierarchyDiagonal.M + 1);
    hierarchyDiagonal.col_ind = AlignedAllocate<I64>(hierarchyDiagonal.nnz);
    hierarchyDiagonal.val = AlignedAllocate<Real>(hierarchyDiagonal.nnz);
    hierarchyDiagonal.row_ptr[0] = 0;
    //
    I64 tempLnnz = 0;
    I64 tempUnnz = 0;
    for(I64 ti = 0; ti < hierarchyStiffness.M; ++ ti){
        I64 tj_start = hierarchyStiffness.row_ptr[ti];
        I64 tj_end = hierarchyStiffness.row_ptr[ti + 1];
        for(I64 tj = tj_start; tj < tj_end; ++ tj){
            I64 col_tj = hierarchyStiffness.col_ind[tj];
            Real val_tj = hierarchyStiffness.val[tj];
            if(ti > col_tj){
                hierStriLowe.col_ind[tempLnnz] = col_tj;
                hierStriLowe.val[tempLnnz] = omega * val_tj;
                ++ tempLnnz;
            }
            else if(ti < col_tj){
                hierStriUppe.col_ind[tempUnnz] = col_tj;
                hierStriUppe.val[tempUnnz] = omega * val_tj;
                ++ tempUnnz;
            }
            else if(ti == col_tj){
                hierarchyDiagonal.row_ptr[ti + 1] = ti + 1;
                hierarchyDiagonal.col_ind[ti] = col_tj;
                hierarchyDiagonal.val[ti] = val_tj;
            }
        }
        hierStriLowe.row_ptr[ti + 1] = tempLnnz;
        hierStriUppe.row_ptr[ti + 1] = tempUnnz;
    }
    p_0.resize(hierarchyStiffness.M);
    bp_0.resize(hierarchyStiffness.M);
    p_1.resize(hierarchyStiffness.M);
    x_1.resize(hierarchyStiffness.M);
}

void RSSOR::Apply(
    const AlignedVectorRx& rhs, 
    AlignedVectorRx& solution, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    I64 numbX = rhs.size();

    //
    std::transform(std::execution::unseq, 
        solution.begin(), solution.end(), 
        hierarchyDiagonal.val, 
        p_0.begin(), 
        [&](const Real x, const Real y) { return (1.0 - omega) * x * y; });
    MV(-1.0, hierStriUppe, solution, 1.0, p_0, p_0, threadTask, nestLevel);
    //
    AXPY(omega, rhs, p_0, bp_0);
    for(I64 ti = 0; ti < numbX; ++ ti){
        I64 tj_start = hierStriLowe.row_ptr[ti];
        I64 tj_end = hierStriLowe.row_ptr[ti + 1];
        alignas(nfsAlign) Real tempSumm = std::transform_reduce(std::execution::unseq, 
            hierStriLowe.col_ind + tj_start, hierStriLowe.col_ind + tj_end, 
            hierStriLowe.val + tj_start, 
            0.0, std::plus<Real>(),
            [this](const I64 c, const Real v) { return v * (this->x_1)[c]; });
        x_1[ti] = (bp_0[ti] - tempSumm) / hierarchyDiagonal.val[ti];
    }
    //
    std::transform(std::execution::unseq, 
        x_1.begin(), x_1.end(), 
        hierarchyDiagonal.val, 
        p_1.begin(), 
        [&](const Real x, const Real y) { return (2.0 - omega) * x * y; });
    AXPY(-1.0, p_0, p_1);
    //
    for(I64 ti = numbX - 1; ti >= 0; -- ti){
        I64 tj_start = hierStriUppe.row_ptr[ti];
        I64 tj_end = hierStriUppe.row_ptr[ti + 1];
        alignas(nfsAlign) Real tempSumm = std::transform_reduce(std::execution::unseq, 
            hierStriUppe.col_ind + tj_start, hierStriUppe.col_ind + tj_end, 
            hierStriUppe.val + tj_start, 
            0.0, std::plus<Real>(),
            [&solution](const I64 c, const Real v) { return v * solution[c]; });
        solution[ti] = (p_1[ti] - tempSumm) / hierarchyDiagonal.val[ti];
    }
}

} // namespace Ddpca