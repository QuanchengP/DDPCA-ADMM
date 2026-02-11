#ifndef _RSYMGS_hpp
#define _RSYMGS_hpp
#if DDPCA_USE_GPU == 1
#include <cusparse.h>         // cusparseSpSV
#endif

#include "../General/SparseMatrix.hpp"

namespace Ddpca{

//relaxation: symmetric Gauss-Seidel
class RSYMGS{
public:
    
    //Cache-friendly: help much
    SparseMatrix hierStriLowe;//1 ~ maxLevel
    SparseMatrix hierStriUppe;//1 ~ maxLevel
    SparseMatrix hierarchyDiagonal;//1 ~ maxLevel

    // forward parallel Gauss-Seidel
    std::vector<std::vector<I64>> forwardLayerRows;
    // backward parallel Gauss-Seidel
    std::vector<std::vector<I64>> backwardLayerRows;
    // Cache-friendly: doesn‘t help much
    std::vector<std::vector<I64>> forwLayeRowsRow_ptr;
    std::vector<std::vector<I64>> forwLayeRowsCol_ind;
    std::vector<std::vector<Real>> forwLayeRowsVal;
    std::vector<std::vector<I64>> backLayeRowsRow_ptr;
    std::vector<std::vector<I64>> backLayeRowsCol_ind;
    std::vector<std::vector<Real>> backLayeRowsVal;

    AlignedVectorRx p0;
    AlignedVectorRx bp0;
    AlignedVectorRx x_1;

public:

    void Establish(const SparseMatrix& hierarchyStiffness);
    void HierarchySplit(const SparseMatrix& hierarchyStiffness); // split rows into layers

    void Apply(
        const AlignedVectorRx& rhs, 
        AlignedVectorRx& solution, 
        AlignedVectorRx& y, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

    void Apply(
        const AlignedVectorRx& rhs, 
        AlignedVectorRx& solution, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

public:

    #if DDPCA_USE_GPU == 1
        //--------------------------------------------------------------------------
        // Device memory management
        I64 A_num_rows;
        I64 A_num_cols;
        I64 A_nnz;
        I64 * dA_csrOffsets;
        I64 * dA_columns;
        Real * dA_values;
        Real * dX;
        Real * dY;
        void* dBufferL;
        void* dBufferU;
        cusparseHandle_t handleL;
        cusparseHandle_t handleU;
        cusparseSpMatDescr_t matAL;
        cusparseSpMatDescr_t matAU;
        cusparseSpSVDescr_t spsvDescrL; 
        cusparseSpSVDescr_t spsvDescrU; 
        double cudaAlpha;
    #endif

    ~RSYMGS(){
        #if DDPCA_USE_GPU == 1
            CHECK_CUDA( cudaFree(dA_csrOffsets) )
            CHECK_CUDA( cudaFree(dA_columns) )
            CHECK_CUDA( cudaFree(dA_values) )
            CHECK_CUDA( cudaFree(dX) )
            CHECK_CUDA( cudaFree(dY) )
            // device memory deallocation
            CHECK_CUDA( cudaFree(dBufferL) )
            CHECK_CUDA( cudaFree(dBufferU) )
            CHECK_CUSPARSE( cusparseDestroy(handleL) )
            CHECK_CUSPARSE( cusparseDestroy(handleU) )
            CHECK_CUSPARSE( cusparseDestroySpMat(matAL) )
            CHECK_CUSPARSE( cusparseDestroySpMat(matAU) )
            CHECK_CUSPARSE( cusparseSpSV_destroyDescr(spsvDescrL) )
            CHECK_CUSPARSE( cusparseSpSV_destroyDescr(spsvDescrU) )
        #endif
    }

}; // class RSYMGS

} // namespace Ddpca

#endif // _RSYMGS_hpp