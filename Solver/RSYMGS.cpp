#include "RSYMGS.hpp"
#include "../Mesh/Coordinate.hpp"

namespace Ddpca{

void RSYMGS::Establish(const SparseMatrix& hierarchyStiffness){
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
                hierStriLowe.val[tempLnnz] = val_tj;
                ++ tempLnnz;
            }
            else if(ti < col_tj){
                hierStriUppe.col_ind[tempUnnz] = col_tj;
                hierStriUppe.val[tempUnnz] = val_tj;
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
    p0.resize(hierarchyStiffness.M);
    bp0.resize(hierarchyStiffness.M);
    x_1.resize(hierarchyStiffness.M);

    #if DDPCA_USE_GPU == 1
        cudaAlpha = 1.0;
        A_num_rows      = hierarchyStiffness.M;
        A_num_cols      = hierarchyStiffness.N;
        A_nnz           = hierarchyStiffness.nnz;
        CHECK_CUDA( cudaMalloc((void**) &(dA_csrOffsets),
                                (A_num_rows + 1) * sizeof(I64)) )
        CHECK_CUDA( cudaMalloc((void**) &(dA_columns), A_nnz * sizeof(I64))        )
        CHECK_CUDA( cudaMalloc((void**) &(dA_values),  A_nnz * sizeof(Real))      )
        CHECK_CUDA( cudaMalloc((void**) &(dX),         A_num_cols * sizeof(Real)) )
        CHECK_CUDA( cudaMalloc((void**) &(dY),         A_num_rows * sizeof(Real)) )

        CHECK_CUDA( cudaMemcpy(dA_csrOffsets, hierarchyStiffness.row_ptr,
                                (A_num_rows + 1) * sizeof(I64),
                                cudaMemcpyHostToDevice) )
        CHECK_CUDA( cudaMemcpy(dA_columns, hierarchyStiffness.col_ind, A_nnz * sizeof(I64),
                                cudaMemcpyHostToDevice) )
        CHECK_CUDA( cudaMemcpy(dA_values, hierarchyStiffness.val, A_nnz * sizeof(Real),
                                cudaMemcpyHostToDevice) )
        dBufferL = NULL;
        dBufferU = NULL;
        CHECK_CUSPARSE( cusparseCreate(&handleL) )
        CHECK_CUSPARSE( cusparseCreate(&handleU) )
        // Create opaque data structure, that holds analysis data between calls.
        CHECK_CUSPARSE( cusparseSpSV_createDescr(&spsvDescrL) )
        CHECK_CUSPARSE( cusparseSpSV_createDescr(&spsvDescrU) )

        // Create sparse matrix A in CSR format
        CHECK_CUSPARSE( cusparseCreateCsr(&matAL, A_num_rows, A_num_cols, A_nnz,
                                            dA_csrOffsets, dA_columns, dA_values,
                                            CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I,
                                            CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )
        // Specify Lower|Upper fill mode.
        cusparseFillMode_t fillmode = CUSPARSE_FILL_MODE_LOWER;
        CHECK_CUSPARSE( cusparseSpMatSetAttribute(matAL, CUSPARSE_SPMAT_FILL_MODE,
                                                    &fillmode, sizeof(fillmode)) )
        // Specify Unit|Non-Unit diagonal type.
        cusparseDiagType_t diagtype = CUSPARSE_DIAG_TYPE_NON_UNIT;
        CHECK_CUSPARSE( cusparseSpMatSetAttribute(matAL, CUSPARSE_SPMAT_DIAG_TYPE,
                                                    &diagtype, sizeof(diagtype)) )
        size_t               bufferSize = 0;
        // allocate an external buffer for analysis
        CHECK_CUSPARSE( cusparseSpSV_bufferSize(
                                    handleL, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    &cudaAlpha, matAL, NULL, NULL, CUDA_R_64F,
                                    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL,
                                    &bufferSize) )
        std::cout << "bufferSize = " << bufferSize << "\n";
        CHECK_CUDA( cudaMalloc(&dBufferL, bufferSize) )
        CHECK_CUSPARSE( cusparseSpSV_analysis(
                                    handleL, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    &cudaAlpha, matAL, NULL, NULL, CUDA_R_64F,
                                    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, dBufferL) )


        // Create sparse matrix A in CSR format
        CHECK_CUSPARSE( cusparseCreateCsr(&matAU, A_num_rows, A_num_cols, A_nnz,
                                            dA_csrOffsets, dA_columns, dA_values,
                                            CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I,
                                            CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )
        // Specify Lower|Upper fill mode.
        fillmode = CUSPARSE_FILL_MODE_UPPER;
        CHECK_CUSPARSE( cusparseSpMatSetAttribute(matAU, CUSPARSE_SPMAT_FILL_MODE,
                                                    &fillmode, sizeof(fillmode)) )
        // Specify Unit|Non-Unit diagonal type.
        diagtype = CUSPARSE_DIAG_TYPE_NON_UNIT;
        CHECK_CUSPARSE( cusparseSpMatSetAttribute(matAU, CUSPARSE_SPMAT_DIAG_TYPE,
                                                    &diagtype, sizeof(diagtype)) )
        bufferSize = 0;
        // allocate an external buffer for analysis
        CHECK_CUSPARSE( cusparseSpSV_bufferSize(
                                    handleU, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    &cudaAlpha, matAU, NULL, NULL, CUDA_R_64F,
                                    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrU,
                                    &bufferSize) )
        std::cout << "bufferSize = " << bufferSize << "\n";
        CHECK_CUDA( cudaMalloc(&dBufferU, bufferSize) )
        CHECK_CUSPARSE( cusparseSpSV_analysis(
                                    handleU, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    &cudaAlpha, matAU, NULL, NULL, CUDA_R_64F,
                                    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrU, dBufferU) )
    #endif
}

void RSYMGS::HierarchySplit(const SparseMatrix& hierarchyStiffness){
    // #pragma omp parallel num_threads(2){#pragma omp single nowait{#pragma omp task
            {
                const I64& M = hierStriLowe.M;
                I64 numbLayer = 0;
                std::vector<I64> rowLayer(M, 0);
                for(I64 ti = 0; ti < M; ++ ti){
                    I64 tempRowLayer = rowLayer[ti];
                    I64 tj_start = hierStriLowe.row_ptr[ti];
                    I64 tj_end = hierStriLowe.row_ptr[ti + 1];
                    for(I64 tj = tj_start; tj < tj_end; ++ tj){
                        I64 col_tj = hierStriLowe.col_ind[tj];
                        tempRowLayer = std::max(tempRowLayer, rowLayer[col_tj] + 1);
                    }
                    rowLayer[ti] = tempRowLayer;
                    numbLayer = std::max(numbLayer, tempRowLayer);
                }
                forwardLayerRows.resize(numbLayer + 1);
                for(I64 ti = 0; ti < M; ++ ti){
                    forwardLayerRows[rowLayer[ti]].emplace_back(ti);
                }
                forwLayeRowsRow_ptr.resize(numbLayer + 1);
                forwLayeRowsCol_ind.resize(numbLayer + 1);
                forwLayeRowsVal.resize(numbLayer + 1);
                for(I64 tla = 0; tla <= numbLayer; ++ tla){
                    const I64 tla_size = forwardLayerRows[tla].size();
                    std::vector<I64>& tempLayeRowsRow_ptr = forwLayeRowsRow_ptr[tla];
                    std::vector<I64>& tempLayeRowsCol_ind = forwLayeRowsCol_ind[tla];
                    std::vector<Real>& tempLayeRowsVal = forwLayeRowsVal[tla];
                    tempLayeRowsRow_ptr.resize(tla_size + 1);
                    tempLayeRowsRow_ptr[0] = 0;
                    for(I64 tia = 0; tia < tla_size; ++ tia){
                        const I64 ti_real = forwardLayerRows[tla][tia];
                        tempLayeRowsRow_ptr[tia + 1] = tempLayeRowsRow_ptr[tia] 
                            + (hierarchyStiffness.row_ptr[ti_real + 1] 
                            - hierarchyStiffness.row_ptr[ti_real]);
                        tempLayeRowsCol_ind.insert(
                            tempLayeRowsCol_ind.end(),
                            hierarchyStiffness.col_ind + hierarchyStiffness.row_ptr[ti_real],
                            hierarchyStiffness.col_ind + hierarchyStiffness.row_ptr[ti_real + 1]);
                        tempLayeRowsVal.insert(
                            tempLayeRowsVal.end(),
                            hierarchyStiffness.val + hierarchyStiffness.row_ptr[ti_real],
                            hierarchyStiffness.val + hierarchyStiffness.row_ptr[ti_real + 1]);
                    }
                }
            }
            // #pragma omp task
            {
                const I64& M = hierStriUppe.M;
                I64 numbLayer = 0;
                std::vector<I64> rowLayer(M, 0);
                for(I64 ti = M - 1; ti >= 0; -- ti){
                    I64 tempRowLayer = rowLayer[ti];
                    I64 tj_start = hierStriUppe.row_ptr[ti];
                    I64 tj_end = hierStriUppe.row_ptr[ti + 1];
                    for(I64 tj = tj_start; tj < tj_end; ++ tj){
                        I64 col_tj = hierStriUppe.col_ind[tj];
                        tempRowLayer = std::max(tempRowLayer, rowLayer[col_tj] + 1);
                    }
                    rowLayer[ti] = tempRowLayer;
                    numbLayer = std::max(numbLayer, tempRowLayer);
                }
                backwardLayerRows.resize(numbLayer + 1);
                for(I64 ti = M - 1; ti >= 0; -- ti){
                    backwardLayerRows[rowLayer[ti]].emplace_back(ti);
                }
                backLayeRowsRow_ptr.resize(numbLayer + 1);
                backLayeRowsCol_ind.resize(numbLayer + 1);
                backLayeRowsVal.resize(numbLayer + 1);
                for(I64 tla = 0; tla <= numbLayer; ++ tla){
                    const I64 tla_size = backwardLayerRows[tla].size();
                    std::vector<I64>& tempLayeRowsRow_ptr = backLayeRowsRow_ptr[tla];
                    std::vector<I64>& tempLayeRowsCol_ind = backLayeRowsCol_ind[tla];
                    std::vector<Real>& tempLayeRowsVal = backLayeRowsVal[tla];
                    tempLayeRowsRow_ptr.resize(tla_size + 1);
                    tempLayeRowsRow_ptr[0] = 0;
                    for(I64 tia = 0; tia < tla_size; ++ tia){
                        const I64 ti_real = backwardLayerRows[tla][tia];
                        tempLayeRowsRow_ptr[tia + 1] = tempLayeRowsRow_ptr[tia] 
                            + (hierarchyStiffness.row_ptr[ti_real + 1] 
                            - hierarchyStiffness.row_ptr[ti_real]);
                        tempLayeRowsCol_ind.insert(
                            tempLayeRowsCol_ind.end(),
                            hierarchyStiffness.col_ind + hierarchyStiffness.row_ptr[ti_real],
                            hierarchyStiffness.col_ind + hierarchyStiffness.row_ptr[ti_real + 1]);
                        tempLayeRowsVal.insert(
                            tempLayeRowsVal.end(),
                            hierarchyStiffness.val + hierarchyStiffness.row_ptr[ti_real],
                            hierarchyStiffness.val + hierarchyStiffness.row_ptr[ti_real + 1]);
                    }
                }
            }
    // }}
}

void RSYMGS::Apply(
    const AlignedVectorRx& rhs, 
    AlignedVectorRx& solution, 
    AlignedVectorRx& y, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
	//Yang, X, Li, S, Yuan, F et al. (2023) Optimizing Multi-grid Computation
	//and Parallelization on Multi-cores.
    I64 numbX = rhs.size();

    //
    MV(-1.0, hierStriUppe, solution, 0.0, p0, p0, threadTask, nestLevel);
    //
    std::copy(std::execution::unseq, rhs.begin(), rhs.end(), bp0.begin());
    XPEY(bp0, p0);
    for(I64 ti = 0; ti < numbX; ++ ti){
        I64 tj_start = hierStriLowe.row_ptr[ti];
        I64 tj_end = hierStriLowe.row_ptr[ti + 1];
        alignas(nfsAlign) Real tempSumm = std::transform_reduce(std::execution::unseq, 
            hierStriLowe.col_ind + tj_start, hierStriLowe.col_ind + tj_end, 
            hierStriLowe.val + tj_start, 
            0.0, std::plus<Real>(),
            [this](const I64 c, const Real v) { return v * (this->x_1)[c]; });
        x_1[ti] = (bp0[ti] - tempSumm) / hierarchyDiagonal.val[ti];
    }
    //y is p_1
    MV(1.0, hierarchyDiagonal, x_1, -1.0, p0, y, threadTask, nestLevel);
    //
    for(I64 ti = numbX - 1; ti >= 0; -- ti){
        I64 tj_start = hierStriUppe.row_ptr[ti];
        I64 tj_end = hierStriUppe.row_ptr[ti + 1];
        alignas(nfsAlign) Real tempSumm = std::transform_reduce(std::execution::unseq, 
            hierStriUppe.col_ind + tj_start, hierStriUppe.col_ind + tj_end, 
            hierStriUppe.val + tj_start, 
            0.0, std::plus<Real>(),
            [&solution](const I64 c, const Real v) { return v * solution[c]; });
        solution[ti] = (y[ti] - tempSumm) / hierarchyDiagonal.val[ti];
    }
    //
    PEMV(hierStriLowe, solution, y, threadTask, nestLevel);
}

void RSYMGS::Apply(
    const AlignedVectorRx& rhs, 
    AlignedVectorRx& solution, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //amgcl/relaxation/gauss_seidel.hpp
    //
    for(I64 fb = 0; fb < 2; ++ fb){
        const std::vector<std::vector<I64>>& layerRows = 
            (fb == 0) ? forwardLayerRows : backwardLayerRows;
        I64 numbLayer = layerRows.size();
        for(I64 tl = 0; tl < numbLayer; ++ tl){
            const std::vector<I64>& tempRows = layerRows[tl];
            const I64 numbRows = tempRows.size();
            const std::vector<I64>& layerRowsRow_ptr = 
                (fb == 0) ? forwLayeRowsRow_ptr[tl] : backLayeRowsRow_ptr[tl];
            const std::vector<I64>& layerRowsCol_ind = 
                (fb == 0) ? forwLayeRowsCol_ind[tl] : backLayeRowsCol_ind[tl];
            const std::vector<Real>& layerRowsVal = 
                (fb == 0) ? forwLayeRowsVal[tl] : backLayeRowsVal[tl];
            //
            I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
            std::vector<I64> startIndex(numbPart), endIndex(numbPart);
            EvenlyDistribute(numbRows, numbPart, startIndex, endIndex);
            std::function<void(I64)> taskFunction = 
                [&](I64 tp){
                    I64 start_tp = startIndex[tp];
                    I64 end_tp = endIndex[tp];
                    for(I64 ti = start_tp; ti < end_tp; ++ ti){
                        double tempSumm = std::transform_reduce(std::execution::unseq,
                            layerRowsCol_ind.begin() + layerRowsRow_ptr[ti], 
                            layerRowsCol_ind.begin() + layerRowsRow_ptr[ti + 1],
                            layerRowsVal.begin() + layerRowsRow_ptr[ti],
                            0.0,
                            std::plus<>(),
                            [&solution](const I64 col, const double val) {return val * solution[col];});
                        // Real tempSumm = 0.0;
                        // for(I64 tj = layerRowsRow_ptr[ti]; tj < layerRowsRow_ptr[ti + 1]; ++ tj){
                        //     tempSumm += layerRowsVal[tj] * solution[layerRowsCol_ind[tj]];
                        // }
                        I64 ti_real = tempRows[ti];
                        solution[ti_real] += (rhs[ti_real] - tempSumm) / hierarchyDiagonal.val[ti_real];
                    }
                };
            switch(nestLevel){
                case 0: case 1: 
                    threadManager.RunTask(nestLevel, threadTask, taskFunction);
                    break;
                default:
                    taskFunction(0);
                    break;
            }
        }
    }
}

} // namespace Ddpca