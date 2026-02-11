#include "SparseMatrix.hpp"
#include "../Mesh/Coordinate.hpp"

#include <cassert>
#include <cstring>
#include <random>

namespace Ddpca {

// 将COO格式转换为CSR格式
void SparseMatrix::Coo2Csr(
    const I64 cooM, 
    const I64 cooN,
    std::vector<Triplet>& coo, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel) {
    //
    const I64 cooNnz = coo.size();
	Coo2Csr(cooM, cooN, coo.data(), cooNnz, threadTask, nestLevel);
}

void SparseMatrix::SortRemove(
    Triplet* coo, 
    const I64* orderCoo, 
    const std::vector<I64> &rowsInde, 
    const I64 maxRowLength, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel){
    //
    I64 rows = rowsInde.size() - 1;
    if(rows > 100000){
        Log("        SparseMatrix::SortRemove");
    }
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(rows, numbPart, startIndex, endIndex);
    //
    row_ptr = AlignedAllocate<I64>(rows + 1);
    row_ptr[0] = 0;
    std::vector<std::vector<I64>> indices_pre(
        numbPart, std::vector<I64>(maxRowLength));
    std::vector<std::vector<I64>> sorted_col_ind_pre(
        numbPart, std::vector<I64>(maxRowLength));
    std::vector<std::vector<Real>> sorted_val_pre(
        numbPart, std::vector<Real>(maxRowLength));
    //
    std::function<void(I64)> taskFunction_0 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                // 获取当前行的起始和结束位置
                I64 row_start = rowsInde[ti];
                I64 row_length = rowsInde[ti + 1] - row_start;
                // if(row_length <= 1) continue;
                // 创建索引数组，用于间接排序
                std::vector<I64>& indices = indices_pre[tp];
                std::iota(indices.begin(), indices.begin() + row_length, 0);
                std::sort(std::execution::unseq, 
                    indices.begin(), indices.begin() + row_length, 
                    [&](I64 a, I64 b) {
                        return coo[orderCoo[row_start + a]].col < coo[orderCoo[row_start + b]].col;
                    }
                );
                std::vector<I64>& sorted_col_ind = sorted_col_ind_pre[tp];
                std::vector<Real>& sorted_val = sorted_val_pre[tp];
                I64 row_size = -1;
                I64 col_real;
                for(I64 tj = 0; tj < row_length; ++ tj){
                    I64 tj_real = orderCoo[row_start + indices[tj]];
                    const Triplet& coo_tj = coo[tj_real];
                    if(tj == 0 || coo_tj.col != col_real){
                        ++ row_size;
                        col_real = coo_tj.col;
                        sorted_col_ind[row_size] = col_real;
                        sorted_val[row_size] = coo_tj.val;
                    }
                    else{
                        sorted_val[row_size] += coo_tj.val;
                    }
                }
                ++ row_size;
                for(I64 tj = 0; tj < row_size; ++ tj){
                    Triplet& coo_tj = coo[orderCoo[row_start + tj]];
                    coo_tj.col = sorted_col_ind[tj];
                    coo_tj.val = sorted_val[tj];
                }
                row_ptr[ti + 1] = row_size;
            }
        };
    switch(nestLevel){
        case 0: case 1: 
            threadManager.RunTask(nestLevel, threadTask, taskFunction_0);
            break;
        default:
            taskFunction_0(0);
            break;
    }
	std::inclusive_scan(std::execution::unseq, 
        row_ptr, row_ptr + rows + 1, row_ptr);
    col_ind = AlignedAllocate<I64>(row_ptr[rows]);
    val = AlignedAllocate<Real>(row_ptr[rows]);
    std::function<void(I64)> taskFunction_1 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                I64 tj_start = rowsInde[ti];
                I64 tj_start_real = row_ptr[ti];
                I64 tj_length = row_ptr[ti + 1] - tj_start_real;
                for(I64 tj = 0; tj < tj_length; ++ tj){
                    const Triplet& coo_tj = coo[orderCoo[tj_start + tj]];
                    const I64 tj_real = tj_start_real + tj;
                    col_ind[tj_real] = coo_tj.col;
                    val[tj_real] = coo_tj.val;
                }

            }
        };
    switch(nestLevel){
        case 0: case 1: 
            threadManager.RunTask(nestLevel, threadTask, taskFunction_1);
            break;
        default:
            taskFunction_1(0);
            break;
    }
}

// 将COO格式转换为CSR格式
void SparseMatrix::Coo2Csr(
    const I64 cooM, 
    const I64 cooN,
    Triplet* coo, 
    const I64 cooNnz, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel) {
    //
    if(cooM > 100000){
        Log("        SparseMatrix::Coo2Csr");
    }

	//with atomic???
	std::vector<I64> rowsInde(1 + cooM, 0);
	for(I64 ti = 0; ti < cooNnz; ++ ti){
		I64 row_ti = 1 + coo[ti].row;
		++ rowsInde[row_ti];
	}
    I64 maxRowLength = *std::max_element(std::execution::unseq, 
        rowsInde.begin(), rowsInde.end());
	std::inclusive_scan(std::execution::unseq, 
        rowsInde.begin(), rowsInde.end(), rowsInde.begin());
    //
    I64* orderCoo = AlignedAllocate<I64>(rowsInde[cooM]);
    std::vector<I64> rowsInde_ = rowsInde;
    for(I64 ti = 0; ti < cooNnz; ++ ti){
        I64& rowsInde_ti = rowsInde_[(coo[ti]).row]; //careful!
        orderCoo[rowsInde_ti] = ti;
        ++ rowsInde_ti;
    }
    SortRemove(coo, orderCoo, rowsInde, maxRowLength, threadTask, nestLevel);
    Deallocate<I64>(orderCoo, rowsInde[cooM]);

    M = cooM;
    N = cooN;
    nnz = row_ptr[M]; // != cooNnz
    sorted = true;
}

void SparseMatrix::Sort(
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    if(sorted) return;
    if(M > 100000) Log("        SparseMatrix::Sort");
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(M, numbPart, startIndex, endIndex);
    //
    // currently not used
    // improvement: move "indices/sorted_col_ind/sorted_val" out of the loop
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                // 获取当前行的起始和结束位置
                I64 row_start = row_ptr[ti];
                I64 row_end = row_ptr[ti + 1];
                I64 row_length = row_end - row_start;
                if(row_length <= 1) continue;
                
                // 创建索引数组，用于间接排序
                std::vector<I64> indices(row_length);
                std::iota(indices.begin(), indices.end(), 0);
                std::sort(indices.begin(), indices.end(), 
                    [&](I64 a, I64 b) {
                        return col_ind[row_start + a] < col_ind[row_start + b];
                    }
                );
                std::vector<I64> sorted_col_ind(row_length);
                std::vector<Real> sorted_val(row_length);
                for(I64 tj = 0; tj < row_length; ++ tj){
                    I64 tj_real = row_start + indices[tj];
                    sorted_col_ind[tj] = col_ind[tj_real];
                    sorted_val[tj] = val[tj_real];
                }
                
                // 将排序后的结果复制回原始数组
                std::copy(std::execution::unseq, 
                    sorted_col_ind.begin(), sorted_col_ind.end(), col_ind + row_start);
                std::copy(std::execution::unseq, 
                    sorted_val.begin(), sorted_val.end(), val + row_start);
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
    sorted = true;
}

// 输出矩阵信息
void SparseMatrix::Output(std::ostream& os) const {
    os << "SparseMatrix dimensions: " << M << "x" << N << ", nnz: " << nnz << std::endl;
    
    os << "CSR format data:" << std::endl;
    os << "row_ptr: ";
    I64 tempM = M;
    if(M > 10) tempM = 11;
    for (I64 i = 0; i <= tempM; ++ i) {
        os << row_ptr[i] << " ";
    }
    os << std::endl;
    
    os << "col_ind: ";
    I64 tempNnz = nnz;
    if(M > 10) tempNnz = row_ptr[11];
    for (I64 i = 0; i < tempNnz; ++ i) {
        os << col_ind[i] << " ";
    }
    os << std::endl;
    
    os << "val: ";
    for (I64 i = 0; i < tempNnz; ++ i) {
        os << val[i] << " ";
    }
    os << std::endl;
}

Real SparseMatrix::SpectralRadius(
    const I64 power_iters, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel) const{
    //amgcl/backend/builtin.hpp
    Real radius = 0.0;
    AlignedVectorRx b0(M), b1(M);
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(M, numbPart, startIndex, endIndex);
    //
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            std::mt19937 rng(tp);
            std::uniform_real_distribution<Real> rnd(-1.0, 1.0);
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                Real v = rnd(rng);
                b0[ti] = v;
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
    Real b0_norm = 1.0 / NRM2(b0);
    SCAL(b0_norm, b0);
    //
    for(I64 iter = 0; iter < power_iters;){
        MV(*this, b0, b1, threadTask, nestLevel);
        if(++iter < power_iters){
            std::copy(std::execution::unseq, b1.begin(), b1.end(), b0.begin());
            Real b1_norm = 1.0 / NRM2(b1);
            SCAL(b1_norm, b0);
        }
        else{
            radius = DOT(b0, b1);
        }
    }
    return radius;
}

CS * SparseMatrix::ToCSparse() const {
	//
	CS *outpCspa;
	I64 nz = nnz, p, *w, *Cp, *Ci;
	Real *Cx;
	//
	outpCspa = CSSpalloc (M, N, nz, 1, 0) ;  /* allocate result */
    w = (I64*)CSCalloc (N, sizeof (I64)) ;      /* get workspace */
    if (!outpCspa || !w){
		CSDone (outpCspa, w, NULL, 0);             /* out of memory */
		std::cout << "ERROR @SparMatr::TO_CSPARSE out of memory" << std::endl;
	}
    Cp = outpCspa->p ; Ci = outpCspa->i ; Cx = outpCspa->x ;
    for(I64 ti = 0; ti < M; ++ ti){
        I64 tj_start = row_ptr[ti];
        I64 tj_end = row_ptr[ti + 1];
        for(I64 tj = tj_start; tj < tj_end; tj ++){
            w[col_ind[tj]] ++;                   /* column counts */
        }
    }
    CSCumsum (Cp, w, N) ;                       /* column pointers */
	for(I64 ti = 0; ti < M; ti ++){
        I64 tj_start = row_ptr[ti];
        I64 tj_end = row_ptr[ti + 1];
        for(I64 tj = tj_start; tj < tj_end; tj ++){
        	Ci [p = w[col_ind[tj]]++] = ti ;    /* A(i,j) is the pth entry in C */
			Cx [p] = val[tj] ;
        }
	}
	CSDone (outpCspa, w, NULL, 1);                 /* success; free w and return C */
	CSDropzeros (outpCspa) ;                       /* drop zero entries */
	return outpCspa;
}

void SparseMatrix::KeepNPAdd(
    const SparseMatrix& in, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    assert(in.M == M && "in.M != M");
    assert(in.N == N && "in.N != N");
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(in.M, numbPart, startIndex, endIndex);
    //
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                I64 tj_start = in.row_ptr[ti];
                I64 tj_end = in.row_ptr[ti + 1];
                for(I64 tj = tj_start; tj < tj_end; ++ tj){
                    I64 col_tj = in.col_ind[tj];
                    Real val_tj = in.val[tj];
                    //
                    I64* itCol = std::lower_bound(
                        col_ind + row_ptr[ti], col_ind + row_ptr[ti + 1], col_tj);
                    assert(*itCol == col_tj && "col_tj not found");
                    I64 indexReal = itCol - col_ind;
                    val[indexReal] += val_tj;
                }
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

void TRANSPOSE(const SparseMatrix& in, SparseMatrix& out){
    if(in.M > 100000) Log("        TRANSPOSE");
    //
    out.M = in.N;
    out.N = in.M;
    out.nnz = in.nnz;
    out.row_ptr = AlignedAllocate<I64>(out.M + 1); // careful!
    out.col_ind = AlignedAllocate<I64>(out.nnz);
    out.val = AlignedAllocate<Real>(out.nnz);
    out.sorted = true;

    std::fill(std::execution::unseq, out.row_ptr, out.row_ptr + out.M + 1, 0);
    for(I64 ti = 0; ti < in.M; ++ ti){
        I64 tj_start = in.row_ptr[ti];
        I64 tj_end = in.row_ptr[ti + 1];
        // col_ind[tj_start~tj_end] has no duplication
        for(I64 tj = tj_start; tj < tj_end; ++ tj){
            ++ out.row_ptr[1 + in.col_ind[tj]];
        }
    }
	std::inclusive_scan(std::execution::unseq, 
        out.row_ptr, out.row_ptr + out.M + 1, out.row_ptr);
    I64* row_ptr_ = AlignedAllocate<I64>(out.M);
    std::copy(std::execution::unseq, out.row_ptr, out.row_ptr + out.M, row_ptr_);
    for(I64 ti = 0; ti < in.M; ++ ti){
        I64 tj_start = in.row_ptr[ti];
        I64 tj_end = in.row_ptr[ti + 1];
        for(I64 tj = tj_start; tj < tj_end; ++ tj){
            I64 row_real = in.col_ind[tj];
            I64& tj_real = row_ptr_[row_real]; // careful!
            out.col_ind[tj_real] = ti;
            out.val[tj_real] = in.val[tj];
            ++ tj_real;
        }
    }
    Deallocate<I64>(row_ptr_, out.M);
}

void MV(
    const SparseMatrix& A, 
    const AlignedVectorRx& x, 
    AlignedVectorRx& y, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel){
    //
    const I64 M = A.M;
    const I64* row_ptr = A.row_ptr;
    const I64* col_ind = A.col_ind;
    const Real* val = A.val;
    //
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(M, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                I64 tj_start = row_ptr[ti];
                I64 tj_end = row_ptr[ti + 1];
                // Real y_local = 0.0;
                // for(I64 tj = tj_start; tj < tj_end; ++ tj){
                //     y_local += val[tj] * x[col_ind[tj]];
                // }
                Real y_local = std::transform_reduce(std::execution::unseq, 
                    col_ind + tj_start, col_ind + tj_end, 
                    val + tj_start, 
                    0.0, std::plus<Real>(),
                    [&x](const I64 c, const Real v) { return v * x[c]; });
                y[ti] = y_local;
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

void PEMV(
    const SparseMatrix& A, 
    const AlignedVectorRx& x, 
    AlignedVectorRx& y, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel){
    //
    const I64 M = A.M;
    const I64* row_ptr = A.row_ptr;
    const I64* col_ind = A.col_ind;
    const Real* val = A.val;
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(M, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                I64 tj_start = row_ptr[ti];
                I64 tj_end = row_ptr[ti + 1];
                Real y_local = y[ti];
                // for(I64 tj = tj_start; tj < tj_end; ++ tj){
                //     y_local += val[tj] * x[col_ind[tj]];
                // }
                y_local += std::transform_reduce(std::execution::unseq, 
                    col_ind + tj_start, col_ind + tj_end, 
                    val + tj_start, 
                    0.0, std::plus<Real>(),
                    [&x](const I64 c, const Real v) { return v * x[c]; });
                y[ti] = y_local;
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

void PEMV(
    const Real alpha, 
    const SparseMatrix& A, 
    const AlignedVectorRx& x, 
    AlignedVectorRx& y, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel){
    //
    const I64 M = A.M;
    const I64* row_ptr = A.row_ptr;
    const I64* col_ind = A.col_ind;
    const Real* val = A.val;
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(M, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                I64 tj_start = row_ptr[ti];
                I64 tj_end = row_ptr[ti + 1];
                Real y_local = y[ti];
                // for(I64 tj = tj_start; tj < tj_end; ++ tj){
                //     y_local += val[tj] * x[col_ind[tj]];
                // }
                y_local += alpha * std::transform_reduce(std::execution::unseq, 
                    col_ind + tj_start, col_ind + tj_end, 
                    val + tj_start, 
                    0.0, std::plus<Real>(),
                    [&x](const I64 c, const Real v) { return v * x[c]; });
                y[ti] = y_local;
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

void MV(
    const Real alpha, 
    const SparseMatrix& A, 
    const AlignedVectorRx& x, 
    const Real beta, 
    const AlignedVectorRx& y, 
    AlignedVectorRx& z, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel){
    //
    const I64 M = A.M;
    const I64* row_ptr = A.row_ptr;
    const I64* col_ind = A.col_ind;
    const Real* val = A.val;
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(M, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                I64 tj_start = row_ptr[ti];
                I64 tj_end = row_ptr[ti + 1];
                Real local_y = y[ti];
                Real local_r = 0.0;
                // __m512d sum = _mm512_setzero_pd();
                // I64 tj_aligned = tj_start + ((tj_end - tj_start) / 8) * 8;
                // for(I64 tj = tj_start; tj < tj_aligned; tj += 8) {
                //     __m512d vals = _mm512_loadu_pd(val + tj);
                //     __m512d xs = _mm512_set_pd(
                //         x[col_ind[tj+7]], x[col_ind[tj+6]], x[col_ind[tj+5]], x[col_ind[tj+4]],
                //         x[col_ind[tj+3]], x[col_ind[tj+2]], x[col_ind[tj+1]], x[col_ind[tj]]
                //     );
                //     __m512d prod = _mm512_mul_pd(vals, xs);
                //     sum = _mm512_add_pd(sum, prod);
                // }
                // local_r = _mm512_reduce_add_pd(sum);
                // for(I64 tj = tj_aligned; tj < tj_end; ++tj) {
                //     local_r += val[tj] * x[col_ind[tj]];
                // }
                for(I64 tj = tj_start; tj < tj_end; ++ tj){
                    local_r += val[tj] * x[col_ind[tj]];
                }
                // Real local_r = std::transform_reduce(std::execution::unseq, 
                //     col_ind + tj_start, col_ind + tj_end, 
                //     val + tj_start, 
                //     0.0, std::plus<Real>(),
                //     [&x](const I64 c, const Real v) { return v * x[c]; });
                Real local_z = beta * local_y + alpha * local_r;
                z[ti] = local_z;
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

void GEMM(
    const SparseMatrix& A, 
    const SparseMatrix& B, 
    SparseMatrix& C, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel){
    //
    if(A.M > 100000) Log("        GEMM");
    //
    C.M = A.M;
    C.N = B.N;
    C.row_ptr = AlignedAllocate<I64>(C.M + 1);
    C.row_ptr[0] = 0;
    std::vector<std::vector<std::pair<I64, Real>>> data(C.M);
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(C.M, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction_0 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 tr = start_tp; tr < end_tp; ++ tr){
                I64 tj_start_A = A.row_ptr[tr];
                I64 tj_end_A = A.row_ptr[tr + 1];
                auto& data_tr = data[tr]; //careful!
                // data_tr.reserve(numbNonzPerRow);
                for(I64 tj = tj_start_A; tj < tj_end_A; ++ tj){
                    I64 col_tj_A = A.col_ind[tj];
                    Real val_tj_A = A.val[tj];
                    I64 tk_start_B = B.row_ptr[col_tj_A];
                    I64 tk_end_B = B.row_ptr[col_tj_A + 1];
                    for(I64 tk = tk_start_B; tk < tk_end_B; ++ tk){
                        I64 col_tk_B = B.col_ind[tk];
                        Real val_tk = val_tj_A * B.val[tk];
                        auto iteratorDt = std::lower_bound(
                            data_tr.begin(), data_tr.end(), 
                            std::make_pair(col_tk_B, Real(0)),
                            [](const auto& a, const auto& b) {
                                return a.first < b.first;
                            });
                        if(iteratorDt != data_tr.end() && iteratorDt->first == col_tk_B){
                            // 找到现有元素，累加值
                            iteratorDt->second += val_tk;
                        } else {
                            // 插入新元素并保持有序
                            data_tr.emplace(iteratorDt, col_tk_B, val_tk);
                        }
                    }
                }
                C.row_ptr[1 + tr] = data_tr.size();
            }
        };
    switch(nestLevel){
        case 0: case 1: 
            threadManager.RunTask(nestLevel, threadTask, taskFunction_0);
            break;
        default:
            taskFunction_0(0);
            break;
    }
	std::inclusive_scan(std::execution::unseq, 
        C.row_ptr, C.row_ptr + C.M + 1, C.row_ptr);
    //
    C.nnz = C.row_ptr[C.M];
    C.col_ind = AlignedAllocate<I64>(C.nnz);
    C.val = AlignedAllocate<Real>(C.nnz);
    C.sorted = true;
    std::function<void(I64)> taskFunction_1 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 tr = start_tp; tr < end_tp; ++ tr){
                I64 tr_start = C.row_ptr[tr];
                I64 tr_length = C.row_ptr[tr + 1] - tr_start;
                auto& data_tr = data[tr]; //careful!
                for(I64 tj = 0; tj < tr_length; ++ tj){
                    I64 tj_real = tr_start + tj;
                    const auto& data_tj = data_tr[tj]; //careful!
                    C.col_ind[tj_real] = data_tj.first;
                    C.val[tj_real] = data_tj.second;
                }
            }
        };
    switch(nestLevel){
        case 0: case 1: 
            threadManager.RunTask(nestLevel, threadTask, taskFunction_1);
            break;
        default:
            taskFunction_1(0);
            break;
    }
}

}// namespace Ddpca