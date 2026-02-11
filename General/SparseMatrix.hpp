#ifndef _SparseMatrix_hpp
#define _SparseMatrix_hpp
#include "../Solver/CSparse/CS.hpp"
#include "Triplet.hpp"
#include "AlignedVector.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <execution>

namespace Ddpca {

class SparseMatrix {
public:
    // 矩阵维度
    I64 M;
    I64 N;
    I64 nnz;
    
    // COO格式数据存储
    I64* row_ptr; // 行索引数组
    I64* col_ind; // sorted/unsorted column index
    Real* val;  // 非零元素值数组

    bool sorted; // the col_ind in the same row are sorted

public:
    // 构造函数
    SparseMatrix() : M(0), N(0), nnz(0), 
        row_ptr(nullptr), col_ind(nullptr), val(nullptr), 
        sorted(false) {}
    
    // 析构函数
    ~SparseMatrix(){
        if(nnz > 0){
            if(row_ptr != nullptr) Deallocate<I64>(row_ptr, M + 1);
            if(col_ind != nullptr) Deallocate<I64>(col_ind, nnz);
            if(val != nullptr) Deallocate<Real>(val, nnz);
            row_ptr = nullptr;
            col_ind = nullptr;
            val = nullptr;
        }
    }
    
    // 拷贝构造函数
    SparseMatrix(const SparseMatrix& other) : 
        M(other.M), 
        N(other.N), 
        nnz(other.nnz),
        row_ptr(nullptr),
        col_ind(nullptr),
        val(nullptr),
        sorted(other.sorted) {
        // 分配内存并复制数据
        if (nnz > 0) {
            row_ptr = AlignedAllocate<I64>(M + 1);  // 假设row_ptr长度为M+1（CSR格式）
            col_ind = AlignedAllocate<I64>(nnz);
            val = AlignedAllocate<Real>(nnz);
            
            // 复制数据
            std::copy(std::execution::unseq, other.row_ptr, other.row_ptr + M + 1, row_ptr);
            std::copy(std::execution::unseq, other.col_ind, other.col_ind + nnz, col_ind);
            std::copy(std::execution::unseq, other.val, other.val + nnz, val);
        }
    }

    // 移动构造函数
    SparseMatrix(SparseMatrix&& other) noexcept : 
        M(other.M), 
        N(other.N), 
        nnz(other.nnz),
        row_ptr(other.row_ptr),
        col_ind(other.col_ind),
        val(other.val),
        sorted(other.sorted) {
        // 将other的指针设为nullptr，避免析构函数释放资源
        other.M = 0;
        other.N = 0;
        other.nnz = 0;
        other.row_ptr = nullptr;
        other.col_ind = nullptr;
        other.val = nullptr;
        other.sorted = false;
    }

    // 拷贝赋值操作符
    SparseMatrix& operator=(const SparseMatrix& other) {
        if (this != &other) {
            if(nnz > 0){
                Deallocate<I64>(row_ptr, M + 1);
                Deallocate<I64>(col_ind, nnz);
                Deallocate<Real>(val, nnz);
                row_ptr = nullptr;
                col_ind = nullptr;
                val = nullptr;
            }
            // 复制数据
            M = other.M;
            N = other.N;
            nnz = other.nnz;
            sorted = other.sorted;
            if (nnz > 0) {
                row_ptr = AlignedAllocate<I64>(M + 1);
                col_ind = AlignedAllocate<I64>(nnz);
                val = AlignedAllocate<Real>(nnz);
                std::copy(std::execution::unseq, other.row_ptr, other.row_ptr + M + 1, row_ptr);
                std::copy(std::execution::unseq, other.col_ind, other.col_ind + nnz, col_ind);
                std::copy(std::execution::unseq, other.val, other.val + nnz, val);
            }
        }
        return *this;
    }

    // 移动赋值操作符
    SparseMatrix& operator=(SparseMatrix&& other) noexcept {
        if (this != &other) {
            if(nnz > 0){
                Deallocate<I64>(row_ptr, M + 1);
                Deallocate<I64>(col_ind, nnz);
                Deallocate<Real>(val, nnz);
                row_ptr = nullptr;
                col_ind = nullptr;
                val = nullptr;
            }
            // 移动数据
            M = other.M;
            N = other.N;
            nnz = other.nnz;
            row_ptr = other.row_ptr;
            col_ind = other.col_ind;
            val = other.val;
            sorted = other.sorted;
            
            // 将other的指针设为nullptr，避免析构函数释放资源
            other.M = 0;
            other.N = 0;
            other.nnz = 0;
            other.row_ptr = nullptr;
            other.col_ind = nullptr;
            other.val = nullptr;
            other.sorted = false;
        }
        return *this;
    }
    
public:

    void Coo2Csr(
        const I64 cooM, 
        const I64 cooN,
        std::vector<Triplet>& coo, 
        const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
        const I64 nestLevel);

    void SortRemove(
        Triplet* coo, 
        const I64* orderCoo, 
	    const std::vector<I64> &rowsInde, 
        const I64 maxRowLength, 
        const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
        const I64 nestLevel);

    void Coo2Csr(
        const I64 cooM, 
        const I64 cooN,
        Triplet* coo, 
        const I64 cooNnz, 
        const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
        const I64 nestLevel);

    void Sort(
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

    void Output(std::ostream& os = std::cout) const;

    Real SpectralRadius(
        const I64 power_iters, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel) const;

    CS * ToCSparse() const;

    //keep the non-zero pattern unchanged, matrix + matrix
    void KeepNPAdd(
        const SparseMatrix& in, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);
};

void TRANSPOSE(const SparseMatrix& in, SparseMatrix& out);

//y = A * x;
void MV(
    const SparseMatrix& A, 
    const AlignedVectorRx& x, 
    AlignedVectorRx& y, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel);

//y += A * x;
void PEMV(const SparseMatrix& A, const AlignedVectorRx& x, 
    AlignedVectorRx& y, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel);

//y += alpha * A * x;
void PEMV(const Real alpha, 
    const SparseMatrix& A, const AlignedVectorRx& x, 
    AlignedVectorRx& y, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel);

//z = alpha * A * x + beta * y;
void MV(
    const Real alpha, 
    const SparseMatrix& A, 
    const AlignedVectorRx& x, 
    const Real beta, 
    const AlignedVectorRx& y, 
    AlignedVectorRx& z, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel);

//C = A * B
void GEMM(
    const SparseMatrix& A, 
    const SparseMatrix& B, 
    SparseMatrix& C, 
    const std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    const I64 nestLevel);

} // namespace Ddpca

#endif // _SparseMatrix_hpp