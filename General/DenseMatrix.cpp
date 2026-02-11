#include "DenseMatrix.hpp"

#include <algorithm>
#include <cassert>
#include <ranges>
#include <numeric>

namespace Ddpca {

// Default constructor
DenseMatrix::DenseMatrix() noexcept : rows(0), cols(0) {
}

// Constructor - specify matrix dimensions
DenseMatrix::DenseMatrix(I64 numRows, I64 numCols)
    : rows(numRows),
      cols(numCols),
      data(numRows * numCols, 0.0) {
    assert(numRows > 0 && numCols > 0 && 
        "Matrix dimensions must be positive");
}

// Constructor - use existing data
DenseMatrix::DenseMatrix(I64 numRows, I64 numCols, std::vector<Real> inputData)
    : rows(numRows),
      cols(numCols),
      data(std::move(inputData)) {
    assert(numRows > 0 && numCols > 0 && 
        "Matrix dimensions must be positive");
    assert((data.size() == static_cast<size_t>(numRows * numCols)) && 
        "Input data size does not match matrix dimensions (should be rows×cols)");
}

// Move constructor
DenseMatrix::DenseMatrix(DenseMatrix&& other) noexcept
    : rows(other.rows),
      cols(other.cols),
      data(std::move(other.data)) {
    other.rows = 0;
    other.cols = 0;
}

// Move assignment operator
DenseMatrix& DenseMatrix::operator=(DenseMatrix&& other) noexcept {
    if (this != &other) {
        rows = other.rows;
        cols = other.cols;
        data = std::move(other.data);
        other.rows = 0;
        other.cols = 0;
    }
    return *this;
}

// Resize matrix dimensions
void DenseMatrix::Resize(I64 newRows, I64 newCols) {
    assert(newRows > 0 && newCols > 0 && 
        "Matrix dimensions must be positive");
    rows = newRows;
    cols = newCols;
    data.resize(newRows * newCols, 0.0);
}

// Create identity matrix (static method) - only for square matrices
DenseMatrix DenseMatrix::CreateIdentity(I64 dimension) {
    assert(dimension > 0 && 
        "Matrix dimension must be positive");
    DenseMatrix matrix(dimension, dimension);
    // Position of diagonal elements in row-major order: [i * dimension + i]
    for (I64 i = 0; i < dimension; ++ i) {
        matrix.data[i * dimension + i] = 1.0;
    }
    return matrix;
}

void DenseMatrix::Output(std::ostream& os) const {
    for (I64 i = 0; i < rows; ++ i) {
        for (I64 j = 0; j < cols; ++ j) {
            os << data[j * rows + i] << "\t";
        }
        os << "\n";
    }
}

Real DenseMatrix::DerterminantAfterLU(const I64* ipiv) const {
    assert(rows == cols && 
        "Matrix must be square to compute determinant");
    I64 sign = 0;
    // 统计行交换次数
    for (I64 i = 0; i < rows; ++ i) {
        if (ipiv[i] != i + 1) ++ sign; // ipiv[i] 是 1-based 索引
    }
    // 计算U的对角线乘积
    Real det = 1.0;
    for (I64 i = 0; i < rows; ++ i) {
        det *= data[i * cols + i];  // 列优先存储的U对角线元素
    }
    // 乘以置换矩阵的符号
    return (sign % 2 == 0) ? det : -det;
}

void GEMM(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix){
    // 确保矩阵维度匹配
    assert(inputMatrix_0.cols == inputMatrix_1.rows 
        && "Matrix dimensions don't match for multiplication");
    assert(outputMatrix.rows == inputMatrix_0.rows 
        && outputMatrix.cols == inputMatrix_1.cols 
        && "Output matrix dimensions incorrect");

    outputMatrix.Fill(0.0);
    I64 outputRows = inputMatrix_0.rows;
    I64 outputCols = inputMatrix_1.cols;
    I64 wCols = inputMatrix_1.rows;
    I64 base_0, base_1, base_2, tc, tw;
    Real weight;
    for(tc = 0; tc < outputCols; ++ tc){
        base_0 = tc * outputRows;
        base_2 = tc * wCols;
        for(tw = 0; tw < wCols; ++ tw){
            weight = inputMatrix_1.data[base_2 + tw];
            if(weight == 0.0) continue;
            base_1 = tw * outputRows;
            std::transform(std::execution::unseq,
                inputMatrix_0.data.begin() + base_1, inputMatrix_0.data.begin() + base_1 + outputRows,
                outputMatrix.data.begin() + base_0,
                outputMatrix.data.begin() + base_0,
                [weight](Real matVal, Real& val){ return val + weight * matVal; });
        }
    }
}

void GEPEMM(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix){
    // 确保矩阵维度匹配
    assert(inputMatrix_0.cols == inputMatrix_1.rows 
        && "Matrix dimensions don't match for multiplication");
    assert(outputMatrix.rows == inputMatrix_0.rows 
        && outputMatrix.cols == inputMatrix_1.cols 
        && "Output matrix dimensions incorrect");

    I64 outputRows = inputMatrix_0.rows;
    I64 outputCols = inputMatrix_1.cols;
    I64 wCols = inputMatrix_1.rows;
    I64 base_0, base_1, base_2, tc, tw;
    Real weight;
    for(tc = 0; tc < outputCols; ++ tc){
        base_0 = tc * outputRows;
        base_2 = tc * wCols;
        for(tw = 0; tw < wCols; ++ tw){
            weight = inputMatrix_1.data[base_2 + tw];
            if(weight == 0.0) continue;
            base_1 = tw * outputRows;
            std::transform(std::execution::unseq,
                inputMatrix_0.data.begin() + base_1, inputMatrix_0.data.begin() + base_1 + outputRows,
                outputMatrix.data.begin() + base_0,
                outputMatrix.data.begin() + base_0,
                [weight](Real matVal, Real& val){ return val + weight * matVal; });
        }
    }
}

void GEPEMM(
    Real alpha, 
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix){
    // 确保矩阵维度匹配
    assert(inputMatrix_0.cols == inputMatrix_1.rows 
        && "Matrix dimensions don't match for multiplication");
    assert(outputMatrix.rows == inputMatrix_0.rows 
        && outputMatrix.cols == inputMatrix_1.cols 
        && "Output matrix dimensions incorrect");

    I64 outputRows = inputMatrix_0.rows;
    I64 outputCols = inputMatrix_1.cols;
    I64 wCols = inputMatrix_1.rows;
    I64 base_0, base_1, base_2, tc, tw;
    Real weight;
    for(tc = 0; tc < outputCols; ++ tc){
        base_0 = tc * outputRows;
        base_2 = tc * wCols;
        for(tw = 0; tw < wCols; ++ tw){
            weight = alpha * inputMatrix_1.data[base_2 + tw];
            if(weight == 0.0) continue;
            base_1 = tw * outputRows;
            std::transform(std::execution::unseq,
                inputMatrix_0.data.begin() + base_1, inputMatrix_0.data.begin() + base_1 + outputRows,
                outputMatrix.data.begin() + base_0,
                outputMatrix.data.begin() + base_0,
                [weight](Real matVal, Real& val){ return val + weight * matVal; });
        }
    }
}

void GEMTM(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix){
    // 确保矩阵维度匹配
    assert(inputMatrix_0.rows == inputMatrix_1.rows 
        && "Matrix dimensions don't match for multiplication");
    assert(outputMatrix.rows == inputMatrix_0.cols 
        && outputMatrix.cols == inputMatrix_1.cols 
        && "Output matrix dimensions incorrect");
    
    I64 outputRows = inputMatrix_0.cols;
    I64 outputCols = inputMatrix_1.cols;
    I64 wCols = inputMatrix_0.rows;
    I64 base_0, base_1, base_2, tc, tr;
    for(tc = 0; tc < outputCols; ++ tc){
        base_1 = tc * outputRows;
        base_2 = tc * wCols;
        for(tr = 0; tr < outputRows; ++ tr){
            base_0 = tr * wCols;
            outputMatrix.data[base_1 + tr] = std::transform_reduce(std::execution::unseq,
                inputMatrix_0.data.data() + base_0, inputMatrix_0.data.data() + base_0 + wCols,
                inputMatrix_1.data.data() + base_2,
                0.0, // 初始累加值
                std::plus<>(), // 累加操作
                std::multiplies<>());
        }
    }
}

//outputMatrix = inputMatrix_0 * inputMatrix_1^T
void GEMMT(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix){
    assert(inputMatrix_0.cols == inputMatrix_1.cols 
        && "Matrix dimensions don't match for multiplication");
    assert(outputMatrix.rows == inputMatrix_0.rows 
        && outputMatrix.cols == inputMatrix_1.rows 
        && "Output matrix dimensions incorrect");

    outputMatrix.Fill(0.0);
    I64 outputRows = inputMatrix_0.rows;
    I64 outputCols = inputMatrix_1.rows;
    I64 wCols = inputMatrix_1.cols;
    I64 base_w, base_0, tc, tw;
    Real weight;
    for(tw = 0; tw < wCols; ++ tw){
        base_w = tw * outputCols;
        base_0 = tw * outputRows;
        for(tc = 0; tc < outputCols; ++ tc){
            weight = inputMatrix_1.data[base_w + tc];
            if(weight == 0.0) continue;
            std::transform(std::execution::unseq,
                inputMatrix_0.data.begin() + base_0, inputMatrix_0.data.begin() + base_0 + outputRows,
                outputMatrix.data.begin() + base_0,
                outputMatrix.data.begin() + base_0,
                [weight](Real matVal, Real& val){ return val += weight * matVal; });
        }
    }
}

//outputMatrix = inputMatrix_0^T * inputMatrix_1^T
void GEMTMT(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix){
    // 确保矩阵维度匹配
    assert(inputMatrix_0.rows == inputMatrix_1.cols 
        && "Matrix dimensions don't match for multiplication");
    assert(outputMatrix.rows == inputMatrix_0.cols 
        && outputMatrix.cols == inputMatrix_1.rows 
        && "Output matrix dimensions incorrect");
    // outputMatrix = (inputMatrix_1 * inputMatrix_0)^T
    // not true outputRows, just duplicate from GEMM
    
    outputMatrix.Fill(0.0);
    I64 outputRows = inputMatrix_1.rows;
    I64 outputCols = inputMatrix_0.cols;
    I64 wCols = inputMatrix_0.rows;
    I64 base_0, base_1, tc, tw;
    Real weight;
    for(tc = 0; tc < outputCols; ++ tc){
        base_0 = tc * wCols;
        for(tw = 0; tw < wCols; ++ tw){
            weight = inputMatrix_0.data[base_0 + tw];
            if(weight == 0.0) continue;
            base_1 = tw * outputRows;
            // Create a view of outputMatrix's j-th column
            auto col_output = std::views::iota(0, outputRows) 
				| std::views::transform([&](int tr) -> Real& { 
                    return outputMatrix.data[tr * outputMatrix.rows + tc]; });
            std::transform(std::execution::unseq,
                inputMatrix_1.data.begin() + base_1, inputMatrix_1.data.begin() + base_1 + outputRows,
                col_output.begin(),
                col_output.begin() ,
                [weight](Real matVal, Real& val){ return val += weight * matVal; });
        }
    }
}

void SCAL(Real alpha, DenseMatrix& inputMatrix){
    std::transform(std::execution::unseq,
        inputMatrix.data.begin(), inputMatrix.data.end(),
        inputMatrix.data.begin(),
        [alpha](Real val){ return val * alpha; });
}

} // namespace Ddpca