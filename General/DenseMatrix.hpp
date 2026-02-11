#ifndef _DenseMatrix_hpp
#define _DenseMatrix_hpp

#include "General.hpp"

#include <stdexcept>
#include <cassert>
#include <execution>

namespace Ddpca {

// DenseMatrix class
class DenseMatrix {
public:
    I64 rows;                   // Number of matrix rows
    I64 cols;                   // Number of matrix columns
    std::vector<Real> data;   // Data storage (must col-major order, total elements: rows×cols)

public:
    // Default constructor
    DenseMatrix() noexcept;

    // Constructor - specify matrix dimensions
    DenseMatrix(I64 numRows, I64 numCols);

    // Constructor - use existing data
    DenseMatrix(I64 numRows, I64 numCols, std::vector<Real> inputData);

    // Copy constructor
    DenseMatrix(const DenseMatrix& other) = default;
    
    // Move constructor
    DenseMatrix(DenseMatrix&& other) noexcept;
    
    // Copy assignment operator
    DenseMatrix& operator=(const DenseMatrix& other) = default;
    
    // Move assignment operator
    DenseMatrix& operator=(DenseMatrix&& other) noexcept;

public:

    // Element access (0-based index, col-major order)
    inline Real& At(I64 row, I64 col) {
        assert(row >= 0 && row < rows && col >= 0 && col < cols && 
            "Matrix indices out of range");
        return data[col * rows + row];
    }
    
    inline const Real& At(I64 row, I64 col) const {
        assert(row >= 0 && row < rows && col >= 0 && col < cols && 
            "Matrix indices out of range");
        return data[col * rows + row];
    };

    // Get total number of elements
    inline I64 GetTotalElements() const noexcept { return rows * cols; }

    // Overload operator() for row-column access
    inline Real& operator()(I64 row, I64 col) { return At(row, col); }
    inline const Real& operator()(I64 row, I64 col) const { return At(row, col); }

    // Resize matrix dimensions
    void Resize(I64 newRows, I64 newCols);

    // Fill matrix with value
    inline void Fill(Real value) {
        std::fill_n(std::execution::unseq, data.data(), data.size(), value);
    }

    template<typename... Args>
    inline void Fill(Args&&... args) {
        assert(sizeof...(args) == data.size() && "Argument count mismatch");
        I64 ti = 0;
        ((data[ti++] = static_cast<Real>(args)), ...);
    }
    
    // Create identity matrix (static method) - only for square matrices
    static DenseMatrix CreateIdentity(I64 dimension);

    void Output(std::ostream& os = std::cout) const;

    //can only be called after LU factorization
    Real DerterminantAfterLU(const I64* ipiv) const;
};

//outputMatrix = inputMatrix_0 * inputMatrix_1
void GEMM(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix);

//outputMatrix += inputMatrix_0 * inputMatrix_1
void GEPEMM(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix);

//outputMatrix += alpha * inputMatrix_0 * inputMatrix_1
void GEPEMM(
    Real alpha, 
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix);

//outputMatrix = inputMatrix_0^T * inputMatrix_1
void GEMTM(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix);

//outputMatrix = inputMatrix_0 * inputMatrix_1^T
void GEMMT(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix);

//outputMatrix = inputMatrix_0^T * inputMatrix_1^T
void GEMTMT(
    const DenseMatrix& inputMatrix_0, 
    const DenseMatrix& inputMatrix_1, 
    DenseMatrix& outputMatrix);

void SCAL(
    Real alpha, DenseMatrix& inputMatrix);

//outputVector = inputMatrix * inputVector
template <typename InputContainer, typename OutputContainer>
void GEMV(
    const DenseMatrix& inputMatrix, 
    const InputContainer& inputVector, 
    OutputContainer& outputVector){
    //
    outputVector.fill(0.0);
    I64 base_0;
    Real weight;
    for(I64 tc = 0; tc < inputMatrix.cols; ++ tc){
        base_0 = tc * inputMatrix.rows;
        weight = inputVector[tc];
        std::transform(std::execution::unseq,
            inputMatrix.data.begin() + base_0, inputMatrix.data.begin() + base_0 + inputMatrix.rows,
            outputVector.begin(),
            outputVector.begin(),
            [weight](Real matVal, Real& val){ return weight * matVal + val; }); // +=?
    }
}

//outputVector = inputMatrix^T * inputVector
template <typename InputContainer, typename OutputContainer>
void GEMTV(
    const DenseMatrix& inputMatrix, 
    const InputContainer& inputVector, 
    OutputContainer& outputVector){
    //
    I64 base_0;
    for(I64 tc = 0; tc < inputMatrix.cols; ++ tc){
        base_0 = tc * inputMatrix.rows;
        outputVector[tc] = std::transform_reduce(std::execution::unseq,
            inputMatrix.data.begin() + base_0, inputMatrix.data.begin() + base_0 + inputMatrix.rows,
            inputVector.begin(),
            0.0, 
            std::plus<>(), 
            [](Real matVal, Real vecVal){ return matVal * vecVal; }); // +=?
    }
}

//outputVector += inputMatrix * inputVector
template <typename InputContainer, typename OutputContainer>
void PEMV(
    const DenseMatrix& inputMatrix, 
    const InputContainer& inputVector, 
    OutputContainer& outputVector){
    //
    I64 base_0;
    Real weight;
    for(I64 tc = 0; tc < inputMatrix.cols; ++ tc){
        base_0 = tc * inputMatrix.rows;
        weight = inputVector[tc];
        std::transform(std::execution::unseq,
            inputMatrix.data.begin() + base_0, inputMatrix.data.begin() + base_0 + inputMatrix.rows,
            outputVector.begin(),
            outputVector.begin(),
            [weight](Real matVal, Real& val){ return val + weight * matVal; }); // +=?
    }
}

//outputVector = inputMatrix^-1 * inputVector
template <typename InputContainer, typename OutputContainer>
void Solve(
    const DenseMatrix& inputMatrix, 
    const InputContainer& inputVector, 
    OutputContainer& outputVector){
    //
    assert((inputMatrix.rows == 2 || inputMatrix.cols == 3) 
        && "Matrix dimensions are not applicable");
    switch(inputMatrix.rows){
        case 2:{
            Real D = inputMatrix(0, 0) * inputMatrix(1, 1) - inputMatrix(0, 1) * inputMatrix(1, 0);
            Real D_1 = inputVector[0] * inputMatrix(1, 1) - inputVector[1] * inputMatrix(0, 1);
            Real D_2 = inputMatrix(0, 0) * inputVector[1] - inputMatrix(1, 0) * inputVector[0];
            outputVector[0] = D_1 / D;
            outputVector[1] = D_2 / D;
            break;
        }
        default:
            Log("Matrix dimensions are not applicable");
    }
}

} // namespace Ddpca

#endif // _DenseMatrix_hpp