#ifndef _PCGeometricMultigrid_hpp
#define _PCGeometricMultigrid_hpp

#include "RSYMGS.hpp"
#include "RChebyshev.hpp"
#include "RJacobian.hpp"
#include "RSSOR.hpp"
#include "DSCholesky.hpp"

namespace Ddpca{

static constexpr I64 RelaxationSSGS = 0; // serial symmetric Gauss-Seidel
static constexpr I64 RelaxationPSGS = 1; // parallel symmetric Gauss-Seidel
static constexpr I64 RelaxationJacobian = 2;
static constexpr I64 RelaxationChebyshev = 3;
static constexpr I64 RelaxationSSOR = 4;

//preconditioner: geometric multigrid
class PCGeometricMultigrid{
public:

    std::vector<SparseMatrix> realProlong;//0 ~ maxLevel-1
    std::vector<SparseMatrix> hierarchyStiffness;//0 ~ maxLevel

public:

    std::vector<SparseMatrix> realProlongT;//0 ~ maxLevel-1
    DSCholesky coarsestDs; // level 0 - direct solver

    I64 relaxation = RelaxationSSGS;

    std::vector<RSYMGS> gaussSeidel;
    std::vector<RJacobian> jacobian;
    std::vector<RChebyshev> chebRela;
    std::vector<RSSOR> ssor;

public:

    void Establish(
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

    void Apply(
        I64 tempLevel, const AlignedVectorRx& rhs, 
        AlignedVectorRx& solution, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

public:

    std::vector<AlignedVectorRx> y;
    std::vector<AlignedVectorRx> residualError;
    std::vector<AlignedVectorRx> coarseResiErro;
    std::vector<AlignedVectorRx> coarseSolution;
    
}; //class PCGeometricMultigrid

} //namespace Ddpca

#endif // _PCGeometricMultigrid_hpp