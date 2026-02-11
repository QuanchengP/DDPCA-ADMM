#ifndef _MPLatin_hpp
#define _MPLatin_hpp

#include "SingleDomain.hpp"
#include "../Contact/ContactInterface.hpp"

namespace Ddpca{

// macroscopic problem in LATIN
class MPLatin{

public:

    const std::vector<SingleDomain>* domains;
    const std::vector<ContactInterface>* interfaces;
    const std::vector<std::array<std::map<I64, I64>, 2>>* nodeId2ContactId;

    // The domain level used to establish global coarse problem
    // is different to the coarsest domain level.
    // The geomMult has level 0 ~ maxLevel+1, but the realLvel
    // can be any in 0 ~ maxLevel.
    std::vector<I64> realDomaLeve;

public:

    std::vector<AlignedVectorRx> inteForc;
    std::vector<AlignedVectorRx> coarDisp;
    AlignedVectorRx globalSolution;
    AlignedVectorRx globalForce;

    // Contact version of GeometricMultigrid.userSolver.
    // Only contact nodes of master side are needed.
    // Auxiliary variable dof is not rotated, so only userSolver and no xyzRtz.
    std::vector<SparseMatrix> contactUserSolver;
    std::vector<SparseMatrix> accuContProl;
    std::vector<SparseMatrix> accuContProlT;
    void SubTransfer();

    // Must have SingleDomain.conjGrad.pcGM.Establish() called first.
    // Need SingleDomain.geomMult.rom/SingleDomain.conjGrad.pcGM.realProlong.
    std::vector<SparseMatrix> accuDomaProl;
    std::vector<SparseMatrix> accuDomaProlT;
    void SubAccuDomaProl();

    // Must have SingleDomain.conjGrad.pcGM.Establish() called first.
    // Need SingleDomain.geomMult.romT/SingleDomain.conjGrad.pcGM.realProlongT.
    std::vector<std::array<SparseMatrix, 2>> domaAuxi;
    void SubDomaAuxi();
    std::vector<SparseMatrix> auxiAuxi;
    void SubAuxiAuxi();

    std::vector<I64> MAccu;
    SparseMatrix globalCouple;
    DSCholesky macroDs;
    void SubGlobalCouple();

    std::vector<std::array<SparseMatrix, 2>> erroMult;
    std::vector<std::array<SparseMatrix, 2>> erroPenaAuxi;
    void SubErrorMA();

    std::vector<std::array<SparseMatrix, 2>> erroPenaDisp;
    void SubErroPenaDisp();

    void Establish();

    void Apply(
        std::vector<AlignedVectorRx>& resuDisp, 
        const std::vector<std::array<AlignedVectorRx, 2>>& inteAuxi, 
        const std::vector<std::array<AlignedVectorRx, 2>>& inteMult);

}; // class MPLatin

} // namespace Ddpca

#endif // _MPLatin_hpp