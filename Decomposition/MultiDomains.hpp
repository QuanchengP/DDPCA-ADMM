#ifndef _MultiDomains_hpp
#define _MultiDomains_hpp

#include "MPLatin.hpp"

namespace Ddpca {

class MultiDomains{

public:

    Real NormalPenaltyCoef = 25.0;
    Real TangentialPenaltyCoef = 25.0;

    // Must be manually given: mesh, boundary
    // No need to manually given: geomMult, stiffness, conjGrad
    std::vector<SingleDomain> domains;
    // Must be manually given: all member variables of ContactInterface
    std::vector<ContactInterface> interfaces;
    // Must be manually given: mpLatin.realDomaLeve
    MPLatin mpLatin;

    std::vector<AlignedVectorRx> resuDisp;
    std::vector<std::array<AlignedVectorRx, 2>> inteAuxi;
    std::vector<std::array<AlignedVectorRx, 2>> inteMult;

public:

    Real CalculateCharacteristicLength();
    //the default calculation method for penalty parameter
    void CalculatePenaltyParameter();

    //domain - domain
    //does not influence the non-zero pattern of subdomain stiffness
    void SubDomainPenalty();

    std::vector<std::array<std::map<I64, I64>, 2>> nodeId2ContactId;
    void SubNi2ci();

    // domain - interface
    // !!!!!new version: displacement dof are not rotated!!!!!
    // domain force = domaInte * interface pressure
    std::vector<std::array<SparseMatrix, 2>> domaInte;
    std::vector<std::array<SparseMatrix, 2>> domaPenaInte;
    std::vector<std::array<SparseMatrix, 2>> domaPenaInte_T;
    void SubDomaInte();

    //interface - interface
    std::vector<std::array<SparseMatrix, 2>> inteInte;
    std::vector<std::array<SparseMatrix, 2>> intePenaInte;
    std::vector<std::array<DSCholesky, 2>> inteInteDs;
    std::vector<std::array<DSCholesky, 2>> intePenaInteDs;
    void SubInterfaceMass();

    // multiplier: integral point - interface
    std::vector<std::array<SparseMatrix, 2>> inpoInteLagr;
    void SubInpoInteLagr();

    // displacement: integral point - domain
    std::vector<std::array<SparseMatrix, 2>> inpoPenaDomaDisp;
    // integral point stress - integral point gap
    std::vector<AlignedVectorRx> inpoNgapStre; // always unchanged
    void SubInpoPenaDomaDisp();

    // interface - integral point
    std::vector<std::array<SparseMatrix, 2>> inteInpo;
    void SubInteInpo();

    void Establish();

    void OutpAuxiMult(
        const std::vector<AlignedVectorRx>& inpoGamm, 
        const std::string& directoryPath);
    void ADMM(const std::string& directoryPath);

}; // class MultiDomains

} // namespace Ddpca

#endif // _MultiDomains_hpp