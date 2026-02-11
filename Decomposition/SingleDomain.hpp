#ifndef _SingleDomain_hpp
#define _SingleDomain_hpp

#include "../Mesh/GeometricMultigrid.hpp"
#include "../Solver/ConjugateGradient.hpp"

namespace Ddpca {

class SingleDomain{

public:

    Mesh mesh;
    BoundaryCondition boundary;
    GeometricMultigrid geomMult;
    SparseMatrix stiffness;
    constexpr static I64 dscoLimi = 50000;
    ConjugateGradient conjGrad;
    DSCholesky cholDs;

    Real elasticity = 210.0E9;
    Real poissonRatio = 0.3;
    
}; // class SingleDomain

} // namespace Ddpca

#endif // _SingleDomain_hpp