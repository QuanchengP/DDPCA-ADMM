#ifndef _BoundaryCondition_hpp
#define _BoundaryCondition_hpp

#include "../General/General.hpp"
#include "../General/DenseMatrix.hpp"

#include <map>

namespace Ddpca {

/****************************************************************************************************/
class BoundaryCondition {
public:
	//only zero displacement constraint is considered
    std::map<I64, Real> constrainedDof; //DOF's displacement is given in rtz > xyz
	//already in urz direction aligned with constrainedDof
    std::map<I64, Real> externalForce; //DOF's external force is given in rtz > xyz
	//(displacement in xyz) = nodeRota * (displacement in rtz)
	//must avoid duplication of node ID
    std::map<I64, DenseMatrix> nodeRotation;//node rotation matrix

    //solver order
    AlignedVectorRx loadVect;

    //0: the i-th DOF (total 3*nodeNumb) is constrained
    std::vector<I64> consFlag;//1 - not constrained, 0 - constrained
    std::vector<I64> consPrefSum;//constraint prefix summation

public:

    //accumulate loadValue to dof
    void LoadAccumulate(I64 dof, Real loadValue){
        //
        assert(constrainedDof.find(dof) == constrainedDof.end() 
            && "dof is already constrained!");
        auto iteratorEf = externalForce.find(dof);
        if(iteratorEf == externalForce.end()){
            externalForce.emplace(dof, loadValue);
        }
        else{
            iteratorEf->second += loadValue;
        }
    }
}; // class BoundaryCondition

} // namespace Ddpca

#endif // _BoundaryCondition_hpp