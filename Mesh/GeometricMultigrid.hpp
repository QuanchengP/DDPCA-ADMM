#ifndef _GeometricMultigrid_hpp
#define _GeometricMultigrid_hpp

#include "Mesh.hpp"
#include "BoundaryCondition.hpp"
#include "FemStiffnessMatrix.hpp"
#include "../Solver/PCGeometricMultigrid.hpp"

namespace Ddpca {

class GeometricMultigrid{
public:
    Mesh* dataMesh;
    BoundaryCondition *dataBoundary;
    SparseMatrix* dataStiffness;
    PCGeometricMultigrid* dataPcg;

	//coarser node to hanging node, may be empty
	std::map<std::vector<I64>,I64> coarser2Hanging;
	//hanging node to coarser node, may be empty
	std::map<I64,std::vector<I64>> hanging2Coarser;

	//nodes of different level: level 0 is the coarsest
	std::vector<std::vector<I64>> levelNode;
	//node is locating at level/position
	std::vector<std::array<I64,2>> nodeLevelPosition;
	std::vector<I64> positionNode;//from position to node

    //xyz = xyzRtz * rtz
    SparseMatrix xyzRtz;
    //user DOF = userSolver * solver DOF
    SparseMatrix userSolver;
    //rotateOrder = xyzRtz * userSolver
    SparseMatrix rotateOrder;

	//prolongation operator: solver layer, but not user layer, subordinate to xyzRtz/userSolver
	std::vector<SparseMatrix> scalarProlong;
    //each level has to consider nodal DOF rotation: add nodal rotation for scalProl
    std::vector<SparseMatrix> rotateProlong;//0 ~ maxiLeve
    //each level has to consider displacement constraint: add constraint for rotaProl
    //it is just a operator in the same level, not a transfer operator between different levels
    std::vector<SparseMatrix> constraintProlong;//0 ~ maxiLeve
    SparseMatrix maxProlong;

    // rotateOrder * maxProlong
    SparseMatrix rom;
    // rom^T
    SparseMatrix romT;

public:
    //
    void Preprocess(
        std::vector<std::map<std::vector<I64>,I64>>& ininTran, 
        std::vector<std::set<I64>>& levelNode_s, 
        bool willMultigrid = true);
    //local refinement: handle the hanging node to be nested
    void NestedHangingNode(
        const std::vector<std::map<std::vector<I64>,I64>>& ininTran);

	//establish scalarProlong
	void Transfer(
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);
    //nodal DOF rotation on dataStiffness
    //user DOF ordering to solver DOF ordering
    void DofRotateOrder(
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);
    void RotateProlong(
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);
    void ConstraintProlong(
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);
    //real transfer operator: consider displacement constraint
    void Establish(
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

    void DispPost(
        const AlignedVectorRx& inpuDisp, 
        AlignedVectorRx& outpDisp, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel) const;

    void Output(
        const std::string& directoryPath, 
        const I64& fileIden, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel) const;
}; //class GeometricMultigrid

} // namespace Ddpca

#endif // _GeometricMultigrid_hpp