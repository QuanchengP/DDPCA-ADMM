#ifndef _Mesh_hpp
#define _Mesh_hpp

#include "../General/General.hpp"
#include "../General/AlignedVector.hpp"
#include "../General/DenseMatrix.hpp"
#include "Coordinate.hpp"
#include "OctreeElement.hpp"

#include <map>
#include <vector>
#include <set>

namespace Ddpca {

/****************************************************************************************************/
class Mesh{
public:
	//from node numbering to coordinate
	std::map<I64, Coordinate> node2Coordinate;//every node must be valid

	//from coordinate to node numbering
	std::map<Coordinate, I64> coordinate2Node;

	//line has been used by which elements
	std::map<std::array<I64, 2>, std::set<I64>> lineUsedByElement;

	//face has been used by which elements
	std::map<std::array<I64, 4>, std::set<I64>> faceUsedByElement;

	//each element is a TREE_ELEM
	std::vector<OctreeElement> elements;

	I64 maxLevel = 0;//the maximum refinement level
	std::vector<I64> effePrefSum;
    
public:
	
	//try to add one node, the added node may already exist
	//WARNNING: it is designed to reuse already existing node,
	//          every added node will be used by final mesh
    I64 TryAddNode(const Coordinate& tempCoordinate);

	//add one element
	I64 AddElement(const OctreeElement& tempElement);

	//refine the elements in spliElem
	//subElements: contains specified sub elements after refinement
	//curvilinearInterpolation: Cartesian - curvilinear - Cartesian
	void Refine(std::set<I64> &elementsToSplit, 
	    const std::map<I64, std::set<I64>> &subElements, 
	    const std::map<std::vector<I64>, Coordinate> &curvilinearInterpolation
	);
	//gradual level check
	void AdjacentLevelCheck(std::set<I64> &elementsToSplit);

	//output element information
	void OutputMesh(const std::string& directoryPath, const I64& fileIden) const;

	//rotation and translation of rigid body
	void RigidRotationTranslation(const DenseMatrix& rotationMatrix, const Coordinate& translationVector);

	//output displacement
	void OutputDisplacement(
	    const std::string& directoryPath, const I64& fileIden, 
    	const AlignedVectorRx& resultDisplacement
	) const;

	I64 EffectiveElements(
		const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
		const I64 nestLevel);

	Real Volume(
		const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
		const I64 nestLevel);

}; // class Mesh

} // namespace Ddpca

#endif // _Mesh_hpp