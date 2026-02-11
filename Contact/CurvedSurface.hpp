#ifndef _CurvedSurface_hpp
#define _CurvedSurface_hpp

#include "../Mesh/Mesh.hpp"

namespace Ddpca{

class CurvedSurface{

public:

    std::vector<std::vector<Coordinate>> indexPoint;
    std::map<Coordinate, std::array<I64, 2>> pointIndex;

public:

    void Resize(I64 rows, I64 cols);

	void Insert(I64 ti, I64 tj, Coordinate inputCoordinate);
    
	//output: true - on the curved surface, outpCoor is valid
	//otherwise - not on the curved surface, outpCoor is not valid
	bool QuadrilateralRefine(
        const std::vector<Coordinate>& inputCoordinate, Coordinate &outputCoordinate);

	//if one face of an element is on CurvedSurface, this element is gonna be refined
	//(put into planSurf)
	void Refine(const Mesh& tempMesh, 
        const std::set<I64>& elementsToSplit, 
		std::map<std::vector<I64>, Coordinate>& curvInte);

	//rotation and translation of rigid body
	void RigidRotationTranslation(
        const DenseMatrix& rotationMatrix, const Coordinate& translationVector);

public:

    I64 ti;
    I64 tj;
    std::array<I64, 4> currentFace;

    void Initialize();

    bool Increment(const Mesh& tempMesh);
}; // class CurvedSurface

} // namespace Ddpca

#endif // _CurvedSurface_hpp