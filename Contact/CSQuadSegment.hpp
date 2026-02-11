#ifndef _CSQuadSegment_hpp
#define _CSQuadSegment_hpp

#include "IntegralPoint.hpp"

namespace Ddpca{

// local contact search between two quadrilateral segments
class CSQuadSegment{

public:

    std::array<std::array<Real, 3>, 4> masterCorners;
    std::array<std::array<Real, 3>, 4> slaveCorners;

    std::vector<IntegralPoint> integralPoints;

    //Projecting 3d faces onto 2d plane: Cartesian xyz TO curvilinear xi,eta.
    //Two faces in 2d are intersecting, the minimum value of the intersecting area.
    //If the area of a triangle is less than miniArea, the triangle is degenerated to a line-segment.
	//Why an absolute value: xi,eta but not xyz
    static constexpr Real minArea = 1.0E-10;
	//If |xi_0 - xi_1| <= minXietaDiff, then it is viewed that xi_0 == xi_1
	static constexpr Real minXietaDiff = 1.0E-6;

public:

	//subroutine of ProjectS2M
	void ProjectS2MSub(
		const std::array<Real, 3>& slavePoint, 
		std::array<std::array<Real, 3>, 2> &PrmaPxie, 
		std::array<Real, 2> &masterXiet, 
		Real &ngap);

	//project slave point slavPoin onto master face mastCorn
	void ProjectS2M(
		const std::array<Real, 3>& slavePoint, 
		std::array<std::array<Real, 3>, 2> &PrmaPxie, 
		std::array<Real, 2> &masterXiet, 
		Real &ngap);

	//project master face point mastXiet/mastPoin onto slave face slavCorn
	void ProjectMB2SSub(
		const std::array<Real, 2>& masterXiet, 
		const std::array<Real, 3>& masterPoint, 
		std::array<std::array<Real, 3>, 2> &PrmaPxie, 
		std::array<Real, 2> &slaveXiet, 
		Real &ngap);

	//project master face point mastXiet onto slave face slavCorn: maslPoin[0] to maslPoin[1]
	void ProjectMB2S(
		const std::array<Real, 2>& masterXiet, 
		std::array<std::array<Real, 4>, 2>& maslShape, 
		std::array<std::array<Real, 3>, 2>& maslPoint, 
		std::array<std::array<Real, 3>, 3>& basisVector, 
		Real &weightFactor);

	//area of 2d triangle
	Real TriangleArea2d(std::array<Real, 2> tempPoint_0, 
		std::array<Real, 2> tempPoint_1, std::array<Real, 2> tempPoint_2);

	//Guass quadrature over a triangle (tempXieta_0~2)
	void TriangleQuadrature(std::array<Real, 2> tempXiet_0, 
		std::array<Real, 2> tempXiet_1, std::array<Real, 2> tempXiet_2, 
		std::vector<std::array<Real, 2>> &listXiet, std::vector<double> &listWeight);

	//sort std::array<Real, 2> by the tempIndex-th component
	void SortBy2d(std::array<Real, 2> &tempPoint_0, 
        std::array<Real, 2> &tempPoint_1, I64 tempIndex);

	//whether line tempPoint_0~1 intersect with line tempPoint_2~3
	bool IsCross2d(
        const std::array<Real, 2>& tempPoint_0, const std::array<Real, 2>& tempPoint_1, 
		const std::array<Real, 2>& tempPoint_2, const std::array<Real, 2>& tempPoint_3);

	//the intersection of two lines
	void LineIntersection2d(std::array<Real, 2> tempPoint_0, std::array<Real, 2> tempPoint_1, 
		std::array<Real, 2> tempPoint_2, std::array<Real, 2> tempPoint_3, 
		std::vector<std::array<Real, 2>> &resultIntersections);

	//whether tempPoin in 2d face tempCorn
	bool InQuadrilateral(const std::array<Real, 2>& tempPoint, 
        const std::array<std::array<Real, 2>, 4>& tempCorners);

	//subroutine of Search
	void SegmentIntersect(
		std::vector<std::array<Real, 2>> &listXiet, std::vector<double> &listWeight);

	//intersection of segment tempMast and segment tempSlav
	void Search(std::array<I64, 4> tempMasterSegment, std::array<I64, 4> tempSlaveSegment);

}; // class CSQuadSegment

} // namespace Ddpca

#endif