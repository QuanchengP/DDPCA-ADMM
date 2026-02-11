#ifndef _IntegralPoint_hpp
#define _IntegralPoint_hpp

#include "../General/General.hpp"

namespace Ddpca{

//Quadrature over master side, one quadrature point = one IntegralPoint.
//Project slave face onto master face (curvilinear plane), a polygon is obtained,
//divide the polygon into sub triangles by its centroid,
//2*2 Guass quadrature over one sub triangle,
//one sub triangle = four quadrature points = four IntegralPoint
struct IntegralPoint{
	std::array<std::array<I64, 4>, 2> node; //four nodes: one element face
	std::array<std::array<Real, 4>, 2> shapeFunction; //four shape/basis functions
	std::array<Real, 4> dualBasis; //four dual basis functions on non-mortar side
	std::array<std::array<Real, 3>, 2> contactPoint; // respectively on master side and slave side
	//basis vector, 0 - normal of master side, 1/2 - tangential
	std::array<std::array<Real, 3>, 3> basisVector;
	Real initialNormalGap;//= normVect * (contactPoint 1 - contactPoint 0)
	Real quadratureWeight;//quadrature weight
};

} // namespace Ddpca

#endif // _IntegralPoint_hpp