#ifndef _ContactInterface_hpp
#define _ContactInterface_hpp

#include "IntegralPoint.hpp"
#include "CurvedSurface.hpp"
#include "../Mesh/Mesh.hpp"
#include "../Mesh/BilinearQuadrature.hpp"

namespace Ddpca{

class ContactInterface{

public:

    //two-body contact
    std::array<I64, 2> domainIndex;
	//fricCoef > 0: frictional contact
	//fricCoef = 0: normal contact
	//fricCoef < 0: perfect interface
	double frictionCoefficient;
	double normPenaPara; // normal penalty parameter
	double tangPenaPara; // tangential penalty parameter

    //n*4, n segments/element faces, 4 nodes per segment
    std::vector<std::array<I64, 4>> masterSegments;
    std::vector<std::array<I64, 4>> slaveSegments;

    std::vector<IntegralPoint> integralPoints;

	//global contact search, bucket search
	//2d local coordinate, 0 - infimum, 1 - supremum, 2 - increment
    std::array<std::array<Real, 3>, 2> bucketLocal;
	//sort master segments into different buckets
    std::vector<std::vector<std::vector<I64>>> bucket;

public:

    void OutputSegments(const std::string& directoryPath, const I64& fileIden) const;
    
	void OutputIntegralPoints(const std::string& directoryPath, const I64& fileIden) const;

    void BucketSort(
        const std::array<std::vector<Real>, 2>& masterLocal, std::array<I64, 2> divisionNumber);

    // MUST: the slave element face is smaller than the master element face!!!
    // slaveLocal must be compatible with masterLocal.
	void LocalSearch(
		const Mesh& masterMesh, 
        const Mesh& slaveMesh, 
		const std::vector<std::array<Real, 8>>& slaveLocal, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
	    const I64 nestLevel, 
        double maxiDist = 1.0E12);

    void InterpolateBasisGap(const Mesh& masterMesh, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
	    const I64 nestLevel);
    
    void AdaptiveRefine(Mesh& tempMast, Mesh& tempSlav, bool &isRefi, 
        CurvedSurface& mastSurf, CurvedSurface& slavSurf, I64 tempLeve, 
        Real distCrit, std::array<I64, 2> divisionNumber, 
        std::function<void(Coordinate, Real &, Real &)> CART_CURV);

}; // class ContactInterface

} // namespace Ddpca

#endif // _ContactInterface_hpp