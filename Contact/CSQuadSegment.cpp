#include "CSQuadSegment.hpp"
#include "../Mesh/BilinearQuadrature.hpp"
#include "../Mesh/TriangularQuadrature.hpp"

namespace Ddpca{

void CSQuadSegment::ProjectS2MSub(
    const std::array<Real, 3>& slavePoint, 
    std::array<std::array<Real, 3>, 2> &PrmaPxie, 
    std::array<Real, 2> &masterXiet, 
    Real &ngap){
    //
    const BilinearQuadrature<2>& biliQuad = GetBilinearQuadrature();
    DenseMatrix factMatr(3, 4);
    for(I64 ti = 0; ti < 3; ++ ti){
        for(I64 tj = 0; tj < 4; ++ tj){
            factMatr(ti, 0) += masterCorners[tj][ti] / 4.0;
            factMatr(ti, 1) += masterCorners[tj][ti] * biliQuad.cornerNodes[tj][0] / 4.0;
            factMatr(ti, 2) += masterCorners[tj][ti] * biliQuad.cornerNodes[tj][1] / 4.0;
            factMatr(ti, 3) += masterCorners[tj][ti] 
                * biliQuad.cornerNodes[tj][0] * biliQuad.cornerNodes[tj][1] / 4.0;
        }
    }
    std::array<std::array<Real, 6>, 2> equaFact;
    equaFact[0][0] = factMatr(0, 0) * factMatr(0, 1) - factMatr(0, 1) * slavePoint[0] 
		+ factMatr(1, 0) * factMatr(1, 1) - factMatr(1, 1) * slavePoint[1] 
		+ factMatr(2, 0) * factMatr(2, 1) - factMatr(2, 1) * slavePoint[2];
    equaFact[0][1] = factMatr(0, 1) * factMatr(0, 1) 
		+ factMatr(1, 1) * factMatr(1, 1) + factMatr(2, 1) * factMatr(2, 1);
    equaFact[0][2] = factMatr(0, 1) * factMatr(0, 2) + factMatr(0, 0) * factMatr(0, 3) 
        - factMatr(0, 3) * slavePoint[0] 
		+ factMatr(1, 1) * factMatr(1, 2) +factMatr(1, 0) * factMatr(1, 3) 
		- factMatr(1, 3) * slavePoint[1] 
		+ factMatr(2, 1) * factMatr(2, 2) +factMatr(2, 0) * factMatr(2, 3) 
		- factMatr(2, 3) * slavePoint[2];
	equaFact[0][3] = 2.0 * factMatr(0, 1) * factMatr(0, 3) 
		+ 2.0 * factMatr(1, 1) * factMatr(1, 3) 
		+ 2.0 * factMatr(2, 1) * factMatr(2, 3);
	equaFact[0][4] = factMatr(0, 2) * factMatr(0, 3) 
		+ factMatr(1, 2) * factMatr(1, 3) + factMatr(2, 2) * factMatr(2, 3);
	equaFact[0][5] = factMatr(0, 3) * factMatr(0, 3) 
		+ factMatr(1, 3) * factMatr(1, 3) + factMatr(2, 3) * factMatr(2, 3);
	equaFact[1][0] = factMatr(0, 0) * factMatr(0, 2) - factMatr(0, 2) * slavePoint[0] 
		+ factMatr(1, 0) * factMatr(1, 2) - factMatr(1, 2) * slavePoint[1] 
		+ factMatr(2, 0) * factMatr(2, 2) - factMatr(2, 2) * slavePoint[2];
	equaFact[1][1] = factMatr(0, 2) * factMatr(0, 2) 
		+ factMatr(1, 2) * factMatr(1, 2) + factMatr(2, 2) * factMatr(2, 2);
	equaFact[1][2] = equaFact[0][2];
	equaFact[1][3] = 2.0 * factMatr(0, 2) * factMatr(0, 3) 
		+ 2.0 * factMatr(1, 2) * factMatr(1, 3) 
		+ 2.0 * factMatr(2, 2) * factMatr(2, 3);
	equaFact[1][4] = factMatr(0, 1) * factMatr(0, 3) 
		+ factMatr(1, 1) * factMatr(1, 3) + factMatr(2, 1) * factMatr(2, 3);
	equaFact[1][5] = equaFact[0][5];
    //
    I64 tj;
    const I64 maxIteration = 60;
    std::array<Real, 2> residual;
    std::array<Real, 2> deltaXiet;
    for(tj = 0; tj < maxIteration; ++ tj){
		residual[0] = equaFact[0][0] + equaFact[0][1] * masterXiet[0] + equaFact[0][2] * masterXiet[1] 
			+ equaFact[0][3] * masterXiet[0] * masterXiet[1] 
			+ equaFact[0][4] * masterXiet[1] * masterXiet[1] 
			+ equaFact[0][5] * masterXiet[0] * masterXiet[1] * masterXiet[1];
		residual[1] = equaFact[1][0] + equaFact[1][1] * masterXiet[1] + equaFact[1][2] * masterXiet[0] 
			+ equaFact[1][3] * masterXiet[0] * masterXiet[1] 
			+ equaFact[1][4] * masterXiet[0] * masterXiet[0] 
			+ equaFact[1][5] * masterXiet[0] * masterXiet[0] * masterXiet[1];
		std::array<std::array<Real, 2>, 2> hessMatr;
		hessMatr[0][0] = equaFact[0][1] + equaFact[0][3] * masterXiet[1] 
			+ equaFact[0][5] * masterXiet[1] * masterXiet[1];
		hessMatr[0][1] = equaFact[0][2] + equaFact[0][3] * masterXiet[0] 
			+ equaFact[0][4] * 2.0 * masterXiet[1] + equaFact[0][5] * masterXiet[0] * 2.0 * masterXiet[1];
		hessMatr[1][0] = equaFact[1][2] + equaFact[1][3] * masterXiet[1] 
			+ equaFact[1][4] * 2.0 * masterXiet[0] + equaFact[1][5] * 2.0 * masterXiet[0] * masterXiet[1];
		hessMatr[1][1] = equaFact[1][1] + equaFact[1][3] * masterXiet[0] 
			+ equaFact[1][5] * masterXiet[0] * masterXiet[0];
        Real detHessMatr = hessMatr[0][0] * hessMatr[1][1] - hessMatr[1][0] * hessMatr[0][1];
        deltaXiet[0] = - (hessMatr[1][1] * residual[0] - hessMatr[0][1] * residual[1]) / detHessMatr;
        deltaXiet[1] = - (hessMatr[0][0] * residual[1] - hessMatr[1][0] * residual[0]) / detHessMatr;
        if(NRM2(deltaXiet) / (NRM2(masterXiet) + 1.0E-16) < 1.0E-16 || NRM2(residual) < 1.0E-16){
            break;
        }
        XPEY(masterXiet, deltaXiet);
    }
    if(tj >= maxIteration){
        Log("            Warning 1 in ContactInterface::ProjectS2MSub! deltaXiet = " 
            + Double2String(deltaXiet[0]) + ", " + Double2String(deltaXiet[1])
            + ", residual = " + Double2String(NRM2(residual)));
    }
    for(I64 ti = 0; ti < 3; ++ ti){
        PrmaPxie[0][ti] = factMatr(ti, 1) + factMatr(ti, 3) * masterXiet[1];
        PrmaPxie[1][ti] = factMatr(ti, 2) + factMatr(ti, 3) * masterXiet[0];
    }
	std::array<Real, 3> tempNorm = Cross(PrmaPxie[0], PrmaPxie[1]);
    SCAL(1.0 / NRM2(tempNorm), tempNorm);
	std::array<Real, 4> tempVect = {1.0, masterXiet[0], masterXiet[1], masterXiet[0] * masterXiet[1]};
    std::array<Real, 3> masterPoint;
    GEMV(factMatr, tempVect, masterPoint);
    AXPY(-1.0, slavePoint, masterPoint);
	ngap = - DOT(tempNorm, masterPoint);
}

void CSQuadSegment::ProjectS2M(
    const std::array<Real, 3>& slavePoint, 
    std::array<std::array<Real, 3>, 2> &PrmaPxie, 
    std::array<Real, 2> &masterXiet, 
    Real &ngap){
    //
    Real elementDiameter = 1.0E15;
    std::array<Real, 3> masterCenter = {0.0, 0.0, 0.0};
    for(I64 ti = 0; ti < 4; ++ ti){
        std::array<Real, 3> tempVector_ti = masterCorners[ti];
        for(I64 tj = ti + 1; tj < 4; ++ tj){
            std::array<Real, 3> tempVector_tj = masterCorners[tj];
            AXPY(-1.0, tempVector_ti, tempVector_tj);
            elementDiameter = std::min(elementDiameter, NRM2(tempVector_tj));
        }
        AXPY(0.25, tempVector_ti, masterCenter);
    }
    //
    masterXiet.fill(0.0);
    std::array<Real, 3> tempDirection = masterCenter;
    AXPY(-1.0, slavePoint, tempDirection);
    Real normDirection = NRM2(tempDirection);
    I64 realTime = normDirection / elementDiameter;
    SCAL(- 1.0 / normDirection, tempDirection);
    std::array<Real, 3> slavePoint_1 = masterCenter;
    //the slavePoint maybe very far away from masterCorners
    for(I64 ti = 0; ti < realTime; ++ ti){
        AXPY(elementDiameter, tempDirection, slavePoint_1);
        ProjectS2MSub(slavePoint_1, PrmaPxie, masterXiet, ngap);
    }
    ProjectS2MSub(slavePoint, PrmaPxie, masterXiet, ngap);
    // if(masterXiet[0] < -1.0 - minXietaDiff || masterXiet[0] > 1.0 + minXietaDiff 
	// 	|| masterXiet[1] < -1.0 - minXietaDiff || masterXiet[1] > 1.0 + minXietaDiff){
    //     Log("            Warning 1 in ContactInterface::ProjectS2M! masterXiet = " 
    //         + Double2String(masterXiet[0]) + ", " + Double2String(masterXiet[1]));
    // }
}

void CSQuadSegment::ProjectMB2SSub(
    const std::array<Real, 2>& masterXiet, 
    const std::array<Real, 3>& masterPoint, 
    std::array<std::array<Real, 3>, 2> &PrmaPxie, 
    std::array<Real, 2> &slaveXiet, 
    Real &ngap){
    //
    const BilinearQuadrature<2>& biliQuad = GetBilinearQuadrature();
    PrmaPxie[0].fill(0.0);
    PrmaPxie[1].fill(0.0);
    for(I64 ti = 0; ti < 3; ++ ti){
        for(I64 tj = 0; tj < 4; ++ tj){
            PrmaPxie[0][ti] += masterCorners[tj][ti] * (biliQuad.cornerNodes[tj][0] / 4.0 
                + biliQuad.cornerNodes[tj][0] * biliQuad.cornerNodes[tj][1] * masterXiet[1] / 4.0);
            PrmaPxie[1][ti] += masterCorners[tj][ti] * (biliQuad.cornerNodes[tj][1] / 4.0 
                + biliQuad.cornerNodes[tj][0] * biliQuad.cornerNodes[tj][1] * masterXiet[0] / 4.0);
        }
    }
    DenseMatrix factMatr(3, 4);
    for(I64 ti = 0; ti < 3; ++ ti){
        for(I64 tj = 0; tj < 4; ++ tj){
            factMatr(ti, 0) += slaveCorners[tj][ti] / 4.0;
            factMatr(ti, 1) += slaveCorners[tj][ti] * biliQuad.cornerNodes[tj][0] / 4.0;
            factMatr(ti, 2) += slaveCorners[tj][ti] * biliQuad.cornerNodes[tj][1] / 4.0;
            factMatr(ti, 3) += slaveCorners[tj][ti] 
                * biliQuad.cornerNodes[tj][0] * biliQuad.cornerNodes[tj][1] / 4.0;
        }
    }
    std::array<std::array<Real, 4>, 2> equaFact;
	equaFact[0][0] = factMatr(0, 0) * PrmaPxie[0][0] + factMatr(1, 0) * PrmaPxie[0][1] 
		+ factMatr(2, 0) * PrmaPxie[0][2] - masterPoint[0] * PrmaPxie[0][0] 
		- masterPoint[1] * PrmaPxie[0][1] - masterPoint[2] * PrmaPxie[0][2];
	equaFact[0][1] = factMatr(0, 1) * PrmaPxie[0][0] + factMatr(1, 1) * PrmaPxie[0][1] 
		+ factMatr(2, 1) * PrmaPxie[0][2];
	equaFact[0][2] = factMatr(0, 2) * PrmaPxie[0][0] + factMatr(1, 2) * PrmaPxie[0][1] 
		+ factMatr(2, 2) * PrmaPxie[0][2];
	equaFact[0][3] = factMatr(0, 3) * PrmaPxie[0][0] + factMatr(1, 3) * PrmaPxie[0][1] 
		+ factMatr(2, 3) * PrmaPxie[0][2];
	equaFact[1][0] = factMatr(0, 0) * PrmaPxie[1][0] + factMatr(1, 0) * PrmaPxie[1][1] 
		+ factMatr(2, 0) * PrmaPxie[1][2] - masterPoint[0] * PrmaPxie[1][0] 
		- masterPoint[1] * PrmaPxie[1][1] - masterPoint[2] * PrmaPxie[1][2];
	equaFact[1][1] = factMatr(0, 1) * PrmaPxie[1][0] + factMatr(1, 1) * PrmaPxie[1][1] 
		+ factMatr(2, 1) * PrmaPxie[1][2];
	equaFact[1][2] = factMatr(0, 2) * PrmaPxie[1][0] + factMatr(1, 2) * PrmaPxie[1][1] 
		+ factMatr(2, 2) * PrmaPxie[1][2];
	equaFact[1][3] = factMatr(0, 3) * PrmaPxie[1][0] + factMatr(1, 3) * PrmaPxie[1][1] 
		+ factMatr(2, 3) * PrmaPxie[1][2];
    slaveXiet.fill(0.0);
    I64 tj;
    const I64 maxIteration = 60;
    std::array<Real, 2> residual;
    std::array<Real, 2> deltaXiet;
    for(tj = 0; tj < maxIteration; ++ tj){
        residual[0] = equaFact[0][0] + equaFact[0][1] * slaveXiet[0] + equaFact[0][2] * slaveXiet[1] 
			+ equaFact[0][3] * slaveXiet[0] * slaveXiet[1];
        residual[1] = equaFact[1][0] + equaFact[1][1] * slaveXiet[0] + equaFact[1][2] * slaveXiet[1] 
			+ equaFact[1][3] * slaveXiet[0] * slaveXiet[1];
        std::array<std::array<Real, 2>, 2> hessMatr;
        hessMatr[0][0] = equaFact[0][1] + equaFact[0][3] * slaveXiet[1];
        hessMatr[0][1] = equaFact[0][2] + equaFact[0][3] * slaveXiet[0];
        hessMatr[1][0] = equaFact[1][1] + equaFact[1][3] * slaveXiet[1];
        hessMatr[1][1] = equaFact[1][2] + equaFact[1][3] * slaveXiet[0];
        Real detHessMatr = hessMatr[0][0] * hessMatr[1][1] - hessMatr[1][0] * hessMatr[0][1];
        deltaXiet[0] = - (hessMatr[1][1] * residual[0] - hessMatr[0][1] * residual[1]) / detHessMatr;
        deltaXiet[1] = - (hessMatr[0][0] * residual[1] - hessMatr[1][0] * residual[0]) / detHessMatr;
        if(NRM2(deltaXiet) / (NRM2(slaveXiet) + 1.0E-16) < 1.0E-16 || NRM2(residual) < 1.0E-16){
            break;
        }
        XPEY(slaveXiet, deltaXiet);
    }
    if(tj >= maxIteration){
        Log("            Warning 1 in ContactInterface::ProjectMB2SSub! deltaXiet = " 
            + Double2String(deltaXiet[0]) + ", " + Double2String(deltaXiet[1])
            + ", residual = " + Double2String(NRM2(residual)));
    }
    //
    std::array<Real, 3> tempNorm = Cross(PrmaPxie[0], PrmaPxie[1]);
    SCAL(1.0 / NRM2(tempNorm), tempNorm);
    std::array<Real, 4> tempVect = {1.0, slaveXiet[0], slaveXiet[1], slaveXiet[0] * slaveXiet[1]};
    std::array<Real, 3> slavePoint;
    GEMV(factMatr, tempVect, slavePoint);
    AXPY(-1.0, masterPoint, slavePoint);
    ngap = DOT(tempNorm, slavePoint);
}

void CSQuadSegment::ProjectMB2S(
    const std::array<Real, 2>& masterXiet, 
    std::array<std::array<Real, 4>, 2>& maslShape, 
    std::array<std::array<Real, 3>, 2>& maslPoint, 
    std::array<std::array<Real, 3>, 3>& basisVector, 
    Real &weightFactor){
    //
    std::array<std::array<Real, 3>, 2> PrmaPxie;
    std::array<Real, 2> slaveXiet;
    Real ngap;
    ProjectMB2SSub(masterXiet, maslPoint[0], PrmaPxie, slaveXiet, ngap);
    //
    basisVector[1] = PrmaPxie[0];
    basisVector[2] = PrmaPxie[1];
    SCAL(1.0 / NRM2(basisVector[1]), basisVector[1]);
    SCAL(1.0 / NRM2(basisVector[2]), basisVector[2]);
    basisVector[0] = Cross(basisVector[1], basisVector[2]);
    SCAL(1.0 / NRM2(basisVector[0]), basisVector[0]);
    //
    weightFactor = std::sqrt(
        std::pow(PrmaPxie[0][1] * PrmaPxie[1][2] - PrmaPxie[1][1] * PrmaPxie[0][2], 2.0) 
		+ std::pow(PrmaPxie[0][2] * PrmaPxie[1][0] - PrmaPxie[1][2] * PrmaPxie[0][0], 2.0) 
		+ std::pow(PrmaPxie[0][0] * PrmaPxie[1][1] - PrmaPxie[1][0] * PrmaPxie[0][1], 2.0));
    //1.0E-6???
    if(slaveXiet[0] < -1.0 - minXietaDiff || slaveXiet[0] > 1.0 + minXietaDiff 
		|| slaveXiet[1] < -1.0 - minXietaDiff || slaveXiet[1] > 1.0 + minXietaDiff){
        Log("            Warning 1 in ContactInterface::ProjectMB2S! slaveXiet = " 
            + Double2String(slaveXiet[0]) + ", " + Double2String(slaveXiet[1]));
    }
    maslPoint[1].fill(0.0);
    const BilinearQuadrature<2>& biliQuad = GetBilinearQuadrature();
    for(I64 ti = 0; ti < 4; ++ ti){
        maslShape[1][ti] = (1.0 + biliQuad.cornerNodes[ti][0] * slaveXiet[0]) 
			* (1.0 + biliQuad.cornerNodes[ti][1] * slaveXiet[1]) / 4.0;
        AXPY(maslShape[1][ti], slaveCorners[ti], maslPoint[1]);
    }
}

Real CSQuadSegment::TriangleArea2d(std::array<Real, 2> tempPoint_0, 
    std::array<Real, 2> tempPoint_1, std::array<Real, 2> tempPoint_2){
    //
    AXPY(-1.0, tempPoint_0, tempPoint_1);
    AXPY(-1.0, tempPoint_0, tempPoint_2);
    return std::abs(tempPoint_1[0] * tempPoint_2[1] - tempPoint_1[1] * tempPoint_2[0]) / 2.0;
}

void CSQuadSegment::TriangleQuadrature(std::array<Real, 2> tempXiet_0, 
    std::array<Real, 2> tempXiet_1, std::array<Real, 2> tempXiet_2, 
    std::vector<std::array<Real, 2>> &listXiet, std::vector<double> &listWeight){
    //
    Real area = TriangleArea2d(tempXiet_0, tempXiet_1, tempXiet_2);
    const TriangularQuadrature<2>& triaQuad = GetTriangularQuadrature();
    for(I64 ti = 0; ti < triaQuad.numbGaussPoints; ++ ti){
        std::array<Real, 2> tempXiet = {0.0, 0.0};
        AXPY(triaQuad.gaussPoints[ti][0], tempXiet_0, tempXiet);
        AXPY(triaQuad.gaussPoints[ti][1], tempXiet_1, tempXiet);
        AXPY(triaQuad.gaussPoints[ti][2], tempXiet_2, tempXiet);
        listXiet.emplace_back(tempXiet);
		//integral variable transformation: from Cartesian coordinate to shape function
		//x:xmin~xmax, y:ymin(x)~ymax(x)
		//A1:0~1, A2:0~1-A1, A3=1-A0-A1
        listWeight.emplace_back(2.0 * area * triaQuad.weights[ti]);
    }
}

void CSQuadSegment::SortBy2d(std::array<Real, 2> &tempPoint_0, 
    std::array<Real, 2> &tempPoint_1, I64 tempIndex){
    //
    if(tempPoint_0[tempIndex] > tempPoint_1[tempIndex]){
        std::swap(tempPoint_0, tempPoint_1);
    }
}

bool CSQuadSegment::IsCross2d(
    const std::array<Real, 2>& tempPoint_0, const std::array<Real, 2>& tempPoint_1, 
    const std::array<Real, 2>& tempPoint_2, const std::array<Real, 2>& tempPoint_3){
    //
    if(std::max(tempPoint_0[0], tempPoint_1[0]) < std::min(tempPoint_2[0], tempPoint_3[0])
		|| std::max(tempPoint_0[1], tempPoint_1[1]) < std::min(tempPoint_2[1], tempPoint_3[1])
		|| std::min(tempPoint_0[0], tempPoint_1[0]) > std::max(tempPoint_2[0], tempPoint_3[0])
		|| std::min(tempPoint_0[1], tempPoint_1[1]) > std::max(tempPoint_2[1], tempPoint_3[1])){
        return false;
    }
	if(((tempPoint_2[0] - tempPoint_0[0]) * (tempPoint_2[1] - tempPoint_3[1]) 
		- (tempPoint_2[1] - tempPoint_0[1]) * (tempPoint_2[0] - tempPoint_3[0])) 
		* ((tempPoint_2[0] - tempPoint_1[0]) * (tempPoint_2[1] - tempPoint_3[1]) 
		- (tempPoint_2[1] - tempPoint_1[1]) * (tempPoint_2[0] - tempPoint_3[0])) <= 0
		&& ((tempPoint_0[0] - tempPoint_2[0]) * (tempPoint_0[1] - tempPoint_1[1]) 
		- (tempPoint_0[1] - tempPoint_2[1]) * (tempPoint_0[0] - tempPoint_1[0])) 
		* ((tempPoint_0[0] - tempPoint_3[0]) * (tempPoint_0[1] - tempPoint_1[1]) 
		- (tempPoint_0[1] - tempPoint_3[1]) * (tempPoint_0[0] - tempPoint_1[0])) <= 0){
		return true;
	}
	else{
		return false;
	}
}

void CSQuadSegment::LineIntersection2d(
    std::array<Real, 2> tempPoint_0, std::array<Real, 2> tempPoint_1, 
    std::array<Real, 2> tempPoint_2, std::array<Real, 2> tempPoint_3, 
    std::vector<std::array<Real, 2>> &resultIntersections){
    //
    if(!IsCross2d(tempPoint_0, tempPoint_1, tempPoint_2, tempPoint_3)){
        return;
    }
    Real area_2 = std::abs(TriangleArea2d(tempPoint_2, tempPoint_0, tempPoint_1));
    Real area_3 = std::abs(TriangleArea2d(tempPoint_3, tempPoint_0, tempPoint_1));
    if(area_2 < minArea && area_3 < minArea){//co-linear
        if(std::abs(tempPoint_0[0] - tempPoint_1[0]) > minXietaDiff){
            SortBy2d(tempPoint_0, tempPoint_1, 0);
            SortBy2d(tempPoint_2, tempPoint_3, 0);
            Real from_x = tempPoint_0[0];
            Real from_y = tempPoint_0[1];
            if(tempPoint_0[0] < tempPoint_2[0]){
                from_x = tempPoint_2[0];
                from_y = tempPoint_2[1];
            }
            Real to_x = tempPoint_1[0];
            Real to_y = tempPoint_1[1];
            if(tempPoint_1[0] > tempPoint_3[0]){
                to_x = tempPoint_3[0];
                to_y = tempPoint_3[1];
            }
            if(std::abs(from_x - to_x) < minXietaDiff){
                resultIntersections.emplace_back(std::array<Real, 2>{from_x, from_y});
            }
            else{
                resultIntersections.emplace_back(std::array<Real, 2>{from_x, from_y});
                resultIntersections.emplace_back(std::array<Real, 2>{to_x, to_y});
            }
        }
        else{
            SortBy2d(tempPoint_0, tempPoint_1, 1);
            SortBy2d(tempPoint_2, tempPoint_3, 1);
            Real from_x = tempPoint_0[0];
            Real from_y = tempPoint_0[1];
            if(tempPoint_0[1] < tempPoint_2[1]){
                from_x = tempPoint_2[0];
                from_y = tempPoint_2[1];
            }
            Real to_x = tempPoint_1[0];
            Real to_y = tempPoint_1[1];
            if(tempPoint_1[1] > tempPoint_3[1]){
                to_x = tempPoint_3[0];
                to_y = tempPoint_3[1];
            }
            if(std::abs(from_x - to_x) < minXietaDiff){
                resultIntersections.emplace_back(std::array<Real, 2>{from_x, from_y});
            }
            else{
                resultIntersections.emplace_back(std::array<Real, 2>{from_x, from_y});
                resultIntersections.emplace_back(std::array<Real, 2>{to_x, to_y});
            }
        }
    }
    else if(area_2 < minArea){//one endpoint lies on the another line-segment
        resultIntersections.emplace_back(tempPoint_2);
    }
    else if(area_3 < minArea){//one endpoint lies on the another line-segment
        resultIntersections.emplace_back(tempPoint_3);
    }
    else{ // true intersect
        Real tempFactor = area_2 / area_3;
        std::array<Real, 2> tempPoint = {
            (tempPoint_2[0] + tempFactor * tempPoint_3[0]) / (1.0 + tempFactor), 
            (tempPoint_2[1] + tempFactor * tempPoint_3[1]) / (1.0 + tempFactor)};
        resultIntersections.emplace_back(tempPoint);
    }
}

bool CSQuadSegment::InQuadrilateral(const std::array<Real, 2>& tempPoint, 
    const std::array<std::array<Real, 2>, 4>& tempCorners){
	Real subqArea = 0.0;
	for(I64 ti = 0; ti < 4; ++ ti){
		subqArea += TriangleArea2d(tempPoint, tempCorners[ti], tempCorners[(ti + 1) % 4]);
	}
	Real totalArea = TriangleArea2d(tempCorners[0], tempCorners[1], tempCorners[2]) + 
		TriangleArea2d(tempCorners[2], tempCorners[3], tempCorners[0]);
	if(subqArea <= (1.0 + 1.0E-12) * totalArea){
		return true;
	}
	else{
		return false;
	}
}

void CSQuadSegment::SegmentIntersect(
    std::vector<std::array<Real, 2>> &listXiet, 
    std::vector<double> &listWeight){
    // projection
    std::array< std::array<Real, 2>, 4> masterProjection = {{
        {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}}};
    std::array< std::array<Real, 2>, 4> slaveProjection;
    std::array<std::array<Real, 3>, 2> PrmaPxie;
    Real ngap;
    for(I64 ti = 0; ti < 4; ++ ti){
        ProjectS2M(slaveCorners[ti], PrmaPxie, slaveProjection[ti], ngap);
    }
    // intersection
    std::vector<std::array<Real, 2>> tempInte_0, tempInte_1, tempInte_2;
    for(I64 ti = 0; ti < 4; ++ ti){
        if(InQuadrilateral(slaveProjection[ti], masterProjection)){
            tempInte_0.emplace_back(slaveProjection[ti]);
        }
        if(InQuadrilateral(masterProjection[ti], slaveProjection)){
            tempInte_0.emplace_back(masterProjection[ti]);
        }
    }
    for(I64 ti = 0; ti < 4; ++ ti){
        for(I64 tj = 0; tj < 4; ++ tj){
            LineIntersection2d(masterProjection[ti], masterProjection[(ti + 1) % 4], 
                slaveProjection[tj], slaveProjection[(tj + 1) % 4], tempInte_0
            );
        }
    }
    I64 teinSize_0 = tempInte_0.size();
    if(teinSize_0 < 3){
		return;
	}
    // no repeat
    std::vector<I64> indices_0(teinSize_0);
    std::iota(indices_0.begin(), indices_0.end(), 0);
    std::sort(std::execution::unseq, indices_0.begin(), indices_0.end(), 
        [&tempInte_0](I64 ti, I64 tj)->bool{
            if(tempInte_0[ti][0] < tempInte_0[tj][0] - minXietaDiff){
                return true;
            }
            else if(tempInte_0[ti][0] <= tempInte_0[tj][0] + minXietaDiff){
                if(tempInte_0[ti][1] < tempInte_0[tj][1] - minXietaDiff){
                    return true;
                }
                else{
                    return false;
                }
            }
            else{
                return false;
            }
        });
    tempInte_1.emplace_back(tempInte_0[indices_0[0]]);
    for(I64 ti = 1; ti < teinSize_0; ++ ti){
        const std::array<Real, 2>& tiInte_0 = tempInte_0[indices_0[ti]];
        const std::array<Real, 2>& timoInte_0 = tempInte_0[indices_0[ti - 1]];
        if(std::abs(tiInte_0[0] - timoInte_0[0]) > minXietaDiff 
            || std::abs(tiInte_0[1] - timoInte_0[1]) > minXietaDiff){
            tempInte_1.emplace_back(tiInte_0);
        }
    }
    // sort
    std::array<Real, 2> tempCenter = {0.0, 0.0};
    for(const auto& iterTein_1 : tempInte_1){
        tempCenter[0] += iterTein_1[0];
        tempCenter[1] += iterTein_1[1];
    }
    I64 teinSize_1 = tempInte_1.size();
    tempCenter[0] /= teinSize_1;
    tempCenter[1] /= teinSize_1;
    std::vector<Real> tempAngle(teinSize_1);
    for(I64 ti = 0; ti < teinSize_1; ++ ti){
        tempAngle[ti] = std::atan2(tempInte_1[ti][1] - tempCenter[1], tempInte_1[ti][0] - tempCenter[0]);
    }
    std::vector<I64> indices_1(teinSize_1);
    std::iota(indices_1.begin(), indices_1.end(), 0);
    std::sort(std::execution::unseq, indices_1.begin(), indices_1.end(), 
        [&tempAngle](I64 ti, I64 tj)->bool{
            if(tempAngle[ti] < tempAngle[tj]){
                return true;
            }
            else{
                return false;
            }
        });
    tempInte_2.resize(teinSize_1);
    std::transform(std::execution::unseq, indices_1.begin(), indices_1.end(), 
        tempInte_2.begin(), 
        [&tempInte_1](I64 ti)->std::array<Real, 2>{
            return tempInte_1[ti];
        });
	//
	//centroid and area
	//R. Nurnberg. Calculating the area and centroid of a polygon in 2d.
	//https://paulbourke.net/geometry/polygonmesh/centroid.pdf
	double area = 0.0;
	tempCenter = {0.0, 0.0};
	for(I64 ti = 0; ti < teinSize_1; ti ++){
		area += tempInte_2[ti][0] * tempInte_2[(ti + 1) % teinSize_1][1] 
			- tempInte_2[(ti + 1) % teinSize_1][0] * tempInte_2[ti][1];
		tempCenter[0] += (tempInte_2[ti][0] + tempInte_2[(ti + 1) % teinSize_1][0]) 
			* (tempInte_2[ti][0] * tempInte_2[(ti + 1) % teinSize_1][1] 
			- tempInte_2[(ti + 1) % teinSize_1][0] * tempInte_2[ti][1]);
		tempCenter[1] += (tempInte_2[ti][1] + tempInte_2[(ti + 1) % teinSize_1][1]) 
			* (tempInte_2[ti][0] * tempInte_2[(ti + 1) % teinSize_1][1] 
			- tempInte_2[(ti + 1) % teinSize_1][0] * tempInte_2[ti][1]);
	}
	area /= 2.0;
	if(std::abs(area) <= minArea){
		return;
	}
    SCAL(0.166666666666666666 / area, tempCenter);
	for(I64 ti = 0; ti < teinSize_1; ti ++){
		TriangleQuadrature(tempCenter, tempInte_2[ti], 
			tempInte_2[(ti + 1) % teinSize_1], listXiet, listWeight);
	}
}

void CSQuadSegment::Search(std::array<I64, 4> tempMasterSegment, std::array<I64, 4> tempSlaveSegment){
    //
    std::vector<std::array<Real, 2>> listXiet;
    std::vector<double> listWeight;
    SegmentIntersect(listXiet, listWeight);
	//shape function, initial gap
    const BilinearQuadrature<2>& biliQuad = GetBilinearQuadrature();
    I64 lixiSize = listXiet.size();
	for(I64 ti = 0; ti < lixiSize; ++ ti){
        IntegralPoint tempInpo;
        tempInpo.node[0] = tempMasterSegment;
        tempInpo.node[1] = tempSlaveSegment;
        tempInpo.contactPoint[0].fill(0.0);
        for(I64 tj = 0; tj < 4; ++ tj){
            tempInpo.shapeFunction[0][tj] = (1.0 + biliQuad.cornerNodes[tj][0] * listXiet[ti][0]) 
				* (1.0 + biliQuad.cornerNodes[tj][1] * listXiet[ti][1]) / 4.0;
            AXPY(tempInpo.shapeFunction[0][tj], masterCorners[tj], tempInpo.contactPoint[0]);
        }
        ProjectMB2S(listXiet[ti], tempInpo.shapeFunction, tempInpo.contactPoint, 
            tempInpo.basisVector, tempInpo.quadratureWeight);
        std::array<Real, 3> master2slave = tempInpo.contactPoint[1];
        AXPY(-1.0, tempInpo.contactPoint[0], master2slave);
        tempInpo.initialNormalGap = DOT(tempInpo.basisVector[0], master2slave);
        tempInpo.quadratureWeight *= listWeight[ti];
        integralPoints.emplace_back(tempInpo);
    }
}

} // namespace Ddpca