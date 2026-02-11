#include "CurvedSurface.hpp"

namespace Ddpca{

void CurvedSurface::Resize(I64 rows, I64 cols){
    indexPoint.resize(rows);
    for(I64 ti = 0; ti < rows; ++ ti){
        indexPoint[ti].resize(cols);
    }
}

void CurvedSurface::Insert(I64 ti, I64 tj, Coordinate inputCoordinate){
    //not considering duplicated key
    indexPoint[ti][tj] = inputCoordinate;
    pointIndex.emplace(inputCoordinate, std::array<I64, 2>{ti, tj});
}

bool CurvedSurface::QuadrilateralRefine(
    const std::vector<Coordinate>& inputCoordinate, Coordinate &outputCoordinate){
    //
    bool tempFlag = true;
    std::array<I64, 2> resultIndex = {0, 0};
    const I64 incoSize = inputCoordinate.size();
    for(I64 ti = 0; ti < incoSize; ti ++){
        auto iterPoin = pointIndex.find(inputCoordinate[ti]);
        if(iterPoin == pointIndex.end()){
            tempFlag = false;
            break;
        }
        resultIndex[0] += (iterPoin->second)[0];
        resultIndex[1] += (iterPoin->second)[1];
    }
    if(tempFlag == false){
        return false;
    }
    resultIndex[0] /= incoSize;
    resultIndex[1] /= incoSize;
    outputCoordinate = indexPoint[resultIndex[0]][resultIndex[1]];
    return true;
}

void CurvedSurface::Refine(
    const Mesh& tempMesh, const std::set<I64>& elementsToSplit, 
    std::map<std::vector<I64>, Coordinate>& curvInte){
    //
    for(const auto currentIndex : elementsToSplit){
        const OctreeElement& tempElement = tempMesh.elements[currentIndex];
        //
        const I64 heliSize = hexaLine.size();
        for(I64 tj = 0; tj < heliSize; ++ tj){
            const I64 tj_size = hexaLine[tj].size();
            std::vector<I64> inputNode(tj_size);
            std::vector<Coordinate> inputCoor(tj_size);
            for(I64 tk = 0; tk < tj_size; ++ tk){
                inputNode[tk] = tempElement.cornerNodes[hexaLine[tj][tk]];
                auto iterNoco = tempMesh.node2Coordinate.find(inputNode[tk]);
                inputCoor[tk] = iterNoco->second;
            }
            Coordinate outputCoor;
            if(QuadrilateralRefine(inputCoor, outputCoor) == false){
                continue;
            }
            std::sort(std::execution::unseq, inputNode.begin(), inputNode.end());
            curvInte.emplace(inputNode, outputCoor);
        }
        //
        const I64 hefaSize = hexaFace.size();
        for(I64 tj = 0; tj < hefaSize; ++ tj){
            const I64 tj_size = hexaFace[tj].size();
            std::vector<I64> inputNode(tj_size);
            std::vector<Coordinate> inputCoor(tj_size);
            for(I64 tk = 0; tk < tj_size; ++ tk){
                inputNode[tk] = tempElement.cornerNodes[hexaFace[tj][tk]];
                auto iterNoco = tempMesh.node2Coordinate.find(inputNode[tk]);
                inputCoor[tk] = iterNoco->second;
            }
            Coordinate outputCoor;
            if(QuadrilateralRefine(inputCoor, outputCoor) == false){
                continue;
            }
            std::sort(std::execution::unseq, inputNode.begin(), inputNode.end());
            curvInte.emplace(inputNode, outputCoor);
        }
    }
}

void CurvedSurface::RigidRotationTranslation(
    const DenseMatrix& rotationMatrix, const Coordinate& translationVector){
    //
    pointIndex.clear();
    const I64 indexSize_0 = indexPoint.size();
    for(I64 ti = 0; ti < indexSize_0; ++ ti){
        const I64 indexSize_1 = indexPoint[ti].size();
        for(I64 tj = 0; tj < indexSize_1; ++ tj){
            Coordinate tempCoor = translationVector;
            PEMV(rotationMatrix, indexPoint[ti][tj].data, tempCoor.data);
            indexPoint[ti][tj] = tempCoor;
            pointIndex.emplace(tempCoor, std::array<I64, 2>{ti, tj});
        }
    }
}

void CurvedSurface::Initialize(){
    ti = 0;
    tj = 0;
}

bool CurvedSurface::Increment(const Mesh& tempMesh){
    //
    const I64 elemSize = tempMesh.elements.size();
    for(; ti < elemSize; ++ ti){
        const OctreeElement& tempElement = tempMesh.elements[ti];
        if(tempElement.children.size() > 0){
            tj = 0;
            continue;
        }
        const I64 hefaSize = hexaFace.size();
        for(; tj < hefaSize; ++ tj){
            const I64 tj_size = hexaFace[tj].size();
            bool tempFlag = true;
            for(I64 tk = 0; tk < tj_size; ++ tk){
                currentFace[tk] = tempElement.cornerNodes[hexaFace[tj][tk]];
                auto iterNoco = tempMesh.node2Coordinate.find(currentFace[tk]);
                auto iterPoin = pointIndex.find(iterNoco->second);
                if(iterPoin == pointIndex.end()){
                    tempFlag = false;
                    break;
                }
            }
            if(tempFlag){
                ++ tj;
                if(tj == hefaSize){
                    tj = 0;
                    ++ ti;
                }
                return true;
            }
        }
        tj = 0;
    }
    return false;
}

} // namespace Ddpca