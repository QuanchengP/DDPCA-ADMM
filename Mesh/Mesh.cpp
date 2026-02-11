#include "Mesh.hpp"
#include "TrilinearQuadrature.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <execution>

namespace Ddpca {

/*************************************************************************************************/
I64 Mesh::TryAddNode(const Coordinate& currentCoordinate) {
	//
	const I64 tempSize = node2Coordinate.size();
    // Use the return value of emplace to get the insertion result directly, avoiding secondary lookup
    auto [iterCono, inserted] = coordinate2Node.try_emplace(
		currentCoordinate, tempSize);
    if (inserted) {
        // Only update node2Coordinate when a new coordinate is successfully inserted
        node2Coordinate.emplace(iterCono->second, currentCoordinate);
    }
    return iterCono->second;
}

I64 Mesh::AddElement(const OctreeElement& currentElement) {
    elements.emplace_back(currentElement);
    I64 elemNumb = elements.size() - 1;
    // Process lines - Use try_emplace to reduce one lookup operation
    I64 heliSize = hexaLine.size();
    for(I64 ti = 0; ti < heliSize; ti++) {
        I64 node1 = currentElement.cornerNodes[hexaLine[ti][0]];
        I64 node2 = currentElement.cornerNodes[hexaLine[ti][1]];
        
        std::array<I64, 2> line = {node1, node2};
		std::sort(line.begin(), line.end());
        
        // Use try_emplace to perform both lookup and insertion simultaneously
        auto [iteratorLube, inserted] = lineUsedByElement.try_emplace(line, std::set<I64>());
        // Need to add elemNumb regardless of insertion, so merge both branches
        iteratorLube->second.emplace(elemNumb);
    }
    // Process faces - Use try_emplace to reduce one lookup operation
    I64 hefaSize = hexaFace.size();
    for(I64 ti = 0; ti < hefaSize; ti++) {
        I64 n0 = currentElement.cornerNodes[hexaFace[ti][0]];
        I64 n1 = currentElement.cornerNodes[hexaFace[ti][1]];
        I64 n2 = currentElement.cornerNodes[hexaFace[ti][2]];
        I64 n3 = currentElement.cornerNodes[hexaFace[ti][3]];
        std::array<I64, 4> face = {n0, n1, n2, n3};
		std::sort(face.begin(), face.end());
        
        auto [iteratorFube, inserted] = faceUsedByElement.try_emplace(face, std::set<I64>());
        iteratorFube->second.emplace(elemNumb);
    }
    // Update maximum level
    if (currentElement.level > maxLevel) {
        maxLevel = currentElement.level;
    }
    return elemNumb;
}

void Mesh::Refine(std::set<I64> &elementsToSplit, 
	    const std::map<I64, std::set<I64>> &subElements, 
	    const std::map<std::vector<I64>, Coordinate> &curvilinearInterpolation){
	//
	AdjacentLevelCheck(elementsToSplit);
	//
	//Refinement template, 8 original nodes, 12 edge centers, 6 face centers, 1 element center
	static const std::vector<std::vector<std::vector<I64>>> refinementTemplate_1 = {{
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, 
		{2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, 
		{0, 3, 7, 4}, {1, 2, 6, 5}, {0, 4, 5, 1}, {3, 7, 6, 2}, {0, 1, 2, 3}, {4, 5, 6, 7}, 
		{0, 1, 2, 3, 4, 5, 6, 7}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 1}, {1, 2}, {2, 3}, {3, 0}, 
		{4, 5}, {5, 6}, {6, 7}, {7, 4}, 
		{0, 1, 2, 3}, {4, 5, 6, 7}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 3}, {3, 7}, {7, 4}, {4, 0}, 
		{1, 2}, {2, 6}, {6, 5}, {5, 1}, 
		{0, 3, 7, 4}, {1, 2, 6, 5}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 4}, {4, 5}, {5, 1}, {1, 0}, 
		{3, 7}, {7, 6}, {6, 2}, {2, 3}, 
		{0, 4, 5, 1}, {3, 7, 6, 2}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 1}, {2, 3}, 
		{4, 5}, {6, 7}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 3}, {1, 2}, 
		{4, 7}, {5, 6}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 4}, {1, 5}, 
		{3, 7}, {2, 6}
	}};
	//8 + 12 + 6 + 1 nodes
	static const std::vector<std::vector<std::vector<I64>>> refinementTemplate_2 = {{
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {1, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 1, 0}, 
		{0, 0, 1}, {2, 0, 1}, {2, 2, 1}, {0, 2, 1}, {1, 0, 2}, {2, 1, 2}, 
		{1, 2, 2}, {0, 1, 2}, {0, 1, 1}, {2, 1, 1}, {1, 0, 1}, {1, 2, 1}, 
		{1, 1, 0}, {1, 1, 2}, {1, 1, 1}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {1, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 1, 0}, 
		{1, 0, 2}, {2, 1, 2}, {1, 2, 2}, {0, 1, 2}, {1, 1, 0}, {1, 1, 2}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {0, 1, 0}, {0, 2, 1}, {0, 1, 2}, {0, 0, 1}, 
		{2, 1, 0}, {2, 2, 1}, {2, 1, 2}, {2, 0, 1}, {0, 1, 1}, {2, 1, 1}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {0, 0, 1}, {1, 0, 2}, {2, 0, 1}, {1, 0, 0}, 
		{0, 2, 1}, {1, 2, 2}, {2, 2, 1}, {1, 2, 0}, {1, 0, 1}, {1, 2, 1}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {1, 0, 0}, {1, 2, 0}, {1, 0, 2}, {1, 2, 2}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {0, 1, 0}, {2, 1, 0}, {0, 1, 2}, {2, 1, 2}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {0, 0, 1}, {2, 0, 1}, {0, 2, 1}, {2, 2, 1}
	}};

	//refined element 8 sub-elements, 8 nodes, 3 coordinates
	static const std::vector<std::vector<std::vector<std::vector<I64>>>> refinedElement = {{
		{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}}, 
		{{1,0,0}, {2,0,0}, {2,1,0}, {1,1,0}, {1,0,1}, {2,0,1}, {2,1,1}, {1,1,1}}, 
		{{0,1,0}, {1,1,0}, {1,2,0}, {0,2,0}, {0,1,1}, {1,1,1}, {1,2,1}, {0,2,1}}, 
		{{1,1,0}, {2,1,0}, {2,2,0}, {1,2,0}, {1,1,1}, {2,1,1}, {2,2,1}, {1,2,1}}, 
		{{0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {0,0,2}, {1,0,2}, {1,1,2}, {0,1,2}}, 
		{{1,0,1}, {2,0,1}, {2,1,1}, {1,1,1}, {1,0,2}, {2,0,2}, {2,1,2}, {1,1,2}}, 
		{{0,1,1}, {1,1,1}, {1,2,1}, {0,2,1}, {0,1,2}, {1,1,2}, {1,2,2}, {0,2,2}}, 
		{{1,1,1}, {2,1,1}, {2,2,1}, {1,2,1}, {1,1,2}, {2,1,2}, {2,2,2}, {1,2,2}}
	}, {
		{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,2}, {1,0,2}, {1,1,2}, {0,1,2}}, 
		{{1,0,0}, {2,0,0}, {2,1,0}, {1,1,0}, {1,0,2}, {2,0,2}, {2,1,2}, {1,1,2}}, 
		{{0,1,0}, {1,1,0}, {1,2,0}, {0,2,0}, {0,1,2}, {1,1,2}, {1,2,2}, {0,2,2}}, 
		{{1,1,0}, {2,1,0}, {2,2,0}, {1,2,0}, {1,1,2}, {2,1,2}, {2,2,2}, {1,2,2}}
	}, {
		{{0,0,0}, {2,0,0}, {2,1,0}, {0,1,0}, {0,0,1}, {2,0,1}, {2,1,1}, {0,1,1}}, 
		{{0,1,0}, {2,1,0}, {2,2,0}, {0,2,0}, {0,1,1}, {2,1,1}, {2,2,1}, {0,2,1}}, 
		{{0,0,1}, {2,0,1}, {2,1,1}, {0,1,1}, {0,0,2}, {2,0,2}, {2,1,2}, {0,1,2}}, 
		{{0,1,1}, {2,1,1}, {2,2,1}, {0,2,1}, {0,1,2}, {2,1,2}, {2,2,2}, {0,2,2}}
	}, {
		{{0,0,0}, {1,0,0}, {1,2,0}, {0,2,0}, {0,0,1}, {1,0,1}, {1,2,1}, {0,2,1}}, 
		{{0,0,1}, {1,0,1}, {1,2,1}, {0,2,1}, {0,0,2}, {1,0,2}, {1,2,2}, {0,2,2}}, 
		{{1,0,0}, {2,0,0}, {2,2,0}, {1,2,0}, {1,0,1}, {2,0,1}, {2,2,1}, {1,2,1}}, 
		{{1,0,1}, {2,0,1}, {2,2,1}, {1,2,1}, {1,0,2}, {2,0,2}, {2,2,2}, {1,2,2}}
	}, {
		{{0,0,0}, {1,0,0}, {1,2,0}, {0,2,0}, {0,0,2}, {1,0,2}, {1,2,2}, {0,2,2}}, 
		{{1,0,0}, {2,0,0}, {2,2,0}, {1,2,0}, {1,0,2}, {2,0,2}, {2,2,2}, {1,2,2}}
	}, {
		{{0,0,0}, {2,0,0}, {2,1,0}, {0,1,0}, {0,0,2}, {2,0,2}, {2,1,2}, {0,1,2}}, 
		{{0,1,0}, {2,1,0}, {2,2,0}, {0,2,0}, {0,1,2}, {2,1,2}, {2,2,2}, {0,2,2}}
	}, {
		{{0,0,0}, {2,0,0}, {2,2,0}, {0,2,0}, {0,0,1}, {2,0,1}, {2,2,1}, {0,2,1}}, 
		{{0,0,1}, {2,0,1}, {2,2,1}, {0,2,1}, {0,0,2}, {2,0,2}, {2,2,2}, {0,2,2}}
	}};
	std::set<I64> elementsToSplit_0;
	for(const auto currentIndex : elementsToSplit){
		const std::array<I64, 8>& beforeCorners = elements[currentIndex].cornerNodes;
		const I64 pattern = elements[currentIndex].refinementPattern;
		//
		std::array<std::array<std::array<I64, 3>, 3>, 3> AfterCorners = {};
		const I64 rt1pSize = refinementTemplate_1[pattern].size();
		for(I64 tk = 0; tk < rt1pSize; tk ++){
			const I64 index_0 = refinementTemplate_2[pattern][tk][0];
			const I64 index_1 = refinementTemplate_2[pattern][tk][1];
			const I64 index_2 = refinementTemplate_2[pattern][tk][2];
			//
			const I64 rt1pkSize = refinementTemplate_1[pattern][tk].size();
			if(rt1pkSize == 1){
				const I64 indexBefore = refinementTemplate_1[pattern][tk][0];
				AfterCorners[index_0][index_1][index_2] = beforeCorners[indexBefore];
				continue;
			}
			//
			std::vector<I64> originalNode(rt1pkSize);
			for(I64 ti = 0; ti < rt1pkSize; ti ++){
				originalNode[ti] = beforeCorners[refinementTemplate_1[pattern][tk][ti]];
			}
			sort(originalNode.begin(), originalNode.end());
			auto iteratorCurvInte = curvilinearInterpolation.find(originalNode);
			if(iteratorCurvInte != curvilinearInterpolation.end()){
				AfterCorners[index_0][index_1][index_2] = TryAddNode(iteratorCurvInte->second);
			}
			else{
				Coordinate newCoordinate(0.0, 0.0, 0.0);
				for(I64 ti = 0; ti < rt1pkSize; ti ++){
					XPEY(newCoordinate.data, node2Coordinate[originalNode[ti]].data);
				}
				SCAL(1.0 / rt1pkSize, newCoordinate.data);
				AfterCorners[index_0][index_1][index_2] = TryAddNode(newCoordinate);
			}
		}
		//
		const I64 reSize = refinedElement[pattern].size();
		elements[currentIndex].children.resize(reSize);
		const I64 newLevel = elements[currentIndex].level + 1;
		for(I64 ti = 0; ti < reSize; ti ++){
			OctreeElement newElement;
			newElement.parent = currentIndex;
			std::array<I64, 8> newCorners;
			for (I64 tj = 0; tj < 8; ++tj) {
				const auto& cornerIndices = refinedElement[pattern][ti][tj];
				newCorners[tj] = AfterCorners[cornerIndices[0]][cornerIndices[1]][cornerIndices[2]];
			}
			newElement.cornerNodes = std::move(newCorners);
			newElement.level = newLevel;
			newElement.refinementPattern = OctreeElement::REFINEMENT_NONE;
			newElement.children.resize(0);
			I64 newIndex = AddElement(newElement);
			elements[currentIndex].children[ti] = newIndex;
		}
		//
		auto iteratorSubElements = subElements.find(currentIndex);
		if(iteratorSubElements == subElements.end()){
			continue;
		}
		for(const auto subElement : (iteratorSubElements->second)){
			elementsToSplit_0.emplace(elements[currentIndex].children[subElement]);
		}
	}
	elementsToSplit.swap(elementsToSplit_0);
}

void Mesh::AdjacentLevelCheck(std::set<I64> &elementsToSplit){
	//[1] The level difference of two neighboring elements is less than or equal to 1:
	//when refine a line L of spliElem: spliElem's parent does not have this line L;
	//denote spliElem's parent as P, the line L is a sub segment of one line in the 12 lines of P;
	//P's neighbor should be refined in the same way of as P;
	//[2] hanging node can only be used under refiPatt = 0, 
	// better to comply, must be very careful when violate
	//
	//when the 0th sub-element of the p-th element is refined, line {0,1/3/4} should be dealed with
	static const std::vector<std::vector<std::vector<std::array<I64, 2>>>> parentLine = {{
		{{0, 1}, {0, 3}, {0, 4}}, {{1, 0}, {1, 2}, {1, 5}}, 
		{{3, 0}, {3, 2}, {3, 7}}, {{2, 1}, {2, 3}, {2, 6}}, 
		{{4, 0}, {4, 5}, {4, 7}}, {{5, 1}, {5, 4}, {5, 6}}, 
		{{7, 3}, {7, 6}, {7, 4}}, {{6, 2}, {6, 5}, {6, 7}}
	}, {
		{{0, 1}, {0, 3}, {4, 5}, {4, 7}}, 
		{{1, 0}, {1, 2}, {5, 4}, {5, 6}}, 
		{{3, 0}, {3, 2}, {7, 4}, {7, 6}}, 
		{{2, 1}, {2, 3}, {6, 5}, {6, 7}}
	}, {
		{{0, 3}, {0, 4}, {1, 2}, {1, 5}}, 
		{{3, 0}, {3, 7}, {2, 1}, {2, 6}}, 
		{{4, 0}, {4, 7}, {5, 1}, {5, 6}}, 
		{{7, 3}, {7, 4}, {6, 2}, {6, 5}}
	}, {
		{{0, 1}, {0, 4}, {3, 2}, {3, 7}}, 
		{{4, 0}, {4, 5}, {7, 3}, {7, 6}}, 
		{{1, 0}, {1, 5}, {2, 3}, {2, 6}}, 
		{{5, 1}, {5, 4}, {6, 2}, {6, 7}}
	}, {
		{{0, 1}, {4, 5}, {2, 3}, {7, 6}}, 
		{{0, 1}, {4, 5}, {2, 3}, {7, 6}}
	}, {
		{{0, 3}, {1, 2}, {4, 7}, {5, 6}}, 
		{{0, 3}, {1, 2}, {4, 7}, {5, 6}}
	}, {
		{{0, 4}, {1, 5}, {3, 7}, {2, 6}}, 
		{{0, 4}, {1, 5}, {3, 7}, {2, 6}}
	}};
	static const std::vector<std::vector<std::vector<std::array<I64, 4>>>> parentFace = {{
		{{0, 1, 2, 3}, {0, 3, 7, 4}, {0, 4, 5, 1}}, 
		{{1, 2, 3, 0}, {1, 2, 6, 5}, {1, 0, 4, 5}}, 
		{{3, 0, 1, 2}, {3, 7, 4, 0}, {3, 7, 6, 2}}, 
		{{2, 3, 0, 1}, {2, 6, 5, 1}, {2, 3, 7, 6}}, 
		{{4, 5, 6, 7}, {4, 0, 3, 7}, {4, 5, 1, 0}}, 
		{{5, 6, 7, 4}, {5, 1, 2, 6}, {5, 1, 0, 4}}, 
		{{7, 4, 5, 6}, {7, 4, 0, 3}, {7, 6, 2, 3}}, 
		{{6, 7, 4, 5}, {6, 5, 1, 2}, {6, 2, 3, 7}}
	}, {
		{{0, 1, 2, 3}, {0, 3, 7, 4}, {0, 4, 5, 1}, {4, 5, 6, 7}}, 
		{{1, 2, 3, 0}, {1, 2, 6, 5}, {1, 0, 4, 5}, {5, 6, 7, 4}}, 
		{{3, 0, 1, 2}, {3, 7, 4, 0}, {3, 7, 6, 2}, {7, 4, 5, 6}}, 
		{{2, 3, 0, 1}, {2, 6, 5, 1}, {2, 3, 7, 6}, {6, 7, 4, 5}}
	}, {
		{{0, 3, 7, 4}, {1, 2, 6, 5}, {0, 4, 5, 1}, {0, 1, 2, 3}}, 
		{{3, 7, 4, 0}, {2, 6, 5, 1}, {3, 7, 6, 2}, {3, 0, 1, 2}}, 
		{{4, 0, 3, 7}, {5, 1, 2, 6}, {4, 5, 6, 7}, {4, 5, 1, 0}}, 
		{{7, 4, 0, 3}, {6, 5, 1, 2}, {7, 6, 2, 3}, {7, 4, 5, 6}}
	}, {
		{{0, 4, 5, 1}, {3, 7, 6, 2}, {0, 3, 7, 4}, {0, 1, 2, 3}}, 
		{{4, 5, 1, 0}, {7, 6, 2, 3}, {4, 5, 6, 7}, {4, 0, 3, 7}}, 
		{{1, 0, 4, 5}, {2, 3, 7, 6}, {1, 2, 3, 0}, {1, 2, 6, 5}}, 
		{{5, 1, 0, 4}, {6, 2, 3, 7}, {5, 6, 7, 4}, {5, 1, 2, 6}}
	}, {
		{{0, 1, 2, 3}, {0, 4, 5, 1}, {4, 5, 6, 7}, {3, 7, 6, 2}}, 
		{{0, 1, 2, 3}, {0, 4, 5, 1}, {4, 5, 6, 7}, {3, 7, 6, 2}}
	}, {
		{{0, 3, 7, 4}, {4, 5, 6, 7}, {1, 2, 6, 5}, {0, 1, 2, 3}}, 
		{{0, 3, 7, 4}, {4, 5, 6, 7}, {1, 2, 6, 5}, {0, 1, 2, 3}}
	}, {
		{{0, 4, 5, 1}, {3, 7, 6, 2}, {0, 3, 7, 4}, {1, 2, 6, 5}}, 
		{{0, 4, 5, 1}, {3, 7, 6, 2}, {0, 3, 7, 4}, {1, 2, 6, 5}}
	}};
	std::set<I64> elementsToSplit_0 = elementsToSplit;
	while(true){
		std::set<I64> elementsToSplit_A;
		for(const auto currentIndex : elementsToSplit_0){
			const I64 parentIndex = elements[currentIndex].parent;
			if(parentIndex != -1){
				const I64 parentPattern = elements[parentIndex].refinementPattern;
				const std::vector<I64>& parentChildren = elements[parentIndex].children;
				const I64 childrenSize = parentChildren.size();
				I64 isub;
				for(isub = 0; isub < childrenSize; isub ++){
					if(parentChildren[isub] == currentIndex){
						break;
					}
				}
				//
				const I64 plpiSize = parentLine[parentPattern][isub].size();
				for(I64 tj = 0; tj < plpiSize; tj ++){
					const std::array<I64, 2>& plppij = parentLine[parentPattern][isub][tj];
					std::array<I64, 2> tempLine = {
						elements[parentIndex].cornerNodes[plppij[0]], 
						elements[parentIndex].cornerNodes[plppij[1]]
					};
					std::sort(tempLine.begin(), tempLine.end());
					auto iteratorLube = lineUsedByElement.find(tempLine);
					//iteratorLube != lineUsed.end()
					for(const auto elementLubeIndex : (iteratorLube->second)){
						if(elements[elementLubeIndex].children.size() == 0){
							auto iteratorEts = elementsToSplit.find(elementLubeIndex);
							if(iteratorEts == elementsToSplit.end()){
								elementsToSplit_A.emplace(elementLubeIndex);
								elements[elementLubeIndex].refinementPattern = 
									OctreeElement::REFINEMENT_FULL;
							}
						}
					}
				}
				//
				const I64 pfpiSize = parentFace[parentPattern][isub].size();
				for(I64 tj = 0; tj < pfpiSize; tj ++){
					const std::array<I64, 4>& pfppij = parentFace[parentPattern][isub][tj];
					std::array<I64, 4> tempFace = {
						elements[parentIndex].cornerNodes[pfppij[0]], 
						elements[parentIndex].cornerNodes[pfppij[1]], 
						elements[parentIndex].cornerNodes[pfppij[2]], 
						elements[parentIndex].cornerNodes[pfppij[3]]
					};
					std::sort(tempFace.begin(), tempFace.end());
					auto iteratorFube = faceUsedByElement.find(tempFace);
					//iteratorFube != faceUsed.end()
					for(const auto elementFubeIndex : (iteratorFube->second)){
						if(elements[elementFubeIndex].children.size() == 0){
							auto iteratorEts = elementsToSplit.find(elementFubeIndex);
							if(iteratorEts == elementsToSplit.end()){
								elementsToSplit_A.emplace(elementFubeIndex);
								elements[elementFubeIndex].refinementPattern = 
									OctreeElement::REFINEMENT_FULL;
							}
						}
					}
				}
			}
		}
		if(elementsToSplit_A.empty()){
			break;
		}
		else{
			elementsToSplit_0.swap(elementsToSplit_A);
			elementsToSplit.insert(std::make_move_iterator(elementsToSplit_0.begin()), 
                              	   std::make_move_iterator(elementsToSplit_0.end()));
		}
	}
}

void Mesh::OutputMesh(const std::string& directoryPath, const I64& fileIden) const {
	Log("Mesh::OutputMesh");
	std::ofstream tempOfst(
		directoryPath + "/resuNode_" + std::to_string(fileIden) + ".txt", 
		std::ios::out);
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(const auto& iterNoco : node2Coordinate){
		const Coordinate& currentCoordinate = iterNoco.second;
		tempOfst << std::setw(30) << currentCoordinate[0]
			<< std::setw(30) << currentCoordinate[1]
			<< std::setw(30) << currentCoordinate[2] << "\n";
	}
	tempOfst.close();
	tempOfst.open(
		directoryPath + "/resuElem_" + std::to_string(fileIden) + ".txt", 
		std::ios::out);
	for(const auto& currentElement : elements){
		if(currentElement.children.size() == 0){
			const std::array<I64, 8>& currentCornerNodes = currentElement.cornerNodes;
			for(const I64& currentNode : currentCornerNodes){
				tempOfst << std::setw(10) << currentNode;
			}
			tempOfst << "\n";
		}
	}
	tempOfst.close();
}

void Mesh::RigidRotationTranslation(
	const DenseMatrix& rotationMatrix, const Coordinate& translationVector){
	coordinate2Node.clear();
	for(auto& currentNodeCoordinate : node2Coordinate){
		Coordinate& resultCoordinate = currentNodeCoordinate.second;
		Coordinate currentCoordinate = resultCoordinate;
		resultCoordinate = translationVector;
		PEMV(rotationMatrix, currentCoordinate.data, resultCoordinate.data);
		coordinate2Node.emplace(currentNodeCoordinate.second, currentNodeCoordinate.first);
	}
}

void Mesh::OutputDisplacement(
    const std::string& directoryPath, const I64& fileIden, 
    const AlignedVectorRx& resultDisplacement) const {
	//
	Log("    Mesh::OutputDisplacement");
	std::ofstream tempOfst(
		directoryPath + "/resuDisp_" + std::to_string(fileIden) + ".txt", std::ios::out);
	const I64 numNode = resultDisplacement.size() / 3;
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(I64 ti = 0; ti < numNode; ti ++){
		I64 temp_ti = 3 * ti;
		tempOfst << std::setw(30) << resultDisplacement[temp_ti + 0]
			<< std::setw(30) << resultDisplacement[temp_ti + 1]
			<< std::setw(30) << resultDisplacement[temp_ti + 2] << "\n";
	}
	tempOfst.close();
}

I64 Mesh::EffectiveElements(
	const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
	const I64 nestLevel){
	//
    if(effePrefSum.empty()){
		//
		I64 elementSize = elements.size();
		std::vector<I64> isEffective(elementSize + 1, 0);
		effePrefSum.resize(elementSize + 1, 0);
		//
		I64 numbPart = (nestLevel < threadManager.maxNestLevel) ? threadTask.size() : 1;
		std::vector<I64> startIndex(numbPart), endIndex(numbPart);
		EvenlyDistribute(elementSize, numbPart, startIndex, endIndex);
		std::function<void(I64)> taskFunction = 
			[&](I64 tp){
				I64 start_tp = startIndex[tp];
				I64 end_tp = endIndex[tp];
				for(I64 ti = start_tp; ti < end_tp; ++ ti){
					if(elements[ti].children.size() == 0){
						isEffective[ti] = 1;
					}
				}
			};
		switch(nestLevel){
			case 0: case 1: 
				threadManager.RunTask(nestLevel, threadTask, taskFunction);
				break;
			default:
				taskFunction(0);
				break;
		}
		//
		std::exclusive_scan(std::execution::unseq, 
			isEffective.begin(), isEffective.end(), 
			effePrefSum.begin(), 0);
		// effePrefSum[numbEs] = effePrefSum[numbEs - 1] + isEffective[numbEs - 1];
	}
	return effePrefSum.back();
}

Real Mesh::Volume(
	const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
	const I64 nestLevel){
	//
    I64 elementSize = elements.size();
	I64 numbPart = (nestLevel < threadManager.maxNestLevel) ? threadTask.size() : 1;
	std::vector<I64> startIndex(numbPart), endIndex(numbPart);
	EvenlyDistribute(elementSize, numbPart, startIndex, endIndex);
	//
    TrilinearQuadrature<3> trilQuad = GetTrilinearQuadrature<3>();
	std::vector<Real> partVolume(numbPart, 0.0);
	std::function<void(I64)> taskFunction = 
		[&](I64 tp){
			I64 start_tp = startIndex[tp];
			I64 end_tp = endIndex[tp];
			DenseMatrix exyz(8, 3);
			DenseMatrix jacobianJ(3, 3);
			Real jacobianDeterminant;
			Real& partVolume_tp = partVolume[tp];
			for(I64 ti = start_tp; ti < end_tp; ti ++){
				const auto& element = elements[ti];
				if(element.children.size() > 0){
					continue;
				}
				for(I64 tj = 0; tj < 8; tj ++){
					auto iteratorImnc = node2Coordinate.find(element.cornerNodes[tj]);
					for(I64 tk = 0; tk < 3; tk ++){
						exyz(tj,tk) = (iteratorImnc->second)[tk];
					}
				}
				for(I64 tj = 0; tj < trilQuad.numbGaussPoints; tj ++){
					GEMM((trilQuad.shapeDerivatives)[tj], exyz, jacobianJ);
					//a11(a22a33-a23a32)-a12(a21a33-a23a31)+a13(a21a32-a22a31)
					jacobianDeterminant = 
						jacobianJ(0,0) * (jacobianJ(1,1) * jacobianJ(2,2) - jacobianJ(1,2) * jacobianJ(2,1)) -
						jacobianJ(0,1) * (jacobianJ(1,0) * jacobianJ(2,2) - jacobianJ(1,2) * jacobianJ(2,0)) +
						jacobianJ(0,2) * (jacobianJ(1,0) * jacobianJ(2,1) - jacobianJ(1,1) * jacobianJ(2,0));
					partVolume_tp += trilQuad.weights[tj] * jacobianDeterminant;
				}
			}
		};
    switch(nestLevel){
        case 0: case 1: 
            threadManager.RunTask(nestLevel, threadTask, taskFunction);
            break;
        default:
            taskFunction(0);
            break;
    }
	//
	return std::reduce(std::execution::par, 
		partVolume.begin(), partVolume.end(), 
		0.0, 
		std::plus<>()
	);
}

} // namespace Ddpca