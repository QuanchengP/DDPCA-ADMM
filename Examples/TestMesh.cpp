#include "../Mesh/Mesh.hpp"

#include <filesystem>
#include <string>

class MeshTest {
private:
    std::string directoryPath;
    
public:
    explicit MeshTest(std::string path) : directoryPath(std::move(path)) {}
    
    void runTests() {
        homogeneous();
        inhomogeneous();
    }
    
private:
    void homogeneous() {
        Ddpca::Mesh currentMesh;
        //
        const Ddpca::Real length = 1.0;
        const Ddpca::I64 numbDivision = 2;
        //nodes
        Ddpca::I64 tempNode[numbDivision + 1][numbDivision + 1][numbDivision + 1];
        for(Ddpca::I64 ti = 0; ti <= numbDivision; ++ ti){
            for(Ddpca::I64 tj = 0; tj <= numbDivision; ++ tj){
                for(Ddpca::I64 tk = 0; tk <= numbDivision; ++ tk){
                    Ddpca::Real tempX = - length / 2.0 + length / numbDivision * (Ddpca::Real)ti;
                    Ddpca::Real tempY = - length / 2.0 + length / numbDivision * (Ddpca::Real)tj;
                    Ddpca::Real tempZ = 0.0 + length / numbDivision * (Ddpca::Real)tk;
                    Ddpca::Coordinate currentCoordinate(tempX, tempY, tempZ);
                    tempNode[ti][tj][tk] = currentMesh.TryAddNode(currentCoordinate);
                }
            }
        }
        //elements
        for(Ddpca::I64 ti = 0; ti < numbDivision; ++ ti){
            for(Ddpca::I64 tj = 0; tj < numbDivision; ++ tj){
                for(Ddpca::I64 tk = 0; tk < numbDivision; ++ tk){
                    Ddpca::OctreeElement tempElement;
                    tempElement.parent = -1;
                    tempElement.cornerNodes = {
                        tempNode[tk][tj][ti], tempNode[tk + 1][tj][ti], 
                        tempNode[tk + 1][tj + 1][ti], tempNode[tk][tj + 1][ti], 
                        tempNode[tk][tj][ti + 1], tempNode[tk + 1][tj][ti + 1], 
                        tempNode[tk + 1][tj + 1][ti + 1], tempNode[tk][tj + 1][ti + 1]
                    };
                    tempElement.level = 0;
                    tempElement.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
                    tempElement.children.resize(0);
                    /*Ddpca::I64 elemNumb = */currentMesh.AddElement(tempElement);
                }
            }
        }
        //refinement level 0
        std::set<Ddpca::I64> elementsToSplit;
	    std::map<Ddpca::I64, std::set<Ddpca::I64>> subElements;
	    std::map<std::vector<Ddpca::I64>, Ddpca::Coordinate> curvilinearInterpolation;
        elementsToSplit.clear();//
        elementsToSplit.emplace(6);
        currentMesh.elements[6].refinementPattern = 
            Ddpca::OctreeElement::REFINEMENT_FULL;
        elementsToSplit.emplace(7);
        currentMesh.elements[7].refinementPattern = 
            Ddpca::OctreeElement::REFINEMENT_FULL;
        subElements.clear();//
        subElements[6] = {5, 7};
        subElements[7] = {4, 6};
        curvilinearInterpolation.clear();//
        currentMesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        //refinement level 1
        for(auto const &iterSpel : elementsToSplit){
        	currentMesh.elements[iterSpel].refinementPattern = 
                Ddpca::OctreeElement::REFINEMENT_FULL;
        }
        subElements.clear();//
        subElements[13] = {5, 7};
        subElements[15] = {5, 7};
        subElements[20] = {4, 6};
        subElements[22] = {4, 6};
        curvilinearInterpolation.clear();//
        currentMesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        //refinement level 2
        for(auto const &iterSpel : elementsToSplit){
        	currentMesh.elements[iterSpel].refinementPattern = 
                Ddpca::OctreeElement::REFINEMENT_FULL;
        }
        subElements.clear();//
        curvilinearInterpolation.clear();//
        currentMesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        //
        currentMesh.OutputMesh(directoryPath, 0);
    }
    
    void inhomogeneous() {
        Ddpca::Mesh currentMesh;
        //
        Ddpca::Real length = 1.0;
        Ddpca::I64 numbDivision = 1;
        //nodes
        Ddpca::I64 tempNode[numbDivision + 1][numbDivision + 1][numbDivision + 1];
        for(Ddpca::I64 ti = 0; ti <= numbDivision; ++ ti){
        	for(Ddpca::I64 tj = 0; tj <= numbDivision; ++ tj){
        		for(Ddpca::I64 tk = 0; tk <= numbDivision; ++ tk){
        			Ddpca::Real tempX = - length / 2.0 + length / numbDivision * (Ddpca::Real)ti;
        			Ddpca::Real tempY = - length / 2.0 + length / numbDivision * (Ddpca::Real)tj;
        			Ddpca::Real tempZ = 0.0 + length / numbDivision * (Ddpca::Real)tk;
        			Ddpca::Coordinate tempCoor(tempX, tempY, tempZ);
        			tempNode[ti][tj][tk] = currentMesh.TryAddNode(tempCoor);
        		}
        	}
        }
        //elements
        for(Ddpca::I64 ti = 0; ti < numbDivision; ++ ti){
        	for(Ddpca::I64 tj = 0; tj < numbDivision; ++ tj){
        		for(Ddpca::I64 tk = 0; tk < numbDivision; ++ tk){
        			Ddpca::OctreeElement tempElement;
        			tempElement.parent = -1;
                    tempElement.cornerNodes = {
                        tempNode[tk][tj][ti], tempNode[tk + 1][tj][ti], 
                        tempNode[tk + 1][tj + 1][ti], tempNode[tk][tj + 1][ti], 
                        tempNode[tk][tj][ti + 1], tempNode[tk + 1][tj][ti + 1], 
                        tempNode[tk + 1][tj + 1][ti + 1], tempNode[tk][tj + 1][ti + 1]
                    };
                    tempElement.level = 0;
                    tempElement.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
                    tempElement.children.resize(0);
                    /*Ddpca::I64 elemNumb = */currentMesh.AddElement(tempElement);
        		}
        	}
        }
        //refinement level 0
        std::set<Ddpca::I64> elementsToSplit;
	    std::map<Ddpca::I64, std::set<Ddpca::I64>> subElements;
	    std::map<std::vector<Ddpca::I64>, Ddpca::Coordinate> curvilinearInterpolation;
        elementsToSplit.clear();//
        elementsToSplit.emplace(0);
        currentMesh.elements[0].refinementPattern = 
            Ddpca::OctreeElement::REFINEMENT_ZETA;
        subElements.clear();//
        subElements[0] = {0, 1};
        curvilinearInterpolation.clear();//
        currentMesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        //refinement level 1
        elementsToSplit.clear();//
        std::vector<Ddpca::OctreeElement>& elements = currentMesh.elements;
        Ddpca::I64 meelSize = currentMesh.elements.size();
        for(Ddpca::I64 ti = 0; ti < meelSize; ++ ti){
        	if(elements[ti].children.size() > 0){
        		continue;
        	}
        	elementsToSplit.emplace(ti);
        	elements[ti].refinementPattern = 
                Ddpca::OctreeElement::REFINEMENT_ZETA_XI;
        }
        subElements.clear();//
        curvilinearInterpolation.clear();//0
        currentMesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        //refinement level 2
        elementsToSplit.clear();//
        meelSize = currentMesh.elements.size();
        for(Ddpca::I64 ti = 0; ti < meelSize; ++ ti){
        	if(elements[ti].children.size() > 0){
        		continue;
        	}
        	elementsToSplit.emplace(ti);
        	elements[ti].refinementPattern = 
                Ddpca::OctreeElement::REFINEMENT_ETA_ZETA;
        }
        subElements.clear();//
        curvilinearInterpolation.clear();//0
        currentMesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        //refinement level 3
        elementsToSplit.clear();//
        meelSize = currentMesh.elements.size();
        for(Ddpca::I64 ti = 0; ti < meelSize; ++ ti){
        	if(elements[ti].children.size() > 0){
        		continue;
        	}
        	elementsToSplit.emplace(ti);
        	elements[ti].refinementPattern = 
                Ddpca::OctreeElement::REFINEMENT_FULL;
        }
        subElements.clear();//
        curvilinearInterpolation.clear();//0
        currentMesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        //
        currentMesh.OutputMesh(directoryPath, 1);
    }
};

int main(/*int argc, char **argv*/){
    std::string directoryPath = "./TestMesh_";
    std::filesystem::create_directory(directoryPath);
    MeshTest test(directoryPath);
    test.runTests();
    return 0;
}