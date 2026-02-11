#include "../Mesh/FemStiffnessMatrix.hpp"

int main(/*int argc, char **argv*/){
    /****************************************************************************************************/
    std::cout << "!********************************************************************************!\n";
    Ddpca::Mesh tempMesh;
    //
    std::array<Ddpca::I64, 8> tempNode;
    tempNode[0] = tempMesh.TryAddNode(Ddpca::Coordinate(0.0, 0.0, 0.0));
    tempNode[1] = tempMesh.TryAddNode(Ddpca::Coordinate(1.0, 0.0, 0.0));
    tempNode[2] = tempMesh.TryAddNode(Ddpca::Coordinate(1.0, 1.0, 0.0));
    tempNode[3] = tempMesh.TryAddNode(Ddpca::Coordinate(0.0, 1.0, 0.0));
    tempNode[4] = tempMesh.TryAddNode(Ddpca::Coordinate(0.0, 0.0, 1.0));
    tempNode[5] = tempMesh.TryAddNode(Ddpca::Coordinate(1.0, 0.0, 1.0));
    tempNode[6] = tempMesh.TryAddNode(Ddpca::Coordinate(1.0, 1.0, 1.0));
    tempNode[7] = tempMesh.TryAddNode(Ddpca::Coordinate(0.0, 1.0, 1.0));
    //
    Ddpca::OctreeElement tempElement;
    tempElement.parent = -1;
    tempElement.cornerNodes[0] = tempNode[0];
    tempElement.cornerNodes[1] = tempNode[1];
    tempElement.cornerNodes[2] = tempNode[2];
    tempElement.cornerNodes[3] = tempNode[3];
    tempElement.cornerNodes[4] = tempNode[4];
    tempElement.cornerNodes[5] = tempNode[5];
    tempElement.cornerNodes[6] = tempNode[6];
    tempElement.cornerNodes[7] = tempNode[7];
    tempElement.level = 0;
    tempElement.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
    tempElement.children.resize(0);
    tempMesh.AddElement(tempElement);
    //
    Ddpca::SparseMatrix tempStiff;
    Ddpca::FemStiffnessMatrix(tempMesh, 210.0E9, 0.3, tempStiff, Ddpca::threadManager.one2oneS2S, 0);
    tempStiff.Output();
    /****************************************************************************************************/
    std::cout << "!********************************************************************************!\n";
    Ddpca::Mesh tempMesh1;
    //
    std::array<Ddpca::I64, 8> tempNode1;
    tempNode1[0] = tempMesh1.TryAddNode(Ddpca::Coordinate(0.0, 0.0, 0.0));
    tempNode1[1] = tempMesh1.TryAddNode(Ddpca::Coordinate(1.0, 0.0, 0.0));
    tempNode1[2] = tempMesh1.TryAddNode(Ddpca::Coordinate(1.0, 0.8, 0.0));
    tempNode1[3] = tempMesh1.TryAddNode(Ddpca::Coordinate(0.0, 1.0, 0.0));
    tempNode1[4] = tempMesh1.TryAddNode(Ddpca::Coordinate(0.0, 0.0, 1.0));
    tempNode1[5] = tempMesh1.TryAddNode(Ddpca::Coordinate(1.0, 0.0, 1.0));
    tempNode1[6] = tempMesh1.TryAddNode(Ddpca::Coordinate(0.8, 0.8, 0.8));
    tempNode1[7] = tempMesh1.TryAddNode(Ddpca::Coordinate(0.0, 1.0, 1.0));
    //
    Ddpca::OctreeElement tempElement1;
    tempElement1.parent = -1;
    tempElement1.cornerNodes[0] = tempNode1[0];
    tempElement1.cornerNodes[1] = tempNode1[1];
    tempElement1.cornerNodes[2] = tempNode1[2];
    tempElement1.cornerNodes[3] = tempNode1[3];
    tempElement1.cornerNodes[4] = tempNode1[4];
    tempElement1.cornerNodes[5] = tempNode1[5];
    tempElement1.cornerNodes[6] = tempNode1[6];
    tempElement1.cornerNodes[7] = tempNode1[7];
    tempElement1.level = 0;
    tempElement1.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
    tempElement1.children.resize(0);
    tempMesh1.AddElement(tempElement1);
    //
    Ddpca::SparseMatrix tempStiff1;
    Ddpca::FemStiffnessMatrix(tempMesh1, 210.0E9, 0.3, tempStiff1, Ddpca::threadManager.one2oneS2S, 0);
    tempStiff1.Output();
    return 0;
}
