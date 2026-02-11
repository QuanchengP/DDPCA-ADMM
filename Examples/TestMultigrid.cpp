#include "../General/General.hpp"
#include "../Mesh/GeometricMultigrid.hpp"
#include "../Solver/ConjugateGradient.hpp"

#include <filesystem>
// #include <string>

class MultigridTest {
private:
    std::string directoryPath;
    
public:
    explicit MultigridTest(std::string path) : directoryPath(std::move(path)) {}
    
    void runTest() {
        Ddpca::threadManager.ThreadDistribute(
            80, 1 + 2, Ddpca::threadManager.levelS2S, Ddpca::threadManager.levelS2M);

        Ddpca::Mesh dataMesh;
        Ddpca::BoundaryCondition dataBoundary;
        Ddpca::SparseMatrix dataStiffness;
        Ddpca::GeometricMultigrid geomMult;
        Ddpca::ConjugateGradient conjGrad;
        conjGrad.preconditioner = Ddpca::PreconditionerMultigrid;
        conjGrad.pcGM.relaxation = Ddpca::RelaxationChebyshev;

        //
        GenerateMesh(dataMesh);
        ConstraintLoad(dataBoundary);

        //
        geomMult.dataMesh = &dataMesh;
        geomMult.dataBoundary = &dataBoundary;
        geomMult.dataStiffness = &dataStiffness;
        geomMult.dataPcg = &(conjGrad.pcGM);
        geomMult.Transfer(Ddpca::threadManager.levelS2S, 0);
        FemStiffnessMatrix(dataMesh, 210.0E9, 0.3, dataStiffness, Ddpca::threadManager.one2oneS2S, 0);
        geomMult.Establish(Ddpca::threadManager.levelS2S, 0);
        geomMult.Output(directoryPath, 0, Ddpca::threadManager.one2oneS2S, 0);

        //
        Ddpca::I64 maxLevel = conjGrad.pcGM.realProlong.size();
        const Ddpca::SparseMatrix& solverStiffness = (conjGrad.pcGM.hierarchyStiffness)[maxLevel];
        conjGrad.Establish(solverStiffness, Ddpca::threadManager.one2oneS2S, 0);
        Ddpca::AlignedVectorRx solvDisp(solverStiffness.M, 0.0);
        conjGrad.Solve(solverStiffness, dataBoundary.loadVect, solvDisp, Ddpca::threadManager.one2oneS2S, 0);

        //
        Ddpca::I64 numbDmnc = 3 * dataMesh.node2Coordinate.size();
        Ddpca::AlignedVectorRx resuDisp(numbDmnc, 0.0);
        geomMult.DispPost(solvDisp, resuDisp, Ddpca::threadManager.one2oneS2S, 0);
        dataMesh.OutputDisplacement(directoryPath, 0, resuDisp);
    }
    
private:

    void GenerateMesh(Ddpca::Mesh &dataMesh){
        //nodes
        dataMesh.TryAddNode(Ddpca::Coordinate(0.0, 0.0, 0.0));
        dataMesh.TryAddNode(Ddpca::Coordinate(1.0, 0.0, 0.0));
        dataMesh.TryAddNode(Ddpca::Coordinate(1.0, 1.0, 0.0));
        dataMesh.TryAddNode(Ddpca::Coordinate(0.0, 1.0, 0.0));
        dataMesh.TryAddNode(Ddpca::Coordinate(0.0, 0.0, 1.0));
        dataMesh.TryAddNode(Ddpca::Coordinate(1.0, 0.0, 1.0));
        dataMesh.TryAddNode(Ddpca::Coordinate(1.0, 1.0, 1.0));
        dataMesh.TryAddNode(Ddpca::Coordinate(0.0, 1.0, 1.0));
        //elements
        Ddpca::OctreeElement tempElement;
        tempElement.parent = -1;
        tempElement.cornerNodes = {0, 1, 2, 3, 4, 5, 6, 7};
        tempElement.level = 0;
        tempElement.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
        tempElement.children.clear();
        dataMesh.AddElement(tempElement);
        //global refinement
        std::set<Ddpca::I64> elementsToSplit;
	    std::map<Ddpca::I64, std::set<Ddpca::I64>> subElements;
	    std::map<std::vector<Ddpca::I64>, Ddpca::Coordinate> curvilinearInterpolation;
        Ddpca::I64 ti = 0;
        elementsToSplit.emplace(ti);
        dataMesh.elements[ti].refinementPattern = 
            Ddpca::OctreeElement::REFINEMENT_FULL;
        subElements.clear();
        curvilinearInterpolation.clear();
        dataMesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        //
        dataMesh.OutputMesh(directoryPath, 0);
    }

    void ConstraintLoad(Ddpca::BoundaryCondition &dataBoundary){
        //displacement constraint
        std::vector<Ddpca::I64> inpuNode = {0, 3, 7, 4, 11, 15, 19, 12, 20};
        for(const auto& node : inpuNode){
            dataBoundary.constrainedDof.emplace(3 * node + 0, 0.0);
            dataBoundary.constrainedDof.emplace(3 * node + 1, 0.0);
            dataBoundary.constrainedDof.emplace(3 * node + 2, 0.0);
        }
        //load
        inpuNode = {1, 2, 6, 5, 9, 14, 17, 13, 21};
        for(Ddpca::I64 ti = 0; ti <= 3; ++ ti){
            dataBoundary.LoadAccumulate(3 * inpuNode[ti] + 0, 100.0);
        }
        for(Ddpca::I64 ti = 4; ti <= 7; ++ ti){
            dataBoundary.LoadAccumulate(3 * inpuNode[ti] + 0, 200.0);
        }
        for(Ddpca::I64 ti = 8; ti <= 8; ++ ti){
            dataBoundary.LoadAccumulate(3 * inpuNode[ti] + 0, 400.0);
        }
    }
};

int main(int argc, char **argv){
    //
    Ddpca::Initialize(argc, argv);
    //
    std::string directoryPath = "./TestMultigrid_";
    std::filesystem::create_directory(directoryPath);
    MultigridTest test(directoryPath);
    test.runTest();
    //
    Ddpca::Finalize();
	return 0;
}
