#include "../Mesh/BilinearQuadrature.hpp"
#include "../Decomposition/SingleDomain.hpp"

#include <filesystem>

//example beam
class Beam : public Ddpca::SingleDomain {
    
public:
	/*********************************************************************************************/
    std::string directoryPath;
	std::array<Ddpca::Real,3> leng = {1.0, 0.12, 0.06}; // geometric length
	Ddpca::Real lengFact = 1.0 / 3.0;//length reduction factor from fixed end to free end
	Ddpca::Real angl = 45.0 * Ddpca::PI / 180.0;;//geometric angle
    //number of elements along different directions
	std::array<Ddpca::I64,3> diviNumb = {64, 4, 2};//must be even
	Ddpca::I64 globLeve = 4;//multigrid level, total levels of global mesh refinement
	Ddpca::Real loadInte = - 8000.0;//centerline load intensity

public:
	/*********************************************************************************************/

    void AdjustCoordinate(){
        mesh.coordinate2Node.clear();
        for(auto& currentNodeCoordinate : mesh.node2Coordinate){
            Ddpca::Coordinate& resultCoordinate = currentNodeCoordinate.second;
            Ddpca::Coordinate currentCoordinate = resultCoordinate;
            Ddpca::Real angl_temp = angl * currentCoordinate[0] / leng[0];
            //col-major
            Ddpca::DenseMatrix rotationMatrix(3, 3, {
                1.0, 0.0, 0.0, 
                0.0, std::cos(angl_temp), std::sin(angl_temp), 
                0.0, -std::sin(angl_temp), std::cos(angl_temp)});
            Ddpca::GEMV(rotationMatrix, currentCoordinate.data, resultCoordinate.data);
            mesh.coordinate2Node.emplace(currentNodeCoordinate.second, currentNodeCoordinate.first);
        }
    }

    void GenerateMesh(){
        //nodes
        std::vector<std::vector<std::vector<Ddpca::I64>>> tempNode(
            diviNumb[0] + 1, 
            std::vector<std::vector<Ddpca::I64>>(
                diviNumb[1] + 1, 
                std::vector<Ddpca::I64>(diviNumb[2] + 1, 0)
            )
        );
        Ddpca::I64 numbPart = Ddpca::threadManager.one2oneS2S.size();
        std::vector<Ddpca::I64> startIndex(numbPart), endIndex(numbPart);
        Ddpca::EvenlyDistribute(diviNumb[0] + 1, numbPart, startIndex, endIndex);
        std::mutex tempMutex;
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.one2oneS2S, [&](Ddpca::I64 tp){
            Ddpca::I64 start_tp = startIndex[tp];
            Ddpca::I64 end_tp = endIndex[tp];
            for(Ddpca::I64 ti = start_tp; ti < end_tp; ++ ti){
                Ddpca::Real heig_ti = leng[1] * (1.0- (Ddpca::Real)ti / diviNumb[0] * lengFact);
                Ddpca::Real widt_ti = leng[2] * (1.0- (Ddpca::Real)ti / diviNumb[0] * lengFact);
                for(Ddpca::I64 tj = 0; tj <= diviNumb[1]; tj ++){
                    Ddpca::Real yCoo = - heig_ti / 2.0 + heig_ti / diviNumb[1] * tj;
                    for(Ddpca::I64 tk = 0; tk <= diviNumb[2]; tk ++){
                        Ddpca::Real zCoo = - widt_ti / 2.0 + widt_ti / diviNumb[2] * tk;
                        Ddpca::Coordinate tempCoor(leng[0] / diviNumb[0] * ti, yCoo, zCoo);
                        {
                            std::lock_guard<std::mutex> lock(tempMutex);
                            tempNode[ti][tj][tk] = mesh.TryAddNode(tempCoor);
                        }
                    }
                }
            }
        });
        //elements
        Ddpca::EvenlyDistribute(diviNumb[0], numbPart, startIndex, endIndex);
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.one2oneS2S, [&](Ddpca::I64 tp){
            Ddpca::I64 start_tp = startIndex[tp];
            Ddpca::I64 end_tp = endIndex[tp];
            for(Ddpca::I64 ti = start_tp; ti < end_tp; ++ ti){
                for(Ddpca::I64 tj = 0; tj < diviNumb[1]; ++ tj){
                    for(Ddpca::I64 tk = 0; tk < diviNumb[2]; ++ tk){
                        Ddpca::OctreeElement tempElem;
                        tempElem.parent = -1;
                        tempElem.cornerNodes[0] = tempNode[ti][tj][tk];
                        tempElem.cornerNodes[1] = tempNode[ti][tj + 1][tk];
                        tempElem.cornerNodes[2] = tempNode[ti][tj + 1][tk + 1];
                        tempElem.cornerNodes[3] = tempNode[ti][tj][tk + 1];
                        tempElem.cornerNodes[4] = tempNode[ti + 1][tj][tk];
                        tempElem.cornerNodes[5] = tempNode[ti + 1][tj + 1][tk];
                        tempElem.cornerNodes[6] = tempNode[ti + 1][tj + 1][tk + 1];
                        tempElem.cornerNodes[7] = tempNode[ti + 1][tj][tk + 1];
                        tempElem.level = 0;
                        tempElem.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
                        tempElem.children.clear();
                        {
                            std::lock_guard<std::mutex> lock(tempMutex);
                            mesh.AddElement(tempElem);
                        }
                    }
                }
            }
        });
        //global refinement
        std::set<Ddpca::I64> elementsToSplit;
	    std::map<Ddpca::I64, std::set<Ddpca::I64>> subElements;
	    std::map<std::vector<Ddpca::I64>, Ddpca::Coordinate> curvilinearInterpolation;
        for(long tr = 0; tr < globLeve; tr ++){
            elementsToSplit.clear();
            Ddpca::I64 numbDael = mesh.elements.size();
            Ddpca::EvenlyDistribute(numbDael, numbPart, startIndex, endIndex);
            Ddpca::threadManager.RunTask(0, Ddpca::threadManager.one2oneS2S, [&](Ddpca::I64 tp){
                Ddpca::I64 start_tp = startIndex[tp];
                Ddpca::I64 end_tp = endIndex[tp];
                for(Ddpca::I64 ti = start_tp; ti < end_tp; ++ ti){
                    if(mesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    mesh.elements[ti].refinementPattern = 
                        Ddpca::OctreeElement::REFINEMENT_FULL;
                    {
                        std::lock_guard<std::mutex> lock(tempMutex);
                        elementsToSplit.emplace(ti);
                    }
                }
            });
            mesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        }
        AdjustCoordinate();
        //
        mesh.OutputMesh(directoryPath, 0);
    }

    void ConstraintLoad(){
        //displacement constraint
        for(const auto& iterNoco: mesh.node2Coordinate){
            if((iterNoco.second)[0] <= 1.0E-10){
                boundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                boundary.constrainedDof.emplace(3 * iterNoco.first + 1, 0.0);
                boundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
            }
        }
        //load must at next
        Ddpca::I64 ti_max = mesh.elements.size();
        Ddpca::I64 numbPart = Ddpca::threadManager.one2oneS2S.size();
        std::vector<Ddpca::I64> startIndex(numbPart), endIndex(numbPart);
        Ddpca::EvenlyDistribute(ti_max, numbPart, startIndex, endIndex);
        std::mutex tempMutex;
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.one2oneS2S, [&](Ddpca::I64 tp){
            Ddpca::I64 start_tp = startIndex[tp];
            Ddpca::I64 end_tp = endIndex[tp];
            const Ddpca::I64 heliSize = Ddpca::hexaLine.size();
            for(Ddpca::I64 ti = start_tp; ti < end_tp; ++ ti){
                const Ddpca::OctreeElement& tempElem = mesh.elements[ti];
                if(tempElem.children.size() > 0){
                    continue;
                }
                for(long tj = 0; tj < heliSize; tj ++){
                    const Ddpca::I64 tk_size = Ddpca::hexaLine[tj].size();
                    std::vector<Ddpca::I64> inpuNode(tk_size);
                    std::vector<Ddpca::Coordinate> inpuCoor(tk_size);
                    bool tempFlag = true;
                    for(long tk = 0; tk < tk_size; tk ++){
                        inpuNode[tk] = tempElem.cornerNodes[Ddpca::hexaLine[tj][tk]];
                        auto iterNoco = mesh.node2Coordinate.find(inpuNode[tk]);
                        inpuCoor[tk] = iterNoco->second;
                        if(std::abs((iterNoco->second)[1]) > 1.0E-10 
                            || std::abs((iterNoco->second)[2]) > 1.0E-10){
                            tempFlag = false;
                            break;
                        }
                    }
                    if(tempFlag == false){
                        continue;
                    }
                    // / 4.0: must (four elements share one edge)
                    Ddpca::Real forcValu = loadInte * std::abs(inpuCoor[0][0] - inpuCoor[1][0]) / 2.0 / 4.0;
                    for(long tk = 0; tk < tk_size; tk ++){
                        if(inpuCoor[tk][0] > 1.0E-10){
                            {
                                std::lock_guard<std::mutex> lock(tempMutex);
                                boundary.LoadAccumulate(3 * inpuNode[tk] + 2, forcValu);
                            }
                        }
                    }
                }
            }
        });
    }

    void Test(){
        //
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.numbThreads / (globLeve + 2), globLeve + 2, 
            Ddpca::threadManager.levelS2S, Ddpca::threadManager.levelS2M);
        //
        conjGrad.preconditioner = Ddpca::PreconditionerMultigrid;
        conjGrad.pcGM.relaxation = Ddpca::RelaxationChebyshev;

        //
        GenerateMesh();
        ConstraintLoad();

        //
        geomMult.dataMesh = &mesh;
        geomMult.dataBoundary = &boundary;
        geomMult.dataStiffness = &stiffness;
        geomMult.dataPcg = &(conjGrad.pcGM);
        geomMult.Transfer(Ddpca::threadManager.levelS2S, 0);
        //must after Transfer
        FemStiffnessMatrix(mesh, 210.0E9, 0.3, stiffness, Ddpca::threadManager.one2oneS2S, 0);
        geomMult.Establish(Ddpca::threadManager.levelS2S, 0);

        //
        Ddpca::I64 maxLevel = conjGrad.pcGM.realProlong.size();
        const Ddpca::SparseMatrix& solverStiffness = (conjGrad.pcGM.hierarchyStiffness)[maxLevel];
        conjGrad.Establish(solverStiffness, Ddpca::threadManager.one2oneS2S, 0);
        Ddpca::AlignedVectorRx solvDisp(solverStiffness.M, 0.0);
        conjGrad.Solve(solverStiffness, boundary.loadVect, solvDisp, Ddpca::threadManager.one2oneS2S, 0);

        //
        Ddpca::I64 numbDmnc = 3 * mesh.node2Coordinate.size();
        Ddpca::AlignedVectorRx resuDisp(numbDmnc, 0.0);
        geomMult.DispPost(solvDisp, resuDisp, Ddpca::threadManager.one2oneS2S, 0);
        mesh.OutputDisplacement(directoryPath, 0, resuDisp);
    }

};

int main(int argc, char **argv){
    //
    Ddpca::Initialize(argc, argv);
	//
    std::string directoryPath = "./TestNDDBeam_";
    std::filesystem::create_directory(directoryPath);
    //
    Beam beam;
    beam.directoryPath = directoryPath;
    beam.Test();
    //
    Ddpca::Finalize();
	return 1;
}
