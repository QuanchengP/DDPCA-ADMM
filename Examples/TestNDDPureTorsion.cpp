#include "../Mesh/BilinearQuadrature.hpp"
#include "../Decomposition/SingleDomain.hpp"

#include <filesystem>

//example pure torsional cylinder
class PureTorsion : public Ddpca::SingleDomain {
public:
	/*********************************************************************************************/
    std::string directoryPath;
	Ddpca::Real axiaLeng = 0.1;
	Ddpca::Real inneRadi = 0.015;
	Ddpca::Real outeRadi = 0.025;//I_p = PI / 32.0 * (D^4 - d^4) = 0.005340707511103
	std::array<Ddpca::I64,3> diviNumb = {2, 16, 16};//8, 32, 64
	Ddpca::I64 globInho = 0;//0
	Ddpca::I64 globHomo = 4;//0
	Ddpca::Real torqLoad = 20.0;//N*m, u = T*l/G/I_p * R = 1.159111630361142e-06

public:
	/*********************************************************************************************/

    Ddpca::Coordinate CoordinateAverage(const std::vector<Ddpca::Coordinate> &inpuCoor){
        Ddpca::Coordinate cyliCoor(0.0, 0.0, 0.0);
        Ddpca::I64 tempFlag_0 = 0;
        Ddpca::I64 tempFlag_1 = 0;
        Ddpca::I64 numbInco = inpuCoor.size();
        for(Ddpca::I64 tk = 0; tk < numbInco; ++ tk){
            cyliCoor[0] += std::sqrt(
                inpuCoor[tk][0] * inpuCoor[tk][0] 
                + inpuCoor[tk][1] * inpuCoor[tk][1]
            );
            Ddpca::Real tempAngl = std::atan2(inpuCoor[tk][1], inpuCoor[tk][0]);
            cyliCoor[1] += tempAngl;
            cyliCoor[2] += inpuCoor[tk][2];
            if(tempAngl > Ddpca::PI / 2.0){
                ++ tempFlag_0;
            }
            if(tempAngl < - Ddpca::PI / 2.0){
                ++ tempFlag_1;
            }
        }
        if(tempFlag_0 > 0 && tempFlag_1 > 0){
            cyliCoor[1] += tempFlag_1 * (Ddpca::PI * 2.0);
        }
        Ddpca::Real averRadi = cyliCoor[0] / inpuCoor.size();
        Ddpca::Real averAngl = cyliCoor[1] / inpuCoor.size();
        Ddpca::Coordinate outpCoor(
            averRadi * std::cos(averAngl), 
            averRadi * std::sin(averAngl), 
            cyliCoor[2] / inpuCoor.size());
        return outpCoor;
    }

    void CurvilinearInterpolation(
        Ddpca::I64 ei, std::map<std::vector<Ddpca::I64>,Ddpca::Coordinate> &planSurf){
        //
        std::vector<Ddpca::I64> inpuNode(8);
        std::vector<Ddpca::Coordinate> inpuCoor(8);
        for(Ddpca::I64 tk = 0; tk < 8; ++ tk){
            inpuNode[tk] = mesh.elements[ei].cornerNodes[tk];
            inpuCoor[tk] = mesh.node2Coordinate[inpuNode[tk]];
        }
        Ddpca::Coordinate outpCoor = CoordinateAverage(inpuCoor);
        std::sort(inpuNode.begin(), inpuNode.end());
        planSurf.emplace(inpuNode, outpCoor);
        //
        const Ddpca::I64 heliSize = Ddpca::hexaLine.size();
        for(Ddpca::I64 tj = 0; tj < heliSize; ++ tj){
            const Ddpca::I64 size_tj = Ddpca::hexaLine[tj].size();
            std::vector<Ddpca::I64> inpuNode(size_tj);
            std::vector<Ddpca::Coordinate> inpuCoor(size_tj);
            for(Ddpca::I64 tk = 0; tk < size_tj; ++ tk){
                inpuNode[tk] = mesh.elements[ei].cornerNodes[Ddpca::hexaLine[tj][tk]];
                inpuCoor[tk] = mesh.node2Coordinate[inpuNode[tk]];
            }
            Ddpca::Coordinate outpCoor = CoordinateAverage(inpuCoor);
            std::sort(inpuNode.begin(), inpuNode.end());
            planSurf.emplace(inpuNode, outpCoor);
        }
        //
        const Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
        for(Ddpca::I64 tj = 0; tj < hefaSize; ++ tj){
            const Ddpca::I64 size_tj = Ddpca::hexaFace[tj].size();
            std::vector<Ddpca::I64> inpuNode(size_tj);
            std::vector<Ddpca::Coordinate> inpuCoor(size_tj);
            for(Ddpca::I64 tk = 0; tk < size_tj; ++ tk){
                inpuNode[tk] = mesh.elements[ei].cornerNodes[Ddpca::hexaFace[tj][tk]];
                inpuCoor[tk] = mesh.node2Coordinate[inpuNode[tk]];
            }
            Ddpca::Coordinate outpCoor = CoordinateAverage(inpuCoor);
            std::sort(inpuNode.begin(), inpuNode.end());
            planSurf.emplace(inpuNode, outpCoor);
        }
    }

    void GenerateMesh(){
        //nodes
        std::vector<std::vector<std::vector<Ddpca::I64>>> blocNode(
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
                Ddpca::Real radi_ti = inneRadi + (outeRadi - inneRadi) / (Ddpca::Real)(diviNumb[0]) * ti;
                for(Ddpca::I64 tj = 0; tj <= diviNumb[1]; ++ tj){
                    Ddpca::Real angl_tj = 0.0 + (2.0 * Ddpca::PI - 0.0) / (Ddpca::Real)(diviNumb[1]) * tj;
                    for(Ddpca::I64 tk =  0; tk <= diviNumb[2]; ++ tk){
                        Ddpca::Coordinate tempCoor(
                            radi_ti * std::cos(angl_tj), radi_ti * std::sin(angl_tj), 
                            axiaLeng / diviNumb[2] * (Ddpca::Real)tk
                        );
                        {
                            std::lock_guard<std::mutex> lock(tempMutex);
                            blocNode[ti][tj][tk] = mesh.TryAddNode(tempCoor);
                        }
                    }
                }
            }
        });
        //elements
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
                        tempElem.cornerNodes[0] = blocNode[ti][tj][tk];
                        tempElem.cornerNodes[1] = blocNode[ti + 1][tj][tk];
                        tempElem.cornerNodes[2] = blocNode[ti + 1][tj + 1][tk];
                        tempElem.cornerNodes[3] = blocNode[ti][tj + 1][tk];
                        tempElem.cornerNodes[4] = blocNode[ti][tj][tk + 1];
                        tempElem.cornerNodes[5] = blocNode[ti + 1][tj][tk + 1];
                        tempElem.cornerNodes[6] = blocNode[ti + 1][tj + 1][tk + 1];
                        tempElem.cornerNodes[7] = blocNode[ti][tj + 1][tk + 1];
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
        for(Ddpca::I64 tr = 0; tr < globInho + globHomo; ++ tr){
            elementsToSplit.clear();
            Ddpca::I64 numbDael = mesh.elements.size();
            for(Ddpca::I64 ti = 0; ti < numbDael; ++ ti){
                if(mesh.elements[ti].children.size() > 0){
                    continue;
                }
                elementsToSplit.emplace(ti);
                if(tr < globInho){
                    mesh.elements[ti].refinementPattern = 
                        Ddpca::OctreeElement::REFINEMENT_ZETA;//_ZETA_XI
                }
                else{
                    mesh.elements[ti].refinementPattern = 
                        Ddpca::OctreeElement::REFINEMENT_FULL;
                }
            }
            curvilinearInterpolation.clear();
            for(const auto &iterSpel : elementsToSplit){
                CurvilinearInterpolation(iterSpel, curvilinearInterpolation);
            }
            mesh.Refine(elementsToSplit, subElements, curvilinearInterpolation);
        }
        //
        mesh.OutputMesh(directoryPath, 0);
    }

    void ConstraintLoad(){
        //displacement constraint
        for(const auto& iterNoco: mesh.node2Coordinate){
            if((iterNoco.second)[2] <= 1.0E-10){
                boundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                boundary.constrainedDof.emplace(3 * iterNoco.first + 1, 0.0);
                boundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
            }
        }
        //torque load
        Ddpca::Real torqScal = 2.0 * torqLoad / (std::pow(outeRadi, 4.0) - std::pow(inneRadi, 4.0)) / Ddpca::PI;
        Ddpca::I64 ti_max = mesh.elements.size();
        Ddpca::I64 numbPart = Ddpca::threadManager.one2oneS2S.size();
        std::vector<Ddpca::I64> startIndex(numbPart), endIndex(numbPart);
        Ddpca::EvenlyDistribute(ti_max, numbPart, startIndex, endIndex);
        std::mutex tempMutex;
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.one2oneS2S, [&](Ddpca::I64 tp){
            Ddpca::I64 start_tp = startIndex[tp];
            Ddpca::I64 end_tp = endIndex[tp];
            for(Ddpca::I64 ti = start_tp; ti < end_tp; ++ ti){
                const Ddpca::OctreeElement& tempElem = mesh.elements[ti];
                if(tempElem.children.size() > 0){
                    continue;
                }
                const Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                for(Ddpca::I64 tj = 0; tj < hefaSize; ++ tj){
                    const Ddpca::I64 tk_size = Ddpca::hexaFace[tj].size();
                    std::vector<Ddpca::I64> inpuNode(tk_size);
                    std::vector<Ddpca::Coordinate> elemCoor(tk_size);
                    bool tempFlag = true;
                    for(Ddpca::I64 tk = 0; tk < tk_size; ++ tk){
                        inpuNode[tk] = tempElem.cornerNodes[Ddpca::hexaFace[tj][tk]];
                        auto iterNoco = mesh.node2Coordinate.find(inpuNode[tk]);
                        if((iterNoco->second)[2] < axiaLeng - 1.0E-10){
                            tempFlag = false;
                            break;
                        }
                        elemCoor[tk][0] = (iterNoco->second)[0];
                        elemCoor[tk][1] = (iterNoco->second)[1];
                        elemCoor[tk][2] = (iterNoco->second)[2];
                    }
                    if(tempFlag == false){
                        continue;
                    }
                    std::vector<Ddpca::Real> tempForc(12, 0.0);
                    const Ddpca::BilinearQuadrature<2>& biliQuad = Ddpca::GetBilinearQuadrature();
                    for(Ddpca::I64 tp = 0; tp < biliQuad.numbGaussPoints; ++ tp){
                        Ddpca::Real jacobian_tp = Ddpca::BilinearQuadratureJacobian<2>(
                            biliQuad.gaussPoints[tp], elemCoor);
                        std::vector<Ddpca::Real> tempCoor(3, 0.0);
                        for(Ddpca::I64 tk = 0; tk < tk_size; ++ tk){
                            for(Ddpca::I64 tn = 0; tn < 3; ++ tn){
                                tempCoor[tn] += biliQuad.shapeFunctions(tp,tk) * elemCoor[tk][tn];
                            }
                        }
                        Ddpca::Real tempAngl = std::atan2(tempCoor[1], tempCoor[0]) + Ddpca::PI / 2.0;
                        Ddpca::Real tempAmpl = torqScal * std::sqrt(std::pow(tempCoor[0], 2.0) + std::pow(tempCoor[1], 2.0));
                        std::vector<Ddpca::Real> tempPres = {
                            tempAmpl * (Ddpca::Real)std::cos(tempAngl), tempAmpl * (Ddpca::Real)std::sin(tempAngl), 0.0};
                        Ddpca::Real tempWeight = biliQuad.weights[tp] * jacobian_tp;
                        for(Ddpca::I64 tk = 0; tk < tk_size; ++ tk){
                            for(Ddpca::I64 tn = 0; tn < 3; ++ tn){
                                tempForc[tk * 3 + tn] += tempWeight * biliQuad.shapeFunctions(tp,tk) * tempPres[tn];
                            }
                        }
                    }
                    {
                        std::lock_guard<std::mutex> lock(tempMutex);
                        for(Ddpca::I64 tk = 0; tk < tk_size; ++ tk){
                            boundary.LoadAccumulate(3 * inpuNode[tk] + 0, tempForc[3 * tk + 0]);
                            boundary.LoadAccumulate(3 * inpuNode[tk] + 1, tempForc[3 * tk + 1]);
                            boundary.LoadAccumulate(3 * inpuNode[tk] + 2, tempForc[3 * tk + 2]);
                        }
                    }
                }
            }
        });
    }

    void Test(){
        //
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.numbThreads / (globHomo + globInho + 2), globHomo + globInho + 2, 
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
    std::string directoryPath = "./TestNDDPureTorsion_";
    std::filesystem::create_directory(directoryPath);
    //
    PureTorsion pureTorsion;
    pureTorsion.directoryPath = directoryPath;
    pureTorsion.Test();
    //
    Ddpca::Finalize();
	return 1;
}
