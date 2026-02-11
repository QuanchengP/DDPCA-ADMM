#include "../Mesh/BilinearQuadrature.hpp"
#include "../Decomposition/MultiDomains.hpp"
#include "../Contact/CurvedSurface.hpp"

#include <filesystem>

//example pure torsional cylinder
class PureTorsion : public Ddpca::MultiDomains {

public:

	/*********************************************************************************************/
    std::string directoryPath;
	Ddpca::Real axiaLeng = 0.1;
	Ddpca::Real inneRadi = 0.015;
	Ddpca::Real outeRadi = 0.025;//I_p = PI / 32.0 * (D^4 - d^4) = 0.005340707511103
	std::array<Ddpca::I64, 3> diviNumb = {2, 16, 16};//8, 32, 64
	Ddpca::I64 globInho = 0;//0
	Ddpca::I64 globHomo = 4;//0
	Ddpca::Real torqLoad = 20.0;//N*m, u = T*l/G/I_p * R = 1.159111630361142e-06

	std::array<Ddpca::I64, 3> domaNumb = {1, 8, 4};//number of decomposed domains

    // intermediate variables: originally can be treated as local variable in function, 
    // but now has to be treated as global variable for the sake of threadManager parallelism.
    std::vector<std::array<Ddpca::CurvedSurface, 6>> blocSurf;
    std::array<Ddpca::I64, 3> xyzN;

public:
	/*********************************************************************************************/

    Ddpca::Coordinate CoordinateAverage(const std::vector<Ddpca::Coordinate> &inpuCoor){
        //
        Ddpca::Coordinate cyliCoor(0.0, 0.0, 0.0);
        Ddpca::I64 tempFlag_0 = 0;
        Ddpca::I64 tempFlag_1 = 0;
        Ddpca::I64 numbInco = inpuCoor.size();
        for(Ddpca::I64 tk = 0; tk < numbInco; ++ tk){
            cyliCoor[0] += std::sqrt(
                inpuCoor[tk][0] * inpuCoor[tk][0] + inpuCoor[tk][1] * inpuCoor[tk][1]);
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
        Ddpca::I64 tg, Ddpca::I64 ei, 
        std::map<std::vector<Ddpca::I64>,Ddpca::Coordinate> &planSurf){
        //
        Ddpca::Mesh& tgMesh = domains[tg].mesh;
        std::vector<Ddpca::I64> inpuNode(8);
        std::vector<Ddpca::Coordinate> inpuCoor(8);
        for(Ddpca::I64 tk = 0; tk < 8; ++ tk){
            inpuNode[tk] = tgMesh.elements[ei].cornerNodes[tk];
            inpuCoor[tk] = tgMesh.node2Coordinate[inpuNode[tk]];
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
                inpuNode[tk] = tgMesh.elements[ei].cornerNodes[Ddpca::hexaLine[tj][tk]];
                inpuCoor[tk] = tgMesh.node2Coordinate[inpuNode[tk]];
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
                inpuNode[tk] = tgMesh.elements[ei].cornerNodes[Ddpca::hexaFace[tj][tk]];
                inpuCoor[tk] = tgMesh.node2Coordinate[inpuNode[tk]];
            }
            Ddpca::Coordinate outpCoor = CoordinateAverage(inpuCoor);
            std::sort(inpuNode.begin(), inpuNode.end());
            planSurf.emplace(inpuNode, outpCoor);
        }
    }

    void GenerateMeshes(){
        //
        domains.resize(domaNumb[0] * domaNumb[1] * domaNumb[2]);
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.threadsPerDomain, domains.size(), 
            Ddpca::threadManager.domainS2S, Ddpca::threadManager.domainS2M);
        //
        std::array<Ddpca::I64, 3> diviReal = {
            diviNumb[0] / domaNumb[0], diviNumb[1] / domaNumb[1], diviNumb[2] / domaNumb[2]};
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            long tg_0 = tg / (domaNumb[1] * domaNumb[2]);
            long tg_1 = (tg % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
            long tg_2 = (tg % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
            //nodes, no need for nested parallelism: diviReal[0] is too small
            std::vector<std::vector<std::vector<Ddpca::I64>>> blocNode(
                diviReal[0] + 1, 
                std::vector<std::vector<Ddpca::I64>>(
                    diviReal[1] + 1, 
                    std::vector<Ddpca::I64>(diviReal[2] + 1, 0)
                )
            );
            Ddpca::Coordinate tempCoor;
            for(long ti = 0; ti <= diviReal[0]; ti ++){
                long ti_real = tg_0 * diviReal[0] + ti;
                double radi_ti = inneRadi + (outeRadi - inneRadi) / (double)(diviNumb[0]) * ti_real;
                for(long tj = 0; tj <= diviReal[1]; tj ++){
                    long tj_real = tg_1 * diviReal[1] + tj;
                    double angl_tj = 0.0 + (2.0 * Ddpca::PI - 0.0) / (double)(diviNumb[1]) * tj_real;
                    for(long tk =  0; tk <= diviReal[2]; tk ++){
                        long tk_real = tg_2 * diviReal[2] + tk;
                        tempCoor[0] = radi_ti * std::cos(angl_tj);
                        tempCoor[1] = radi_ti * std::sin(angl_tj);
                        tempCoor[2] = axiaLeng / diviNumb[2] * (double)tk_real;
                        blocNode[ti][tj][tk] = tgMesh.TryAddNode(tempCoor);
                    }
                }
            }
            //elements
            Ddpca::OctreeElement tempElem;
            for(Ddpca::I64 ti = 0; ti < diviReal[0]; ++ ti){
                for(Ddpca::I64 tj = 0; tj < diviReal[1]; ++ tj){
                    for(Ddpca::I64 tk = 0; tk < diviReal[2]; ++ tk){
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
                        tgMesh.AddElement(tempElem);
                    }
                }
            }
            //global refinement
            std::set<Ddpca::I64> elementsToSplit;
            std::map<Ddpca::I64, std::set<Ddpca::I64>> subElements;
            std::map<std::vector<Ddpca::I64>, Ddpca::Coordinate> curvInte;
            for(Ddpca::I64 tr = 0; tr < globInho + globHomo; ++ tr){
                elementsToSplit.clear();
                Ddpca::I64 numbDael = tgMesh.elements.size();
                for(Ddpca::I64 ti = 0; ti < numbDael; ++ ti){
                    if(tgMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    elementsToSplit.emplace(ti);
                    if(tr < globInho){
                        tgMesh.elements[ti].refinementPattern = 
                            Ddpca::OctreeElement::REFINEMENT_ZETA;//_ZETA_XI
                    }
                    else{
                        tgMesh.elements[ti].refinementPattern = 
                            Ddpca::OctreeElement::REFINEMENT_FULL;
                    }
                }
                curvInte.clear();
                for(const auto &iterSpel : elementsToSplit){
                    CurvilinearInterpolation(tg, iterSpel, curvInte);
                }
                tgMesh.Refine(elementsToSplit, subElements, curvInte);
            }
            //
            tgMesh.OutputMesh(directoryPath, tg);
        });
    }

    void ConstraintLoad(){
        //
        Ddpca::Real torqScal = 2.0 * torqLoad / (std::pow(outeRadi, 4.0) - std::pow(inneRadi, 4.0)) / Ddpca::PI;
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
            // displacement constraint
            for(const auto& iterNoco: tgMesh.node2Coordinate){
                if((iterNoco.second)[2] <= 1.0E-10){
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 1, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                }
            }
            //torque load
            Ddpca::I64 ti_max = tgMesh.elements.size();
            const Ddpca::BilinearQuadrature<2>& biliQuad = Ddpca::GetBilinearQuadrature();
            std::array<Ddpca::Real, 12> tempForc;
            std::array<Ddpca::Real, 3> tempCoor;
            std::array<Ddpca::Real, 3> tempPres;
            for(Ddpca::I64 ti = 0; ti < ti_max; ++ ti){
                const Ddpca::OctreeElement& tempElem = tgMesh.elements[ti];
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
                        auto iterNoco = tgMesh.node2Coordinate.find(inpuNode[tk]);
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
                    std::fill(std::execution::unseq, 
                        tempForc.begin(), tempForc.end(), 0.0);
                    for(Ddpca::I64 tp = 0; tp < biliQuad.numbGaussPoints; ++ tp){
                        Ddpca::Real jacobian_tp = Ddpca::BilinearQuadratureJacobian<2>(
                            biliQuad.gaussPoints[tp], elemCoor);
                        std::fill(std::execution::unseq, 
                            tempCoor.begin(), tempCoor.end(), 0.0);
                        for(Ddpca::I64 tk = 0; tk < tk_size; ++ tk){
                            for(Ddpca::I64 tn = 0; tn < 3; ++ tn){
                                tempCoor[tn] += biliQuad.shapeFunctions(tp,tk) * elemCoor[tk][tn];
                            }
                        }
                        Ddpca::Real tempAngl = std::atan2(tempCoor[1], tempCoor[0]) + Ddpca::PI / 2.0;
                        Ddpca::Real tempAmpl = torqScal * std::sqrt(std::pow(tempCoor[0], 2.0) + std::pow(tempCoor[1], 2.0));
                        tempPres = {
                            tempAmpl * (Ddpca::Real)std::cos(tempAngl), tempAmpl * (Ddpca::Real)std::sin(tempAngl), 0.0};
                        Ddpca::Real tempWeight = biliQuad.weights[tp] * jacobian_tp;
                        for(Ddpca::I64 tk = 0; tk < tk_size; ++ tk){
                            for(Ddpca::I64 tn = 0; tn < 3; ++ tn){
                                tempForc[tk * 3 + tn] += tempWeight 
                                    * biliQuad.shapeFunctions(tp,tk) * tempPres[tn];
                            }
                        }
                    }
                    for(Ddpca::I64 tk = 0; tk < tk_size; ++ tk){
                        tgBoundary.LoadAccumulate(3 * inpuNode[tk] + 0, tempForc[3 * tk + 0]);
                        tgBoundary.LoadAccumulate(3 * inpuNode[tk] + 1, tempForc[3 * tk + 1]);
                        tgBoundary.LoadAccumulate(3 * inpuNode[tk] + 2, tempForc[3 * tk + 2]);
                    }
                }
            }
        });
    }

    void GenerateSurfaces(){
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tb){
            //
            std::array<Ddpca::I64, 3> diviNumb_es = {
                diviNumb[0] * (1 << (globInho + globHomo)), 
                diviNumb[1] * (1 << (globInho + globHomo)), 
                diviNumb[2] * (1 << (globInho + globHomo))};
            std::array<Ddpca::I64, 3> diviReal = {
                diviNumb_es[0] / domaNumb[0], 
                diviNumb_es[1] / domaNumb[1], 
                diviNumb_es[2] / domaNumb[2]};
            //
            blocSurf[tb][0].Resize(diviReal[1] + 1, diviReal[2] + 1);
            blocSurf[tb][1].Resize(diviReal[1] + 1, diviReal[2] + 1);
            blocSurf[tb][2].Resize(diviReal[2] + 1, diviReal[0] + 1);
            blocSurf[tb][3].Resize(diviReal[2] + 1, diviReal[0] + 1);
            blocSurf[tb][4].Resize(diviReal[0] + 1, diviReal[1] + 1);
            blocSurf[tb][5].Resize(diviReal[0] + 1, diviReal[1] + 1);
            Ddpca::I64 tb_0 = tb / (domaNumb[1] * domaNumb[2]);
            Ddpca::I64 tb_1 = (tb % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
            Ddpca::I64 tb_2 = (tb % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
            //nodes
            for(Ddpca::I64 ti = 0; ti <= diviReal[0]; ti += diviReal[0]){
                Ddpca::I64 ti_real = tb_0 * diviReal[0] + ti;
                Ddpca::Real radi_ti = inneRadi + (outeRadi - inneRadi) 
                    / (Ddpca::Real)(diviNumb_es[0]) * (Ddpca::Real)ti_real;
                for(Ddpca::I64 tj = 0; tj <= diviReal[1]; ++ tj){
                    Ddpca::I64 tj_real = tb_1 * diviReal[1] + tj;
                    Ddpca::Real angl_tj = 0.0 + (2.0 * Ddpca::PI - 0.0) 
                        / (Ddpca::Real)(diviNumb_es[1]) * tj_real;
                    for(Ddpca::I64 tk = 0; tk <= diviReal[2]; ++ tk){
                        Ddpca::I64 tk_real = tb_2 * diviReal[2] + tk;
                        Ddpca::Coordinate tempCoor(radi_ti * std::cos(angl_tj), radi_ti * std::sin(angl_tj), 
                            axiaLeng / diviNumb_es[2] * (Ddpca::Real)tk_real
                        );
                        if(ti == 0){
                            blocSurf[tb][0].Insert(tj, tk, tempCoor);
                        }
                        else{
                            blocSurf[tb][1].Insert(tj, tk, tempCoor);
                        }
                    }
                }
            }
            for(Ddpca::I64 ti = 0; ti <= diviReal[0]; ++ ti){
                Ddpca::I64 ti_real = tb_0 * diviReal[0] + ti;
                Ddpca::Real radi_ti = inneRadi + (outeRadi - inneRadi) 
                    / (Ddpca::Real)(diviNumb_es[0]) * (Ddpca::Real)ti_real;
                for(Ddpca::I64 tj = 0; tj <= diviReal[1]; tj += diviReal[1]){
                    Ddpca::I64 tj_real = tb_1 * diviReal[1] + tj;
                    Ddpca::Real angl_tj = 0.0 + (2.0 * Ddpca::PI - 0.0) 
                        / (Ddpca::Real)(diviNumb_es[1]) * tj_real;
                    for(Ddpca::I64 tk = 0; tk <= diviReal[2]; ++ tk){
                        Ddpca::I64 tk_real = tb_2 * diviReal[2] + tk;
                        Ddpca::Coordinate tempCoor(radi_ti * std::cos(angl_tj), radi_ti * std::sin(angl_tj), 
                            axiaLeng / diviNumb_es[2] * (Ddpca::Real)tk_real
                        );
                        if(tj == 0){
                            blocSurf[tb][2].Insert(tk, ti, tempCoor);
                        }
                        else{
                            blocSurf[tb][3].Insert(tk, ti, tempCoor);
                        }
                    }
                }
            }
            for(Ddpca::I64 ti = 0; ti <= diviReal[0]; ++ ti){
                Ddpca::I64 ti_real = tb_0 * diviReal[0] + ti;
                Ddpca::Real radi_ti = inneRadi + (outeRadi - inneRadi) 
                    / (Ddpca::Real)(diviNumb_es[0]) * (Ddpca::Real)ti_real;
                for(Ddpca::I64 tj = 0; tj <= diviReal[1]; ++ tj){
                    Ddpca::I64 tj_real = tb_1 * diviReal[1] + tj;
                    Ddpca::Real angl_tj = 0.0 + (2.0 * Ddpca::PI - 0.0) 
                        / (Ddpca::Real)(diviNumb_es[1]) * tj_real;
                    for(Ddpca::I64 tk = 0; tk <= diviReal[2]; tk += diviReal[2]){
                        Ddpca::I64 tk_real = tb_2 * diviReal[2] + tk;
                        Ddpca::Coordinate tempCoor(radi_ti * std::cos(angl_tj), radi_ti * std::sin(angl_tj), 
                            axiaLeng / diviNumb_es[2] * (Ddpca::Real)tk_real
                        );
                        if(tk == 0){
                            blocSurf[tb][4].Insert(ti, tj, tempCoor);
                        }
                        else{
                            blocSurf[tb][5].Insert(ti, tj, tempCoor);
                        }
                    }
                }
            }
        });
    }

    void GenerateInterfaces(){
        //
        const Ddpca::I64 numbDomains = domains.size();
        blocSurf.resize(numbDomains);
        GenerateSurfaces();
        //
        xyzN = {
            (domaNumb[0] - 1) * domaNumb[1] * domaNumb[2], 
            (domaNumb[1] - 0) * domaNumb[0] * domaNumb[2], 
            (domaNumb[2] - 1) * domaNumb[0] * domaNumb[1]};
        interfaces.resize(xyzN[0] + xyzN[1] + xyzN[2]);
        Ddpca::I64 inteSize = interfaces.size();
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.threadsPerInterface, inteSize, 
            Ddpca::threadManager.interfaceS2S, Ddpca::threadManager.interfaceS2M);
        //
        for(Ddpca::I64 tg_0 = 0; tg_0 < domaNumb[0]; ++ tg_0){
            for(Ddpca::I64 tg_1 = 0; tg_1 < domaNumb[1]; ++ tg_1){
                for(Ddpca::I64 tg_2 = 0; tg_2 < domaNumb[2]; ++ tg_2){
                    Ddpca::I64 tg_m = tg_0 * domaNumb[1] * domaNumb[2] + tg_1 * domaNumb[2] + tg_2;
                    if(tg_0 <= domaNumb[0] - 2){
                        Ddpca::I64 tg_s = tg_m + domaNumb[1] * domaNumb[2];
                        Ddpca::I64 ts = tg_m;
                        interfaces[ts].domainIndex = {tg_m, tg_s};
                    }
                    if(tg_1 <= domaNumb[1] - 1){
                        Ddpca::I64 tg_s = tg_m + domaNumb[2];
                        if(tg_1 == domaNumb[1] - 1){
                            tg_s -= domaNumb[1] * domaNumb[2];
                        }
                        Ddpca::I64 ts = xyzN[0] + tg_1 * domaNumb[0] * domaNumb[2] 
                            + tg_0 * domaNumb[2] + tg_2;
                        interfaces[ts].domainIndex = {tg_m, tg_s};
                    }
                    if(tg_2 <= domaNumb[2] - 2){
                        Ddpca::I64 tg_s = tg_m + 1;
                        Ddpca::I64 ts = xyzN[0] + xyzN[1] 
                            + tg_2 * domaNumb[0] * domaNumb[1] + tg_0 * domaNumb[1] + tg_1;
                        interfaces[ts].domainIndex = {tg_m, tg_s};
                    }
                    
                }
            }
        }
        //
        for(Ddpca::I64 ts = 0; ts < inteSize; ++ ts){
            interfaces[ts].frictionCoefficient = -1.0;
        }
        CalculatePenaltyParameter();
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.interfaceS2S, [&](Ddpca::I64 ts){
            //
            Ddpca::I64 masterIndex, slaveIndex;
            if(ts < xyzN[0]){
                masterIndex = 1;
                slaveIndex = 0;
            }
            else if(ts < xyzN[0] + xyzN[1]){
                masterIndex = 3;
                slaveIndex = 2;
            }
            else{
                masterIndex = 5;
                slaveIndex = 4;
            }
            Ddpca::CurvedSurface& mastSurf = 
                blocSurf[interfaces[ts].domainIndex[0]][masterIndex];
            Ddpca::CurvedSurface& slavSurf = 
                blocSurf[interfaces[ts].domainIndex[1]][slaveIndex];
            //
            mastSurf.Initialize();
            const Ddpca::Mesh& masterMesh = domains[interfaces[ts].domainIndex[0]].mesh;
            while(mastSurf.Increment(masterMesh)){
                interfaces[ts].masterSegments.emplace_back(mastSurf.currentFace);
            }
            slavSurf.Initialize();
            const Ddpca::Mesh& slaveMesh = domains[interfaces[ts].domainIndex[1]].mesh;
            while(slavSurf.Increment(slaveMesh)){
                interfaces[ts].slaveSegments.emplace_back(slavSurf.currentFace);
            }
            interfaces[ts].OutputSegments(directoryPath, ts);
            //
            Ddpca::I64 inmaSize = interfaces[ts].masterSegments.size();
            std::array<std::vector<Ddpca::Real>, 2> masterLocal;
            masterLocal[0].resize(inmaSize);
            masterLocal[1].resize(inmaSize);
            Ddpca::I64 node_tj;
            Ddpca::Real tempXi, tempEta, radius_tj, angle_tj;
            Ddpca::Coordinate tempCoor;
            for(Ddpca::I64 ti = 0; ti < inmaSize; ++ ti){
                tempXi = 0.0;
                tempEta = 0.0;
                for(Ddpca::I64 tj = 0; tj < 4; ++ tj){
                    node_tj = interfaces[ts].masterSegments[ti][tj];
                    auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                    tempCoor = iterNoco->second;
                    radius_tj = std::sqrt(std::pow(tempCoor[0], 2.0) 
                        + std::pow(tempCoor[1], 2.0));
                    angle_tj = std::atan2(tempCoor[1], tempCoor[0]);
                    if(angle_tj < 0.0){
                        angle_tj += 2.0 * Ddpca::PI;
                    }
                    if(tempCoor[0] > 0.0 && std::abs(tempCoor[1]) < 1.0E-10){
                        angle_tj = 0.0;
                    }
                    if(ts < xyzN[0]){
                        tempXi += angle_tj;
                        tempEta += tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1]){
                        tempXi += radius_tj;
                        tempEta += tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
                        tempXi += radius_tj;
                        tempEta += angle_tj;
                    }
                }
                masterLocal[0][ti] = tempXi / 4.0;
                masterLocal[1][ti] = tempEta / 4.0;
            }
            if(ts < xyzN[0]){
                interfaces[ts].BucketSort(masterLocal, {
                    diviNumb[1] / domaNumb[1] * (1 << globHomo), 
                    diviNumb[2] / domaNumb[2] * (1 << globHomo)});
            }
            else if(ts < xyzN[0] + xyzN[1]){
                interfaces[ts].BucketSort(masterLocal, {
                    diviNumb[0] / domaNumb[0] * (1 << globHomo), 
                    diviNumb[2] / domaNumb[2] * (1 << globHomo)});
            }
            else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
                interfaces[ts].BucketSort(masterLocal, {
                    diviNumb[0] / domaNumb[0] * (1 << globHomo), 
                    diviNumb[1] / domaNumb[1] * (1 << globHomo)});
            }
            //
            Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
            std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
            for(Ddpca::I64 ti = 0; ti < inslSize; ++ ti){
                for(Ddpca::I64 tj = 0; tj < 4; ++ tj){
                    node_tj = interfaces[ts].slaveSegments[ti][tj];
                    auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                    tempCoor = iterNoco->second;
                    radius_tj = std::sqrt(std::pow(tempCoor[0], 2.0) 
                        + std::pow(tempCoor[1], 2.0));
                    angle_tj = std::atan2(tempCoor[1], tempCoor[0]);
                    if(angle_tj < 0.0){
                        angle_tj += 2.0 * Ddpca::PI;
                    }
                    if(tempCoor[0] > 0.0 && std::abs(tempCoor[1]) < 1.0E-10){
                        angle_tj = 0.0;
                    }
                    if(ts < xyzN[0]){
                        slaveLocal[ti][tj * 2 + 0] = angle_tj;
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1]){
                        slaveLocal[ti][tj * 2 + 0] = radius_tj;
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
                        slaveLocal[ti][tj * 2 + 0] = radius_tj;
                        slaveLocal[ti][tj * 2 + 1] = angle_tj;
                    }
                }
            }
            interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                Ddpca::threadManager.interfaceS2M[ts], 1, 1.0E12);
            interfaces[ts].OutputIntegralPoints(directoryPath, ts);
        });
    }

    void Test(){
        //
        GenerateMeshes();
        ConstraintLoad();
        GenerateInterfaces();

        mpLatin.realDomaLeve.assign(domains.size(), 1);
        Establish();
        ADMM(directoryPath);
    }
};

int main(int argc, char **argv){
    //
    Ddpca::Initialize(argc, argv);
	//
    std::string directoryPath = "./TestPureTorsion_";
    std::filesystem::create_directory(directoryPath);
    //
    PureTorsion pureTorsion;
    pureTorsion.directoryPath = directoryPath;
    pureTorsion.Test();
    //
	Ddpca::Finalize();
	return 1;
}
