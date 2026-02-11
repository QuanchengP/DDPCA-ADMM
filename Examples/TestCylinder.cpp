#include "../Mesh/BilinearQuadrature.hpp"
#include "../Decomposition/MultiDomains.hpp"
#include "../Contact/CSQuadSegment.hpp"
#include "../Contact/CurvedSurface.hpp"

#include <filesystem>

//example Hertz contact
class Cylinder : public Ddpca::MultiDomains {

public:

	/*********************************************************************************************/
    std::string directoryPath;
	Ddpca::I64 copyNumb = 8;//axial number of domains for each cylinder
	//0 - lower cylinder, 1/2 - middle cylinder, 3 - upper cylinder
	std::array<Ddpca::Real,4> radi = {0.02, 0.022, 0.022, 0.02};//radius of cylinder
	std::array<Ddpca::Real,4> leng = {0.02 / copyNumb, 0.02 / copyNumb, 
		0.02 / copyNumb, 0.02 / copyNumb};//length of cylinder
	std::array<Ddpca::Real,2> diviAngl 
        = {- 5.0 / 8.0 * Ddpca::PI, - 3.0 / 8.0 * Ddpca::PI};;//relative to - PI / 2.0
	std::array<std::array<std::array<Ddpca::Real,2>,3>,4> diviPoin;//auxiliary point for meshing
    std::array<std::array<Ddpca::I64,4>,4> diviNumb = {{ {{2, 2, 1, 16 / copyNumb}}, 
        {{2, 2, 1, 16 / copyNumb}}, {{2, 2, 1, 16 / copyNumb}}, {{2, 2, 1, 16 / copyNumb}} }};
	Ddpca::I64 globInho = 3;//level of inhomogeneous global refinement
	Ddpca::I64 globHomo = 0;//level of homogeneous global refinement
	Ddpca::I64 locaLeve = 7;//level of local refinement
	Ddpca::Real bandWidt = 100.0E-6;//predicted contact band width > real value
	Ddpca::Real loadInte = -50.0E3;//load intensity

	//cylindrical surface: 0 - lower cylinder, 1/2 - middle cylinder, 3 - upper cylinder
	std::array<Ddpca::CurvedSurface,4> cyliSurf;
	std::array<Ddpca::CurvedSurface,4> cyliSurf_1;//left
	std::array<Ddpca::CurvedSurface,4> cyliSurf_2;//right
	std::array<Ddpca::CurvedSurface,3> inteSurf;//interface between two middle cylinders

public:
	/*********************************************************************************************/

    void CoordinateTransform_0(Ddpca::Coordinate &tempCoor){
        tempCoor[0] = - tempCoor[0];
        tempCoor[1] = - tempCoor[1] - radi[0] - radi[1] - radi[2] - radi[3];
    }

    void CoordinateTransform_2(Ddpca::Coordinate &tempCoor){
        tempCoor[0] = - tempCoor[0];
        tempCoor[1] = - tempCoor[1] - radi[2] - radi[3];
    }

    void GenerateSurfaces(){
        //
        Ddpca::Coordinate tempCoor;
        for(Ddpca::I64 tg = 0; tg < 4; tg ++){
            Ddpca::I64 totaDivi_0 = diviNumb[tg][1] * (1 << (globInho + globHomo + locaLeve));
            Ddpca::I64 totaDivi_1 = diviNumb[tg][3] * (1 << (globHomo + locaLeve));
            cyliSurf[tg].Resize(totaDivi_0 + 1, totaDivi_1 + 1);
            for(Ddpca::I64 ti = 0; ti <= totaDivi_0; ti ++){
                Ddpca::Real angl_ti = diviAngl[0] + (diviAngl[1] - diviAngl[0]) / (Ddpca::Real)totaDivi_0 * ti;
                for(Ddpca::I64 tj = 0; tj <= totaDivi_1; tj ++){
                    tempCoor[0] = radi[tg] * std::cos(angl_ti);
                    tempCoor[1] = radi[tg] * std::sin(angl_ti);
                    tempCoor[2] = leng[tg] / (Ddpca::Real)totaDivi_1 * tj;
                    if(tg == 0){
                        CoordinateTransform_0(tempCoor);
                    }
                    else if(tg == 1){
                        tempCoor[1] = tempCoor[1] - radi[0] - radi[1];
                    }
                    else if(tg == 2){
                        CoordinateTransform_2(tempCoor);
                    }
                    cyliSurf[tg].Insert(ti, tj, tempCoor);
                }
            }
        }
        //
        for(Ddpca::I64 tg = 0; tg < 4; tg ++){
            Ddpca::I64 totaDivi_0 = diviNumb[tg][0] * (1 << (globInho + globHomo));
            Ddpca::I64 totaDivi_1 = diviNumb[tg][3] * (1 << (globHomo));
            cyliSurf_1[tg].Resize(totaDivi_0 + 1, totaDivi_1 + 1);
            cyliSurf_2[tg].Resize(totaDivi_0 + 1, totaDivi_1 + 1);
            for(Ddpca::I64 ti = 0; ti <= totaDivi_0; ti ++){
                Ddpca::Real angl_ti = - Ddpca::PI + (diviAngl[0] + Ddpca::PI) / (Ddpca::Real)totaDivi_0 * ti;
                for(Ddpca::I64 tj = 0; tj <= totaDivi_1; tj ++){
                    tempCoor[0] = radi[tg] * std::cos(angl_ti);
                    tempCoor[1] = radi[tg] * std::sin(angl_ti);
                    tempCoor[2] = leng[tg] / (Ddpca::Real)totaDivi_1 * tj;
                    if(tg == 0){
                        CoordinateTransform_0(tempCoor);
                    }
                    else if(tg == 1){
                        tempCoor[1] = tempCoor[1] - radi[0] - radi[1];
                    }
                    else if(tg == 2){
                        CoordinateTransform_2(tempCoor);
                    }
                    cyliSurf_1[tg].Insert(ti, tj, tempCoor);
                    tempCoor[0] = - tempCoor[0];
                    cyliSurf_2[tg].Insert(totaDivi_0 - ti, tj, tempCoor);
                }
            }
        }
        //
        for(Ddpca::I64 ta = 0; ta < 3; ta ++){
            Ddpca::I64 totaDivi_0;
            Ddpca::I64 totaDivi_1 = diviNumb[1][3] * (1 << globHomo);
            if(ta == 1){
                totaDivi_0 = diviNumb[1][1] * (1 << (globInho + globHomo));
            }
            else{
                totaDivi_0 = diviNumb[1][2] * (1 << (globInho + globHomo));
            }
            inteSurf[ta].Resize(totaDivi_0 + 1, totaDivi_1 + 1);
            for(Ddpca::I64 ti = 0; ti <= totaDivi_0; ti ++){
                Ddpca::Real xcoo_ti;
                if(ta == 0){
                    xcoo_ti = - radi[1] + (diviPoin[1][0][0] + radi[1]) / (Ddpca::Real)totaDivi_0 * ti;
                }
                else if(ta == 1){
                    xcoo_ti = diviPoin[1][0][0] 
                        + (- diviPoin[1][0][0] - diviPoin[1][0][0]) / (Ddpca::Real)totaDivi_0 * ti;
                }
                else if(ta == 2){
                    xcoo_ti = - diviPoin[1][0][0] 
                        + (radi[1] + diviPoin[1][0][0]) / (Ddpca::Real)totaDivi_0 * ti;
                }
                for(Ddpca::I64 tj = 0; tj <= totaDivi_1; tj ++){
                    tempCoor[0] = xcoo_ti;
                    tempCoor[1] = - radi[3] - radi[2];
                    tempCoor[2] = leng[1] / (Ddpca::Real)totaDivi_1 * tj;
                    inteSurf[ta].Insert(ti, tj, tempCoor);
                }
            }
        }
    }

    void MeshConstraintLoad(){
        //
        for(Ddpca::I64 tg = 0; tg < 4; tg ++){
            diviPoin[tg][0] = {- radi[tg] / 3.0, 0.0};
            diviPoin[tg][1] = {- radi[tg] / 5.0, - radi[tg] / 2.0};
            diviPoin[tg][2] = {radi[tg] / 5.0, - radi[tg] / 2.0};
        }
        std::array<Ddpca::Real,4> refeRadi = {- radi[3] - radi[2] - radi[1], 
            - radi[3] - radi[2] - radi[1], - radi[3], - radi[3]};
        GenerateSurfaces();
        //
        domains.resize(4 * 2 * copyNumb);
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.threadsPerDomain, domains.size(), 
            Ddpca::threadManager.domainS2S, Ddpca::threadManager.domainS2M);
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            if(tg >= 4){
                return;
            }
            //boundary nodes of block
            std::vector<std::array<Ddpca::Real,2>> uppeLine_0, uppeLine_1, uppeLine_2, downLine_0, downLine_1;
            uppeLine_0.resize(diviNumb[tg][0] + 1);
            downLine_0.resize(diviNumb[tg][0] + 1);
            uppeLine_1.resize(diviNumb[tg][1] / 2 + 1);
            downLine_1.resize(diviNumb[tg][1] / 2 + 1);
            uppeLine_2.resize(diviNumb[tg][1] / 2 + 1);
            for(Ddpca::I64 ti = 0; ti <= diviNumb[tg][0]; ti ++){
                uppeLine_0[ti][0] = (1.0 - (Ddpca::Real)ti / diviNumb[tg][0]) * diviPoin[tg][0][0] 
                    + (Ddpca::Real)ti / diviNumb[tg][0] * diviPoin[tg][1][0];
                uppeLine_0[ti][1] = (1.0 - (Ddpca::Real)ti / diviNumb[tg][0]) * diviPoin[tg][0][1] 
                    + (Ddpca::Real)ti / diviNumb[tg][0] * diviPoin[tg][1][1];
                Ddpca::Real tempAngl = - Ddpca::PI + (diviAngl[0] + Ddpca::PI) / diviNumb[tg][0] * (Ddpca::Real)ti;
                downLine_0[ti][0] = radi[tg] * std::cos(tempAngl);
                downLine_0[ti][1] = radi[tg] * std::sin(tempAngl);
            }
            for(Ddpca::I64 ti = 0; ti <= diviNumb[tg][1] / 2; ti ++){
                uppeLine_1[ti][0] = (1.0 - (Ddpca::Real)ti / diviNumb[tg][1]) * diviPoin[tg][1][0] 
                    + (Ddpca::Real)ti / diviNumb[tg][1] * diviPoin[tg][2][0];
                uppeLine_1[ti][1] = (1.0 - (Ddpca::Real)ti / diviNumb[tg][1]) * diviPoin[tg][1][1] 
                    + (Ddpca::Real)ti / diviNumb[tg][1] * diviPoin[tg][2][1];
                Ddpca::Real tempAngl = diviAngl[0] 
                    + (diviAngl[1] - diviAngl[0]) / diviNumb[tg][1] * (Ddpca::Real)ti;
                downLine_1[ti][0] = radi[tg] * std::cos(tempAngl);
                downLine_1[ti][1] = radi[tg] * std::sin(tempAngl);
                uppeLine_2[ti][0] = (1.0 - (Ddpca::Real)ti / diviNumb[tg][1]) * (-radi[tg] / 3.0)
                    + (Ddpca::Real)ti / diviNumb[tg][1] * (radi[tg] / 3.0);
                uppeLine_2[ti][1] = 0.0;
            }
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            //node of block 0~2
            Ddpca::Coordinate tempCoor;
            std::vector<std::vector<std::vector<std::vector<Ddpca::I64>>>> blocNode(3);
            blocNode[0].resize(diviNumb[tg][0] + 1);
            for(Ddpca::I64 ti = 0; ti <= diviNumb[tg][0]; ti ++){
                blocNode[0][ti].resize(diviNumb[tg][2] + 1);
                for(Ddpca::I64 tj = 0; tj <= diviNumb[tg][2]; tj ++){
                    blocNode[0][ti][tj].resize(diviNumb[tg][3] + 1);
                    Ddpca::Real tempX = (1.0 - (Ddpca::Real)tj / diviNumb[tg][2]) * downLine_0[ti][0] 
                        + (Ddpca::Real)tj / diviNumb[tg][2] * uppeLine_0[ti][0];
                    Ddpca::Real tempY = (1.0 - (Ddpca::Real)tj / diviNumb[tg][2]) * downLine_0[ti][1] 
                        + (Ddpca::Real)tj / diviNumb[tg][2] * uppeLine_0[ti][1];
                    for(Ddpca::I64 tk =  0; tk <= diviNumb[tg][3]; tk ++){
                        tempCoor[0] = tempX;
                        tempCoor[1] = tempY;
                        tempCoor[2] = leng[tg] / diviNumb[tg][3] * (Ddpca::Real)tk;
                        if(tg == 0){
                            CoordinateTransform_0(tempCoor);
                        }
                        else if(tg == 1){
                            tempCoor[1] = tempCoor[1] - radi[3] - radi[2];
                        }
                        else if(tg == 2){
                            CoordinateTransform_2(tempCoor);
                        }
                        blocNode[0][ti][tj][tk] = tgMesh.TryAddNode(tempCoor);
                    }
                }
            }
            blocNode[1].resize(diviNumb[tg][1] / 2 + 1);
            for(Ddpca::I64 ti = 0; ti <= diviNumb[tg][1] / 2; ti ++){
                blocNode[1][ti].resize(diviNumb[tg][2] + 1);
                for(Ddpca::I64 tj = 0; tj <= diviNumb[tg][2]; tj ++){
                    blocNode[1][ti][tj].resize(diviNumb[tg][3] + 1);
                    Ddpca::Real tempX = (1.0 - (Ddpca::Real)tj / diviNumb[tg][2]) * downLine_1[ti][0] 
                        + (Ddpca::Real)tj / diviNumb[tg][2] * uppeLine_1[ti][0];
                    Ddpca::Real tempY = (1.0 - (Ddpca::Real)tj / diviNumb[tg][2]) * downLine_1[ti][1] 
                        + (Ddpca::Real)tj / diviNumb[tg][2] * uppeLine_1[ti][1];
                    for(Ddpca::I64 tk = 0; tk <= diviNumb[tg][3]; tk ++){
                        tempCoor[0] = tempX;
                        tempCoor[1] = tempY;
                        tempCoor[2] = leng[tg] / diviNumb[tg][3] * (Ddpca::Real)tk;
                        if(tg == 0){
                            CoordinateTransform_0(tempCoor);
                        }
                        else if(tg == 1){
                            tempCoor[1] = tempCoor[1] - radi[3] - radi[2];
                        }
                        else if(tg == 2){
                            CoordinateTransform_2(tempCoor);
                        }
                        blocNode[1][ti][tj][tk] = tgMesh.TryAddNode(tempCoor);
                    }
                }
            }
            blocNode[2].resize(diviNumb[tg][1] / 2 + 1);
            for(Ddpca::I64 ti = 0; ti <= diviNumb[tg][1] / 2; ti ++){
                blocNode[2][ti].resize(diviNumb[tg][0] + 1);
                for(Ddpca::I64 tj = 0; tj <= diviNumb[tg][0]; tj ++){
                    blocNode[2][ti][tj].resize(diviNumb[tg][3] + 1);
                    Ddpca::Real tempX = (1.0 - (Ddpca::Real)tj / diviNumb[tg][0]) * uppeLine_1[ti][0] 
                        + (Ddpca::Real)tj / diviNumb[tg][0] * uppeLine_2[ti][0];
                    Ddpca::Real tempY = (1.0 - (Ddpca::Real)tj / diviNumb[tg][0]) * uppeLine_1[ti][1] 
                        + (Ddpca::Real)tj / diviNumb[tg][0] * uppeLine_2[ti][1];
                    for(Ddpca::I64 tk = 0; tk <= diviNumb[tg][3]; tk ++){
                        tempCoor[0] = tempX;
                        tempCoor[1] = tempY;
                        tempCoor[2] = leng[tg] / diviNumb[tg][3] * (Ddpca::Real)tk;
                        if(tg == 0){
                            CoordinateTransform_0(tempCoor);
                        }
                        else if(tg == 1){
                            tempCoor[1] = tempCoor[1] - radi[3] - radi[2];
                        }
                        else if(tg == 2){
                            CoordinateTransform_2(tempCoor);
                        }
                        blocNode[2][ti][tj][tk] = tgMesh.TryAddNode(tempCoor);
                    }
                }
            }
            //elements
            Ddpca::OctreeElement tempElem;
            const Ddpca::I64 blnoSize_0 = blocNode.size();
            for(Ddpca::I64 ti = 0; ti < blnoSize_0; ti ++){
                const Ddpca::I64 blnoSize_1 = blocNode[ti].size();
                for(Ddpca::I64 tj = 0; tj < blnoSize_1 - 1; tj ++){
                    const Ddpca::I64 blnoSize_2 = blocNode[ti][tj].size();
                    for(Ddpca::I64 tk = 0; tk < blnoSize_2 - 1; tk ++){
                        const Ddpca::I64 blnoSize_3 = blocNode[ti][tj][tk].size();
                        for(Ddpca::I64 tm = 0; tm < blnoSize_3 - 1; tm ++){
                            tempElem.parent = -1;
                            tempElem.cornerNodes[0] = blocNode[ti][tj][tk][tm];
                            tempElem.cornerNodes[1] = blocNode[ti][tj + 1][tk][tm];
                            tempElem.cornerNodes[2] = blocNode[ti][tj + 1][tk + 1][tm];
                            tempElem.cornerNodes[3] = blocNode[ti][tj][tk + 1][tm];
                            tempElem.cornerNodes[4] = blocNode[ti][tj][tk][tm + 1];
                            tempElem.cornerNodes[5] = blocNode[ti][tj + 1][tk][tm + 1];
                            tempElem.cornerNodes[6] = blocNode[ti][tj + 1][tk + 1][tm + 1];
                            tempElem.cornerNodes[7] = blocNode[ti][tj][tk + 1][tm + 1];
                            tempElem.level = 0;
                            tempElem.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
                            tempElem.children.resize(0);
                            tgMesh.AddElement(tempElem);
                        }
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
                    tgMesh.elements[ti].refinementPattern = (tr < globInho) ? 
                        Ddpca::OctreeElement::REFINEMENT_XI_ETA : Ddpca::OctreeElement::REFINEMENT_FULL;
                }
                curvInte.clear();
                cyliSurf[tg].Refine(tgMesh, elementsToSplit, curvInte);
                cyliSurf_1[tg].Refine(tgMesh, elementsToSplit, curvInte);
                cyliSurf_2[tg].Refine(tgMesh, elementsToSplit, curvInte);
                if(tg == 1 || tg == 2){
                    const Ddpca::I64 insuSize = inteSurf.size();
                    for(Ddpca::I64 ti = 0; ti < insuSize; ti ++){
                        inteSurf[ti].Refine(tgMesh, elementsToSplit, curvInte);
                    }
                }
                tgMesh.Refine(elementsToSplit, subElements, curvInte);
            }
            //local refinement
            for(Ddpca::I64 tr = 0; tr < locaLeve; tr ++){
                elementsToSplit.clear();
                //***********************************************************************************
                Ddpca::I64 numbDael = tgMesh.elements.size();
                for(Ddpca::I64 ti = 0; ti < numbDael; ti ++){
                    if(tgMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    const auto &tempElem = tgMesh.elements[ti];
                    bool tempFlag = false;
                    for(Ddpca::I64 tj = 0; tj < 8; tj ++){
                        auto iterNoco = tgMesh.node2Coordinate.find(tempElem.cornerNodes[tj]);
                        if(tg == 0 && std::abs((iterNoco->second)[0]) <= bandWidt
                            && (iterNoco->second)[1] >= refeRadi[tg] - 2.0 * bandWidt){
                            tempFlag = true;
                            break;
                        }
                        else if(tg == 1 && std::abs((iterNoco->second)[0]) <= bandWidt
                            && (iterNoco->second)[1] <= refeRadi[tg] + 2.0 * bandWidt){
                            tempFlag = true;
                            break;
                        }
                        else if(tg == 2 && std::abs((iterNoco->second)[0]) <= bandWidt
                            && (iterNoco->second)[1] >= refeRadi[tg] - 2.0 * bandWidt){
                            tempFlag = true;
                            break;
                        }
                        else if(tg == 3 && std::abs((iterNoco->second)[0]) <= bandWidt
                            && (iterNoco->second)[1] <= refeRadi[tg] + 2.0 * bandWidt){
                            tempFlag = true;
                            break;
                        }
                    }
                    if(tempFlag == true){
                        elementsToSplit.emplace(ti);
                        tgMesh.elements[ti].refinementPattern = Ddpca::OctreeElement::REFINEMENT_FULL;
                    }
                }
                //***********************************************************************************
                curvInte.clear();
                cyliSurf[tg].Refine(tgMesh, elementsToSplit, curvInte);
                tgMesh.Refine(elementsToSplit, subElements, curvInte);
            }
            //
            tgMesh.OutputMesh(directoryPath, tg);
            //displacement constraint must at first
            Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
            for(const auto& iterNoco : tgMesh.node2Coordinate){
                if(tg == 0 && 
                    (iterNoco.second)[1] <= - radi[0] - radi[1] - radi[2] - radi[3] + 1.0E-10){
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 1, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                }
                else if((tg == 1 || tg == 2) 
                    && std::abs((iterNoco.second)[1] + radi[2] + radi[3]) <= 1.0E-10){
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                }
                else if(tg == 3 && (iterNoco.second)[1] >= -1.0E-10){
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                }
            }
            //load must at next
            if(tg == 3){
                //the normal value of one node: no need for * 0.5!
                Ddpca::Real loadIncr = loadInte * leng[tg] / (diviNumb[tg][3] * (1 << (globHomo)));
                for(const auto& iterNoco : tgMesh.node2Coordinate){
                    if((iterNoco.second)[1] >= -1.0E-10 && std::abs((iterNoco.second)[0]) <= 1.0E-10){
                        Ddpca::Real tempFact = 0.5; // two domain share the load
                        if((iterNoco.second)[2] <= 1.0E-10 
                            || (iterNoco.second)[2] >= leng[tg] - 1.0E-10){
                            tempFact = 0.25; // endpoint
                        }
                        Ddpca::I64 free_tm = 3 * iterNoco.first + 1;
                        tgBoundary.LoadAccumulate(free_tm, tempFact * loadIncr);
                    }
                }
            }
        });
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            if(tg < 4 || tg >= 8){
                return;
            }
            domains[tg] = domains[tg - 4];
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::DenseMatrix rotationMatrix(3,3, {-1.0, 0.0, 0.0, 
                0.0, 1.0, 0.0, 
                0.0, 0.0, -1.0});
            Ddpca::Coordinate translationVector(0.0, 0.0, leng[tg - 4]);
            tgMesh.RigidRotationTranslation(rotationMatrix, translationVector);
            tgMesh.OutputMesh(directoryPath, tg);
        });
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            if(tg < 8){
                return;
            }
            domains[tg] = domains[tg % 8];
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::DenseMatrix rotationMatrix(3,3, {1.0, 0.0, 0.0, 
                0.0, 1.0, 0.0, 
                0.0, 0.0, 1.0});
            Ddpca::Coordinate translationVector(0.0, 0.0, (tg / 8) * leng[(tg % 8) % 4]);
            tgMesh.RigidRotationTranslation(rotationMatrix, translationVector);
            tgMesh.OutputMesh(directoryPath, tg);
        });
    }

    void GenerateInterfaces(){
        //
	    Ddpca::I64 inteSize = 6 * copyNumb + 8 * (copyNumb - 1) + 4 * copyNumb;
        interfaces.resize(inteSize);
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.threadsPerInterface, inteSize, 
            Ddpca::threadManager.interfaceS2S, Ddpca::threadManager.interfaceS2M);
        //
        for(Ddpca::I64 ta = 0; ta < copyNumb; ta ++){
            interfaces[6 * ta + 0].domainIndex = {ta * 8 + 0, ta * 8 + 5};
            interfaces[6 * ta + 1].domainIndex = {ta * 8 + 4, ta * 8 + 1};
            interfaces[6 * ta + 2].domainIndex = {ta * 8 + 2, ta * 8 + 7};
            interfaces[6 * ta + 3].domainIndex = {ta * 8 + 6, ta * 8 + 3};
            interfaces[6 * ta + 4].domainIndex = {ta * 8 + 5, ta * 8 + 2};
            interfaces[6 * ta + 5].domainIndex = {ta * 8 + 1, ta * 8 + 6};
            //
            for(Ddpca::I64 ti = 0; ti < 6; ti ++){
                interfaces[6 * ta + ti].frictionCoefficient = 0.0;
            }
        }
        for(Ddpca::I64 ta = 0; ta < copyNumb - 1; ta ++){
            for(Ddpca::I64 tb = 0; tb < 8; tb ++){
                Ddpca::I64 ts = 6 * copyNumb + 8 * ta + tb;
                interfaces[ts].domainIndex = {ta * 8 + tb, (ta + 1) * 8 + tb};
                interfaces[ts].frictionCoefficient = -1.0;
            }
        }
        for(Ddpca::I64 ta = 0; ta < copyNumb; ta ++){
            for(Ddpca::I64 tb = 0; tb < 4; tb ++){
                Ddpca::I64 ts = 6 * copyNumb + 8 * (copyNumb - 1) + 4 * ta + tb;
                interfaces[ts].domainIndex = {ta * 8 + tb, ta * 8 + tb + 4};
                interfaces[ts].frictionCoefficient = -1.0;
            }
        }
        CalculatePenaltyParameter();
        //
        std::array<std::array<Ddpca::I64,2>,11> buckNumb;
        //buckNumb can not be too large, the smaller is the safer (although lower efficiency)
        buckNumb[0] = {1 << (locaLeve - 4), 
            std::min(diviNumb[0][3], diviNumb[1][3]) * (1 << (globHomo + locaLeve - 1))
        };
        buckNumb[1] = {1 << (locaLeve - 4), 
            std::min(diviNumb[2][3], diviNumb[3][3]) * (1 << (globHomo + locaLeve - 1))
        };
        buckNumb[2] = {(diviNumb[1][2] + diviNumb[1][1] / 2) * (1 << (globInho + globHomo - 1)), 
            diviNumb[1][3] * (1 << std::max((Ddpca::I64)0, globHomo - 1))
        };
        for(Ddpca::I64 ts = 0; ts < 4; ts ++){
            buckNumb[3 + ts] = {
                (diviNumb[ts][2] + diviNumb[ts][1] / 2) * (1 << (globInho + globHomo - 1)), 
                (diviNumb[ts][0] + diviNumb[ts][1]) * (1 << (globInho + globHomo - 1))
            };
        }
        for(Ddpca::I64 ts = 0; ts < 4; ts ++){
            buckNumb[7 + ts] = {
                (diviNumb[ts][0] + diviNumb[ts][1]) * (1 << (globInho + globHomo - 1)), 
                1 << (globHomo + locaLeve - 1)//matching meshes
            };
        }
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.interfaceS2S, [&](Ddpca::I64 ts){
            //
            if(!((ts < 6) 
                || (0 <= ts - 6 * copyNumb && ts - 6 * copyNumb < 8) 
                || (0 <= ts - (6 * copyNumb + 8 * (copyNumb - 1)) && ts - (6 * copyNumb + 8 * (copyNumb - 1)) < 4))){
                return;
            }
            //
            Ddpca::I64 tg_mast = interfaces[ts].domainIndex[0];
            Ddpca::I64 tg_slav = interfaces[ts].domainIndex[1];
            const Ddpca::Mesh& masterMesh = domains[tg_mast].mesh;
            const Ddpca::Mesh& slaveMesh = domains[tg_slav].mesh;
            if(ts < 6){
                //
                if(ts <= 3){
                    Ddpca::CurvedSurface& mastSurf = cyliSurf[tg_mast % 4];
                    Ddpca::CurvedSurface& slavSurf = cyliSurf[tg_slav % 4];
                    mastSurf.Initialize();
                    while(mastSurf.Increment(masterMesh)){
                        std::array<Ddpca::I64, 4> tempFace = mastSurf.currentFace;
                        bool tempFlag = true;
                        for(Ddpca::I64 ti = 0; ti < 4; ti ++){
                            auto iterNoco = masterMesh.node2Coordinate.find(tempFace[ti]);
                            if(std::abs((iterNoco->second)[0]) > bandWidt){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == true){
                            interfaces[ts].masterSegments.emplace_back(mastSurf.currentFace);
                        }
                    }
                    slavSurf.Initialize();
                    while(slavSurf.Increment(slaveMesh)){
                        std::array<Ddpca::I64, 4> tempFace = slavSurf.currentFace;
                        bool tempFlag = true;
                        for(Ddpca::I64 ti = 0; ti < 4; ti ++){
                            auto iterNoco = slaveMesh.node2Coordinate.find(tempFace[ti]);
                            if(std::abs((iterNoco->second)[0]) > bandWidt){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == true){
                            interfaces[ts].slaveSegments.emplace_back(slavSurf.currentFace);
                        }
                    }
                }
                else{
                    const Ddpca::I64 insuSize = inteSurf.size();
                    for(Ddpca::I64 ta = 0; ta < insuSize; ta ++){
                        Ddpca::CurvedSurface& mastSurf = inteSurf[ta];
                        Ddpca::CurvedSurface& slavSurf = inteSurf[ta];
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
                    }
                }
                interfaces[ts].OutputSegments(directoryPath, ts);
                //
                Ddpca::I64 inmaSize = interfaces[ts].masterSegments.size();
                std::array<std::vector<Ddpca::Real>, 2> masterLocal;
                masterLocal[0].resize(inmaSize);
                masterLocal[1].resize(inmaSize);
                Ddpca::I64 node_tj;
                Ddpca::Real tempXi, tempEta;
                Ddpca::Coordinate tempCoor;
                for(Ddpca::I64 ti = 0; ti < inmaSize; ++ ti){
                    tempXi = 0.0;
                    tempEta = 0.0;
                    for(Ddpca::I64 tj = 0; tj < 4; ++ tj){
                        node_tj = interfaces[ts].masterSegments[ti][tj];
                        auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                        tempCoor = iterNoco->second;
                        tempXi += tempCoor[0];
                        tempEta += tempCoor[2];
                    }
                    masterLocal[0][ti] = tempXi / 4.0;
                    masterLocal[1][ti] = tempEta / 4.0;
                }
                interfaces[ts].BucketSort(masterLocal, {buckNumb[ts / 2][0], buckNumb[ts / 2][1]});
                Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
                std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
                for(Ddpca::I64 ti = 0; ti < inslSize; ++ ti){
                    for(Ddpca::I64 tj = 0; tj < 4; ++ tj){
                        node_tj = interfaces[ts].slaveSegments[ti][tj];
                        auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                        tempCoor = iterNoco->second;
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[0];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[2];
                    }
                }
                interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                    Ddpca::threadManager.interfaceS2M[ts], 1, 1.0E12);
                interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            }
            else if(0 <= ts - 6 * copyNumb && ts - 6 * copyNumb < 8){
                //
                const Ddpca::I64 tb = ts - 6 * copyNumb;
                Ddpca::I64 numbDael = masterMesh.elements.size();
                std::array<Ddpca::I64, 4> tempNode;
                const Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                for(Ddpca::I64 ti = 0; ti < numbDael; ++ ti){
                    if(masterMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    for(Ddpca::I64 tj = 0; tj < hefaSize; tj ++){
                        bool tempFlag = true;
                        for(Ddpca::I64 tk = 0; tk < 4; tk ++){
                            tempNode[tk] = masterMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                            auto iterNoco = masterMesh.node2Coordinate.find(tempNode[tk]);
                            if(std::abs((iterNoco->second)[2] - leng[tb % 4]) > 1.0E-10){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == true){
                            interfaces[ts].masterSegments.emplace_back(tempNode);
                        }
                    }
                }
                Ddpca::I64 numbDael_1 = slaveMesh.elements.size();
                for(Ddpca::I64 ti = 0; ti < numbDael_1; ++ ti){
                    if(slaveMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    for(Ddpca::I64 tj = 0; tj < hefaSize; tj ++){
                        bool tempFlag = true;
                        for(Ddpca::I64 tk = 0; tk < 4; tk ++){
                            tempNode[tk] = slaveMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                            auto iterNoco = slaveMesh.node2Coordinate.find(tempNode[tk]);
                            if(std::abs((iterNoco->second)[2] - leng[tb % 4]) > 1.0E-10){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == true){
                            interfaces[ts].slaveSegments.emplace_back(tempNode);
                        }
                    }
                }
                interfaces[ts].OutputSegments(directoryPath, ts);
                //
                Ddpca::I64 inmaSize = interfaces[ts].masterSegments.size();
                std::array<std::vector<Ddpca::Real>, 2> masterLocal;
                masterLocal[0].resize(inmaSize);
                masterLocal[1].resize(inmaSize);
                Ddpca::I64 node_tj;
                Ddpca::Real tempXi, tempEta;
                Ddpca::Coordinate tempCoor;
                for(Ddpca::I64 ti = 0; ti < inmaSize; ti ++){
                    tempXi = 0.0;
                    tempEta = 0.0;
                    for(Ddpca::I64 tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].masterSegments[ti][tj];
                        auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                        tempCoor = iterNoco->second;
                        tempXi += tempCoor[0];
                        tempEta += tempCoor[1];
                    }
                    masterLocal[0][ti] = tempXi / 4.0;
                    masterLocal[1][ti] = tempEta / 4.0;
                }
                interfaces[ts].BucketSort(masterLocal, {buckNumb[3 + tb%4][0], buckNumb[3 + tb%4][1]});
                Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
                std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
                for(Ddpca::I64 ti = 0; ti < inslSize; ti ++){
                    for(Ddpca::I64 tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].slaveSegments[ti][tj];
                        auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                        tempCoor = iterNoco->second;
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[0];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[1];
                    }
                }
                interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                    Ddpca::threadManager.interfaceS2M[ts], 1, 1.0E12);
                interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            }
            else{
                Ddpca::I64 tb = ts - (6 * copyNumb + 8 * (copyNumb - 1));
                Ddpca::I64 numbDael = masterMesh.elements.size();
                std::array<Ddpca::I64, 4> tempNode;
                const Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                for(Ddpca::I64 ti = 0; ti < numbDael; ti ++){
                    if(masterMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    for(Ddpca::I64 tj = 0; tj < hefaSize; tj ++){
                        bool tempFlag = true;
                        for(Ddpca::I64 tk = 0; tk < 4; tk ++){
                            tempNode[tk] = masterMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                            auto iterNoco = masterMesh.node2Coordinate.find(tempNode[tk]);
                            if(std::abs((iterNoco->second)[0]) > 1.0E-10){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == true){
                            interfaces[ts].masterSegments.emplace_back(tempNode);
                        }
                    }
                }
                Ddpca::I64 numbDael_1 = slaveMesh.elements.size();
                for(Ddpca::I64 ti = 0; ti < numbDael_1; ti ++){
                    if(slaveMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    for(Ddpca::I64 tj = 0; tj < hefaSize; tj ++){
                        bool tempFlag = true;
                        for(Ddpca::I64 tk = 0; tk < 4; tk ++){
                            tempNode[tk] = slaveMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                            auto iterNoco = slaveMesh.node2Coordinate.find(tempNode[tk]);
                            if(std::abs((iterNoco->second)[0]) > 1.0E-10){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == true){
                            interfaces[ts].slaveSegments.emplace_back(tempNode);
                        }
                    }
                }
                interfaces[ts].OutputSegments(directoryPath, ts);
                //
                Ddpca::I64 inmaSize = interfaces[ts].masterSegments.size();
                std::array<std::vector<Ddpca::Real>, 2> masterLocal;
                masterLocal[0].resize(inmaSize);
                masterLocal[1].resize(inmaSize);
                Ddpca::I64 node_tj;
                Ddpca::Real tempXi, tempEta;
                Ddpca::Coordinate tempCoor;
                for(Ddpca::I64 ti = 0; ti < inmaSize; ti ++){
                    tempXi = 0.0;
                    tempEta = 0.0;
                    for(Ddpca::I64 tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].masterSegments[ti][tj];
                        auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                        tempCoor = iterNoco->second;
                        tempXi += tempCoor[1];
                        tempEta += tempCoor[2];
                    }
                    masterLocal[0][ti] = tempXi / 4.0;
                    masterLocal[1][ti] = tempEta / 4.0;
                }
                interfaces[ts].BucketSort(masterLocal, {buckNumb[7 + tb][0], buckNumb[7 + tb][1]});
                Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
                std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
                for(Ddpca::I64 ti = 0; ti < inslSize; ti ++){
                    for(Ddpca::I64 tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].slaveSegments[ti][tj];
                        auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                        tempCoor = iterNoco->second;
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[1];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[2];
                    }
                }
                interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                    Ddpca::threadManager.interfaceS2M[ts], 1, 1.0E12);
                interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            }
        });
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.interfaceS2S, [&](Ddpca::I64 ts){
            //
            if((ts < 6) 
                || (0 <= ts - 6 * copyNumb && ts - 6 * copyNumb < 8) 
                || (0 <= ts - (6 * copyNumb + 8 * (copyNumb - 1)) && ts - (6 * copyNumb + 8 * (copyNumb - 1)) < 4)){
                return;
            }
            Ddpca::ContactInterface& tempInterface = interfaces[ts];
            Ddpca::Real deltZcoo;
            Ddpca::I64 interfaceId;
            if(ts < 6 * copyNumb){
                interfaceId = ts % 6;
                deltZcoo = (ts / 6) * leng[0];
            }
            else if(ts - 6 * copyNumb < 8 * (copyNumb - 1)){
                interfaceId = 6 * copyNumb + (ts - 6 * copyNumb) % 8;
                deltZcoo = ((ts - 6 * copyNumb) / 8) * leng[0];
            }
            else{
                interfaceId = 6 * copyNumb + 8 * (copyNumb - 1) + (ts - (6 * copyNumb + 8 * (copyNumb - 1))) % 4;
                deltZcoo = ((ts - (6 * copyNumb + 8 * (copyNumb - 1))) / 4) * leng[0];
            }
            tempInterface.masterSegments = interfaces[interfaceId].masterSegments;
            tempInterface.slaveSegments = interfaces[interfaceId].slaveSegments;
            tempInterface.integralPoints = interfaces[interfaceId].integralPoints;
            for(Ddpca::IntegralPoint& tempInpo : tempInterface.integralPoints){
                tempInpo.contactPoint[0][2] += deltZcoo;
                tempInpo.contactPoint[1][2] += deltZcoo;
            }
            tempInterface.OutputSegments(directoryPath, ts);
            tempInterface.OutputIntegralPoints(directoryPath, ts);
        });
    }

    void Test(){
        //
        MeshConstraintLoad();
        GenerateInterfaces();

        mpLatin.realDomaLeve.assign(domains.size(), 2);
        Establish();
        ADMM(directoryPath);
    }
};

int main(int argc, char **argv){
    //
    Ddpca::Initialize(argc, argv);
	//
    std::string directoryPath = "./TestCylinder_";
    std::filesystem::create_directory(directoryPath);
    //
    Cylinder cylinder;
    cylinder.directoryPath = directoryPath;
    cylinder.Test();
    //
	Ddpca::Finalize();
	return 1;
}