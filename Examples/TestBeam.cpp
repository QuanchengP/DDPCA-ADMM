#include "../Mesh/BilinearQuadrature.hpp"
#include "../Decomposition/MultiDomains.hpp"
#include "../Contact/CurvedSurface.hpp"

#include <filesystem>

//example pure torsional cylinder
class Beam : public Ddpca::MultiDomains {

public:

	/*********************************************************************************************/
    std::string directoryPath;
    std::array<Ddpca::Real,3> leng = {1.0, 0.12, 0.06};//length of the beam
	Ddpca::Real lengFact = 1.0 / 3.0;//length reduction factor from fixed end to free end
	Ddpca::Real angl = 45.0 * Ddpca::PI / 180.0;//geometric angle
	std::array<Ddpca::I64,3> diviNumb = {64, 4, 2};//number of elements along different directions
	Ddpca::I64 globLeve = 4;//multigrid level, total levels of global mesh refinement
	std::array<Ddpca::I64,3> domaNumb = {16, 2, 1};//number of decomposed domains
	Ddpca::Real loadInte = - 8000.0;//centerline load intensity

    // intermediate variables: originally can be treated as local variable in function, 
    // but now has to be treated as global variable for the sake of threadManager parallelism.
    std::vector<std::array<Ddpca::CurvedSurface, 6>> blocSurf;
    std::array<Ddpca::I64, 3> xyzN;

public:
	/*********************************************************************************************/

    void CoordinateRotate(Ddpca::Coordinate& resultCoordinate, Ddpca::Real dire){
        Ddpca::Coordinate currentCoordinate = resultCoordinate;
        Ddpca::Real angl_temp = dire * angl * currentCoordinate[0] / leng[0];
        //col-major
        Ddpca::DenseMatrix rotationMatrix(3, 3, {
            1.0, 0.0, 0.0, 
            0.0, std::cos(angl_temp), std::sin(angl_temp), 
            0.0, -std::sin(angl_temp), std::cos(angl_temp)});
        Ddpca::GEMV(rotationMatrix, currentCoordinate.data, resultCoordinate.data);
    }

    void CoordinateAdjustment(Ddpca::I64 tg){
        Ddpca::Mesh& tgMesh = domains[tg].mesh;
        tgMesh.coordinate2Node.clear();
        for(auto& currentNodeCoordinate : tgMesh.node2Coordinate){
            CoordinateRotate(currentNodeCoordinate.second, 1.0);
            tgMesh.coordinate2Node.emplace(currentNodeCoordinate.second, currentNodeCoordinate.first);
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
            Ddpca::I64 tg_0 = tg / (domaNumb[1] * domaNumb[2]);
            Ddpca::I64 tg_1 = (tg % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
            Ddpca::I64 tg_2 = (tg % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
            //nodes, no need for nested parallelism: diviReal[0] is too small
            std::vector<std::vector<std::vector<Ddpca::I64>>> blocNode(
                diviReal[0] + 1, 
                std::vector<std::vector<Ddpca::I64>>(
                    diviReal[1] + 1, 
                    std::vector<Ddpca::I64>(diviReal[2] + 1, 0)
                )
            );
            Ddpca::Coordinate tempCoor;
            for(Ddpca::I64 ti = 0; ti <= diviReal[0]; ti ++){
                Ddpca::I64 ti_real = tg_0 * diviReal[0] + ti;
                tempCoor[0] = leng[0] / diviNumb[0] * ti_real;
                Ddpca::Real heig_ti = leng[1] * (1.0- (Ddpca::Real)ti_real / diviNumb[0] * lengFact);
                Ddpca::Real widt_ti = leng[2] * (1.0- (Ddpca::Real)ti_real / diviNumb[0] * lengFact);
                for(Ddpca::I64 tj = 0; tj <= diviReal[1]; tj ++){
                    Ddpca::I64 tj_real = tg_1 * diviReal[1] + tj;
                    tempCoor[1] = - heig_ti / 2.0 + heig_ti / diviNumb[1] * (Ddpca::Real)tj_real;
                    for(Ddpca::I64 tk =  0; tk <= diviReal[2]; tk ++){
                        Ddpca::I64 tk_real = tg_2 * diviReal[2] + tk;
                        tempCoor[2] = - widt_ti / 2.0 + widt_ti / diviNumb[2] * (Ddpca::Real)tk_real;
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
            for(Ddpca::I64 tr = 0; tr < globLeve; ++ tr){
                elementsToSplit.clear();
                Ddpca::I64 numbDael = tgMesh.elements.size();
                for(Ddpca::I64 ti = 0; ti < numbDael; ++ ti){
                    if(tgMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    elementsToSplit.emplace(ti);
                    tgMesh.elements[ti].refinementPattern = 
                        Ddpca::OctreeElement::REFINEMENT_FULL;
                }
                curvInte.clear();
                tgMesh.Refine(elementsToSplit, subElements, curvInte);
            }
            //
            CoordinateAdjustment(tg);
            tgMesh.OutputMesh(directoryPath, tg);
        });
    }

    void ConstraintLoad(){
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
            // displacement constraint
            for(const auto& iterNoco: tgMesh.node2Coordinate){
                if((iterNoco.second)[0] <= 1.0E-10){
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 1, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                }
            }
            //centerline load
            Ddpca::I64 ti_max = tgMesh.elements.size();
            const Ddpca::I64 heliSize = Ddpca::hexaLine.size();
            for(Ddpca::I64 ti = 0; ti < ti_max; ++ ti){
                const Ddpca::OctreeElement& tempElem = tgMesh.elements[ti];
                if(tempElem.children.size() > 0){
                    continue;
                }
                for(Ddpca::I64 tj = 0; tj < heliSize; tj ++){
                    const Ddpca::I64 tk_size = Ddpca::hexaLine[tj].size();
                    std::vector<Ddpca::I64> inpuNode(tk_size);
                    std::vector<Ddpca::Coordinate> inpuCoor(tk_size);
                    bool tempFlag = true;
                    for(Ddpca::I64 tk = 0; tk < tk_size; tk ++){
                        inpuNode[tk] = tempElem.cornerNodes[Ddpca::hexaLine[tj][tk]];
                        auto iterNoco = tgMesh.node2Coordinate.find(inpuNode[tk]);
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
                    for(Ddpca::I64 tk = 0; tk < tk_size; tk ++){
                        if(inpuCoor[tk][0] > 1.0E-10){
                            tgBoundary.LoadAccumulate(3 * inpuNode[tk] + 2, forcValu);
                        }
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
                diviNumb[0] * (1 << (globLeve)), 
                diviNumb[1] * (1 << (globLeve)), 
                diviNumb[2] * (1 << (globLeve))};
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
            Ddpca::Coordinate tempCoor;
            for(Ddpca::I64 ti = 0; ti <= diviReal[0]; ti += diviReal[0]){
                Ddpca::I64 ti_real = tb_0 * diviReal[0] + ti;
                Ddpca::Real tempX = leng[0] / diviNumb_es[0] * ti_real;
                Ddpca::Real heig_ti = leng[1] * (1.0- (Ddpca::Real)ti_real / diviNumb_es[0] * lengFact);
                Ddpca::Real widt_ti = leng[2] * (1.0- (Ddpca::Real)ti_real / diviNumb_es[0] * lengFact);
                for(Ddpca::I64 tj = 0; tj <= diviReal[1]; tj ++){
                    Ddpca::I64 tj_real = tb_1 * diviReal[1] + tj;
                    Ddpca::Real tempY = - heig_ti / 2.0 + heig_ti / diviNumb_es[1] * (Ddpca::Real)tj_real;
                    for(Ddpca::I64 tk =  0; tk <= diviReal[2]; tk ++){
                        Ddpca::I64 tk_real = tb_2 * diviReal[2] + tk;
                        Ddpca::Real tempZ = - widt_ti / 2.0 + widt_ti / diviNumb_es[2] * (Ddpca::Real)tk_real;
                        tempCoor = {tempX, tempY, tempZ};
                        CoordinateRotate(tempCoor, 1);
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
                Ddpca::Real tempX = leng[0] / diviNumb_es[0] * ti_real;
                Ddpca::Real heig_ti = leng[1] * (1.0- (Ddpca::Real)ti_real / diviNumb_es[0] * lengFact);
                Ddpca::Real widt_ti = leng[2] * (1.0- (Ddpca::Real)ti_real / diviNumb_es[0] * lengFact);
                for(Ddpca::I64 tj = 0; tj <= diviReal[1]; tj += diviReal[1]){
                    Ddpca::I64 tj_real = tb_1 * diviReal[1] + tj;
                    Ddpca::Real tempY = - heig_ti / 2.0 + heig_ti / diviNumb_es[1] * (Ddpca::Real)tj_real;
                    for(Ddpca::I64 tk =  0; tk <= diviReal[2]; tk ++){
                        Ddpca::I64 tk_real = tb_2 * diviReal[2] + tk;
                        Ddpca::Real tempZ = - widt_ti / 2.0 + widt_ti / diviNumb_es[2] * (Ddpca::Real)tk_real;
                        tempCoor = {tempX, tempY, tempZ};
                        CoordinateRotate(tempCoor, 1);
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
                Ddpca::Real tempX = leng[0] / diviNumb_es[0] * ti_real;
                Ddpca::Real heig_ti = leng[1] * (1.0- (Ddpca::Real)ti_real / diviNumb_es[0] * lengFact);
                Ddpca::Real widt_ti = leng[2] * (1.0- (Ddpca::Real)ti_real / diviNumb_es[0] * lengFact);
                for(Ddpca::I64 tj = 0; tj <= diviReal[1]; ++ tj){
                    Ddpca::I64 tj_real = tb_1 * diviReal[1] + tj;
                    Ddpca::Real tempY = - heig_ti / 2.0 + heig_ti / diviNumb_es[1] * (Ddpca::Real)tj_real;
                    for(Ddpca::I64 tk =  0; tk <= diviReal[2]; tk += diviReal[2]){
                        Ddpca::I64 tk_real = tb_2 * diviReal[2] + tk;
                        Ddpca::Real tempZ = - widt_ti / 2.0 + widt_ti / diviNumb_es[2] * (Ddpca::Real)tk_real;
                        tempCoor = {tempX, tempY, tempZ};
                        CoordinateRotate(tempCoor, 1);
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
            (domaNumb[1] - 1) * domaNumb[0] * domaNumb[2], 
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
                    if(tg_1 <= domaNumb[1] - 2){
                        Ddpca::I64 tg_s = tg_m + domaNumb[2];
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
            Ddpca::Real tempXi, tempEta;
            Ddpca::Coordinate tempCoor;
            for(Ddpca::I64 ti = 0; ti < inmaSize; ++ ti){
                tempXi = 0.0;
                tempEta = 0.0;
                for(Ddpca::I64 tj = 0; tj < 4; ++ tj){
                    node_tj = interfaces[ts].masterSegments[ti][tj];
                    auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                    tempCoor = iterNoco->second;
                    CoordinateRotate(tempCoor, -1.0);
                    if(ts < xyzN[0]){
                        tempXi += tempCoor[1];
                        tempEta += tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1]){
                        tempXi += tempCoor[0];
                        tempEta += tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
                        tempXi += tempCoor[0];
                        tempEta += tempCoor[1];
                    }
                }
                masterLocal[0][ti] = tempXi / 4.0;
                masterLocal[1][ti] = tempEta / 4.0;
            }
            if(ts < xyzN[0]){
                interfaces[ts].BucketSort(masterLocal, {
                    diviNumb[1] / domaNumb[1] * (1 << globLeve), 
                    diviNumb[2] / domaNumb[2] * (1 << globLeve)});
            }
            else if(ts < xyzN[0] + xyzN[1]){
                interfaces[ts].BucketSort(masterLocal, {
                    diviNumb[0] / domaNumb[0] * (1 << globLeve), 
                    diviNumb[2] / domaNumb[2] * (1 << globLeve)});
            }
            else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
                interfaces[ts].BucketSort(masterLocal, {
                    diviNumb[0] / domaNumb[0] * (1 << globLeve), 
                    diviNumb[1] / domaNumb[1] * (1 << globLeve)});
            }
            //
            Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
            std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
            for(Ddpca::I64 ti = 0; ti < inslSize; ++ ti){
                for(Ddpca::I64 tj = 0; tj < 4; ++ tj){
                    node_tj = interfaces[ts].slaveSegments[ti][tj];
                    auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                    tempCoor = iterNoco->second;
                    CoordinateRotate(tempCoor, -1.0);
                    if(ts < xyzN[0]){
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[1];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1]){
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[0];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[0];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[1];
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
    std::string directoryPath = "./TestBeam_";
    std::filesystem::create_directory(directoryPath);
    //
    Beam beam;
    beam.directoryPath = directoryPath;
    beam.Test();
    //
	Ddpca::Finalize();
	return 1;
}
