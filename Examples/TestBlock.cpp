#include "../Mesh/BilinearQuadrature.hpp"
#include "../Decomposition/MultiDomains.hpp"
#include "../Contact/CSQuadSegment.hpp"
#include "../Contact/CurvedSurface.hpp"

#include <filesystem>

//example patch test
class Block : public Ddpca::MultiDomains {

public:

	/*********************************************************************************************/
    std::string directoryPath;
    std::array<Ddpca::Real,3> leng = {0.03, 0.025, 0.02};//edge length of each block
	std::array<Ddpca::I64,3> diviNumb = {6, 6, 6};//number of elements along the edge of each block
	Ddpca::I64 globLeve = 4;//multigrid level, total levels of global mesh refinement
	//to avoid cross corner problem, each block has one complete upper block
	//and one complete lower block
	std::array<Ddpca::Real, 3> deltZlen;//height of the complete upper/lower block
	std::array<Ddpca::Real, 3> uppeZlen;//z coordinate of the complete upper block
	std::array<Ddpca::I64,3> domaNumb = {2, 2, 2};//number of decomposed domains
	Ddpca::DenseMatrix loadPres = Ddpca::DenseMatrix(3, 1, {0.0, 0.0, -1.0E7});

	//0 - upper surface of bottom block, 1 - bottom surface of middle block
	//2 - upper surface of middle block, 3 - bottom surface of upper block
	//4 - upper surface of upper block
	std::vector<Ddpca::CurvedSurface> bmupSurf;
	std::vector<Ddpca::CurvedSurface> bmupSurf_t;//rigid translation of bmupSurf
    // intermediate variables: originally can be treated as local variable in function, 
    // but now has to be treated as global variable for the sake of threadManager parallelism.
    std::vector<std::array<Ddpca::CurvedSurface, 6>> blocSurf;
    std::array<Ddpca::I64, 5> xyzN;

public:
	/*********************************************************************************************/

    void GenerateSurfaces(){
        //
        bmupSurf.resize(5);
        bmupSurf_t.resize(6);
        //
        Ddpca::I64 totaDivi = diviNumb[0] * (1 << (globLeve));
        bmupSurf[0].Resize(totaDivi + 1, totaDivi + 1);
        Ddpca::Coordinate tempCoor;
        for(Ddpca::I64 ti = 0; ti <= totaDivi; ti ++){
            for(Ddpca::I64 tj = 0; tj <= totaDivi; tj ++){
                tempCoor = Ddpca::Coordinate(- leng[0] / 2.0 + leng[0] / totaDivi * (Ddpca::Real)ti,
                    - leng[0] / 2.0 + leng[0] / totaDivi * (Ddpca::Real)tj, leng[0]);
                bmupSurf[0].Insert(ti, tj, tempCoor);
            }
        }
        bmupSurf_t[0] = bmupSurf[0];
        Ddpca::Coordinate tempTran(0.0, 0.0, - uppeZlen[0]);
        bmupSurf_t[0].RigidRotationTranslation(Ddpca::DenseMatrix::CreateIdentity(3), tempTran);
        bmupSurf_t[1] = bmupSurf[0];
        tempTran = Ddpca::Coordinate(0.0, 0.0, - deltZlen[0]);
        bmupSurf_t[1].RigidRotationTranslation(Ddpca::DenseMatrix::CreateIdentity(3), tempTran);
        //
        totaDivi = diviNumb[1] * (1 << (globLeve));
        bmupSurf[1].Resize(totaDivi + 1, totaDivi + 1);
        bmupSurf[2].Resize(totaDivi + 1, totaDivi + 1);
        for(Ddpca::I64 ti = 0; ti <= totaDivi; ti ++){
            for(Ddpca::I64 tj = 0; tj <= totaDivi; tj ++){
                //
                tempCoor = Ddpca::Coordinate(- leng[1] / 2.0 + leng[1] / totaDivi * (Ddpca::Real)ti,
                    - leng[1] / 2.0 + leng[1] / totaDivi * (Ddpca::Real)tj, leng[0]);
                bmupSurf[1].Insert(ti, tj, tempCoor);
                //
                tempCoor[2] += leng[1];
                bmupSurf[2].Insert(ti, tj, tempCoor);
            }
        }
        bmupSurf_t[2] = bmupSurf[2];
        tempTran = Ddpca::Coordinate(0.0, 0.0, - uppeZlen[1]);
        bmupSurf_t[2].RigidRotationTranslation(Ddpca::DenseMatrix::CreateIdentity(3), tempTran);
        bmupSurf_t[3] = bmupSurf[2];
        tempTran = Ddpca::Coordinate(0.0, 0.0, - deltZlen[1]);
        bmupSurf_t[3].RigidRotationTranslation(Ddpca::DenseMatrix::CreateIdentity(3), tempTran);
        //
        totaDivi = diviNumb[2] * (1 << (globLeve));
        bmupSurf[3].Resize(totaDivi + 1, totaDivi + 1);
        bmupSurf[4].Resize(totaDivi + 1, totaDivi + 1);
        for(Ddpca::I64 ti = 0; ti <= totaDivi; ti ++){
            for(Ddpca::I64 tj = 0; tj <= totaDivi; tj ++){
                //
                tempCoor = Ddpca::Coordinate(- leng[2] / 2.0 + leng[2] / totaDivi * (Ddpca::Real)ti,
                    - leng[2] / 2.0 + leng[2] / totaDivi * (Ddpca::Real)tj, leng[0] + leng[1]);
                bmupSurf[3].Insert(ti, tj, tempCoor);
                //
                tempCoor[2] += leng[2];
                bmupSurf[4].Insert(ti, tj, tempCoor);
            }
        }
        bmupSurf_t[4] = bmupSurf[4];
        tempTran = Ddpca::Coordinate(0.0, 0.0, - uppeZlen[2]);
        bmupSurf_t[4].RigidRotationTranslation(Ddpca::DenseMatrix::CreateIdentity(3), tempTran);
        bmupSurf_t[5] = bmupSurf[4];
        tempTran = Ddpca::Coordinate(0.0, 0.0, - deltZlen[2]);
        bmupSurf_t[5].RigidRotationTranslation(Ddpca::DenseMatrix::CreateIdentity(3), tempTran);
        //
        std::array<Ddpca::I64, 3> diviNumb_es = {
            diviNumb[0] * (1 << (globLeve)), 
            diviNumb[1] * (1 << (globLeve)), 
            diviNumb[2] * (1 << (globLeve))};
        blocSurf.resize(3 * domaNumb[0] * domaNumb[1] * domaNumb[2]);
        Ddpca::I64 blsuSize = blocSurf.size();
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            //
            if(tg >= blsuSize){
                return;
            }
            Ddpca::I64 tb = tg / (domaNumb[0] * domaNumb[1] * domaNumb[2]);
            Ddpca::I64 tg_b = tg % (domaNumb[0] * domaNumb[1] * domaNumb[2]);
            Ddpca::I64 tg_0 = tg_b / (domaNumb[1] * domaNumb[2]);
            Ddpca::I64 tg_1 = (tg_b % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
            Ddpca::I64 tg_2 = (tg_b % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
            std::array<Ddpca::I64,3> diviReal = {diviNumb_es[tb] / domaNumb[0], 
                diviNumb_es[tb] / domaNumb[1], diviNumb_es[tb] / domaNumb[2]};
            blocSurf[tg][0].Resize(diviReal[1] + 1, diviReal[2] + 1);
            blocSurf[tg][1].Resize(diviReal[1] + 1, diviReal[2] + 1);
            blocSurf[tg][2].Resize(diviReal[2] + 1, diviReal[0] + 1);
            blocSurf[tg][3].Resize(diviReal[2] + 1, diviReal[0] + 1);
            blocSurf[tg][4].Resize(diviReal[0] + 1, diviReal[1] + 1);
            blocSurf[tg][5].Resize(diviReal[0] + 1, diviReal[1] + 1);
            Ddpca::Real totaZlen = uppeZlen[tb] - deltZlen[tb];
            //nodes
            Ddpca::Coordinate tempCoor;
            for(Ddpca::I64 ti = 0; ti <= diviReal[0]; ti ++){
                Ddpca::I64 ti_real = tg_0 * diviReal[0] + ti;
                Ddpca::Real xCoo = - leng[tb] / 2.0 + leng[tb] / diviNumb_es[tb] * (Ddpca::Real)ti_real;
                for(Ddpca::I64 tj = 0; tj <= diviReal[1]; tj ++){
                    Ddpca::I64 tj_real = tg_1 * diviReal[1] + tj;
                    Ddpca::Real yCoo = - leng[tb] / 2.0 + leng[tb] / diviNumb_es[tb] * (Ddpca::Real)tj_real;
                    for(Ddpca::I64 tk = 0; tk <= diviReal[2]; tk ++){
                        Ddpca::I64 tk_real = tg_2 * diviReal[2] + tk;
                        Ddpca::Real zCoo = deltZlen[tb] + totaZlen / diviNumb_es[tb] * (Ddpca::Real)tk_real;
                        if(tb == 1){
                            zCoo += leng[0];
                        }
                        else if(tb == 2){
                            zCoo += leng[0] + leng[1];
                        }
                        tempCoor = Ddpca::Coordinate(xCoo, yCoo, zCoo);
                        if(ti == 0){
                            blocSurf[tg][0].Insert(tj, tk, tempCoor);
                        }
                        if(ti == diviReal[0]){
                            blocSurf[tg][1].Insert(tj, tk, tempCoor);
                        }
                        if(tj == 0){
                            blocSurf[tg][2].Insert(tk, ti, tempCoor);
                        }
                        if(tj == diviReal[1]){
                            blocSurf[tg][3].Insert(tk, ti, tempCoor);
                        }
                        if(tk == 0){
                            blocSurf[tg][4].Insert(ti, tj, tempCoor);
                        }
                        if(tk == diviReal[2]){
                            blocSurf[tg][5].Insert(ti, tj, tempCoor);
                        }
                    }
                }
            }
        });
    }

    void GenerateMeshes(){
        //
        for(Ddpca::I64 tb = 0; tb < 3; tb ++){
            deltZlen[tb] = leng[tb] / (diviNumb[tb] * (1 << globLeve));
            uppeZlen[tb] = leng[tb] - deltZlen[tb];
        }
        //
        domains.resize(3 * domaNumb[0] * domaNumb[1] * domaNumb[2] + 3 * 2);
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.threadsPerDomain, domains.size(), 
            Ddpca::threadManager.domainS2S, Ddpca::threadManager.domainS2M);
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::I64 tb = tg / (domaNumb[0] * domaNumb[1] * domaNumb[2]);
            if(tb <= 2){
                Ddpca::Real totaZlen = uppeZlen[tb] - deltZlen[tb];
                //
                Ddpca::I64 tg_b = tg % (domaNumb[0] * domaNumb[1] * domaNumb[2]);
                Ddpca::I64 tg_0 = tg_b / (domaNumb[1] * domaNumb[2]);
                Ddpca::I64 tg_1 = (tg_b % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
                Ddpca::I64 tg_2 = (tg_b % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
                //nodes
                std::array<Ddpca::I64, 3> diviReal = {
                    diviNumb[0] / domaNumb[0], diviNumb[1] / domaNumb[1], diviNumb[2] / domaNumb[2]};
                std::vector<std::vector<std::vector<Ddpca::I64>>> tempNode(
                    diviReal[0] + 1, 
                    std::vector<std::vector<Ddpca::I64>>(
                        diviReal[1] + 1, 
                        std::vector<Ddpca::I64>(diviReal[2] + 1, 0)));
                Ddpca::Coordinate tempCoor;
                for(Ddpca::I64 ti = 0; ti <= diviReal[0]; ++ ti){
                    Ddpca::I64 ti_real = tg_0 * diviReal[0] + ti;
                    tempCoor[0] = - leng[tb] / 2.0 + leng[tb] / diviNumb[tb] * (Ddpca::Real)ti_real;
                    for(Ddpca::I64 tj = 0; tj <= diviReal[1]; ++ tj){
                        Ddpca::I64 tj_real = tg_1 * diviReal[1] + tj;
                        tempCoor[1] = - leng[tb] / 2.0 + leng[tb] / diviNumb[tb] * (Ddpca::Real)tj_real;
                        for(Ddpca::I64 tk = 0; tk <= diviReal[2]; ++ tk){
                            Ddpca::I64 tk_real = tg_2 * diviReal[2] + tk;
                            tempCoor[2] = deltZlen[tb] + totaZlen / diviNumb[tb] * (Ddpca::Real)tk_real;
                            if(tb == 1){
                                tempCoor[2] += leng[0];
                            }
                            else if(tb == 2){
                                tempCoor[2] += leng[0] + leng[1];
                            }
                            tempNode[ti][tj][tk] = tgMesh.TryAddNode(tempCoor);
                        }
                    }
                }
                //elements
                Ddpca::OctreeElement tempElem;
                for(Ddpca::I64 ti = 0; ti < diviReal[0]; ++ ti){
                    for(Ddpca::I64 tj = 0; tj < diviReal[1]; ++ tj){
                        for(Ddpca::I64 tk = 0; tk < diviReal[2]; ++ tk){
                            tempElem.parent = -1;
                            tempElem.cornerNodes[0] = tempNode[ti][tj][tk];
                            tempElem.cornerNodes[1] = tempNode[ti + 1][tj][tk];
                            tempElem.cornerNodes[2] = tempNode[ti + 1][tj + 1][tk];
                            tempElem.cornerNodes[3] = tempNode[ti][tj + 1][tk];
                            tempElem.cornerNodes[4] = tempNode[ti][tj][tk + 1];
                            tempElem.cornerNodes[5] = tempNode[ti + 1][tj][tk + 1];
                            tempElem.cornerNodes[6] = tempNode[ti + 1][tj + 1][tk + 1];
                            tempElem.cornerNodes[7] = tempNode[ti][tj + 1][tk + 1];
                            tempElem.level = 0;
                            tempElem.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
                            tempElem.children.resize(0);
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
                tgMesh.OutputMesh(directoryPath, tg);
            }
            else{
                Ddpca::I64 tg_temp = tg - 3 * (domaNumb[0] * domaNumb[1] * domaNumb[2]);
                tb = tg_temp / 2;
                Ddpca::I64 tg_b = tg_temp % 2;
                //nodes
                Ddpca::I64 diviReal = 1;
                std::vector<std::vector<std::vector<Ddpca::I64>>> tempNode(
                    diviNumb[0] + 1, 
                    std::vector<std::vector<Ddpca::I64>>(
                        diviNumb[1] + 1, 
                        std::vector<Ddpca::I64>(diviReal + 1, 0)));
                Ddpca::Coordinate tempCoor;
                for(Ddpca::I64 ti = 0; ti <= diviNumb[0]; ti ++){
                    tempCoor[0] = - leng[tb] / 2.0 + leng[tb] / diviNumb[tb] * (Ddpca::Real)ti;
                    for(Ddpca::I64 tj = 0; tj <= diviNumb[1]; tj ++){
                        tempCoor[1] = - leng[tb] / 2.0 + leng[tb] / diviNumb[tb] * (Ddpca::Real)tj;
                        for(Ddpca::I64 tk = 0; tk <= diviReal; tk ++){
                            tempCoor[2] = deltZlen[tb] / diviReal * (Ddpca::Real)tk;
                            if(tg_b == 1){
                                tempCoor[2] += uppeZlen[tb];
                            }
                            if(tb == 1){
                                tempCoor[2] += leng[0];
                            }
                            else if(tb == 2){
                                tempCoor[2] += leng[0] + leng[1];
                            }
                            tempNode[ti][tj][tk] = tgMesh.TryAddNode(tempCoor);
                        }
                    }
                }
                //elements
                Ddpca::OctreeElement tempElem;
                for(Ddpca::I64 ti = 0; ti < diviNumb[0]; ti ++){
                    for(Ddpca::I64 tj = 0; tj < diviNumb[1]; tj ++){
                        for(Ddpca::I64 tk = 0; tk < diviReal; tk ++){
                            tempElem.parent = -1;
                            tempElem.cornerNodes[0] = tempNode[ti][tj][tk];
                            tempElem.cornerNodes[1] = tempNode[ti + 1][tj][tk];
                            tempElem.cornerNodes[2] = tempNode[ti + 1][tj + 1][tk];
                            tempElem.cornerNodes[3] = tempNode[ti][tj + 1][tk];
                            tempElem.cornerNodes[4] = tempNode[ti][tj][tk + 1];
                            tempElem.cornerNodes[5] = tempNode[ti + 1][tj][tk + 1];
                            tempElem.cornerNodes[6] = tempNode[ti + 1][tj + 1][tk + 1];
                            tempElem.cornerNodes[7] = tempNode[ti][tj + 1][tk + 1];
                            tempElem.level = 0;
                            tempElem.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
                            tempElem.children.resize(0);
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
                            Ddpca::OctreeElement::REFINEMENT_XI_ETA;
                    }
                    curvInte.clear();
                    tgMesh.Refine(elementsToSplit, subElements, curvInte);
                }
                //
                tgMesh.OutputMesh(directoryPath, tg);
            }
        });
        GenerateSurfaces();
    }

    void SubLoad_0(Ddpca::I64 tb, Ddpca::I64 tg){
        //
        Ddpca::Real zCoo = (tb == 0) ? leng[0] : leng[0] + leng[1];
        std::array<std::array<Ddpca::Coordinate,4>,4> pseuElem_0;
        pseuElem_0[0][0] = Ddpca::Coordinate(- leng[0 + tb] / 2.0, - leng[0 + tb] / 2.0, zCoo);
        pseuElem_0[0][1] = Ddpca::Coordinate(- leng[0 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[0][2] = Ddpca::Coordinate(leng[0 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[0][3] = Ddpca::Coordinate(leng[0 + tb] / 2.0, - leng[0 + tb] / 2.0, zCoo);
        pseuElem_0[1][0] = Ddpca::Coordinate(- leng[0 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[1][1] = Ddpca::Coordinate(- leng[0 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[1][2] = Ddpca::Coordinate(- leng[1 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[1][3] = Ddpca::Coordinate(- leng[1 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[2][0] = Ddpca::Coordinate(leng[1 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[2][1] = Ddpca::Coordinate(leng[1 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[2][2] = Ddpca::Coordinate(leng[0 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[2][3] = Ddpca::Coordinate(leng[0 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[3][0] = Ddpca::Coordinate(- leng[0 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo);
        pseuElem_0[3][1] = Ddpca::Coordinate(- leng[0 + tb] / 2.0, leng[0 + tb] / 2.0, zCoo);
        pseuElem_0[3][2] = Ddpca::Coordinate(leng[0 + tb] / 2.0, leng[0 + tb] / 2.0, zCoo);
        pseuElem_0[3][3] = Ddpca::Coordinate(leng[0 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo);
        Ddpca::I64 refiNumb = 1;
        std::vector<std::array<Ddpca::Coordinate,4>> pseuElem(4 * refiNumb * refiNumb);
        Ddpca::I64 pselSize = pseuElem.size();
        for(Ddpca::I64 ti = 0; ti < 4; ti ++){
            for(Ddpca::I64 tj = 0; tj < refiNumb; tj ++){
                for(Ddpca::I64 tk = 0; tk < refiNumb; tk ++){
                    Ddpca::I64 id_tijk = ti * refiNumb * refiNumb + tj * refiNumb + tk;
                    pseuElem[id_tijk][0][0] = pseuElem_0[ti][0][0] 
                        + (pseuElem_0[ti][3][0] - pseuElem_0[ti][0][0]) / refiNumb * (Ddpca::Real)tj;
                    pseuElem[id_tijk][0][1] = pseuElem_0[ti][0][1] 
                        + (pseuElem_0[ti][1][1] - pseuElem_0[ti][0][1]) / refiNumb * (Ddpca::Real)tk;
                    pseuElem[id_tijk][0][2] = zCoo;
                    pseuElem[id_tijk][1][0] = pseuElem_0[ti][0][0] 
                        + (pseuElem_0[ti][3][0] - pseuElem_0[ti][0][0]) / refiNumb * (Ddpca::Real)tj;
                    pseuElem[id_tijk][1][1] = pseuElem_0[ti][0][1] 
                        + (pseuElem_0[ti][1][1] - pseuElem_0[ti][0][1]) / refiNumb * (Ddpca::Real)(tk + 1);
                    pseuElem[id_tijk][1][2] = zCoo;
                    pseuElem[id_tijk][2][0] = pseuElem_0[ti][0][0] 
                        + (pseuElem_0[ti][3][0] - pseuElem_0[ti][0][0]) / refiNumb * (Ddpca::Real)(tj + 1);
                    pseuElem[id_tijk][2][1] = pseuElem_0[ti][0][1] 
                        + (pseuElem_0[ti][1][1] - pseuElem_0[ti][0][1]) / refiNumb * (Ddpca::Real)(tk + 1);
                    pseuElem[id_tijk][2][2] = zCoo;
                    pseuElem[id_tijk][3][0] = pseuElem_0[ti][0][0] 
                        + (pseuElem_0[ti][3][0] - pseuElem_0[ti][0][0]) / refiNumb * (Ddpca::Real)(tj + 1);
                    pseuElem[id_tijk][3][1] = pseuElem_0[ti][0][1] 
                        + (pseuElem_0[ti][1][1] - pseuElem_0[ti][0][1]) / refiNumb * (Ddpca::Real)tk;
                    pseuElem[id_tijk][3][2] = zCoo;
                }
            }
        }
        //
        const Ddpca::BilinearQuadrature<2>& biliQuad = Ddpca::GetBilinearQuadrature();
        Ddpca::CurvedSurface& tgSurf = bmupSurf[(tb == 0) ? 0 : 2];
        const Ddpca::Mesh& tgMesh = domains[tg].mesh;
        Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
        std::array<Ddpca::I64,4> tempNode;
        std::array<std::array<Ddpca::Real, 3>, 4> mastCorn, slavCorn;
        std::vector<std::array<Ddpca::Real,2>> listXiet; 
        std::vector<Ddpca::Real> listWeig;
        Ddpca::CSQuadSegment tempCsea;
        std::array<Ddpca::Real,4> N_e;
        std::vector<Ddpca::Coordinate> elementCoordinates(4);
        Ddpca::DenseMatrix N_e_0(3,12), F_e(12,1), tkF_e(12,1);
        tgSurf.Initialize();
        while(tgSurf.Increment(tgMesh)){
            tempNode = tgSurf.currentFace;
            Ddpca::Real minX = 1.0E20, maxX = -1.0E20, minY = 1.0E20, maxY = -1.0E20;
            for(Ddpca::I64 ti = 0; ti < 4; ti ++){
                auto iterNoco = tgMesh.node2Coordinate.find(tempNode[ti]);
                mastCorn[ti] = {(iterNoco->second)[0], (iterNoco->second)[1], (iterNoco->second)[2]};
                elementCoordinates[ti] = Ddpca::Coordinate(
                    mastCorn[ti][0], mastCorn[ti][1], mastCorn[ti][2]);
                minX = std::min(minX, mastCorn[ti][0]);
                maxX = std::max(maxX, mastCorn[ti][0]);
                minY = std::min(minY, mastCorn[ti][1]);
                maxY = std::max(maxY, mastCorn[ti][1]);
            }
            for(Ddpca::I64 ti = 0; ti < pselSize; ti ++){
                Ddpca::Real minX_1 = 1.0E20, maxX_1 = -1.0E20, minY_1 = 1.0E20, maxY_1 = -1.0E20;
                for(Ddpca::I64 tj = 0; tj < 4; ++ tj){
                    slavCorn[tj] = {(pseuElem[ti][tj][0]), (pseuElem[ti][tj][1]), (pseuElem[ti][tj][2])};
                    minX_1 = std::min(minX_1, slavCorn[tj][0]);
                    maxX_1 = std::max(maxX_1, slavCorn[tj][0]);
                    minY_1 = std::min(minY_1, slavCorn[tj][1]);
                    maxY_1 = std::max(maxY_1, slavCorn[tj][1]);
                }
                if(minX_1 + Ddpca::Coordinate::coorError >= maxX 
                    || maxX_1 - Ddpca::Coordinate::coorError <= minX 
                    || minY_1 + Ddpca::Coordinate::coorError >= maxY 
                    || maxY_1 - Ddpca::Coordinate::coorError <= minY) continue;
                if(minX_1 - Ddpca::Coordinate::coorError <= minX 
                    && maxX <= maxX_1 + Ddpca::Coordinate::coorError 
                    && minY_1 - Ddpca::Coordinate::coorError <= minY 
                    && maxY <= maxY_1 + Ddpca::Coordinate::coorError){
                    F_e.Fill(0.0);
                    for(Ddpca::I64 tk = 0; tk < biliQuad.numbGaussPoints; tk ++){
                        Ddpca::Real tkJacobian = Ddpca::BilinearQuadratureJacobian<2>(
                            biliQuad.gaussPoints[tk], elementCoordinates);
                        // col-major
                        N_e_0.Fill(
                            biliQuad.shapeFunctions(tk,0), 0.0, 0.0, 
                            0.0, biliQuad.shapeFunctions(tk,0), 0.0, 
                            0.0, 0.0, biliQuad.shapeFunctions(tk,0), 
                            biliQuad.shapeFunctions(tk,1), 0.0, 0.0, 
                            0.0, biliQuad.shapeFunctions(tk,1), 0.0, 
                            0.0, 0.0, biliQuad.shapeFunctions(tk,1), 
                            biliQuad.shapeFunctions(tk,2), 0.0, 0.0, 
                            0.0, biliQuad.shapeFunctions(tk,2), 0.0, 
                            0.0, 0.0, biliQuad.shapeFunctions(tk,2), 
                            biliQuad.shapeFunctions(tk,3), 0.0, 0.0, 
                            0.0, biliQuad.shapeFunctions(tk,3), 0.0, 
                            0.0, 0.0, biliQuad.shapeFunctions(tk,3));
                        Ddpca::GEMTM(N_e_0, loadPres, tkF_e);
                        Ddpca::AXPY(biliQuad.weights[tk] * tkJacobian, tkF_e.data, F_e.data);
                    }
                }
                else{
                    listXiet.clear();
                    listWeig.clear();
                    tempCsea.integralPoints.clear();
                    tempCsea.masterCorners = mastCorn;
                    tempCsea.slaveCorners = slavCorn;
                    tempCsea.SegmentIntersect(listXiet, listWeig);
                    Ddpca::I64 lixiSize = listXiet.size();
                    F_e.Fill(0.0);
                    for(Ddpca::I64 tj = 0; tj < lixiSize; tj ++){
                        for(Ddpca::I64 tk = 0; tk < 4; tk ++){
                            N_e[tk] = (1.0 + biliQuad.cornerNodes[tk][0] * listXiet[tj][0]) 
                                * (1.0 + biliQuad.cornerNodes[tk][1] * listXiet[tj][1]) / 4.0;
                        }
                        Ddpca::Real tjJacobian = Ddpca::BilinearQuadratureJacobian<2>(
                            biliQuad.gaussPoints[tj], elementCoordinates);
                        // col-major
                        N_e_0.Fill(
                            N_e[0], 0.0, 0.0, 
                            0.0, N_e[0], 0.0, 
                            0.0, 0.0, N_e[0], 
                            N_e[1], 0.0, 0.0, 
                            0.0, N_e[1], 0.0, 
                            0.0, 0.0, N_e[1], 
                            N_e[2], 0.0, 0.0, 
                            0.0, N_e[2], 0.0, 
                            0.0, 0.0, N_e[2], 
                            N_e[3], 0.0, 0.0, 
                            0.0, N_e[3], 0.0, 
                            0.0, 0.0, N_e[3]);
                        Ddpca::GEMTM(N_e_0, loadPres, tkF_e);
                        Ddpca::AXPY(listWeig[tj] * tjJacobian, tkF_e.data, F_e.data);
                    }
                }
                for(Ddpca::I64 tk = 0; tk < 4; tk ++){
                    for(Ddpca::I64 tm = 0; tm < 3; tm ++){
                        Ddpca::I64 free_tm = 3 * tempNode[tk] + tm;
                        tgBoundary.LoadAccumulate(free_tm, F_e(3 * tk + tm, 0));
                    }
                }
            }
        }
    }

    void SubLoad_1(Ddpca::I64 tg){
        Ddpca::CurvedSurface& tgSurf = bmupSurf[4];
        const Ddpca::Mesh& tgMesh = domains[tg].mesh;
        Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
        const Ddpca::BilinearQuadrature<2>& biliQuad = Ddpca::GetBilinearQuadrature();
        std::array<Ddpca::I64,4> tempNode;
        std::vector<Ddpca::Coordinate> elementCoordinates(4);
        Ddpca::DenseMatrix F_e(12,1), tkF_e(12,1), N_e_0(3,12);
        tgSurf.Initialize();
        while(tgSurf.Increment(tgMesh)){
            tempNode = tgSurf.currentFace;
            for(Ddpca::I64 tk = 0; tk < 4; tk ++){
                auto iterNoco = tgMesh.node2Coordinate.find(tempNode[tk]);
                elementCoordinates[tk] = Ddpca::Coordinate((iterNoco->second)[0], 
                    (iterNoco->second)[1], (iterNoco->second)[2]);
            }
            F_e.Fill(0.0);
            for(Ddpca::I64 tk = 0; tk < biliQuad.numbGaussPoints; tk ++){
                Ddpca::Real tkJacobian = Ddpca::BilinearQuadratureJacobian<2>(
                    biliQuad.gaussPoints[tk], elementCoordinates);
                // col-major
                N_e_0.Fill(
                    biliQuad.shapeFunctions(tk,0), 0.0, 0.0, 
                    0.0, biliQuad.shapeFunctions(tk,0), 0.0, 
                    0.0, 0.0, biliQuad.shapeFunctions(tk,0), 
                    biliQuad.shapeFunctions(tk,1), 0.0, 0.0, 
                    0.0, biliQuad.shapeFunctions(tk,1), 0.0, 
                    0.0, 0.0, biliQuad.shapeFunctions(tk,1), 
                    biliQuad.shapeFunctions(tk,2), 0.0, 0.0, 
                    0.0, biliQuad.shapeFunctions(tk,2), 0.0, 
                    0.0, 0.0, biliQuad.shapeFunctions(tk,2), 
                    biliQuad.shapeFunctions(tk,3), 0.0, 0.0, 
                    0.0, biliQuad.shapeFunctions(tk,3), 0.0, 
                    0.0, 0.0, biliQuad.shapeFunctions(tk,3));
                Ddpca::GEMTM(N_e_0, loadPres, tkF_e);
                Ddpca::AXPY(biliQuad.weights[tk] * tkJacobian, tkF_e.data, F_e.data);
            }
            for(Ddpca::I64 tk = 0; tk < 4; tk ++){
                for(Ddpca::I64 tm = 0; tm < 3; tm ++){
                    Ddpca::I64 free_tk = 3 * tempNode[tk] + tm;
                    tgBoundary.LoadAccumulate(free_tk, F_e(3 * tk + tm, 0));
                }
            }
        }
    }

    void ConstraintLoad(){
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tg){
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
            Ddpca::I64 tb = tg / (domaNumb[0] * domaNumb[1] * domaNumb[2]);
            // displacement constraint
            if(tb <= 2){
                Ddpca::Real zlenMini = deltZlen[tb];
                Ddpca::Real zlenMaxi = uppeZlen[tb];
                if(tb == 1){
                    zlenMini += leng[0];
                    zlenMaxi += leng[0];
                }
                else if(tb == 2){
                    zlenMini += leng[0] + leng[1];
                    zlenMaxi += leng[0] + leng[1];
                }
                for(const auto& iterNoco : tgMesh.node2Coordinate){
                    if((iterNoco.second)[2] <= zlenMini + 1.0E-12 
                        || (iterNoco.second)[2] >= zlenMaxi - 1.0E-12){
                        continue;
                    }
                    if((iterNoco.second)[0] <= - leng[tb] / 2.0 + 1.0E-12){
                        tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    }
                    if((iterNoco.second)[1] <= - leng[tb] / 2.0 + 1.0E-12){
                        tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 1, 0.0);
                    }
                }
            }
            else{
                Ddpca::I64 tg_temp = tg - 3 * (domaNumb[0] * domaNumb[1] * domaNumb[2]);
                tb = tg_temp / 2;
                Ddpca::I64 tg_b = tg_temp % 2;
                //displacement constraint must at first
                for(const auto& iterNoco : tgMesh.node2Coordinate){
                    if((iterNoco.second)[2] <= 1.0E-10){//displacement constraint
                        tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                    }
                    if((iterNoco.second)[0] <= - leng[tb] / 2.0 + 1.0E-12){
                        tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    }
                    if((iterNoco.second)[1] <= - leng[tb] / 2.0 + 1.0E-12){
                        tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 1, 0.0);
                    }
                }
                //load must at next
                if(tg_b == 1){
                    if(tb == 2){
                        SubLoad_1(tg);
                    }
                    else{
                        SubLoad_0(tb, tg);
                    }
                }
            }
        });
    }

    void GenerateInterfaces(){
        //
        xyzN = {
            3 * (domaNumb[0] - 1) * domaNumb[1] * domaNumb[2], 
            3 * (domaNumb[1] - 1) * domaNumb[0] * domaNumb[2], 
            3 * (domaNumb[2] - 1) * domaNumb[0] * domaNumb[1], 
            3 * 2 * domaNumb[0] * domaNumb[1], 2};
        interfaces.resize(xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3] + xyzN[4]);
        Ddpca::I64 inteSize = interfaces.size();
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.threadsPerInterface, inteSize, 
            Ddpca::threadManager.interfaceS2S, Ddpca::threadManager.interfaceS2M);
        //
        for(Ddpca::I64 tb = 0; tb < 3; tb ++){
            for(Ddpca::I64 tg_0 = 0; tg_0 < domaNumb[0]; tg_0 ++){
                for(Ddpca::I64 tg_1 = 0; tg_1 < domaNumb[1]; tg_1 ++){
                    for(Ddpca::I64 tg_2 = 0; tg_2 < domaNumb[2]; tg_2 ++){
                        Ddpca::I64 tg_m = tb * domaNumb[0] * domaNumb[1] * domaNumb[2]
                            + tg_0 * domaNumb[1] * domaNumb[2] + tg_1 * domaNumb[2] + tg_2;
                        if(tg_0 <= domaNumb[0] - 2){
                            Ddpca::I64 tg_s = tg_m + domaNumb[1] * domaNumb[2];
                            Ddpca::I64 ts = tg_0 * domaNumb[1] * domaNumb[2] + tg_1 * domaNumb[2] + tg_2 
                                + tb * (domaNumb[0] - 1) * domaNumb[1] * domaNumb[2];
                            interfaces[ts].domainIndex = {tg_m, tg_s};
                        }
                        if(tg_1 <= domaNumb[1] - 2){
                            Ddpca::I64 tg_s = tg_m + domaNumb[2];
                            Ddpca::I64 ts = xyzN[0] 
                                + tg_1 * domaNumb[0] * domaNumb[2] + tg_0 * domaNumb[2] + tg_2 
                                + tb * (domaNumb[1] - 1) * domaNumb[0] * domaNumb[2];
                            interfaces[ts].domainIndex = {tg_m, tg_s};
                        }
                        if(tg_2 <= domaNumb[2] - 2){
                            Ddpca::I64 tg_s = tg_m + 1;
                            Ddpca::I64 ts = xyzN[0] + xyzN[1] 
                                + tg_2 * domaNumb[0] * domaNumb[1] + tg_0 * domaNumb[1] + tg_1 
                                + tb * (domaNumb[2] - 1) * domaNumb[0] * domaNumb[1];
                            interfaces[ts].domainIndex = {tg_m, tg_s};
                        }
                    }
                }
            }
        }
        Ddpca::I64 contBase = xyzN[0] + xyzN[1] + xyzN[2];
        Ddpca::I64 domaBase = 3 * domaNumb[0] * domaNumb[1] * domaNumb[2];
        for(Ddpca::I64 tb = 0; tb < 3; tb ++){
            for(Ddpca::I64 bu = 0; bu < 2; bu ++){
                for(Ddpca::I64 tg_0 = 0; tg_0 < domaNumb[0]; tg_0 ++){
                    for(Ddpca::I64 tg_1 = 0; tg_1 < domaNumb[1]; tg_1 ++){
                        Ddpca::I64 tg_m = tb * domaNumb[0] * domaNumb[1] * domaNumb[2] 
                            + tg_0 * domaNumb[1] * domaNumb[2] + tg_1 * domaNumb[2] 
                            + (bu == 0 ? 0 : domaNumb[2] - 1);
                        interfaces[contBase].domainIndex = {tg_m, domaBase + tb * 2 + bu};
                        contBase ++;
                    }
                }
            }
        }
        interfaces[contBase + 0].domainIndex = {domaBase + 1, domaBase + 2};
        interfaces[contBase + 1].domainIndex = {domaBase + 3, domaBase + 4};
        //
        for(Ddpca::I64 ts = 0; ts < inteSize; ++ ts){
            if(ts < xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]){
                interfaces[ts].frictionCoefficient = -1.0;
            }
            else{
                interfaces[ts].frictionCoefficient = 0.0;
            }
        }
        CalculatePenaltyParameter();
        //
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.interfaceS2S, [&](Ddpca::I64 ts){
            //
            Ddpca::CurvedSurface* mastSurf;
            Ddpca::CurvedSurface* slavSurf;
            if(ts < xyzN[0]){
                mastSurf = &(blocSurf[interfaces[ts].domainIndex[0]][1]);
                slavSurf = &(blocSurf[interfaces[ts].domainIndex[1]][0]);
            }
            else if(ts < xyzN[0] + xyzN[1]){
                mastSurf = &(blocSurf[interfaces[ts].domainIndex[0]][3]);
                slavSurf = &(blocSurf[interfaces[ts].domainIndex[1]][2]);
            }
            else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
                mastSurf = &(blocSurf[interfaces[ts].domainIndex[0]][5]);
                slavSurf = &(blocSurf[interfaces[ts].domainIndex[1]][4]);
            }
            else if(ts < xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]){
                long tempSunu = (ts - (xyzN[0] + xyzN[1] + xyzN[2])) 
                    / (domaNumb[0] * domaNumb[1]);
                mastSurf = new Ddpca::CurvedSurface;
                *mastSurf = (blocSurf[interfaces[ts].domainIndex[0]][4 + tempSunu % 2]);
                slavSurf = new Ddpca::CurvedSurface;
                *slavSurf = (bmupSurf_t[tempSunu]);
            }
            else if(ts < xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3] + 1){
                mastSurf = &(bmupSurf[0]);
                slavSurf = &(bmupSurf[1]);
            }
            else{
                mastSurf = &(bmupSurf[2]);
                slavSurf = &(bmupSurf[3]);
            }
            //
            mastSurf->Initialize();
            const Ddpca::Mesh& masterMesh = domains[interfaces[ts].domainIndex[0]].mesh;
            while(mastSurf->Increment(masterMesh)){
                interfaces[ts].masterSegments.emplace_back(mastSurf->currentFace);
            }
            slavSurf->Initialize();
            const Ddpca::Mesh& slaveMesh = domains[interfaces[ts].domainIndex[1]].mesh;
            while(slavSurf->Increment(slaveMesh)){
                interfaces[ts].slaveSegments.emplace_back(slavSurf->currentFace);
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
                    if(ts < xyzN[0]){
                        tempXi += tempCoor[1];
                        tempEta += tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1]){
                        tempXi += tempCoor[0];
                        tempEta += tempCoor[2];
                    }
                    else{
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
            else if(ts < xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]){
                interfaces[ts].BucketSort(masterLocal, {
                    diviNumb[0] / domaNumb[0] * (1 << globLeve), 
                    diviNumb[1] / domaNumb[1] * (1 << globLeve)});
            }
            else{
                interfaces[ts].BucketSort(masterLocal, {
                    diviNumb[0] * (1 << globLeve), 
				    diviNumb[1] * (1 << globLeve)});
            }
            //
            Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
            std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
            for(Ddpca::I64 ti = 0; ti < inslSize; ++ ti){
                for(Ddpca::I64 tj = 0; tj < 4; ++ tj){
                    node_tj = interfaces[ts].slaveSegments[ti][tj];
                    auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                    tempCoor = iterNoco->second;
                    if(ts < xyzN[0]){
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[1];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[2];
                    }
                    else if(ts < xyzN[0] + xyzN[1]){
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[0];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[2];
                    }
                    else{
                        slaveLocal[ti][tj * 2 + 0] = tempCoor[0];
                        slaveLocal[ti][tj * 2 + 1] = tempCoor[1];
                    }
                }
            }
            interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                Ddpca::threadManager.interfaceS2M[ts], 1, 1.0E12);
            interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            if(xyzN[0] + xyzN[1] + xyzN[2] <= ts && ts < xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]){
                delete mastSurf;
                delete slavSurf;
            }
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
    std::string directoryPath = "./TestBlock_";
    std::filesystem::create_directory(directoryPath);
    //
    Block block;
    block.directoryPath = directoryPath;
    block.Test();
    //
	Ddpca::Finalize();
	return 1;
}
