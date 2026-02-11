#include "DehwSurf.hpp"

#include "../Mesh/BilinearQuadrature.hpp"
#include "../Decomposition/MultiDomains.hpp"
#include "../Contact/CSQuadSegment.hpp"
#include "../Contact/CurvedSurface.hpp"

#include <filesystem>

//example double enveloping hourglass worm drive
class Dehw : public Ddpca::MultiDomains {

public:

	/*********************************************************************************************/
    std::string directoryPath;
    DehwSurf dehwSurf;
    bool isSelf; // is self locking analysis
    double tempTangPenaCoef = 0.25;
	Ddpca::Real centErro;//center distance error
	std::array<Ddpca::Real,2> analAngl;//rotating angle of worm and worm wheel
	std::array<Ddpca::Real,3> distCrit;//critical gap used for adaptive mesh refinement
	std::array<Ddpca::Real,2> loadIncr;//load intensity of the worm inner hub

	std::vector<double> wodeAuan;//worm decomposition auxiliary angle
	std::vector<std::vector<Ddpca::CurvedSurface>> wodeAucu;//worm decomposition auxiliary CURVEDS
	Ddpca::CurvedSurface whdeAucu;//wheel decomposition auxiliary CURVEDS
	Ddpca::CurvedSurface whdeAucu_midd;//wheel decomposition auxiliary CURVEDS

public:
	/*********************************************************************************************/

    void COOR_AVER(const std::vector<Ddpca::Coordinate> &inpuCoor, Ddpca::Coordinate &outpCoor){
        Ddpca::Coordinate cyliCoor(0.0, 0.0, 0.0);
        long tempFlag_0 = 0;
        long tempFlag_1 = 0;
        const Ddpca::I64 incoSize = inpuCoor.size();
        for(long tk = 0; tk < incoSize; tk ++){
            cyliCoor[0] += std::sqrt(inpuCoor[tk][0] * inpuCoor[tk][0] + 
                inpuCoor[tk][1] * inpuCoor[tk][1]);
            double tempAngl = std::atan2(inpuCoor[tk][1], inpuCoor[tk][0]);
            cyliCoor[1] += tempAngl;
            cyliCoor[2] += inpuCoor[tk][2];
            if(tempAngl > Ddpca::PI / 2.0){
                tempFlag_0 ++;
            }
            if(tempAngl < - Ddpca::PI / 2.0){
                tempFlag_1 ++;
            }
        }
        if(tempFlag_0 > 0 && tempFlag_1 > 0){
            cyliCoor[1] += tempFlag_1 * (Ddpca::PI * 2.0);
        }
        outpCoor[0] = cyliCoor[0] / incoSize * std::cos(cyliCoor[1] / incoSize);
        outpCoor[1] = cyliCoor[0] / incoSize * std::sin(cyliCoor[1] / incoSize);
        outpCoor[2] = cyliCoor[2] / incoSize;
    }

    void COOR_AVER_1(const std::vector<Ddpca::Coordinate> &inpuCoor, Ddpca::Coordinate &outpCoor){
        Ddpca::Coordinate cyliCoor(0.0, 0.0, 0.0);
        long tempFlag_0 = 0;
        long tempFlag_1 = 0;
        bool inneHubf = true;
        const Ddpca::I64 incoSize = inpuCoor.size();
        for(long tk = 0; tk < incoSize; tk ++){
            double tempRadi = std::sqrt(inpuCoor[tk][0] * inpuCoor[tk][0] + 
                inpuCoor[tk][1] * inpuCoor[tk][1]);
            if(std::abs(tempRadi - dehwSurf.inneRadi[1]) > 1.0E-10){
                inneHubf = false;
            }
            double tempAngl = std::atan2(inpuCoor[tk][1], inpuCoor[tk][0]);
            cyliCoor[1] += tempAngl;
            if(tempAngl > Ddpca::PI / 2.0){
                tempFlag_0 ++;
            }
            if(tempAngl < - Ddpca::PI / 2.0){
                tempFlag_1 ++;
            }
        }
        if(inneHubf == true){
            COOR_AVER(inpuCoor, outpCoor);
            return;
        }
        if(tempFlag_0 > 0 && tempFlag_1 > 0){
            cyliCoor[1] += tempFlag_1 * (Ddpca::PI * 2.0);
        }
        //
        double toruRadi = 0.0;
        double toruAngl = 0.0;
        for(long tk = 0; tk < incoSize; tk ++){
            double tempRadi = dehwSurf.a_h2 + centErro - 
                std::sqrt(inpuCoor[tk][0] * inpuCoor[tk][0] + inpuCoor[tk][1] * inpuCoor[tk][1]);
            toruRadi += std::sqrt(std::pow(tempRadi, 2.0) + std::pow(inpuCoor[tk][2], 2.0));
            toruAngl += std::atan2(inpuCoor[tk][2], tempRadi);
        }
        toruRadi /= incoSize;
        toruAngl /= incoSize;
        //
        cyliCoor[0] = dehwSurf.a_h2 + centErro - toruRadi * std::cos(toruAngl);
        cyliCoor[1] /= incoSize;
        cyliCoor[2] = toruRadi * std::sin(toruAngl);
        outpCoor[0] = cyliCoor[0] * std::cos(cyliCoor[1]);
        outpCoor[1] = cyliCoor[0] * std::sin(cyliCoor[1]);
        outpCoor[2] = cyliCoor[2];
    }

    void UPDA_WODE(const std::vector<long> &inpuNode, 
        const Ddpca::Coordinate &outpCoor, long face_tw, long face_sub, long worm_tv){
        //
        Ddpca::Mesh& twMesh = domains[worm_tv].mesh;
        long ti_aver = 0;
        long tj_aver = 0;
        bool tempFlag = true;
        const Ddpca::I64 innoSize = inpuNode.size();
        for(long ic = 0; ic < innoSize; ic ++){
            auto iterNoco = twMesh.node2Coordinate.find(inpuNode[ic]);
            auto iterPoin = wodeAucu[face_tw][face_sub].pointIndex.find(iterNoco->second);
            if(iterPoin == wodeAucu[face_tw][face_sub].pointIndex.end()){
                tempFlag = false;
                break;
            }
            ti_aver += (iterPoin->second)[0];
            tj_aver += (iterPoin->second)[1];
        }
        if(tempFlag == false){
            return;
        }
        ti_aver /= innoSize;
        tj_aver /= innoSize;
        wodeAucu[face_tw][face_sub].Insert(ti_aver, tj_aver, outpCoor);
    }

    long UPDA_WHDE(const std::vector<long> &inpuNode, const Ddpca::Coordinate &outpCoor, long whee_tv){
        //
        Ddpca::Mesh& whMesh = domains[whee_tv].mesh;
        long ti_aver = 0;
        long tj_aver = 0;
        bool tempFlag = true;
        const Ddpca::I64 innoSize = inpuNode.size();
        for(long ic = 0; ic < innoSize; ic ++){
            auto iterNoco = whMesh.node2Coordinate.find(inpuNode[ic]);
            auto iterPoin = whdeAucu.pointIndex.find(iterNoco->second);
            if(iterPoin == whdeAucu.pointIndex.end()){
                tempFlag = false;
                break;
            }
            ti_aver += (iterPoin->second)[0];
            tj_aver += (iterPoin->second)[1];
        }
        if(tempFlag == false){
            return 1;
        }
        ti_aver /= innoSize;
        tj_aver /= innoSize;
        whdeAucu.Insert(ti_aver, tj_aver, outpCoor);
        return 1;
    }

    long UPDA_WHDE_MIDD(const std::vector<long> &inpuNode, const Ddpca::Coordinate &outpCoor, long whee_tv){
        //
        Ddpca::Mesh& whMesh = domains[whee_tv].mesh;
        long ti_aver = 0;
        long tj_aver = 0;
        bool tempFlag = true;
        const Ddpca::I64 innoSize = inpuNode.size();
        for(long ic = 0; ic < innoSize; ic ++){
            auto iterNoco = whMesh.node2Coordinate.find(inpuNode[ic]);
            auto iterPoin = whdeAucu_midd.pointIndex.find(iterNoco->second);
            if(iterPoin == whdeAucu_midd.pointIndex.end()){
                tempFlag = false;
                break;
            }
            ti_aver += (iterPoin->second)[0];
            tj_aver += (iterPoin->second)[1];
        }
        if(tempFlag == false){
            return 1;
        }
        ti_aver /= innoSize;
        tj_aver /= innoSize;
        whdeAucu_midd.Insert(ti_aver, tj_aver, outpCoor);
        return 1;
    }

    void GenerateMesh(){
        //
        Ddpca::Log("TestDehw::GenerateMesh");
        //worm
        //column-major order
        Ddpca::DenseMatrix wormRota_0(3,3,{
            std::cos(analAngl[0]),std::sin(analAngl[0]),0.0,
            -std::sin(analAngl[0]),std::cos(analAngl[0]),0.0,
            0.0,0.0,1.0
        });
        Ddpca::DenseMatrix wormRota_1(3,3,{
            1.0,0.0,0.0,
            0.0,0.0,-1.0,
            0.0,1.0,0.0
        });
        Ddpca::DenseMatrix tempMatr = wormRota_1;
        Ddpca::GEMM(tempMatr, wormRota_0, wormRota_1);
        Ddpca::Coordinate wormTran(- (dehwSurf.a_h2 + centErro), 0.0, 0.0);
        wodeAuan.resize(dehwSurf.gridNumb[0][6] - 1);
        wodeAucu.resize(dehwSurf.gridNumb[0][6]);
        for(Ddpca::I64 ti = 0; ti < dehwSurf.gridNumb[0][6]; ++ ti){
            wodeAucu[ti].resize(2);
        }
        Ddpca::I64 wodeFact_0 = (1 << (dehwSurf.globHomo));
        Ddpca::I64 wodeFact_1 = (1 << (dehwSurf.globInho + dehwSurf.globHomo));
        //wheel
        double wheeAngl_1 = analAngl[1] - 2.0 * Ddpca::PI / dehwSurf.z[1] * 2.0;
        Ddpca::DenseMatrix wheeRota(3,3,{
            std::cos(wheeAngl_1),std::sin(wheeAngl_1),0.0,
            -std::sin(wheeAngl_1),std::cos(wheeAngl_1),0.0,
            0.0,0.0,1.0
        });
        Ddpca::Coordinate wheeTran(0.0, 0.0, 0.0);
        long numbFace = dehwSurf.gridNumb[1][4];
        long whdeFact_0 = (1 << (dehwSurf.globHomo));
        long whdeFact_1 = (1 << (dehwSurf.globInho + dehwSurf.globHomo));
        whdeAucu.Resize(
            (dehwSurf.gridNumb[1][1] + dehwSurf.gridNumb[1][2]) * whdeFact_0 + 1, 
            numbFace * whdeFact_1 + 1);
        whdeAucu_midd.Resize(
            (dehwSurf.gridNumb[1][1] + dehwSurf.gridNumb[1][3]) * whdeFact_0 + 1, 
            numbFace * whdeFact_1 + 1);
        //points
        Ddpca::I64 tw_endi = dehwSurf.gridNumb[0][6];
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tw){
            Ddpca::Mesh& twMesh = domains[tw].mesh;
            if(tw < tw_endi){
                long numb_tw, numbStar;
                if(tw == 0 || tw == dehwSurf.gridNumb[0][6] - 1){
                    numb_tw = dehwSurf.gridNumb[0][5];
                    if(tw == 0){
                        numbStar = 0;
                    }
                    else{
                        numbStar = dehwSurf.gridNumb[0][5] + (tw - 1) * dehwSurf.gridNumb[0][4];
                    }
                }
                else{
                    numb_tw = dehwSurf.gridNumb[0][4];
                    numbStar = dehwSurf.gridNumb[0][5] + (tw - 1) * dehwSurf.gridNumb[0][4];
                }
                for(long ti = 0; ti < 2; ti ++){
                    wodeAucu[tw][ti].Resize(
                        (dehwSurf.gridNumb[0][1] + dehwSurf.gridNumb[0][2]) * wodeFact_0 + 1, 
                        numb_tw * wodeFact_1 + 1);
                }
                if(tw >= 1){
                    long ti_real = numbStar 
                        * (1 << (dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve));
                    double tempAuan = dehwSurf.curvCoor[0][ti_real][0][0];
                    while(tempAuan > Ddpca::PI){
                        tempAuan -= 2.0 * Ddpca::PI;
                    }
                    wodeAuan[tw - 1] = tempAuan;
                }
                std::vector<std::vector<std::vector<std::vector<Ddpca::I64>>>> blocPoin;
                blocPoin.resize(4);
                for(Ddpca::I64 ti = 0; ti < 4; ++ ti){
                    blocPoin[ti].resize(numb_tw + 1);
                }
                for(long ti = 0; ti <= numb_tw; ti ++){
                    blocPoin[0][ti].resize(dehwSurf.gridNumb[0][1] + 1);
                    for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[0][1] + 1; ++ tj){
                        blocPoin[0][ti][tj].resize(dehwSurf.gridNumb[0][0] + 1);
                    }
                    blocPoin[1][ti].resize(dehwSurf.gridNumb[0][2] + 1);
                    for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[0][2] + 1; ++ tj){
                        blocPoin[1][ti][tj].resize(dehwSurf.gridNumb[0][0] / 2 + 1);
                    }
                    blocPoin[2][ti].resize(dehwSurf.gridNumb[0][3] + 1);
                    for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[0][3] + 1; ++ tj){
                        blocPoin[2][ti][tj].resize(2 * dehwSurf.gridNumb[0][2] + 1);
                    }
                    blocPoin[3][ti].resize(dehwSurf.gridNumb[0][2] + 1);
                    for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[0][2] + 1; ++ tj){
                        blocPoin[3][ti][tj].resize(dehwSurf.gridNumb[0][0] / 2 + 1);
                    }
                    Ddpca::I64 ti_real = (numbStar + ti) 
                        * (1 << (dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve));
                    std::vector<Ddpca::DenseMatrix> wormProf(dehwSurf.gridNumb[0][3] + 1);
                    for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[0][3] + 1; ++ tj){
                        wormProf[tj].Resize(3,2);
                    }
                    for(long tj = 0; tj <= dehwSurf.gridNumb[0][3]; tj ++){
                        long tj_real = tj * (1 << (dehwSurf.globHomo + dehwSurf.locaLeve));
                        wormProf[tj](0,0) = dehwSurf.wormTosu.indexPoint[ti_real][tj_real][0];
                        wormProf[tj](1,0) = dehwSurf.wormTosu.indexPoint[ti_real][tj_real][1];
                        wormProf[tj](2,0) = dehwSurf.wormTosu.indexPoint[ti_real][tj_real][2];
                        wormProf[tj](0,1) = dehwSurf.wormToba.indexPoint[ti_real][tj_real][0];
                        wormProf[tj](1,1) = dehwSurf.wormToba.indexPoint[ti_real][tj_real][1];
                        wormProf[tj](2,1) = dehwSurf.wormToba.indexPoint[ti_real][tj_real][2];
                    }
                    ti_real = (numbStar + ti) 
                        * (1 << (dehwSurf.globInho + dehwSurf.globHomo));
                    Ddpca::DenseMatrix rootProf_1(3, dehwSurf.gridNumb[0][0] / 2 + 1);
                    Ddpca::DenseMatrix rootProf_2(3, dehwSurf.gridNumb[0][0] / 2 + 1);
                    for(long tj = 0; tj <= dehwSurf.gridNumb[0][0] / 2; tj ++){
                        long tj_real = tj * (1 << (dehwSurf.globHomo));
                        rootProf_1(0,tj) = dehwSurf.wormRtsu.indexPoint[ti_real][tj_real][0];
                        rootProf_1(1,tj) = dehwSurf.wormRtsu.indexPoint[ti_real][tj_real][1];
                        rootProf_1(2,tj) = dehwSurf.wormRtsu.indexPoint[ti_real][tj_real][2];
                        rootProf_2(0,tj) = dehwSurf.wormRtba.indexPoint[ti_real][tj_real][0];
                        rootProf_2(1,tj) = dehwSurf.wormRtba.indexPoint[ti_real][tj_real][1];
                        rootProf_2(2,tj) = dehwSurf.wormRtba.indexPoint[ti_real][tj_real][2];
                    }
                    std::array<Ddpca::Real,2> profRadi;
                    profRadi = {
                        std::sqrt(rootProf_1(0,0) * rootProf_1(0,0) + rootProf_1(1,0) * rootProf_1(1,0)), 
                        std::sqrt(rootProf_2(0,0) * rootProf_2(0,0) + rootProf_2(1,0) * rootProf_2(1,0))};
                    std::array<Ddpca::Real,2> tranRadi;
                    tranRadi[0] = profRadi[0] - Ddpca::PI / 4.0 * dehwSurf.m_t;
                    tranRadi[1] = profRadi[1] - Ddpca::PI / 4.0 * dehwSurf.m_t;
                    Ddpca::DenseMatrix blocCoor(3,4,{
                        rootProf_1(0,0) * dehwSurf.inneRadi[0] / profRadi[0],rootProf_1(1,0) * dehwSurf.inneRadi[0] / profRadi[0],rootProf_1(2,0),
                        rootProf_2(0,0) * dehwSurf.inneRadi[0] / profRadi[1],rootProf_2(1,0) * dehwSurf.inneRadi[0] / profRadi[1],rootProf_2(2,0),
                        rootProf_2(0,0) * tranRadi[1] / profRadi[1],rootProf_2(1,0) * tranRadi[1] / profRadi[1],rootProf_2(2,0),
                        rootProf_1(0,0) * tranRadi[0] / profRadi[0],rootProf_1(1,0) * tranRadi[0] / profRadi[0],rootProf_1(2,0)
                    });
                    //0
                    for(long tj = 0; tj <= dehwSurf.gridNumb[0][1]; tj ++){
                        for(long tk = 0; tk <= dehwSurf.gridNumb[0][0]; tk ++){
                            Ddpca::DenseMatrix downUp(3,2);
                            std::array<Ddpca::Real,3> tempCoor;
                            for(Ddpca::I64 tl = 0; tl < 3; ++ tl){
                                downUp(tl,0) = 
                                    (1.0 - (double)tk / dehwSurf.gridNumb[0][0]) * blocCoor(tl,0) 
                                    + (double)tk / dehwSurf.gridNumb[0][0] * blocCoor(tl,1);
                                downUp(tl,1) = 
                                    (1.0 - (double)tk / dehwSurf.gridNumb[0][0]) * blocCoor(tl,3) 
                                    + (double)tk / dehwSurf.gridNumb[0][0] * blocCoor(tl,2);
                                tempCoor[tl] = 
                                    (1.0 - (double)tj / dehwSurf.gridNumb[0][1]) * downUp(tl,0) 
                                    + (double)tj / dehwSurf.gridNumb[0][1] * downUp(tl,1);
                            }
                            Ddpca::Coordinate inpuCoor = {tempCoor[0], tempCoor[1], tempCoor[2]};
                            blocPoin[0][ti][tj][tk] = twMesh.TryAddNode(inpuCoor);
                            if(tk == 0){
                                wodeAucu[tw][0].Insert(tj * wodeFact_0, ti * wodeFact_1, inpuCoor);
                            }
                            else if(tk == dehwSurf.gridNumb[0][0]){
                                wodeAucu[tw][1].Insert(tj * wodeFact_0, ti * wodeFact_1, inpuCoor);
                            }
                        }
                    }
                    //1
                    Ddpca::DenseMatrix lineCoor(3,2);
                    for(Ddpca::I64 tj = 0; tj < 3; ++ tj){
                        lineCoor(tj,0) = blocCoor(tj,3);
                        lineCoor(tj,1) = 0.5 * blocCoor(tj,3) + 0.5 * blocCoor(tj,2);
                    }
                    for(long tj = 0; tj <= dehwSurf.gridNumb[0][2]; tj ++){
                        for(long tk = 0; tk <= dehwSurf.gridNumb[0][0] / 2; tk ++){
                            Ddpca::DenseMatrix downUp(3,2);
                            std::array<Ddpca::Real,3> tempCoor;
                            for(Ddpca::I64 tl = 0; tl < 3; ++ tl){
                                downUp(tl,0) = 
                                    (1.0 - (double)tk / (dehwSurf.gridNumb[0][0] / 2.0)) 
                                    * lineCoor(tl,0) 
                                    + (double)tk / (dehwSurf.gridNumb[0][0]/2.0) 
                                    * lineCoor(tl,1);
                                downUp(tl,1) = rootProf_1(tl,tk);
                                tempCoor[tl] = 
                                    (1.0 - (double)tj / dehwSurf.gridNumb[0][2]) 
                                    * downUp(tl,0) 
                                    + (double)tj / dehwSurf.gridNumb[0][2] * downUp(tl,1);
                            }
                            Ddpca::Coordinate inpuCoor = {tempCoor[0], tempCoor[1], tempCoor[2]};
                            blocPoin[1][ti][tj][tk] = twMesh.TryAddNode(inpuCoor);
                            if(tk == 0){
                                wodeAucu[tw][0].Insert(
                                    (dehwSurf.gridNumb[0][1] + tj) * wodeFact_0, 
                                    ti * wodeFact_1, inpuCoor
                                );
                            }
                        }
                    }
                    //2
                    for(Ddpca::I64 tj = 0; tj < 3; ++ tj){
                        lineCoor(tj,0) = lineCoor(tj,1);
                        lineCoor(tj,1) = 0.5 * wormProf[dehwSurf.gridNumb[0][3]](tj,0) 
                            + 0.5 * wormProf[dehwSurf.gridNumb[0][3]](tj,1);
                    }
                    for(long tj = 0; tj <= dehwSurf.gridNumb[0][3]; tj ++){
                        for(long tk = 0; tk <= dehwSurf.gridNumb[0][2]; tk ++){
                            Ddpca::DenseMatrix downUp(3,2);
                            std::array<Ddpca::Real,3> tempCoor;
                            for(Ddpca::I64 tl = 0; tl < 3; ++ tl){
                                downUp(tl,0) = wormProf[tj](tl,0);
                                downUp(tl,1) = 
                                    (1.0 - (double)tj / dehwSurf.gridNumb[0][3]) * lineCoor(tl,0) 
                                    + (double)tj / dehwSurf.gridNumb[0][3] * lineCoor(tl,1);
                                tempCoor[tl] = 
                                    (1.0 - (double)tk / dehwSurf.gridNumb[0][2]) 
                                    * downUp(tl,0) 
                                    + (double)tk / dehwSurf.gridNumb[0][2] * downUp(tl,1);
                            }
                            Ddpca::Coordinate inpuCoor = {tempCoor[0], tempCoor[1], tempCoor[2]};
                            blocPoin[2][ti][tj][tk] = twMesh.TryAddNode(inpuCoor);
                        }
                        for(long tk = 1; tk <= dehwSurf.gridNumb[0][2]; tk ++){
                            Ddpca::DenseMatrix downUp(3,2);
                            std::array<Ddpca::Real,3> tempCoor;
                            for(Ddpca::I64 tl = 0; tl < 3; ++ tl){
                                downUp(tl,0) = 
                                    (1.0 - (double)tj / dehwSurf.gridNumb[0][3]) 
                                    * lineCoor(tl,0) 
                                    + (double)tj / dehwSurf.gridNumb[0][3] * lineCoor(tl,1);
                                downUp(tl,1) = wormProf[tj](tl,1);
                                tempCoor[tl] = 
                                    (1.0 - (double)tk / dehwSurf.gridNumb[0][2]) * downUp(tl,0) 
                                    + (double)tk / dehwSurf.gridNumb[0][2] * downUp(tl,1);
                            }
                            Ddpca::Coordinate inpuCoor = {tempCoor[0], tempCoor[1], tempCoor[2]};
                            blocPoin[2][ti][tj][dehwSurf.gridNumb[0][2] + tk] = twMesh.TryAddNode(inpuCoor);
                        }
                    }
                    //3
                    for(Ddpca::I64 tj = 0; tj < 3; ++ tj){
                        lineCoor(tj,1) = blocCoor(tj,2);
                    }
                    for(long tj = 0; tj <= dehwSurf.gridNumb[0][2]; tj ++){
                        for(long tk = 0; tk <= dehwSurf.gridNumb[0][0] / 2; tk ++){
                            Ddpca::DenseMatrix downUp(3,2);
                            std::array<Ddpca::Real,3> tempCoor;
                            for(Ddpca::I64 tl = 0; tl < 3; ++ tl){
                                downUp(tl,0) = 
                                    (1.0 - (double)tk / (dehwSurf.gridNumb[0][0] / 2.0)) 
                                    * lineCoor(tl,0) 
                                    + (double)tk / (dehwSurf.gridNumb[0][0]/2.0) 
                                    * lineCoor(tl,1);
                                downUp(tl,1) = rootProf_2(tl,dehwSurf.gridNumb[0][0] / 2 - tk);
                                tempCoor[tl] = 
                                    (1.0 - (double)tj / dehwSurf.gridNumb[0][2]) * downUp(tl,0) 
                                    + (double)tj / dehwSurf.gridNumb[0][2] * downUp(tl,1);
                            }
                            Ddpca::Coordinate inpuCoor = {tempCoor[0], tempCoor[1], tempCoor[2]};
                            blocPoin[3][ti][tj][tk] = twMesh.TryAddNode(inpuCoor);
                            if(tk == dehwSurf.gridNumb[0][0] / 2){
                                wodeAucu[tw][1].Insert(
                                    (dehwSurf.gridNumb[0][1] + tj) * wodeFact_0, 
                                    ti * wodeFact_1, inpuCoor
                                );
                            }
                        }
                    }
                }
                //elements
                Ddpca::OctreeElement tempElem;
                const Ddpca::I64 blpoSize_0 = blocPoin.size();
                for(long tj = 0; tj < blpoSize_0; tj ++){
                    const Ddpca::I64 blpoSize_1 = blocPoin[tj].size();
                    for(long tk = 0; tk < blpoSize_1 - 1; tk ++){
                        const Ddpca::I64 blpoSize_2 = blocPoin[tj][tk].size();
                        for(long tm = 0; tm < blpoSize_2 - 1; tm ++){
                            const Ddpca::I64 blpoSize_3 = blocPoin[tj][tk][tm].size();
                            for(long tn = 0; tn < blpoSize_3 - 1; tn ++){
                                tempElem.parent = -1;
                                tempElem.cornerNodes[0] = blocPoin[tj][tk][tm][tn];
                                tempElem.cornerNodes[1] = blocPoin[tj][tk][tm + 1][tn];
                                tempElem.cornerNodes[2] = blocPoin[tj][tk][tm + 1][tn + 1];
                                tempElem.cornerNodes[3] = blocPoin[tj][tk][tm][tn + 1];
                                tempElem.cornerNodes[4] = blocPoin[tj][tk + 1][tm][tn];
                                tempElem.cornerNodes[5] = blocPoin[tj][tk + 1][tm + 1][tn];
                                tempElem.cornerNodes[6] = blocPoin[tj][tk + 1][tm + 1][tn + 1];
                                tempElem.cornerNodes[7] = blocPoin[tj][tk + 1][tm][tn + 1];
                                tempElem.level = 0;
                                tempElem.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
                                tempElem.children.clear();
                                twMesh.AddElement(tempElem);
                            }
                        }
                    }
                }
                //global refinement
                std::set<Ddpca::I64> elementsToSplit;
                std::map<Ddpca::I64, std::set<Ddpca::I64>> subElements;
                std::map<std::vector<Ddpca::I64>, Ddpca::Coordinate> curvInte;
                for(long tr = 0; tr < dehwSurf.globInho + dehwSurf.globHomo; tr ++){
                    elementsToSplit.clear();
                    const Ddpca::I64 tmelSize = twMesh.elements.size();
                    for(long ti = 0; ti < tmelSize; ti ++){
                        if(twMesh.elements[ti].children.size() > 0){
                            continue;
                        }
                        elementsToSplit.insert(ti);
                        if(tr < dehwSurf.globInho){
                            twMesh.elements[ti].refinementPattern = Ddpca::OctreeElement::REFINEMENT_ZETA;
                        }
                        else{
                            twMesh.elements[ti].refinementPattern = Ddpca::OctreeElement::REFINEMENT_FULL;
                        }
                    }
                    curvInte.clear();
                    std::array<Ddpca::CurvedSurface*,4> wormSurf = {
                        &(dehwSurf.wormTosu), &(dehwSurf.wormToba), 
                        &(dehwSurf.wormRtsu), &(dehwSurf.wormRtba)
                    };
                    for(long ts = 0; ts < 4; ts ++){
                        (* wormSurf[ts]).Refine(twMesh, elementsToSplit, curvInte);
                    }
                    for(const auto &iterSpel : elementsToSplit){
                        //
                        std::vector<long> inpuNode(8);
                        std::vector<Ddpca::Coordinate> inpuCoor(8);
                        for(long tk = 0; tk < 8; tk ++){
                            inpuNode[tk] = twMesh.elements[iterSpel].cornerNodes[tk];
                            auto iterNoco = twMesh.node2Coordinate.find(inpuNode[tk]);
                            inpuCoor[tk] = iterNoco->second;
                        }
                        Ddpca::Coordinate outpCoor;
                        COOR_AVER(inpuCoor, outpCoor);
                        std::sort(inpuNode.begin(), inpuNode.end());
                        curvInte.emplace(inpuNode, outpCoor);
                        //
                        const Ddpca::I64 hsliSize = Ddpca::hexaLine.size();
                        for(long tj = 0; tj < hsliSize; tj ++){
                            const Ddpca::I64 hltjSize = Ddpca::hexaLine[tj].size();
                            std::vector<long> inpuNode(hltjSize);
                            std::vector<Ddpca::Coordinate> inpuCoor(hltjSize);
                            for(long tk = 0; tk < hltjSize; tk ++){
                                inpuNode[tk] = twMesh.elements[iterSpel].cornerNodes[Ddpca::hexaLine[tj][tk]];
                                auto iterNoco = twMesh.node2Coordinate.find(inpuNode[tk]);
                                inpuCoor[tk] = iterNoco->second;
                            }
                            Ddpca::Coordinate outpCoor;
                            COOR_AVER(inpuCoor, outpCoor);
                            std::sort(inpuNode.begin(), inpuNode.end());
                            //the same key: new key will be abandoned
                            curvInte.emplace(inpuNode, outpCoor);
                        }
                        //
                        Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                        for(long tj = 0; tj < hefaSize; tj ++){
                            Ddpca::I64 hftjSize = Ddpca::hexaFace[tj].size();
                            std::vector<long> inpuNode(hftjSize);
                            std::vector<Ddpca::Coordinate> inpuCoor(hftjSize);
                            for(long tk = 0; tk < hftjSize; tk ++){
                                inpuNode[tk] = twMesh.elements[iterSpel].cornerNodes[Ddpca::hexaFace[tj][tk]];
                                auto iterNoco = twMesh.node2Coordinate.find(inpuNode[tk]);
                                inpuCoor[tk] = iterNoco->second;
                            }
                            Ddpca::Coordinate outpCoor;
                            COOR_AVER(inpuCoor, outpCoor);
                            std::sort(inpuNode.begin(), inpuNode.end());
                            //the same key: new key will be abandoned
                            curvInte.emplace(inpuNode, outpCoor);
                        }
                    }
                    twMesh.Refine(elementsToSplit, subElements, curvInte);
                    for(const auto &iterPlsu : curvInte){
                        UPDA_WODE(iterPlsu.first, iterPlsu.second, tw, 0, tw);
                        UPDA_WODE(iterPlsu.first, iterPlsu.second, tw, 1, tw);
                    }
                }
                //
                twMesh.RigidRotationTranslation(wormRota_1, wormTran);
                wodeAucu[tw][0].RigidRotationTranslation(wormRota_1, wormTran);
                wodeAucu[tw][1].RigidRotationTranslation(wormRota_1, wormTran);
            }
            else{
                domains[tw].elasticity = 110.0E9;
                long toot_tw = (tw - dehwSurf.gridNumb[0][6]) / dehwSurf.gridNumb[1][6];
                long leri_tw = (tw - dehwSurf.gridNumb[0][6]) % dehwSurf.gridNumb[1][6];
                //point
                std::vector<std::vector<std::vector<std::vector<Ddpca::I64>>>> blocPoin;
                blocPoin.resize(4);
                for(Ddpca::I64 ti = 0; ti < 4; ++ ti){
                    blocPoin[ti].resize(numbFace + 1);
                }
                for(long ti = 0; ti <= numbFace; ti ++){
                    blocPoin[0][ti].resize(dehwSurf.gridNumb[1][1] + 1);
                    for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[1][1] + 1; ++ tj){
                        blocPoin[0][ti][tj].resize(dehwSurf.gridNumb[1][0] / 2 + 1);
                    }
                    if(leri_tw == 0){
                        blocPoin[1][ti].resize(dehwSurf.gridNumb[1][2] + 1);
                        for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[1][2] + 1; ++ tj){
                            blocPoin[1][ti][tj].resize(dehwSurf.gridNumb[1][0] / 2 + 1);
                        }
                    }
                    else{
                        blocPoin[1][ti].clear();
                    }
                    blocPoin[2][ti].resize(dehwSurf.gridNumb[1][3] + 1);
                    for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[1][3] + 1; ++ tj){
                        blocPoin[2][ti][tj].resize(dehwSurf.gridNumb[1][2] + 1);
                    }
                    if(leri_tw == 1){
                        blocPoin[3][ti].resize(dehwSurf.gridNumb[1][2] + 1);
                        for(Ddpca::I64 tj = 0; tj < dehwSurf.gridNumb[1][2] + 1; ++ tj){
                            blocPoin[3][ti][tj].resize(dehwSurf.gridNumb[1][0] / 2 + 1);
                        }
                    }
                    else{
                        blocPoin[3][ti].clear();
                    }
                    long ti_real = ti;
                    ti_real *= (1 << (dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve));
                    std::vector<Ddpca::DenseMatrix> wheeProf(dehwSurf.gridNumb[1][3] + 1);
                    for(Ddpca::I64 tj = 0; tj < (dehwSurf.gridNumb[1][3] + 1); ++ tj){
                        wheeProf[tj].Resize(3, 2);
                    }
                    for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
                        long tj_real = tj * (1 << (dehwSurf.globHomo + dehwSurf.locaLeve));
                        wheeProf[tj](0,0) = dehwSurf.wheeTosu.indexPoint[ti_real][tj_real][0];
                        wheeProf[tj](1,0) = dehwSurf.wheeTosu.indexPoint[ti_real][tj_real][1];
                        wheeProf[tj](2,0) = dehwSurf.wheeTosu.indexPoint[ti_real][tj_real][2];
                        wheeProf[tj](0,1) = dehwSurf.wheeToba.indexPoint[ti_real][tj_real][0];
                        wheeProf[tj](1,1) = dehwSurf.wheeToba.indexPoint[ti_real][tj_real][1];
                        wheeProf[tj](2,1) = dehwSurf.wheeToba.indexPoint[ti_real][tj_real][2];
                    }
                    //transition into "unfolded cone surface"
                    double tempAlph_3 = - dehwSurf.curvCoor[1][ti_real][0][0];
                    std::vector<std::array<std::array<Ddpca::Real,2>,2>> wheeProf_(dehwSurf.gridNumb[1][3] + 1);
                    for(long tk = 0; tk <= 1; tk ++){
                        for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
                            long tj_real = tj * (1 << (dehwSurf.globHomo + dehwSurf.locaLeve));
                            std::array<Ddpca::Real,3> tempCoor;
                            if(tk == 0){
                                tempCoor[0] = dehwSurf.wheeTosu.indexPoint[ti_real][tj_real][0]; 
                                tempCoor[1] = dehwSurf.wheeTosu.indexPoint[ti_real][tj_real][1]; 
                                tempCoor[2] = dehwSurf.wheeTosu.indexPoint[ti_real][tj_real][2];
                            }
                            else{
                                tempCoor[0] = dehwSurf.wheeToba.indexPoint[ti_real][tj_real][0]; 
                                tempCoor[1] = dehwSurf.wheeToba.indexPoint[ti_real][tj_real][1]; 
                                tempCoor[2] = dehwSurf.wheeToba.indexPoint[ti_real][tj_real][2];
                            }
                            wheeProf_[tj][tk]  = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
                        }
                    }
                    ti_real = ti;
                    ti_real *= (1 << (dehwSurf.globInho + dehwSurf.globHomo));
                    double r_1f = dehwSurf.a_h2 / std::cos(tempAlph_3) 
                        - (dehwSurf.a_h2 - dehwSurf.d_f[1] / 2.0);
                    std::array<std::vector<Ddpca::Real>,2> rootProf_0, rootProf_1;
                    rootProf_0[0].resize(dehwSurf.gridNumb[1][0] / 2 + 1);
                    rootProf_0[1].resize(dehwSurf.gridNumb[1][0] / 2 + 1);
                    rootProf_1[0].resize(dehwSurf.gridNumb[1][0] / 2 + 1);
                    rootProf_1[1].resize(dehwSurf.gridNumb[1][0] / 2 + 1);
                    for(long tj = 0; tj <= dehwSurf.gridNumb[1][0] / 2; tj ++){
                        long tj_real = tj * (1 << (dehwSurf.globHomo));
                        std::array<Ddpca::Real,3> tempCoor;
                        tempCoor[0] = dehwSurf.wheeRtsu.indexPoint[ti_real][tj_real][0];
                        tempCoor[1] = dehwSurf.wheeRtsu.indexPoint[ti_real][tj_real][1];
                        tempCoor[2] = dehwSurf.wheeRtsu.indexPoint[ti_real][tj_real][2];
                        std::array<Ddpca::Real,2> tempResu = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
                        rootProf_0[0][tj] = tempResu[0];
                        rootProf_0[1][tj] = tempResu[1];
                        tempCoor[0] = dehwSurf.wheeRtba.indexPoint[ti_real][tj_real][0];
                        tempCoor[1] = dehwSurf.wheeRtba.indexPoint[ti_real][tj_real][1];
                        tempCoor[2] = dehwSurf.wheeRtba.indexPoint[ti_real][tj_real][2];
                        tempResu = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
                        rootProf_1[0][tj] = tempResu[0];
                        rootProf_1[1][tj] = tempResu[1];
                    }
                    double tranRadi_0 = r_1f - M_PI / 4.0 * dehwSurf.m_t;
                    std::array<Ddpca::Real,2> tempAngl_0;
                    tempAngl_0[0] = std::atan2(rootProf_0[1][0], rootProf_0[0][0]);
                    tempAngl_0[1] = std::atan2(rootProf_1[1][0], rootProf_1[0][0]);
                    std::vector<std::vector<Ddpca::Real>> tranProf_0(2), tranProf_1(3), inneProf(3);
                    tranProf_0[0].resize(dehwSurf.gridNumb[1][0] + 1);
                    tranProf_0[1].resize(dehwSurf.gridNumb[1][0] + 1);
                    tranProf_1[0].resize(dehwSurf.gridNumb[1][0] + 1);
                    tranProf_1[1].resize(dehwSurf.gridNumb[1][0] + 1);
                    tranProf_1[2].resize(dehwSurf.gridNumb[1][0] + 1);
                    inneProf[0].resize(dehwSurf.gridNumb[1][0] + 1);
                    inneProf[1].resize(dehwSurf.gridNumb[1][0] + 1);
                    inneProf[2].resize(dehwSurf.gridNumb[1][0] + 1);
                    for(long tj = 0; tj <= dehwSurf.gridNumb[1][0]; tj ++){
                        double tempAngl_j = tempAngl_0[0] 
                            + (tempAngl_0[1] - tempAngl_0[0]) / dehwSurf.gridNumb[1][0] * tj;
                        tranProf_0[0][tj] = tranRadi_0 * std::cos(tempAngl_j);
                        tranProf_0[1][tj] = tranRadi_0 * std::sin(tempAngl_j);
                        std::array<Ddpca::Real,2> tempTemp = {tranProf_0[0][tj], tranProf_0[1][tj]};
                        std::array<Ddpca::Real,3> tempResu = dehwSurf.WHEE_CONE(tempTemp, tempAlph_3);
                        tranProf_1[0][tj] = tempResu[0];
                        tranProf_1[1][tj] = tempResu[1];
                        tranProf_1[2][tj] = tempResu[2];
                        double tempRadi = std::sqrt(std::pow(tranProf_1[0][tj], 2.0) + std::pow(tranProf_1[1][tj], 2.0));
                        inneProf[0][tj] = dehwSurf.inneRadi[1] / tempRadi * tranProf_1[0][tj];
                        inneProf[1][tj] = dehwSurf.inneRadi[1] / tempRadi * tranProf_1[1][tj];
                        inneProf[2][tj] = tranProf_1[2][tj];
                    }
                    //0
                    long tk_star = leri_tw * dehwSurf.gridNumb[1][0] / 2;
                    long tk_endi = tk_star + dehwSurf.gridNumb[1][0] / 2;
                    for(long tj = 0; tj <= dehwSurf.gridNumb[1][1]; tj ++){
                        for(long tk = tk_star; tk <= tk_endi; tk ++){
                            Ddpca::Coordinate tempCoor;
                            for(Ddpca::I64 tl = 0; tl < 3; ++ tl){
                                tempCoor[tl] = 
                                    (1.0 - (double)tj / dehwSurf.gridNumb[1][1]) 
                                    * inneProf[tl][tk] 
                                    + (double)tj / dehwSurf.gridNumb[1][1] 
                                    * tranProf_1[tl][tk];
                            }
                            blocPoin[0][ti][tj][tk - tk_star] = twMesh.TryAddNode(tempCoor);
                            if(toot_tw == 0 && leri_tw == 1 && tk == tk_endi){
                                whdeAucu.Insert(tj * whdeFact_0, ti * whdeFact_1, tempCoor);
                            }
                            if(toot_tw == 0 && leri_tw == 0 && tk == tk_endi){
                                whdeAucu_midd.Insert(tj * whdeFact_0, ti * whdeFact_1, tempCoor);
                            }
                        }
                    }
                    //1
                    if(leri_tw == 0){
                        for(long tj = 0; tj <= dehwSurf.gridNumb[1][2]; tj ++){
                            for(long tk = 0; tk <= dehwSurf.gridNumb[1][0] / 2; tk ++){
                                std::array<Ddpca::Real,2> tempPoin;
                                for(Ddpca::I64 tl = 0; tl < 2; ++ tl){
                                    tempPoin[tl] = 
                                        (1.0 - (double)tj / dehwSurf.gridNumb[1][2]) 
                                        * tranProf_0[tl][tk] 
                                        + (double)tj / dehwSurf.gridNumb[1][2] 
                                        * rootProf_0[tl][tk];
                                }
                                std::array<Ddpca::Real,3> tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
                                blocPoin[1][ti][tj][tk] = twMesh.TryAddNode(Ddpca::Coordinate(
                                    tempCoor[0], tempCoor[1], tempCoor[2]));
                            }
                        }
                    }
                    //2
                    std::array<std::vector<Ddpca::Real>,2> lineCoor;
                    lineCoor[0].resize(dehwSurf.gridNumb[1][3] + 1);
                    lineCoor[1].resize(dehwSurf.gridNumb[1][3] + 1);
                    lineCoor[0][0] = tranProf_0[0][dehwSurf.gridNumb[1][0] / 2];
                    lineCoor[1][0] = tranProf_0[1][dehwSurf.gridNumb[1][0] / 2];
                    for(Ddpca::I64 tj = 0; tj < 2; ++ tj){
                        lineCoor[tj][dehwSurf.gridNumb[1][3]] = 
                            0.5 * (wheeProf_[dehwSurf.gridNumb[1][3]][0][tj] 
                            + wheeProf_[dehwSurf.gridNumb[1][3]][1][tj]);
                    }
                    for(long tj = 1; tj < dehwSurf.gridNumb[1][3]; tj ++){
                        for(Ddpca::I64 tk = 0; tk < 2; ++ tk){
                            lineCoor[tk][tj] = lineCoor[tk][0] 
                                + (lineCoor[tk][dehwSurf.gridNumb[1][3]] 
                                - lineCoor[tk][0]) / dehwSurf.gridNumb[1][3] * tj;
                        }
                    }
                    for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
                        for(long tk = 0; tk <= dehwSurf.gridNumb[1][2]; tk ++){
                            if(leri_tw == 0){
                                std::array<Ddpca::Real,2> tempPoin;
                                for(Ddpca::I64 tl = 0; tl < 2; ++ tl){
                                    tempPoin[tl] = 
                                        (1.0 - (double)tk / dehwSurf.gridNumb[1][2]) * wheeProf_[tj][0][tl] 
                                        + (double)tk / dehwSurf.gridNumb[1][2] 
                                        * lineCoor[tl][tj];
                                }
                                std::array<Ddpca::Real,3> tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
                                blocPoin[2][ti][tj][tk] = twMesh.TryAddNode(Ddpca::Coordinate(
                                    tempCoor[0], tempCoor[1], tempCoor[2]));
                                if(toot_tw == 0 && tk == dehwSurf.gridNumb[1][2]){
                                    Ddpca::Coordinate inpuCoor = {tempCoor[0], tempCoor[1], tempCoor[2]};
                                    whdeAucu_midd.Insert((dehwSurf.gridNumb[1][1] + tj) * whdeFact_0, 
                                        ti * whdeFact_1, inpuCoor);
                                }
                            }
                            else{
                                std::array<Ddpca::Real,2> tempPoin;
                                for(Ddpca::I64 tl = 0; tl < 2; ++ tl){
                                    tempPoin[tl] = 
                                        (1.0 - (double)tk / dehwSurf.gridNumb[1][2]) 
                                        * lineCoor[tl][tj] 
                                        + (double)tk / dehwSurf.gridNumb[1][2] * wheeProf_[tj][1][tl];
                                }
                                std::array<Ddpca::Real,3> tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
                                blocPoin[2][ti][tj][tk] = twMesh.TryAddNode(Ddpca::Coordinate(
                                    tempCoor[0], tempCoor[1], tempCoor[2]));
                            }
                        }
                    }
                    //3
                    if(leri_tw == 1){
                        for(long tj = 0; tj <= dehwSurf.gridNumb[1][2]; tj ++){
                            for(long tk = 0; tk <= dehwSurf.gridNumb[1][0] / 2; tk ++){
                                std::array<Ddpca::Real,2> tempPoin;
                                for(Ddpca::I64 tl = 0; tl < 2; ++ tl){
                                    tempPoin[tl] = (1.0 - (double)tj / dehwSurf.gridNumb[1][2]) 
                                        * tranProf_0[tl][dehwSurf.gridNumb[1][0] / 2 + tk] 
                                        + (double)tj / dehwSurf.gridNumb[1][2] 
                                        * rootProf_1[tl][dehwSurf.gridNumb[1][0] / 2 - tk];
                                }
                                std::array<Ddpca::Real,3> tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
                                blocPoin[3][ti][tj][tk] = twMesh.TryAddNode(Ddpca::Coordinate(
                                        tempCoor[0], tempCoor[1], tempCoor[2]));
                                if(toot_tw == 0 && tk == dehwSurf.gridNumb[1][0] / 2){
                                    Ddpca::Coordinate inpuCoor = {tempCoor[0], tempCoor[1], tempCoor[2]};
                                    whdeAucu.Insert(
                                        (dehwSurf.gridNumb[1][1] + tj) * whdeFact_0, 
                                        ti * whdeFact_1, inpuCoor);
                                }
                            }
                        }
                    }
                }
                //volume
                Ddpca::OctreeElement tempElem;
                long tempSize_tj = blocPoin.size();
                for(long tj = 0; tj < tempSize_tj; tj ++){
                    long tempSize_tk = blocPoin[tj].size() - 1;
                    for(long tk = 0; tk < tempSize_tk; tk ++){
                        long tempSize_tm = blocPoin[tj][tk].size() - 1;
                        for(long tm = 0; tm < tempSize_tm; tm ++){
                            long tempSize_tn = blocPoin[tj][tk][tm].size() - 1;
                            for(long tn = 0; tn < tempSize_tn; tn ++){
                                tempElem.parent = -1;
                                tempElem.cornerNodes[0] = blocPoin[tj][tk][tm][tn];
                                tempElem.cornerNodes[1] = blocPoin[tj][tk][tm + 1][tn];
                                tempElem.cornerNodes[2] = blocPoin[tj][tk][tm + 1][tn + 1];
                                tempElem.cornerNodes[3] = blocPoin[tj][tk][tm][tn + 1];
                                tempElem.cornerNodes[4] = blocPoin[tj][tk + 1][tm][tn];
                                tempElem.cornerNodes[5] = blocPoin[tj][tk + 1][tm + 1][tn];
                                tempElem.cornerNodes[6] = blocPoin[tj][tk + 1][tm + 1][tn + 1];
                                tempElem.cornerNodes[7] = blocPoin[tj][tk + 1][tm][tn + 1];
                                tempElem.level = 0;
                                tempElem.refinementPattern = Ddpca::OctreeElement::REFINEMENT_NONE;
                                tempElem.children.clear();
                                twMesh.AddElement(tempElem);
                            }
                        }
                    }
                }
                //global refinement
                std::set<Ddpca::I64> elementsToSplit;
                std::map<Ddpca::I64, std::set<Ddpca::I64>> subElements;
                std::map<std::vector<Ddpca::I64>, Ddpca::Coordinate> curvInte;
                for(long tr = 0; tr < dehwSurf.globInho + dehwSurf.globHomo; tr ++){
                    elementsToSplit.clear();
                    const Ddpca::I64 twelSize = twMesh.elements.size();
                    for(long ti = 0; ti < twelSize; ti ++){
                        if(twMesh.elements[ti].children.size() > 0){
                            continue;
                        }
                        elementsToSplit.insert(ti);
                        if(tr < dehwSurf.globInho){
                            twMesh.elements[ti].refinementPattern = Ddpca::OctreeElement::REFINEMENT_ZETA;
                        }
                        else{
                            twMesh.elements[ti].refinementPattern = Ddpca::OctreeElement::REFINEMENT_FULL;
                        }
                    }
                    curvInte.clear();
                    std::array<Ddpca::CurvedSurface*,4> wheeSurf = {
                        &(dehwSurf.wheeTosu), &(dehwSurf.wheeToba), 
                        &(dehwSurf.wheeRtsu), &(dehwSurf.wheeRtba)};
                    for(long ts = 0; ts < 4; ts ++){
                        (* wheeSurf[ts]).Refine(twMesh, elementsToSplit, curvInte);
                    }
                    for(const auto &iterSpel : elementsToSplit){
                        //
                        std::vector<long> inpuNode(8);
                        std::vector<Ddpca::Coordinate> inpuCoor(8);
                        for(long tk = 0; tk < 8; tk ++){
                            inpuNode[tk] = twMesh.elements[iterSpel].cornerNodes[tk];
                            auto iterNoco = twMesh.node2Coordinate.find(inpuNode[tk]);
                            inpuCoor[tk] = iterNoco->second;
                        }
                        Ddpca::Coordinate outpCoor;
                        COOR_AVER_1(inpuCoor, outpCoor);
                        std::sort(inpuNode.begin(), inpuNode.end());
                        //the same key: new key will be abandoned
                        curvInte.emplace(inpuNode, outpCoor);
                        //
                        const Ddpca::I64 heliSize = Ddpca::hexaLine.size();
                        for(long tj = 0; tj < heliSize; tj ++){
                            const Ddpca::I64 hltjSize = Ddpca::hexaLine[tj].size();
                            std::vector<long> inpuNode(hltjSize);
                            std::vector<Ddpca::Coordinate> inpuCoor(hltjSize);
                            for(long tk = 0; tk < hltjSize; tk ++){
                                inpuNode[tk] = twMesh.elements[iterSpel].cornerNodes[Ddpca::hexaLine[tj][tk]];
                                auto iterNoco = twMesh.node2Coordinate.find(inpuNode[tk]);
                                inpuCoor[tk] = iterNoco->second;
                            }
                            Ddpca::Coordinate outpCoor;
                            COOR_AVER_1(inpuCoor, outpCoor);
                            std::sort(inpuNode.begin(), inpuNode.end());
                            //the same key: new key will be abandoned
                            curvInte.emplace(inpuNode, outpCoor);
                        }
                        //
                        const Ddpca::I64 hefiSize = Ddpca::hexaFace.size();
                        for(long tj = 0; tj < hefiSize; tj ++){
                            const Ddpca::I64 hftjSize = Ddpca::hexaFace[tj].size();
                            std::vector<long> inpuNode(hftjSize);
                            std::vector<Ddpca::Coordinate> inpuCoor(hftjSize);
                            for(long tk = 0; tk < hftjSize; tk ++){
                                inpuNode[tk] = twMesh.elements[iterSpel].cornerNodes[Ddpca::hexaFace[tj][tk]];
                                auto iterNoco = twMesh.node2Coordinate.find(inpuNode[tk]);
                                inpuCoor[tk] = iterNoco->second;
                            }
                            Ddpca::Coordinate outpCoor;
                            COOR_AVER_1(inpuCoor, outpCoor);
                            std::sort(inpuNode.begin(), inpuNode.end());
                            //the same key: new key will be abandoned
                            curvInte.emplace(inpuNode, outpCoor);
                        }
                    }
                    twMesh.Refine(elementsToSplit, subElements, curvInte);
                    if(toot_tw  == 0 && leri_tw == 1){
                        for(const auto &iterPlsu : curvInte){
                            UPDA_WHDE(iterPlsu.first, iterPlsu.second, tw);
                        }
                    }
                    if(toot_tw  == 0 && leri_tw == 0){
                        for(const auto &iterPlsu : curvInte){
                            UPDA_WHDE_MIDD(iterPlsu.first, iterPlsu.second, tw);
                        }
                    }
                }
                //
                Ddpca::DenseMatrix rotaMatr(3,3);
                double tempAngl = 2.0 * Ddpca::PI / dehwSurf.z[1] * (double)toot_tw;
                rotaMatr.Fill(
                    std::cos(tempAngl),std::sin(tempAngl),0.0,
                    -std::sin(tempAngl),std::cos(tempAngl),0.0,
                    0.0,0.0,1.0
                );
                Ddpca::Coordinate tranVect(0.0,0.0,0.0);
                twMesh.RigidRotationTranslation(rotaMatr, tranVect);
                //
                twMesh.RigidRotationTranslation(wheeRota, wheeTran);
                if(toot_tw  == 0 && leri_tw == 1){
                    whdeAucu.RigidRotationTranslation(rotaMatr, tranVect);
                    whdeAucu.RigidRotationTranslation(wheeRota, wheeTran);
                }
                if(toot_tw  == 0 && leri_tw == 0){
                    whdeAucu_midd.RigidRotationTranslation(rotaMatr, tranVect);
                    whdeAucu_midd.RigidRotationTranslation(wheeRota, wheeTran);
                }
            }
        });
    }

    void SUBR_COLO_WORM(long tg){
        if(tg == -1){
            long tg_endi = dehwSurf.gridNumb[0][6];
            const Ddpca::BilinearQuadrature<2>& biliQuad = Ddpca::GetBilinearQuadrature();
            //caculate the total area
            double totaArea = 0.0;
            for(long tw = 0; tw < tg_endi; tw ++){
                Ddpca::Mesh& twMesh = domains[tw].mesh;
                Ddpca::I64 twelSize = twMesh.elements.size();
                for(long ti = 0; ti < twelSize; ti ++){
                    if(twMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                    for(long tj = 0; tj < hefaSize; tj ++){
                        bool tempFlag = true;
                        std::vector<Ddpca::Coordinate> elemCoor(4);
                        for(long tk = 0; tk < 4; tk ++){
                            long tempNode = twMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                            auto iterNoco = twMesh.node2Coordinate.find(tempNode);
                            elemCoor[tk] = Ddpca::Coordinate((iterNoco->second)[0], 
                                (iterNoco->second)[1], (iterNoco->second)[2]);
                            double x_loca = (iterNoco->second)[0] + (dehwSurf.a_h2 + centErro);
                            double y_loca = - (iterNoco->second)[2];
                            double tempRadi = std::sqrt(std::pow(x_loca, 2.0) + std::pow(y_loca, 2.0));
                            if(std::abs(tempRadi - dehwSurf.inneRadi[0]) > 1.0E-10){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == false){
                            continue;
                        }
                        double tempArea = 0.0;
                        for(long tk = 0; tk < biliQuad.numbGaussPoints; tk ++){
                            Ddpca::Real tkJacobian = Ddpca::BilinearQuadratureJacobian<2>(
                                biliQuad.gaussPoints[tk], elemCoor);
                            tempArea += biliQuad.weights[tk] * tkJacobian;
                        }
                        totaArea += tempArea;
                    }
                }
            }
            Ddpca::Log("The total area of worm inner hub: " + Ddpca::Double2String(totaArea));
            loadIncr[0] = dehwSurf.inpuTorq / dehwSurf.inneRadi[0] / totaArea;
        }
        else if(!isSelf){
            //rotate the nodal coordinate system
            //NO: node coupling!
            //displacement constraint must at first
            const Ddpca::BilinearQuadrature<2>& biliQuad = Ddpca::GetBilinearQuadrature();
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
            for(const auto& iterNoco : tgMesh.node2Coordinate){
                double x_loca = (iterNoco.second)[0] + (dehwSurf.a_h2 + centErro);
                double y_loca = - (iterNoco.second)[2];
                double tempRadi = std::sqrt(std::pow(x_loca, 2.0) + std::pow(y_loca, 2.0));
                if(std::abs(tempRadi - dehwSurf.inneRadi[0]) <= 1.0E-10){
                    double tempAngl = std::atan2(y_loca, x_loca);
                    Ddpca::DenseMatrix tempRota(3,3,{
                        std::cos(tempAngl),0.0,-std::sin(tempAngl),
                        -std::sin(tempAngl),0.0,-std::cos(tempAngl),
                        0.0,1.0,0.0
                    });
                    tgBoundary.nodeRotation.emplace(iterNoco.first, tempRota);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                }
            }
            //load must at next
            double totaForc = 0.0;
            Ddpca::I64 tgelSize = tgMesh.elements.size();
            for(long ti = 0; ti < tgelSize; ti ++){
                if(tgMesh.elements[ti].children.size() > 0){
                    continue;
                }
                Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                for(long tj = 0; tj < hefaSize; tj ++){
                    bool tempFlag = true;
                    std::vector<long> inpuNode(4);
                    std::vector<Ddpca::Coordinate> elemCoor(4);
                    for(long tk = 0; tk < 4; tk ++){
                        inpuNode[tk] = tgMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                        auto iterNoco = tgMesh.node2Coordinate.find(inpuNode[tk]);
                        elemCoor[tk] = Ddpca::Coordinate((iterNoco->second)[0], 
                            (iterNoco->second)[1], (iterNoco->second)[2]);
                        double x_loca = (iterNoco->second)[0] + (dehwSurf.a_h2 + centErro);
                        double y_loca = - (iterNoco->second)[2];
                        double tempRadi = std::sqrt(std::pow(x_loca, 2.0) + std::pow(y_loca, 2.0));
                        if(std::abs(tempRadi - dehwSurf.inneRadi[0]) > 1.0E-10){
                            tempFlag = false;
                            break;
                        }
                    }
                    if(tempFlag == false){
                        continue;
                    }
                    Ddpca::DenseMatrix N_e(3,12);
                    std::array<Ddpca::Real,12> tempForc, tkF_e;
                    tempForc.fill(0.0);
                    for(long tk = 0; tk < biliQuad.numbGaussPoints; tk ++){
                        Ddpca::Real tkJacobian = Ddpca::BilinearQuadratureJacobian<2>(
                            biliQuad.gaussPoints[tk], elemCoor);
                        std::array<Ddpca::Real,3> tempLoad;
                        tempLoad = {0.0, loadIncr[0], 0.0};
                        // col-major
                        N_e.Fill(
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
                        Ddpca::GEMTV(N_e, tempLoad, tkF_e);
                        Ddpca::AXPY(biliQuad.weights[tk] * tkJacobian, tkF_e, tempForc);
                    }
                    for(long tk = 0; tk < 4; tk ++){
                        tgBoundary.LoadAccumulate(3 * inpuNode[tk] + 1, tempForc[3 * tk + 1]);
                        totaForc += tempForc[3 * tk + 1];
                    }
                }
            }
            Ddpca::Log("The total force of worm inner hub: " + Ddpca::Double2String(totaForc));
        }
        else{
            //rotate the nodal coordinate system
            //NO: node coupling!
            //displacement constraint must at first
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
            for(const auto& iterNoco : tgMesh.node2Coordinate){
                double x_loca = (iterNoco.second)[0] + (dehwSurf.a_h2 + centErro);
                double y_loca = - (iterNoco.second)[2];
                double tempRadi = std::sqrt(std::pow(x_loca, 2.0) + std::pow(y_loca, 2.0));
                if(std::abs(tempRadi - dehwSurf.inneRadi[0]) <= 1.0E-10){
                    double tempAngl = std::atan2(y_loca, x_loca);
                    Ddpca::DenseMatrix tempRota(3,3,{
                        std::cos(tempAngl),0.0,-std::sin(tempAngl),
                        -std::sin(tempAngl),0.0,-std::cos(tempAngl),
                        0.0,1.0,0.0
                    });
                    tgBoundary.nodeRotation.emplace(iterNoco.first, tempRota);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                    tgBoundary.LoadAccumulate(3 * iterNoco.first + 1, 1.0E-10);
                }
            }
        }
    }

    void SUBR_COLO_WHEE(long tg){
        if(tg == -1){
            long tg_star = dehwSurf.gridNumb[0][6];
            Ddpca::I64 domaSize = domains.size();
            const Ddpca::BilinearQuadrature<2>& biliQuad = Ddpca::GetBilinearQuadrature();
            //caculate the total area
            double totaArea = 0.0;
            for(long tw = tg_star; tw < domaSize; tw ++){
                Ddpca::Mesh& twMesh = domains[tw].mesh;
                Ddpca::I64 twelSize = twMesh.elements.size();
                for(long ti = 0; ti < twelSize; ti ++){
                    if(twMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                    for(long tj = 0; tj < hefaSize; tj ++){
                        bool tempFlag = true;
                        std::vector<Ddpca::Coordinate> elemCoor(4);
                        for(long tk = 0; tk < 4; tk ++){
                            long tempNode = twMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                            auto iterNoco = twMesh.node2Coordinate.find(tempNode);
                            elemCoor[tk] = Ddpca::Coordinate((iterNoco->second)[0], 
                                (iterNoco->second)[1], (iterNoco->second)[2]);
                            double x_loca = (iterNoco->second)[0];
                            double y_loca = (iterNoco->second)[1];
                            double tempRadi = std::sqrt(std::pow(x_loca, 2.0) + std::pow(y_loca, 2.0));
                            if(std::abs(tempRadi - dehwSurf.inneRadi[1]) > 1.0E-10){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == false){
                            continue;
                        }
                        double tempArea = 0.0;
                        for(long tk = 0; tk < biliQuad.numbGaussPoints; tk ++){
                            Ddpca::Real tkJacobian = Ddpca::BilinearQuadratureJacobian<2>(
                                biliQuad.gaussPoints[tk], elemCoor);
                            tempArea += biliQuad.weights[tk] * tkJacobian;
                        }
                        totaArea += tempArea;
                    }
                }
            }
            Ddpca::Log("The total area of worm wheel inner hub: " + Ddpca::Double2String(totaArea));
            loadIncr[1] = - dehwSurf.inpuTorq * dehwSurf.i_h2 / dehwSurf.inneRadi[1] / totaArea;
        }
        else if(!isSelf){
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
            //displacement constraint must at first
            for(const auto& iterNoco : tgMesh.node2Coordinate){
                double tempRadi = std::sqrt(std::pow((iterNoco.second)[0], 2.0) 
                    + std::pow((iterNoco.second)[1], 2.0)
                );
                if(std::abs(tempRadi - dehwSurf.inneRadi[1]) <= 1.0E-10){
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 1, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                }
            }
        }
        else{
            //rotate the nodal coordinate system
            //NO: node coupling!
            //displacement constraint must at first
            const Ddpca::BilinearQuadrature<2>& biliQuad = Ddpca::GetBilinearQuadrature();
            Ddpca::Mesh& tgMesh = domains[tg].mesh;
            Ddpca::BoundaryCondition& tgBoundary = domains[tg].boundary;
            for(const auto& iterNoco : tgMesh.node2Coordinate){
                double x_loca = (iterNoco.second)[0];
                double y_loca = (iterNoco.second)[1];
                double tempRadi = std::sqrt(std::pow(x_loca, 2.0) + std::pow(y_loca, 2.0));
                if(std::abs(tempRadi - dehwSurf.inneRadi[1]) <= 1.0E-10){
                    double tempAngl = std::atan2(y_loca, x_loca);
                    Ddpca::DenseMatrix tempRota(3, 3, {
                        std::cos(tempAngl),std::sin(tempAngl),0.0,
                        -std::sin(tempAngl),std::cos(tempAngl),0.0,
                        0.0,0.0,1.0
                    });
                    tgBoundary.nodeRotation.emplace(iterNoco.first, tempRota);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 0, 0.0);
                    tgBoundary.constrainedDof.emplace(3 * iterNoco.first + 2, 0.0);
                }
            }
            //load must at next
            double totaForc = 0.0;
            Ddpca::I64 tgelSize = tgMesh.elements.size();
            for(long ti = 0; ti < tgelSize; ti ++){
                if(tgMesh.elements[ti].children.size() > 0){
                    continue;
                }
                Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                for(long tj = 0; tj < hefaSize; tj ++){
                    bool tempFlag = true;
                    std::vector<long> inpuNode(4);
                    std::vector<Ddpca::Coordinate> elemCoor(4);
                    for(long tk = 0; tk < 4; tk ++){
                        inpuNode[tk] = tgMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                        auto iterNoco = tgMesh.node2Coordinate.find(inpuNode[tk]);
                        elemCoor[tk] = Ddpca::Coordinate((iterNoco->second)[0], 
                            (iterNoco->second)[1], (iterNoco->second)[2]);
                        double x_loca = (iterNoco->second)[0];
                        double y_loca = (iterNoco->second)[1];
                        double tempRadi = std::sqrt(std::pow(x_loca, 2.0) + std::pow(y_loca, 2.0));
                        if(std::abs(tempRadi - dehwSurf.inneRadi[1]) > 1.0E-10){
                            tempFlag = false;
                            break;
                        }
                    }
                    if(tempFlag == false){
                        continue;
                    }
                    Ddpca::DenseMatrix N_e(3,12);
                    std::array<Ddpca::Real,12> tempForc, tkF_e;
                    tempForc.fill(0.0);
                    for(long tk = 0; tk < biliQuad.numbGaussPoints; tk ++){
                        Ddpca::Real tkJacobian = Ddpca::BilinearQuadratureJacobian<2>(
                            biliQuad.gaussPoints[tk], elemCoor);
                        std::array<Ddpca::Real,3> tempLoad;
                        tempLoad = {0.0, loadIncr[1], 0.0};
                        // col-major
                        N_e.Fill(
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
                        Ddpca::GEMTV(N_e, tempLoad, tkF_e);
                        Ddpca::AXPY(biliQuad.weights[tk] * tkJacobian, tkF_e, tempForc);
                    }
                    for(long tk = 0; tk < 4; tk ++){
                        tgBoundary.LoadAccumulate(3 * inpuNode[tk] + 1, tempForc[3 * tk + 1]);
                        totaForc += tempForc[3 * tk + 1];
                    }
                }
            }
            Ddpca::Log("The total force of wheel inner hub: " + Ddpca::Double2String(totaForc));
        }
    }

    void GenerateInterfaces(){
        //
        Ddpca::Real charLeng = CalculateCharacteristicLength();
        //****************************************************************************************
        Ddpca::Log("TestDehw::GenerateInterfaces local mesh refinement");
        auto CART_CURV = [&](Ddpca::Coordinate tempCoor, double &tempXico, double &tempEtac){
            tempXico += tempCoor[1];
            tempEtac += std::sqrt(std::pow(tempCoor[0] + (dehwSurf.a_h2 + centErro), 2.0)
                + std::pow(tempCoor[2], 2.0));
        };
        std::array<std::array<bool,3>,4> isnoRefi;
        std::array<std::array<Ddpca::CurvedSurface,2>,4> cusuTabl_0;
        std::array<std::vector<std::array<Ddpca::I64,2>>,4> contBody_0;
        for(Ddpca::I64 ti = 0; ti < 4; ++ ti){
            contBody_0[ti].resize(isnoRefi[0].size());
        }
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.one2oneS2S, [&](Ddpca::I64 tt){
            if(tt >= 4){
                return;
            }
            contBody_0[tt][0] = {2 + 8 * tt, dehwSurf.gridNumb[0][6] + 6 + 2 * tt};
            contBody_0[tt][1] = {3 + 8 * tt, dehwSurf.gridNumb[0][6] + 6 + 2 * tt};
            contBody_0[tt][2] = {4 + 8 * tt, dehwSurf.gridNumb[0][6] + 6 + 2 * tt};
            //
            Ddpca::DenseMatrix wormRota_0(3,3,{
                std::cos(analAngl[0]),std::sin(analAngl[0]),0.0,
                -std::sin(analAngl[0]),std::cos(analAngl[0]),0.0,
                0.0,0.0,1.0
            });
            Ddpca::DenseMatrix wormRota_1(3,3,{
                1.0,0.0,0.0,
                0.0,0.0,-1.0,
                0.0,1.0,0.0
            });
            Ddpca::DenseMatrix tempMatr = wormRota_1;
            Ddpca::GEMM(tempMatr, wormRota_0, wormRota_1);
            Ddpca::Coordinate wormTran(- (dehwSurf.a_h2 + centErro), 0.0, 0.0);
            Ddpca::CurvedSurface mastSurf = dehwSurf.wormTosu;
            mastSurf.RigidRotationTranslation(wormRota_1, wormTran);
            double wheeAngl_1 = analAngl[1] + 2.0 * Ddpca::PI / dehwSurf.z[1] * (1.0 + tt);
            Ddpca::DenseMatrix wheeRota(3,3,{
                std::cos(wheeAngl_1),std::sin(wheeAngl_1),0.0,
                -std::sin(wheeAngl_1),std::cos(wheeAngl_1),0.0,
                0.0,0.0,1.0
            });
            Ddpca::Coordinate wheeTran(0.0, 0.0, 0.0);
            Ddpca::CurvedSurface slavSurf = dehwSurf.wheeTosu;
            slavSurf.RigidRotationTranslation(wheeRota, wheeTran);
            cusuTabl_0[tt][0] = mastSurf;
            cusuTabl_0[tt][1] = slavSurf;
            //
            for(long tr = 0; tr < dehwSurf.locaLeve; tr ++){
                long buckFact = (1 << 
                    (dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve - 1 
                    - (dehwSurf.locaLeve - tr)));
                std::array<Ddpca::I64, 2> buckNumb = {dehwSurf.gridNumb[0][4] * buckFact, 
                    dehwSurf.gridNumb[0][3] * buckFact};
                const Ddpca::I64 isreSize = isnoRefi[tt].size();
                for(long tc = 0; tc < isreSize; tc ++){
                    Ddpca::ContactInterface searCont_0;
                    searCont_0.AdaptiveRefine(
                        domains[contBody_0[tt][tc][0]].mesh, domains[contBody_0[tt][tc][1]].mesh, 
                        isnoRefi[tt][tc], mastSurf, slavSurf, 
                        dehwSurf.globInho + dehwSurf.globHomo + tr, 
                        distCrit[tr], buckNumb, CART_CURV);
                }
            }
        });
        SUBR_COLO_WORM(-1);
        SUBR_COLO_WHEE(-1);
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tw){
            if(tw < dehwSurf.gridNumb[0][6]){
                SUBR_COLO_WORM(tw);
            }
            else{
                SUBR_COLO_WHEE(tw);
            }
        });
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.domainS2S, [&](Ddpca::I64 tw){
            domains[tw].mesh.OutputMesh(directoryPath, tw);
        });
        //****************************************************************************************
        long totaSear = 0;
        for(long tt = 0; tt < 4; tt ++){
            const Ddpca::I64 isreSize = isnoRefi[tt].size();
            for(long tc = 0; tc < isreSize; tc ++){
                if(isnoRefi[tt][tc] == true){
                    totaSear ++;
                }
            }
        }
        std::vector<std::array<Ddpca::CurvedSurface,2>> cusuTabl(totaSear);
        totaSear += dehwSurf.gridNumb[0][6] - 1;
        totaSear += dehwSurf.gridNumb[0][6] - dehwSurf.circNumb;
        totaSear += dehwSurf.gridNumb[1][5] * (dehwSurf.gridNumb[1][6] - 1);
        totaSear += dehwSurf.gridNumb[1][5] - 1;
        interfaces.resize(totaSear);
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.threadsPerInterface, totaSear, 
            Ddpca::threadManager.interfaceS2S, Ddpca::threadManager.interfaceS2M);
        long searCoun = 0;
        for(long tt = 0; tt < 4; tt ++){
            const Ddpca::I64 isreSize = isnoRefi[tt].size();
            for(long tc = 0; tc < isreSize; tc ++){
                if(isnoRefi[tt][tc] == true){
                    interfaces[searCoun].domainIndex = contBody_0[tt][tc];
                    Ddpca::Real averageElasticity = (domains[interfaces[searCoun].domainIndex[0]].elasticity 
                        + domains[interfaces[searCoun].domainIndex[1]].elasticity) / 2.0;
                    interfaces[searCoun].normPenaPara = averageElasticity / charLeng 
                        * NormalPenaltyCoef;
                    interfaces[searCoun].tangPenaPara = averageElasticity / charLeng * tempTangPenaCoef;
                    interfaces[searCoun].frictionCoefficient = (!isSelf) ? 0.08 : 0.2;
                    cusuTabl[searCoun][0] = cusuTabl_0[tt][0];
                    cusuTabl[searCoun][1] = cusuTabl_0[tt][1];
                    searCoun ++;
                }
            }
        }
        // cusuTabl_0.clear();
        for(long tv = 0; tv < dehwSurf.gridNumb[0][6] - 1; tv ++){
            long ts = searCoun + tv;
            interfaces[ts].domainIndex = {tv, tv + 1};
            Ddpca::Real averageElasticity = (domains[interfaces[ts].domainIndex[0]].elasticity 
                + domains[interfaces[ts].domainIndex[1]].elasticity) / 2.0;
            interfaces[ts].normPenaPara = averageElasticity / charLeng * NormalPenaltyCoef;
            interfaces[ts].tangPenaPara = averageElasticity / charLeng * TangentialPenaltyCoef;
            interfaces[ts].frictionCoefficient = -1.0;
        }
        for(long tv = 0; tv < dehwSurf.gridNumb[0][6] - dehwSurf.circNumb; tv ++){
            long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 + tv;
            interfaces[ts].domainIndex = {tv, tv + dehwSurf.circNumb};
            Ddpca::Real averageElasticity = (domains[interfaces[ts].domainIndex[0]].elasticity 
                + domains[interfaces[ts].domainIndex[1]].elasticity) / 2.0;
            interfaces[ts].normPenaPara = averageElasticity / charLeng * NormalPenaltyCoef;
            interfaces[ts].tangPenaPara = averageElasticity / charLeng * TangentialPenaltyCoef;
            interfaces[ts].frictionCoefficient = -1.0;
        }
        for(long ti = 0; ti < dehwSurf.gridNumb[1][5]; ti ++){
            for(long tj = 0; tj < dehwSurf.gridNumb[1][6] - 1; tj ++){
                long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 
                    + dehwSurf.gridNumb[0][6] - dehwSurf.circNumb
                    + ti * (dehwSurf.gridNumb[1][6] - 1) + tj;
                long tv_0 = dehwSurf.gridNumb[0][6] + ti * dehwSurf.gridNumb[1][6] + tj;
                long tv_1 = tv_0 + 1;
                interfaces[ts].domainIndex = {tv_0, tv_1};
                Ddpca::Real averageElasticity = (domains[interfaces[ts].domainIndex[0]].elasticity 
                    + domains[interfaces[ts].domainIndex[1]].elasticity) / 2.0;
                interfaces[ts].normPenaPara = averageElasticity / charLeng * NormalPenaltyCoef;
                interfaces[ts].tangPenaPara = averageElasticity / charLeng * TangentialPenaltyCoef;
                interfaces[ts].frictionCoefficient = -1.0;
            }
        }
        for(long ti = 0; ti < dehwSurf.gridNumb[1][5] - 1; ti ++){
            long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 
                + dehwSurf.gridNumb[0][6] - dehwSurf.circNumb 
                + dehwSurf.gridNumb[1][5] * (dehwSurf.gridNumb[1][6] - 1) 
                + ti;
            long tv_0 = dehwSurf.gridNumb[0][6] + ti * dehwSurf.gridNumb[1][6] 
                + dehwSurf.gridNumb[1][6] - 1;
            long tv_1 = tv_0 + 1;
            interfaces[ts].domainIndex = {tv_0, tv_1};
            Ddpca::Real averageElasticity = (domains[interfaces[ts].domainIndex[0]].elasticity 
                + domains[interfaces[ts].domainIndex[1]].elasticity) / 2.0;
            interfaces[ts].normPenaPara = averageElasticity / charLeng * NormalPenaltyCoef;
            interfaces[ts].tangPenaPara = averageElasticity / charLeng * TangentialPenaltyCoef;
            interfaces[ts].frictionCoefficient = -1.0;
        }
        //****************************************************************************************
        Ddpca::Log("TestDehw::GenerateInterfaces contact surface");
        long totaLeve = dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve;
        std::array<Ddpca::I64,2> buckNumb_0 = {dehwSurf.gridNumb[0][4] * (1 << (totaLeve - 1)), 
            dehwSurf.gridNumb[0][3] * (1 << (totaLeve - 1))
        };
        //column-major order
        Ddpca::DenseMatrix wormRota_0(3,3,{
            std::cos(analAngl[0]),-std::sin(analAngl[0]),0.0,
            std::sin(analAngl[0]),std::cos(analAngl[0]),0.0,
            0.0,0.0,1.0
        });
        Ddpca::DenseMatrix wormRota_1(3,3,{
            1.0,0.0,0.0,
            0.0,0.0,1.0,
            0.0,-1.0,0.0
        });
        Ddpca::DenseMatrix tempMatr = wormRota_0;
        Ddpca::GEMM(tempMatr, wormRota_1, wormRota_0);
        std::array<Ddpca::Real,3> wormTran;
        wormTran = {+ (dehwSurf.a_h2 + centErro), 0.0, 0.0};
        auto CART_CURV_1 = [&](Ddpca::Coordinate tempCoor, double tempAngl) -> bool {
            std::array<Ddpca::Real,3> tempVect, coorVect;
            tempVect = {tempCoor[0] + wormTran[0], tempCoor[1] + wormTran[1], tempCoor[2] + wormTran[2]};
            Ddpca::GEMV(wormRota_0, tempVect, coorVect);
            double coorAngl = std::atan2(coorVect[1], coorVect[0]);
            //!!!plus, but not minus!!!
            if(std::abs(coorAngl + tempAngl) < 1.0E-10){
                return true;
            }
            else{
                return false;
            }
        };
        std::array<Ddpca::I64,2> buckNumb_1 = {dehwSurf.gridNumb[0][0] * (1 << (dehwSurf.globHomo - 1)), 
            dehwSurf.gridNumb[0][1] * (1 << (dehwSurf.globHomo))
        };
        Ddpca::threadManager.RunTask(0, Ddpca::threadManager.interfaceS2S, [&](Ddpca::I64 ts){
            if(ts < searCoun){
                //contact surface
                Ddpca::CurvedSurface& mastSurf = cusuTabl[ts][0];
                const Ddpca::Mesh& masterMesh = domains[interfaces[ts].domainIndex[0]].mesh;
                mastSurf.Initialize();
                while(mastSurf.Increment(masterMesh)){
                    //NO NEED: if(multGrid[contBody[ts][0]].elemVect[iterEfsu_0.eid].level == totaLeve){
                    interfaces[ts].masterSegments.emplace_back(mastSurf.currentFace);
                }
                Ddpca::CurvedSurface& slavSurf = cusuTabl[ts][1];
                const Ddpca::Mesh& slaveMesh = domains[interfaces[ts].domainIndex[1]].mesh;
                slavSurf.Initialize();
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
                Ddpca::Real tempXi = 0.0, tempEta = 0.0;
                Ddpca::Coordinate tempCoor;
                for(long ti = 0; ti < inmaSize; ti ++){
                    tempXi = 0.0;
                    tempEta = 0.0;
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].masterSegments[ti][tj];
                        auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                        CART_CURV(iterNoco->second, tempXi, tempEta);
                    }
                    masterLocal[0][ti] = tempXi / 4.0;
                    masterLocal[1][ti] = tempEta / 4.0;
                }
                interfaces[ts].BucketSort(masterLocal, buckNumb_0);
                Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
                std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
                for(long ti = 0; ti < inslSize; ti ++){
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].slaveSegments[ti][tj];
                        auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                        tempXi = 0.0;
                        tempEta = 0.0;
                        CART_CURV(iterNoco->second, tempXi, tempEta);
                        slaveLocal[ti][tj * 2 + 0] = tempXi;
                        slaveLocal[ti][tj * 2 + 1] = tempEta;
                    }
                }
                interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                    Ddpca::threadManager.interfaceS2M[ts], 1, distCrit[dehwSurf.locaLeve - 1]);
                interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            }
            else if(ts < searCoun + dehwSurf.gridNumb[0][6] - 1){
                //surfaces between domains of worm tooth
                long tv = ts - searCoun;
                long tg_mast = interfaces[ts].domainIndex[0];
                long tg_slav = interfaces[ts].domainIndex[1];
                const Ddpca::Mesh& masterMesh = domains[tg_mast].mesh;
                const Ddpca::Mesh& slaveMesh = domains[tg_slav].mesh;
                Ddpca::I64 maelSize = masterMesh.elements.size();
                Ddpca::I64 slelSize = slaveMesh.elements.size();
                std::array<Ddpca::I64, 4> tempNode;
                const Ddpca::I64 hefaSize = Ddpca::hexaFace.size();
                for(long ti = 0; ti < maelSize; ti ++){
                    if(masterMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    for(long tj = 0; tj < hefaSize; tj ++){
                        bool tempFlag = true;
                        for(long tk = 0; tk < 4; tk ++){
                            tempNode[tk] = masterMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                            auto iterNoco = masterMesh.node2Coordinate.find(tempNode[tk]);
                            if(CART_CURV_1(iterNoco->second, wodeAuan[tv]) == false){
                                tempFlag = false;
                                break;
                            }
                        }
                        if(tempFlag == true){
                            interfaces[ts].masterSegments.emplace_back(tempNode);
                        }
                    }
                }
                for(long ti = 0; ti < slelSize; ti ++){
                    if(slaveMesh.elements[ti].children.size() > 0){
                        continue;
                    }
                    for(long tj = 0; tj < hefaSize; tj ++){
                        bool tempFlag = true;
                        for(long tk = 0; tk < 4; tk ++){
                            tempNode[tk] = slaveMesh.elements[ti].cornerNodes[Ddpca::hexaFace[tj][tk]];
                            auto iterNoco = slaveMesh.node2Coordinate.find(tempNode[tk]);
                            if(CART_CURV_1(iterNoco->second, wodeAuan[tv]) == false){
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
                for(long ti = 0; ti < inmaSize; ti ++){
                    tempXi = 0.0;
                    tempEta = 0.0;
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].masterSegments[ti][tj];
                        auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                        CART_CURV(iterNoco->second, tempXi, tempEta);
                    }
                    masterLocal[0][ti] = tempXi / 4.0;
                    masterLocal[1][ti] = tempEta / 4.0;
                }
                interfaces[ts].BucketSort(masterLocal, buckNumb_1);
                Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
                std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
                for(long ti = 0; ti < inslSize; ti ++){
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].slaveSegments[ti][tj];
                        auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                        tempXi = 0.0;
                        tempEta = 0.0;
                        CART_CURV(iterNoco->second, tempXi, tempEta);
                        slaveLocal[ti][tj * 2 + 0] = tempXi;
                        slaveLocal[ti][tj * 2 + 1] = tempEta;
                    }
                }
                interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                    Ddpca::threadManager.interfaceS2M[ts], 1, 1.0E12);
                interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            }
            else if(ts < searCoun + dehwSurf.gridNumb[0][6] - 1 + dehwSurf.gridNumb[0][6] - dehwSurf.circNumb){
                //surfaces between worm tooth
                long tv = ts - (searCoun + dehwSurf.gridNumb[0][6] - 1);
                std::array<Ddpca::I64,2> buckNumb;
                if(tv == 0){
                    buckNumb = {dehwSurf.gridNumb[0][1] * (1 << (dehwSurf.globHomo)), 
                        dehwSurf.gridNumb[0][5] * (1 << (dehwSurf.globInho + dehwSurf.globHomo - 1))
                    };
                }
                else{
                    buckNumb = {dehwSurf.gridNumb[0][1] * (1 << (dehwSurf.globHomo)), 
                        dehwSurf.gridNumb[0][4] * (1 << (dehwSurf.globInho + dehwSurf.globHomo - 1))
                    };
                }
                //
                long tg_mast = interfaces[ts].domainIndex[0];
                long tg_slav = interfaces[ts].domainIndex[1];
                Ddpca::CurvedSurface& mastSurf = wodeAucu[tg_mast][0];
                const Ddpca::Mesh& masterMesh = domains[tg_mast].mesh;
                mastSurf.Initialize();
                while(mastSurf.Increment(masterMesh)){
                    interfaces[ts].masterSegments.emplace_back(mastSurf.currentFace);
                }
                Ddpca::CurvedSurface& slavSurf = wodeAucu[tg_slav][1];
                const Ddpca::Mesh& slaveMesh = domains[tg_slav].mesh;
                slavSurf.Initialize();
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
                Ddpca::Real tempXi = 0.0, tempEta = 0.0;
                Ddpca::Coordinate tempCoor;
                for(long ti = 0; ti < inmaSize; ti ++){
                    tempXi = 0.0;
                    tempEta = 0.0;
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].masterSegments[ti][tj];
                        auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                        double tempX = (iterNoco->second)[0] + (dehwSurf.a_h2 + centErro);
                        double tempY = - (iterNoco->second)[2];
                        tempXi += std::sqrt(std::pow(tempX, 2.0) + std::pow(tempY, 2.0));
                        tempEta += (iterNoco->second)[1];
                    }
                    masterLocal[0][ti] = tempXi / 4.0;
                    masterLocal[1][ti] = tempEta / 4.0;
                }
                interfaces[ts].BucketSort(masterLocal, buckNumb);
                Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
                std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
                for(long ti = 0; ti < inslSize; ti ++){
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].slaveSegments[ti][tj];
                        auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                        double tempX = (iterNoco->second)[0] + (dehwSurf.a_h2 + centErro);
                        double tempY = - (iterNoco->second)[2];
                        slaveLocal[ti][tj * 2 + 0] = std::sqrt(std::pow(tempX, 2.0) + std::pow(tempY, 2.0));
                        slaveLocal[ti][tj * 2 + 1] = (iterNoco->second)[1];
                    }
                }
                interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                    Ddpca::threadManager.interfaceS2M[ts], 1, distCrit[dehwSurf.locaLeve - 1]);
                interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            }
            else if(ts < searCoun + dehwSurf.gridNumb[0][6] - 1 
                + dehwSurf.gridNumb[0][6] - dehwSurf.circNumb 
                + dehwSurf.gridNumb[1][5] * (dehwSurf.gridNumb[1][6] - 1)){
                //surfaces between domains of wheel tooth
                Ddpca::I64 ts_ = ts - (searCoun + dehwSurf.gridNumb[0][6] - 1 
                    + dehwSurf.gridNumb[0][6] - dehwSurf.circNumb);
                Ddpca::I64 tv_0 = ts_ / (dehwSurf.gridNumb[1][6] - 1);
                // Ddpca::I64 tv_1 = ts_ % (dehwSurf.gridNumb[1][6] - 1);
                std::array<Ddpca::I64,2> buckNumb = {(dehwSurf.gridNumb[1][1] + dehwSurf.gridNumb[1][3]) 
                    * (1 << (dehwSurf.globHomo - 1)), 
                    dehwSurf.gridNumb[1][4] * (1 << (dehwSurf.globInho + dehwSurf.globHomo - 1))};
                //
                Ddpca::CurvedSurface tempCusu = whdeAucu_midd;
                double tempAngl = 2.0 * Ddpca::PI / dehwSurf.z[1] * (double)tv_0;
                Ddpca::DenseMatrix rotaMatr(3,3,{
                    std::cos(tempAngl),std::sin(tempAngl),0.0,
                    -std::sin(tempAngl),std::cos(tempAngl),0.0,
                    0.0,0.0,1.0
                });
                Ddpca::Coordinate tranVect(0.0,0.0,0.0);
                tempCusu.RigidRotationTranslation(rotaMatr, tranVect);
                //
                long tg_mast = interfaces[ts].domainIndex[0];
                long tg_slav = interfaces[ts].domainIndex[1];
                const Ddpca::Mesh& masterMesh = domains[tg_mast].mesh;
                tempCusu.Initialize();
                while(tempCusu.Increment(masterMesh)){
                    interfaces[ts].masterSegments.emplace_back(tempCusu.currentFace);
                }
                const Ddpca::Mesh& slaveMesh = domains[tg_slav].mesh;
                tempCusu.Initialize();
                while(tempCusu.Increment(slaveMesh)){
                    interfaces[ts].slaveSegments.emplace_back(tempCusu.currentFace);
                }
                interfaces[ts].OutputSegments(directoryPath, ts);
                //
                Ddpca::I64 inmaSize = interfaces[ts].masterSegments.size();
                std::array<std::vector<Ddpca::Real>, 2> masterLocal;
                masterLocal[0].resize(inmaSize);
                masterLocal[1].resize(inmaSize);
                Ddpca::I64 node_tj;
                Ddpca::Real tempXi = 0.0, tempEta = 0.0;
                Ddpca::Coordinate tempCoor;
                for(long ti = 0; ti < inmaSize; ti ++){
                    tempXi = 0.0;
                    tempEta = 0.0;
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].masterSegments[ti][tj];
                        auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                        tempXi += std::sqrt(std::pow((iterNoco->second)[0], 2.0) 
                            + std::pow((iterNoco->second)[1], 2.0));
                        tempEta += (iterNoco->second)[2];
                    }
                    masterLocal[0][ti] = tempXi / 4.0;
                    masterLocal[1][ti] = tempEta / 4.0;
                }
                interfaces[ts].BucketSort(masterLocal, buckNumb);
                Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
                std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
                for(long ti = 0; ti < inslSize; ti ++){
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].slaveSegments[ti][tj];
                        auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                        slaveLocal[ti][tj * 2 + 0] = std::sqrt(std::pow((iterNoco->second)[0], 2.0) 
                            + std::pow((iterNoco->second)[1], 2.0));
                        slaveLocal[ti][tj * 2 + 1] = (iterNoco->second)[2];
                    }
                }
                interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                    Ddpca::threadManager.interfaceS2M[ts], 1, distCrit[dehwSurf.locaLeve - 1]);
                interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            }
            else{
                //surfaces between wheel teeth
                std::array<Ddpca::I64, 2> buckNumb = {dehwSurf.gridNumb[1][1] * (1 << (dehwSurf.globHomo)), 
                    dehwSurf.gridNumb[1][4] * (1 << (dehwSurf.globInho + dehwSurf.globHomo - 1))};
                Ddpca::I64 tv_0 = ts - (searCoun + dehwSurf.gridNumb[0][6] - 1 
                    + dehwSurf.gridNumb[0][6] - dehwSurf.circNumb 
                    + dehwSurf.gridNumb[1][5] * (dehwSurf.gridNumb[1][6] - 1));
                //
                Ddpca::CurvedSurface tempCusu = whdeAucu;
                double tempAngl = 2.0 * Ddpca::PI / dehwSurf.z[1] * (double)tv_0;
                Ddpca::DenseMatrix rotaMatr(3,3,{
                    std::cos(tempAngl),std::sin(tempAngl),0.0,
                    -std::sin(tempAngl),std::cos(tempAngl),0.0,
                    0.0,0.0,1.0
                });
                Ddpca::Coordinate tranVect(0.0,0.0,0.0);
                tempCusu.RigidRotationTranslation(rotaMatr, tranVect);
                //
                long tg_mast = interfaces[ts].domainIndex[0];
                long tg_slav = interfaces[ts].domainIndex[1];
                const Ddpca::Mesh& masterMesh = domains[tg_mast].mesh;
                tempCusu.Initialize();
                while(tempCusu.Increment(masterMesh)){
                    interfaces[ts].masterSegments.emplace_back(tempCusu.currentFace);
                }
                const Ddpca::Mesh& slaveMesh = domains[tg_slav].mesh;
                tempCusu.Initialize();
                while(tempCusu.Increment(slaveMesh)){
                    interfaces[ts].slaveSegments.emplace_back(tempCusu.currentFace);
                }
                interfaces[ts].OutputSegments(directoryPath, ts);
                //
                Ddpca::I64 inmaSize = interfaces[ts].masterSegments.size();
                std::array<std::vector<Ddpca::Real>, 2> masterLocal;
                masterLocal[0].resize(inmaSize);
                masterLocal[1].resize(inmaSize);
                Ddpca::I64 node_tj;
                Ddpca::Real tempXi = 0.0, tempEta = 0.0;
                Ddpca::Coordinate tempCoor;
                for(long ti = 0; ti < inmaSize; ti ++){
                    tempXi = 0.0;
                    tempEta = 0.0;
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].masterSegments[ti][tj];
                        auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
                        tempXi += std::sqrt(std::pow((iterNoco->second)[0], 2.0) 
                            + std::pow((iterNoco->second)[1], 2.0));
                        tempEta += (iterNoco->second)[2];
                    }
                    masterLocal[0][ti] = tempXi / 4.0;
                    masterLocal[1][ti] = tempEta / 4.0;
                }
                interfaces[ts].BucketSort(masterLocal, buckNumb);
                Ddpca::I64 inslSize = interfaces[ts].slaveSegments.size();
                std::vector<std::array<Ddpca::Real, 8>> slaveLocal(inslSize);
                for(long ti = 0; ti < inslSize; ti ++){
                    for(long tj = 0; tj < 4; tj ++){
                        node_tj = interfaces[ts].slaveSegments[ti][tj];
                        auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
                        slaveLocal[ti][tj * 2 + 0] = std::sqrt(std::pow((iterNoco->second)[0], 2.0) 
                            + std::pow((iterNoco->second)[1], 2.0));
                        slaveLocal[ti][tj * 2 + 1] = (iterNoco->second)[2];
                    }
                }
                interfaces[ts].LocalSearch(masterMesh, slaveMesh, slaveLocal, 
                    Ddpca::threadManager.interfaceS2M[ts], 1, distCrit[dehwSurf.locaLeve - 1]);
                interfaces[ts].OutputIntegralPoints(directoryPath, ts);
            }
        });
    }

    void Test(){
        //
	    dehwSurf.ESTABLISH(directoryPath);
        isSelf = false;
        tempTangPenaCoef = isSelf ? 2.5 : 0.25;
        centErro = 0.0E-6;
        analAngl = {0.0, 0.0};//must be zero
        if(!isSelf){
            distCrit = {55.0E-6, 35.0E-6, 15.0E-6};
        }
        else{
            distCrit = {65.0E-6, 45.0E-6, 25.0E-6};
        }
        //
        const Ddpca::I64 numbDomains = dehwSurf.gridNumb[0][6] + dehwSurf.gridNumb[1][5] * dehwSurf.gridNumb[1][6];
        domains.resize(numbDomains);
        Ddpca::threadManager.ThreadDistribute(
            Ddpca::threadManager.threadsPerDomain, domains.size(), 
            Ddpca::threadManager.domainS2S, Ddpca::threadManager.domainS2M);
        GenerateMesh();
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
    Dehw dehw;
    dehw.directoryPath = "./TestDehw_";
    std::filesystem::create_directory(dehw.directoryPath);
    dehw.Test();
    //
	Ddpca::Finalize();
	return 1;
}