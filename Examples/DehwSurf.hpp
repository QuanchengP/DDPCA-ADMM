#ifndef _DehwSurf_hpp
#define _DehwSurf_hpp

#include "../Contact/CurvedSurface.hpp"

#include <list>

class DehwSurf{

public:

	//****************************************Input-1****************************************
	DehwSurf();
	//****************************************Input-2****************************************
	Ddpca::I64 z[2];//teeth number
	Ddpca::Real a_h2;//working center distance
	Ddpca::Real modiTran;//modification of transmission ratio
	Ddpca::Real modiCent;//modification of center distance
	Ddpca::Real r_b2;//radius of base circle of worm wheel
	Ddpca::Real beta_c;//inclination angle of generating plane
	//****************************************Input-3****************************************
	Ddpca::Real z_k;//encircled teeth
	Ddpca::Real d[2];//reference circle diameter at throat
	Ddpca::Real h_a_s[2];//addendum coefficient
	Ddpca::Real h_f_s[2];//dedendum coefficient
	Ddpca::Real R_a[2];//tip arc
	Ddpca::Real offsR_a;//offset of tip arc center of worm wheel
	Ddpca::Real wheeWidt;//width of worm wheel
	Ddpca::Real inneRadi[2];//inner rdius of hub
	Ddpca::Real inpuTorq;//input torque
	//****************************************Input-4****************************************
	Ddpca::I64 globInho;//global inhomogeneous level
	Ddpca::I64 globHomo;//global homogeneous level
	Ddpca::I64 locaLeve;//local refinement level
	Ddpca::I64 reliSwit;//0 - no tooth flank relif, 1 - tooth flank relief
	std::array<std::array<Ddpca::I64,7>,2> gridNumb;//number of grid division
    
public:

	//***************************************************************************************
	Ddpca::Real a_1c;//center distance of primary enveloping
	Ddpca::Real i_1c;//transmission ratio of primary enveloping
	Ddpca::Real i_c1;
	Ddpca::Real i_h2;//transmission ratio of second enveloping
	Ddpca::Real i_2h;
	Ddpca::Real m_t;//transverse module
	Ddpca::Real h_a[2];//addendum
	Ddpca::Real h_f[2];//dedendum
	Ddpca::Real d_f[2];//root circle
	Ddpca::Real d_a[2];//tip circle
	Ddpca::Real R_f[2];//root arc
	Ddpca::Real R_t[2];//transition arc
	Ddpca::Real alph;//nominal pressure angle
	Ddpca::Real leadAngl;//nominal lead angle
	Ddpca::Real pitcAngl;//pitch angle
	Ddpca::Real tootThicCoef[2];//tooth thickness coefficient
	Ddpca::Real halfAngl;//half working angle
	Ddpca::Real starAngl;//starting angle
	Ddpca::Real termAngl;//terminating angle
	Ddpca::Real wormCurv[3];//curvilinear coordinate of worm
	Ddpca::Real widtAngl;//face width angle of worm wheel
	Ddpca::Real backlash;//backlash
	//axial tooth thickness of worm, transverse tooth thickness of worm wheel
	Ddpca::Real tootThic[2];
	Ddpca::Real tootThicAngl[2];//tooth thickness angle
	//simultaneous envelope of tooth surface and tooth back
	//the angle between the two rigidly connected coordinates
	Ddpca::Real backAngl[2];
	void BASIC_PARAMETER();//calculate basic parameter
    
public:

	//***************************************************************************************
	//singular thet_c to thet_h
	void SINGULAR_C2H(Ddpca::Real thet_c, Ddpca::Real &thet_hs, Ddpca::Real &thet_hm);
	//first and second meshing equations
	void FSME(Ddpca::Real thet_1, Ddpca::Real thet_h, Ddpca::Real &x_d, Ddpca::Real &y_d);
	//partial derivative of meshing equations
	void PD_FSME(Ddpca::Real thet_1, Ddpca::Real thet_h, 
		Ddpca::Real &x_d, Ddpca::Real &y_d, Ddpca::DenseMatrix &Pxy_d);
	//x_d, y_d, thet_c to r_1
	void WORM_DC2R(Ddpca::Real x_d, Ddpca::Real y_d, Ddpca::Real thet_c, std::array<Ddpca::Real,3> &r_1_1);
	//x_d, y_d, thet_1, thet_h to r_2
	void WHEE_1H2R(Ddpca::Real x_d, Ddpca::Real y_d, Ddpca::Real thet_1, Ddpca::Real thet_h, 
		std::array<Ddpca::Real,3> &r_2_2);
	//partial derivative of r_2_2 relative to thet_1, thet_h
	void PD_WHEE_1H2R(Ddpca::Real x_d, Ddpca::Real y_d, Ddpca::Real thet_1, Ddpca::Real thet_h, 
		std::array<Ddpca::Real,3> &r_2_2, Ddpca::DenseMatrix &Dr_2_2);
	//curvature interference limit function of first envelope
	void CILFOFE(Ddpca::Real thet_1, Ddpca::Real x_d, Ddpca::Real y_d,
		Ddpca::Real &Psi_1, Ddpca::Real &kapp_cxd, Ddpca::Real &kapp_cyd, Ddpca::Real &tau_cxd);
	//curvature interference limit function of second envelope, non-interval version
	Ddpca::Real CILFOSE_NI(Ddpca::Real thet_1, Ddpca::Real thet_h, Ddpca::Real &kapp_h2N);
	//***************************************************************************************
	//worm tooth surface: curvilinear coordinate to Cartesian coordinate
	void WORM_CURV_2_CART(Ddpca::Real xi_11, Ddpca::Real xi_12, std::array<Ddpca::Real,3> &r_1_1, Ddpca::Real &thet_c);
	//global coordinate r_2_2 to local coordinate of worm wheel
	void WHEE_G2L(std::array<Ddpca::Real,3> r_2_2, Ddpca::Real &angl_f, Ddpca::Real &radi_f, 
		Ddpca::Real &R_fmini, Ddpca::Real &R_fmaxi);
	//worm wheel tooth surface: curvilinear coordinate to Cartesian coordinate
	//thet_c, thet_h are initial values, new contact zone
	void WHEE_CURV_2_CART_1(Ddpca::Real xi_21, Ddpca::Real xi_22, std::array<Ddpca::Real,3> &r_2_2, 
		Ddpca::Real &thet_c, Ddpca::Real &thet_h, Ddpca::I64 f_lr, Ddpca::Real &x_d, Ddpca::Real &y_d);
	//x_d, y_d are initial values, former contact zone
	void WHEE_CURV_2_CART_2(Ddpca::Real xi_21, Ddpca::Real xi_22, std::array<Ddpca::Real,3> &r_c_c, 
		Ddpca::Real &thet_c, Ddpca::Real &x_d, Ddpca::Real &y_d);
	//transition zone
	void WHEE_CURV_2_CART_3(Ddpca::Real xi_21, Ddpca::Real xi_22, std::array<Ddpca::Real,3> &r_2_2, 
		Ddpca::Real &thet_c, Ddpca::Real &thet_h, Ddpca::Real xi_11);
	//transition zone of worm wheel, xi_11 - head transition zone, rear transition zone
	void WHEE_TRAN(Ddpca::Real thet_c, Ddpca::Real thet_h, Ddpca::Real xi_11, 
		std::array<Ddpca::Real,3> &r_2_2, Ddpca::DenseMatrix &Dr_2_2);
	//phase analysis of worm wheel tooth surface
	void WHEE_PHAS(Ddpca::I64 ti, Ddpca::I64 tj, Ddpca::I64 f_ij, std::array<Ddpca::Real,3> r_2_2);
	//new contact zone, 1 - left, 2 - right
	void NEW_CONT_ZONE(Ddpca::I64 f_lr);
	void FORMER_CONT_ZONE();//former contact zone
	void TRANSITION_ZONE(Ddpca::I64 f_hr);//transition zone, 1 - head, 2 - rear
	//tooth flank relief of worm
	void WORM_RELI(std::array<Ddpca::Real,3> &tempXYZ, Ddpca::I64 ti, Ddpca::I64 tj);
	//tooth flank relief of worm wheel
	void WHEE_RELI(std::array<Ddpca::Real,3> &tempXYZ, Ddpca::I64 ti, Ddpca::I64 tj);
	// //***************************************************************************************
	//curvilinear coordinate of tooth profile
	std::vector<std::vector<std::vector<std::array<Ddpca::Real,2>>>> curvCoor;
	//Cartesian coordinate of tooth profile
	std::vector<std::vector<std::vector<std::array<Ddpca::Real,3>>>> cartCoor;
	//flag of worm wheel tooth surface
	//1 - left new contact zone, 2 - right new contact zone, 3 - former contact zone
	//4 - head transition zone,  5 - rear transition zone,   0 - fail to solve
	std::vector<std::vector<Ddpca::I64>> fpha;
	Ddpca::CurvedSurface wormTosu;//worm tooth surface
	Ddpca::CurvedSurface wormToba;//worm tooth back surface
	Ddpca::CurvedSurface wormRtsu;//worm root transition surface
	Ddpca::CurvedSurface wormRtba;//worm root transition back
	Ddpca::CurvedSurface wheeTosu;//worm wheel tooth surface
	Ddpca::CurvedSurface wheeToba;//worm wheel tooth back surface
	Ddpca::CurvedSurface wheeRtsu;//worm wheel root transition surface
	Ddpca::CurvedSurface wheeRtba;//worm wheel root transition back
	Ddpca::I64 circNumb;//number of worm domains in 2 * PI
	void WORM_TS_GRID();//grid discretization of worm tooth surface
	void WHEE_TS_GRID();//grid discretization of worm wheel tooth surface
	void TOOT_SURF_GRID();//grid discretization of tooth surface
	void ROOT_TRAN_GRID();//grid discretization of root transition surface
	// //***************************************************************************************
	//calculate the radius of root transition arc
	void WORM_ROOT_RADIUS(Ddpca::I64 flag, Ddpca::DenseMatrix tempPoin, 
		std::array<Ddpca::Real,2> &tempCent, Ddpca::Real &tempRadi, std::array<Ddpca::Real,2> &tempAngl
	);
	//calculate nodes on root transition arc
	void WORM_ROOT(Ddpca::I64 indi, Ddpca::I64 flag, Ddpca::DenseMatrix &rootProf);
	std::array<Ddpca::Real,2> WHEE_UNCONE(std::array<Ddpca::Real,3> tempXYZ, Ddpca::Real tempAlph_3);
	//transfer from "in unfolded cone surface" to r_2_2
	std::array<Ddpca::Real,3> WHEE_CONE(std::array<Ddpca::Real,2> tempXY, Ddpca::Real tempAlph_3);
	//calculate nodes on root transition arc "in unfolded cone surface"
	void WHEE_ROOT(Ddpca::I64 indi, Ddpca::I64 flag, Ddpca::DenseMatrix &rootProf);
	void OUTPUT(std::string directoryPath);//output wormTosu~wheeRtba to files
	void ESTABLISH(std::string directoryPath);

}; //class DehwSurf

class INDE_INIT{
public:
	Ddpca::I64 ti;
	Ddpca::I64 tj;
	Ddpca::Real init_1;
	Ddpca::Real init_2;
	INDE_INIT(Ddpca::I64 inde_1, Ddpca::I64 inde_2, Ddpca::Real valu_1, Ddpca::Real valu_2)
		: ti(inde_1), tj(inde_2), init_1(valu_1), init_2(valu_2){}
};

DehwSurf::DehwSurf(){
	//
	z[0] = 1;//must be 1
	z[1] = 40;
	a_h2 = 0.25;//m
	modiTran = 0.0;
	modiCent = 0.0;
	r_b2 = 0.158/2.0;//m
	beta_c = 11.0 * Ddpca::PI / 180.0;//rad
	//
	z_k = 4.2;
	d[0] = 0.082;//m
	h_a_s[0] = 0.6;//smaller
	h_a_s[1] = 0.7;//unused
	h_f_s[0] = 0.95;//larger
	h_f_s[1] = 1.05;//larger
	R_a[1] = 0.0385;//m
	offsR_a = 0.003;//m
	wheeWidt = 0.06;//m
	inneRadi[0] = 0.018;//m
	inneRadi[1] = 0.15;//m
	inpuTorq = 180.0;//N*m
	//
	//0 - hub width(even), 1 - hub height, 
	//2 - half tooth width, 3 - tooth height, 
	//4 - number of xi_11 per block (=PI/2)/facewidth, 
	//5 - number of xi_11 for two ending blocks(<PI/2)/tooth number, 
	//6 - block number/block number along facewidth
	//gridNumb[1][4] is divisible by gridNumb[1][6]; gridNumb[1][6] must be 2
	gridNumb = {{ {{4, 2, 2, 4, 4, 0, 0}}, {{4, 4, 2, 4, 8, 8 + z[0], 2}} }};
	globInho = 1;//direction of xi_11/facewidth
	globHomo = 2;//must >= 1
	locaLeve = 3;
	reliSwit = 1;
	circNumb = 8;
}

void DehwSurf::BASIC_PARAMETER(){
	//
	a_1c = a_h2 + modiCent;
	i_h2 = (Ddpca::Real)(z[1]) / z[0];
	i_1c = i_h2 + modiTran;
	i_c1 = 1.0 / i_1c;
	i_2h = 1.0 / i_h2;
	//
	d[1] = 2.0 * a_h2 - d[0];
	m_t = d[1] / z[1];
	h_a[0] = h_a_s[0] * m_t;
	h_a[1] = h_a_s[1] * m_t;
	h_f[0] = h_f_s[0] * m_t;
	h_f[1] = h_f_s[1] * m_t;
	d_f[0] = d[0] - 2.0 * h_f[0];
	d_f[1] = d[1] - 2.0 * h_f[1];
	d_a[0] = d[0] + 2.0 * h_a[0];
	d_a[1] = d[1] + 2.0 * h_a[1];
	R_a[0] = a_h2 - 0.5 * d_a[0];
	R_f[0] = a_h2 - 0.5 * d_f[0];
	R_f[1] = a_h2 - 0.5 * d_f[1];
	R_t[0] = a_h2 - 0.5 * d[0] + 0.8 * m_t;//larger
	R_t[1] = a_h2 - 0.5 * d[1] + 0.9 * m_t;//larger
	//
	alph = std::asin(2.0 * r_b2 / d[1]);
	leadAngl = std::atan(d[1] / i_h2 / d[0]);
	pitcAngl = 2.0 * Ddpca::PI / z[1];
	tootThicCoef[0] = 0.45;
	tootThicCoef[1] = 0.55;
	halfAngl = 0.5 * (z_k - tootThicCoef[0]) * pitcAngl;
	starAngl = alph - halfAngl;
	termAngl = starAngl + z_k * pitcAngl;
	wormCurv[0] = i_h2 * starAngl;
	wormCurv[2] = i_h2 * termAngl;
	wormCurv[1] = (wormCurv[0] + wormCurv[2]) / 2.0;
	while(wormCurv[1] - 2.0 * Ddpca::PI >= wormCurv[0]){
		wormCurv[1] = wormCurv[1] - 2.0 * Ddpca::PI;
	}
	//
	widtAngl = std::asin(wheeWidt / 2.0 / R_f[1]);
	backlash = 0.0;
	tootThic[0] = tootThicCoef[0] * Ddpca::PI * m_t - backlash;
	tootThic[1] = tootThicCoef[1] * Ddpca::PI * m_t;
	tootThicAngl[0] = tootThic[0] / (d[1] / 2.0);
	tootThicAngl[1] = tootThic[1] / (d[1] / 2.0);
	//figure 3.11 in [Zhou, L. Modification principle and manufacturing technology for hourglass 
	//worm drives (National University of Defense Technology Press, 2005)].
	backAngl[0] = 2.0 * alph + tootThicAngl[0];
	backAngl[1] = 2.0 * alph - tootThicAngl[1];
}

void DehwSurf::SINGULAR_C2H(Ddpca::Real thet_c, Ddpca::Real &thet_hs, Ddpca::Real &thet_hm){
	//0 < thet_c < PI / 2.0
	Ddpca::Real thet_1 = i_1c * thet_c;
	Ddpca::Real C_m11 = - i_2h * std::cos(beta_c) * std::sin(thet_c);
	Ddpca::Real C_m12 = i_c1 * i_2h * std::cos(beta_c) * std::cos(thet_c)
		+ i_2h * std::sin(beta_c);
	Ddpca::Real C_m13 = i_c1 * std::cos(beta_c) * std::sin(thet_c);
	Ddpca::Real a2CC = std::atan2(C_m11, C_m12);
	if(C_m13 > std::sqrt(C_m11 * C_m11 + C_m12 * C_m12)){
		thet_hs = thet_1 - a2CC - Ddpca::PI / 2.0;
		thet_hm = thet_hs;
	}
	else{
		thet_hs = thet_1 - Ddpca::PI - a2CC + std::asin(C_m13 / std::sqrt(C_m11 * C_m11 + C_m12 * C_m12));
		thet_hm = thet_1 - a2CC - std::asin(C_m13 / std::sqrt(C_m11 * C_m11 + C_m12 * C_m12));
	}
}

void DehwSurf::FSME(Ddpca::Real thet_1, Ddpca::Real thet_h, Ddpca::Real &x_d, Ddpca::Real &y_d){
	//
	Ddpca::Real thet_c = i_c1 * thet_1;
	//
    //column-major order
	Ddpca::DenseMatrix coefA(2,2,{
        - std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c), std::sin(beta_c) * std::cos(thet_c)
		+ i_2h * std::cos(beta_c) * std::cos(thet_h - thet_1) 
		- i_2h * std::sin(beta_c) * std::sin(thet_c) * std::sin(thet_h - thet_1),
        std::sin(thet_c), -std::sin(thet_c) - i_2h * std::cos(thet_c) * std::sin(thet_h -thet_1)
    });
	std::array<Ddpca::Real,2> coefB = {
        - r_b2 * std::sin(beta_c) * std::sin(thet_c) + a_1c * std::sin(beta_c), 
        + r_b2 * std::sin(beta_c) * std::sin(thet_c) 
		+ i_2h * r_b2 * std::sin(beta_c) * std::cos(thet_c) * std::sin(thet_h - thet_1)
		- a_1c * std::sin(beta_c)
		- i_2h * a_1c * std::cos(beta_c) * std::cos(thet_c) * std::cos(thet_h - thet_1)
		+ i_2h * a_h2 * std::cos(beta_c) * std::cos(thet_c)};
    std::array<Ddpca::Real,2> xy_d;
    Ddpca::SCAL(-1.0, coefB);
    Ddpca::Solve(coefA, coefB, xy_d);
	x_d = xy_d[0];
	y_d = xy_d[1];
}

void DehwSurf::PD_FSME(Ddpca::Real thet_1, Ddpca::Real thet_h, Ddpca::Real &x_d, Ddpca::Real &y_d, 
	Ddpca::DenseMatrix &Pxy_d){
	//
	Ddpca::Real thet_c = i_c1 * thet_1;
	//
	Ddpca::DenseMatrix coefA(2,2,{
        - std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c), std::sin(beta_c) * std::cos(thet_c)
		+ i_2h * std::cos(beta_c) * std::cos(thet_h - thet_1) 
		- i_2h * std::sin(beta_c) * std::sin(thet_c) * std::sin(thet_h - thet_1),
        std::sin(thet_c), -std::sin(thet_c) - i_2h * std::cos(thet_c) * std::sin(thet_h -thet_1)
    });
	std::array<Ddpca::Real,2> coefB = {
        - r_b2 * std::sin(beta_c) * std::sin(thet_c) + a_1c * std::sin(beta_c),
		+ r_b2 * std::sin(beta_c) * std::sin(thet_c) 
		+ i_2h * r_b2 * std::sin(beta_c) * std::cos(thet_c) * std::sin(thet_h - thet_1)
		- a_1c * std::sin(beta_c)
		- i_2h * a_1c * std::cos(beta_c) * std::cos(thet_c) * std::cos(thet_h - thet_1)
		+ i_2h * a_h2 * std::cos(beta_c) * std::cos(thet_c)
    };
	std::array<Ddpca::Real,2> xy_d;
    Ddpca::SCAL(-1.0, coefB);
    Ddpca::Solve(coefA, coefB, xy_d);
	x_d = xy_d[0];
	y_d = xy_d[1];
	//derivative relative to thet_1
	Ddpca::DenseMatrix PcoefA(2,2,{
        std::sin(beta_c) * std::sin(thet_c) * i_c1, - std::sin(beta_c) * std::sin(thet_c) * i_c1
		- i_2h * std::cos(beta_c) * std::sin(thet_h - thet_1) * -1.0
		- i_2h * std::sin(beta_c) * std::cos(thet_c) * i_c1 * std::sin(thet_h - thet_1)
		- i_2h * std::sin(beta_c) * std::sin(thet_c) * std::cos(thet_h - thet_1) * -1.0,
        std::cos(thet_c) * i_c1, - std::cos(thet_c) * i_c1 
		+ i_2h * std::sin(thet_c) * i_c1 * std::sin(thet_h -thet_1) 
		- i_2h * std::cos(thet_c) * std::cos(thet_h -thet_1) * -1.0
    });
	std::array<Ddpca::Real,2> PcoefB = {
        - r_b2 * std::sin(beta_c) * std::cos(thet_c) * i_c1,
		r_b2 * std::sin(beta_c) * std::cos(thet_c) * i_c1
		- i_2h * r_b2 * std::sin(beta_c) * std::sin(thet_c) * i_c1 * std::sin(thet_h - thet_1)
		+ i_2h * r_b2 * std::sin(beta_c) * std::cos(thet_c) * std::cos(thet_h - thet_1) * -1.0
		+ i_2h * a_1c * std::cos(beta_c) * std::sin(thet_c) * i_c1 * std::cos(thet_h - thet_1)
		+ i_2h * a_1c * std::cos(beta_c) * std::cos(thet_c) * std::sin(thet_h - thet_1) * -1.0
		- i_2h * a_h2 * std::cos(beta_c) * std::sin(thet_c) * i_c1
    };
    std::array<Ddpca::Real,2> tempResu;
    Ddpca::GEMV(PcoefA, xy_d, tempResu);
    Ddpca::SCAL(-1.0, PcoefB);
    Ddpca::AXPY(-1.0, tempResu, PcoefB);
    Ddpca::Solve(coefA, PcoefB, tempResu);
    Pxy_d(0,0) = tempResu[0];
	Pxy_d(1,0) = tempResu[1];
	//derivative relative to thet_h
    PcoefA.Fill(
        0.0, 0.0 
		- i_2h * std::cos(beta_c) * std::sin(thet_h - thet_1) 
		- i_2h * std::sin(beta_c) * std::sin(thet_c) * std::cos(thet_h - thet_1),
        0.0, - 0.0 - i_2h * std::cos(thet_c) * std::cos(thet_h -thet_1)
    );
    PcoefB = {
        0.0 + 0.0,
		+ 0.0 
		+ i_2h * r_b2 * std::sin(beta_c) * std::cos(thet_c) * std::cos(thet_h - thet_1)
		- 0.0
		+ i_2h * a_1c * std::cos(beta_c) * std::cos(thet_c) * std::sin(thet_h - thet_1)
		+ 0.0
    };
    Ddpca::GEMV(PcoefA, xy_d, tempResu);
    Ddpca::SCAL(-1.0, PcoefB);
    Ddpca::AXPY(-1.0, tempResu, PcoefB);
    Ddpca::Solve(coefA, PcoefB, tempResu);
    Pxy_d(0,1) = tempResu[0];
	Pxy_d(1,1) = tempResu[1];
}

void DehwSurf::WORM_DC2R(Ddpca::Real x_d, Ddpca::Real y_d, Ddpca::Real thet_c, std::array<Ddpca::Real,3> &r_1_1){
	Ddpca::Real thet_1 = i_1c * thet_c;
	std::array<Ddpca::Real,3> xyz = {
        - x_d,
        r_b2 - y_d * std::sin(beta_c),
        y_d * std::cos(beta_c),
    };
	//R_oc,c, column-major order
	Ddpca::DenseMatrix rota(3,3,{
        std::cos(thet_c),std::sin(thet_c),0.0,
        -std::sin(thet_c),std::cos(thet_c),0.0,
        0.0,0.0,1.0
    });
    std::array<Ddpca::Real,3> tempXyz = xyz;
    Ddpca::GEMV(rota, tempXyz, xyz);
	//R_o1,oc
    rota.Fill(
        1.0,0.0,0.0,
        0.0,0.0,1.0,
        0.0,-1.0,0.0
    );
    tempXyz = xyz;
    Ddpca::GEMV(rota, tempXyz, xyz);
	//T_o1,oc
	xyz[0] = xyz[0] + a_1c;
	//R_1,o1
    rota.Fill(
        std::cos(thet_1),-std::sin(thet_1),0.0,
        std::sin(thet_1),std::cos(thet_1),0.0,
        0.0,0.0,1.0
    );
    Ddpca::GEMV(rota,xyz,r_1_1);
}

void DehwSurf::WHEE_1H2R(Ddpca::Real x_d, Ddpca::Real y_d, Ddpca::Real thet_1, Ddpca::Real thet_h, 
	std::array<Ddpca::Real,3> &r_2_2){
	//
	Ddpca::Real thet_c = i_c1 * thet_1;
	Ddpca::Real thet_2 = i_2h * thet_h;
	//
	std::array<Ddpca::Real,3> xyz;
	Ddpca::DenseMatrix rota(3,3);
	WORM_DC2R(x_d, y_d, thet_c, xyz);
	//R_oh,h
    rota.Fill(
        std::cos(thet_h),std::sin(thet_h),0.0,
        -std::sin(thet_h),std::cos(thet_h),0.0,
        0.0,0.0,1.0
    );
    std::array<Ddpca::Real,3> tempXyz = xyz;
    Ddpca::GEMV(rota, tempXyz, xyz);
	//R_o2,oh
    rota.Fill(
        1.0,0.0,0.0,
        0.0,0.0,-1.0,
        0.0,1.0,0.0
    );
    tempXyz = xyz;
    Ddpca::GEMV(rota, tempXyz, xyz);
	//T_o2,oh
	xyz[0] = xyz[0] - a_h2;
	//R_2,o2
    rota.Fill(
        std::cos(thet_2),-std::sin(thet_2),0.0,
        std::sin(thet_2),std::cos(thet_2),0.0,
        0.0,0.0,1.0
    );
    Ddpca::GEMV(rota, xyz, r_2_2);
}

void DehwSurf::PD_WHEE_1H2R(Ddpca::Real x_d, Ddpca::Real y_d, Ddpca::Real thet_1, Ddpca::Real thet_h, 
	std::array<Ddpca::Real,3> &r_2_2, Ddpca::DenseMatrix &Dr_2_2){
	//
	Ddpca::Real thet_c = i_c1 * thet_1;
	Ddpca::Real thet_2 = i_2h * thet_h;
	//relative to thet_1
	Ddpca::DenseMatrix Dxy_d(2,2);
	PD_FSME(thet_1, thet_h, x_d, y_d, Dxy_d);
	Ddpca::Real Dx_d = Dxy_d(0,0);
	Ddpca::Real Dy_d = Dxy_d(1,0);
	std::array<Ddpca::Real,3> r_c_c, Dr_c_c;
    r_c_c = {- x_d, r_b2 - y_d * std::sin(beta_c), y_d * std::cos(beta_c)};
    Dr_c_c = {- Dx_d, - Dy_d * std::sin(beta_c), Dy_d * std::cos(beta_c)};
	Ddpca::DenseMatrix R_oc_c(3,3), DR_oc_c(3,3);
    //column-major order
    R_oc_c.Fill(
        std::cos(thet_c),std::sin(thet_c),0.0,
        - std::sin(thet_c),std::cos(thet_c),0.0,
        0.0,0.0,1.0
    );
    DR_oc_c.Fill(
        - i_c1 * std::sin(thet_c),i_c1 * std::cos(thet_c),0.0,
        - i_c1 * std::cos(thet_c),- i_c1 * std::sin(thet_c),0.0,
        0.0,0.0,0.0
    );
	std::array<Ddpca::Real,3> r_c_oc;
    Ddpca::GEMV(R_oc_c, r_c_c, r_c_oc);
	std::array<Ddpca::Real,3> Dr_c_oc, tempVect;
    Ddpca::GEMV(DR_oc_c, r_c_c, tempVect);
    Ddpca::GEMV(R_oc_c, Dr_c_c, Dr_c_oc);
    Ddpca::XPEY(Dr_c_oc, tempVect);
	Ddpca::DenseMatrix R_o1_oc(3,3);
	std::array<Ddpca::Real,3> T_o1_oc;
    R_o1_oc.Fill(
        1.0,0.0,0.0,
        0.0,0.0,1.0,
        0.0,-1.0,0.0
    );
    T_o1_oc = {a_1c, 0.0, 0.0};
    std::array<Ddpca::Real,3> r_c_o1, Dr_c_o1;
    Ddpca::GEMV(R_o1_oc, r_c_oc, r_c_o1);
    Ddpca::XPEY(r_c_o1, T_o1_oc);
    Ddpca::GEMV(R_o1_oc, Dr_c_oc, Dr_c_o1);
    Ddpca::DenseMatrix R_1_o1(3,3), DR_1_o1(3,3);
    R_1_o1.Fill(
        std::cos(thet_1),-std::sin(thet_1),0.0,
        std::sin(thet_1),std::cos(thet_1),0.0,
        0.0,0.0,1.0
    );
    DR_1_o1.Fill(
        -std::sin(thet_1),-std::cos(thet_1),0.0,
        std::cos(thet_1),-std::sin(thet_1),0.0,
        0.0,0.0,0.0
    );
    std::array<Ddpca::Real,3> r_1_1, Dr_1_1;
    Ddpca::GEMV(R_1_o1, r_c_o1, r_1_1);
    Ddpca::GEMV(DR_1_o1, r_c_o1, tempVect);
    Ddpca::GEMV(R_1_o1, Dr_c_o1, Dr_1_1);
    Ddpca::XPEY(Dr_1_1, tempVect);
	Ddpca::DenseMatrix R_oh_h(3,3);
    R_oh_h.Fill(
        std::cos(thet_h),std::sin(thet_h),0.0,
        - std::sin(thet_h),std::cos(thet_h),0.0,
        0.0,0.0,1.0
    );
    std::array<Ddpca::Real,3> r_h_oh, Dr_h_oh;
    Ddpca::GEMV(R_oh_h, r_1_1, r_h_oh);
    Ddpca::GEMV(R_oh_h, Dr_1_1, Dr_h_oh);
    Ddpca::DenseMatrix R_o2_oh(3,3);
    std::array<Ddpca::Real,3> T_o2_oh;
    R_o2_oh.Fill(
        1.0,0.0,0.0,
        0.0,0.0,-1.0,
        0.0,1.0,0.0
    );
    T_o2_oh = {- a_h2, 0.0, 0.0};
    std::array<Ddpca::Real,3> r_h_o2, Dr_h_o2;
    Ddpca::GEMV(R_o2_oh, r_h_oh, r_h_o2);
    Ddpca::XPEY(r_h_o2, T_o2_oh);
    Ddpca::GEMV(R_o2_oh, Dr_h_oh, Dr_h_o2);
    Ddpca::DenseMatrix R_2_o2(3,3);
    R_2_o2.Fill(
        std::cos(thet_2),-std::sin(thet_2),0.0,
        std::sin(thet_2),std::cos(thet_2),0.0,
        0.0,0.0,1.0
    );
    Ddpca::GEMV(R_2_o2, r_h_o2, r_2_2);
    Ddpca::GEMV(R_2_o2, Dr_h_o2, tempVect);
    Dr_2_2(0,0) = tempVect[0];
    Dr_2_2(1,0) = tempVect[1];
    Dr_2_2(2,0) = tempVect[2];
	//relative to thet_h
	Dx_d = Dxy_d(0,1);
	Dy_d = Dxy_d(1,1);
    Dr_c_c = {- Dx_d, - Dy_d * std::sin(beta_c), Dy_d * std::cos(beta_c)};
    Ddpca::GEMV(R_oc_c, Dr_c_c, Dr_c_oc);
    Ddpca::GEMV(R_o1_oc, Dr_c_oc, Dr_c_o1);
    Ddpca::GEMV(R_1_o1, Dr_c_o1, Dr_1_1);
    Ddpca::DenseMatrix DR_oh_h(3,3);
    DR_oh_h.Fill(
        -std::sin(thet_h),std::cos(thet_h),0.0,
        -std::cos(thet_h),-std::sin(thet_h),0.0,
        0.0,0.0,0.0
    );
    Ddpca::GEMV(DR_oh_h, r_1_1, tempVect);
    Ddpca::GEMV(R_oh_h, Dr_1_1, Dr_h_oh);
    Ddpca::XPEY(Dr_h_oh, tempVect);
    Ddpca::GEMV(R_o2_oh, Dr_h_oh, Dr_h_o2);
    Ddpca::DenseMatrix DR_2_o2(3,3);
    DR_2_o2.Fill(
        -i_2h * std::sin(thet_2),-i_2h * std::cos(thet_2),0.0,
        i_2h * std::cos(thet_2),-i_2h * std::sin(thet_2),0.0,
        0.0,0.0,0.0
    );
    std::array<Ddpca::Real,3> tempVect_1;
    Ddpca::GEMV(DR_2_o2, r_h_o2, tempVect);
    Ddpca::GEMV(R_2_o2, Dr_h_o2, tempVect_1);
    Ddpca::XPEY(tempVect_1, tempVect);
    Dr_2_2(0,1) = tempVect_1[0];
    Dr_2_2(1,1) = tempVect_1[1];
    Dr_2_2(2,1) = tempVect_1[2];
}

void DehwSurf::CILFOFE(Ddpca::Real thet_1, Ddpca::Real x_d, Ddpca::Real y_d, 
	Ddpca::Real &Psi_1, Ddpca::Real &kapp_1xd, Ddpca::Real &kapp_1yd, Ddpca::Real &tau_1xd){
	//
	Ddpca::Real thet_c = thet_1 / i_1c;
	//
	Ddpca::Real kapp_cxd = 0.0;
	Ddpca::Real kapp_cyd = 0.0;
	Ddpca::Real tau_cxd = 0.0;
	std::array<Ddpca::Real,3> i_d_c, j_d_c, i_d_oc, j_d_oc, omeg_c1_oc, v_c1_oc;
	i_d_c = {-1.0, 0.0, 0.0};
	j_d_c = {0.0, - std::sin(beta_c), std::cos(beta_c)};
	//R_oc,c
	Ddpca::DenseMatrix rota(3,3,{
		std::cos(thet_c),std::sin(thet_c),0.0,
		-std::sin(thet_c),std::cos(thet_c),0.0,
		0.0,0.0,1.0
	});
	Ddpca::GEMV(rota, i_d_c, i_d_oc);
	Ddpca::GEMV(rota, j_d_c, j_d_oc);
	omeg_c1_oc = {0.0, -1.0, i_c1};
	v_c1_oc = {- y_d * std::cos(beta_c) 
		- i_c1 * (- x_d * std::sin(thet_c) + std::cos(thet_c) * (r_b2 - y_d * std::sin(beta_c))),
		i_c1 * (- x_d * std::cos(thet_c) - std::sin(thet_c) * (r_b2 - y_d * std::sin(beta_c))),
		- x_d * std::cos(thet_c) - std::sin(thet_c) * (r_b2 - y_d * std::sin(beta_c)) + a_1c};
	Ddpca::Real N_1xd = kapp_cxd * Ddpca::DOT(v_c1_oc, i_d_oc) + tau_cxd * Ddpca::DOT(v_c1_oc, j_d_oc)
		+ Ddpca::DOT(omeg_c1_oc, j_d_oc);
	Ddpca::Real N_1yd = tau_cxd * Ddpca::DOT(v_c1_oc, i_d_oc) + kapp_cyd * Ddpca::DOT(v_c1_oc, j_d_oc)
		- Ddpca::DOT(omeg_c1_oc, i_d_oc);
	std::array<Ddpca::Real,3> N_1_oc = i_d_oc;
	Ddpca::SCAL(N_1xd, N_1_oc);
	Ddpca::AXPY(N_1yd, j_d_oc, N_1_oc);
	Ddpca::Real PPhi_1Pthet_1 = x_d * std::sin(beta_c) * std::sin(thet_c) / i_1c
		+ y_d * std::cos(thet_c) / i_1c
		- r_b2 * std::sin(beta_c) * std::cos(thet_c) / i_1c;
	Psi_1 = Ddpca::DOT(N_1_oc, v_c1_oc) + PPhi_1Pthet_1;
	Ddpca::Real kapp_c1xd = N_1xd * N_1xd / Psi_1;
	Ddpca::Real kapp_c1yd = N_1yd * N_1yd / Psi_1;	
	Ddpca::Real tau_c1xd = N_1xd * N_1yd / Psi_1;
	kapp_1xd = kapp_cxd - kapp_c1xd;
	kapp_1yd = kapp_cyd - kapp_c1yd;
	tau_1xd = tau_cxd - tau_c1xd;
}

Ddpca::Real DehwSurf::CILFOSE_NI(Ddpca::Real thet_1, Ddpca::Real thet_h, Ddpca::Real &kapp_h2N){
	//
	Ddpca::Real thet_c = thet_1 / i_1c;
	// Ddpca::Real thet_2 = thet_h / i_h2;
	//
	Ddpca::Real x_d, y_d;
	FSME(thet_1, thet_h, x_d, y_d);
	Ddpca::Real Psi_1, kapp_1xd, kapp_1yd, tau_1xd;
	CILFOFE(thet_1, x_d, y_d, Psi_1, kapp_1xd, kapp_1yd, tau_1xd);
	Ddpca::Real kapp_hxd = kapp_1xd;
	Ddpca::Real kapp_hyd = kapp_1yd;
	Ddpca::Real tau_hxd = tau_1xd;
	//
	std::array<Ddpca::Real,3> i_d_c, j_d_c, i_d_oh, j_d_oh;
	i_d_c = {-1.0, 0.0, 0.0};
	j_d_c = {0.0, - std::sin(beta_c), std::cos(beta_c)};
	std::array<Ddpca::Real,3> r_c_c, r_h_oh;
	r_c_c = {- x_d, r_b2 - y_d * std::sin(beta_c), y_d * std::cos(beta_c)};
	//R_oc,c
	Ddpca::DenseMatrix rota(3,3,{
		std::cos(thet_c),std::sin(thet_c),0.0,
		-std::sin(thet_c),std::cos(thet_c),0.0,
		0.0,0.0,1.0
	});
	Ddpca::GEMV(rota, i_d_c, i_d_oh);
	Ddpca::GEMV(rota, j_d_c, j_d_oh);
	Ddpca::GEMV(rota, r_c_c, r_h_oh);
	//R_o1,oc
	rota.Fill(
		1.0,0.0,0.0,
		0.0,0.0,1.0,
		0.0,-1.0,0.0
	);
	std::array<Ddpca::Real,3> tempVect = i_d_oh;
	Ddpca::GEMV(rota, tempVect, i_d_oh);
	tempVect = j_d_oh;
	Ddpca::GEMV(rota, tempVect, j_d_oh);
	tempVect = r_h_oh;
	Ddpca::GEMV(rota, tempVect, r_h_oh);
	//T_o1,oc
	r_h_oh[0] = r_h_oh[0] + a_1c;
	//R_1,o1
	rota.Fill(
		std::cos(thet_1),-std::sin(thet_1),0.0,
		std::sin(thet_1),std::cos(thet_1),0.0,
		0.0,0.0,1.0
	);
	tempVect = i_d_oh;
	Ddpca::GEMV(rota, tempVect, i_d_oh);
	tempVect = j_d_oh;
	Ddpca::GEMV(rota, tempVect, j_d_oh);
	tempVect = r_h_oh;
	Ddpca::GEMV(rota, tempVect, r_h_oh);
	//R_oh,h
	rota.Fill(
		std::cos(thet_h),std::sin(thet_h),0.0,
		-std::sin(thet_h),std::cos(thet_h),0.0,
		0.0,0.0,1.0
	);
	tempVect = i_d_oh;
	Ddpca::GEMV(rota, tempVect, i_d_oh);
	tempVect = j_d_oh;
	Ddpca::GEMV(rota, tempVect, j_d_oh);
	tempVect = r_h_oh;
	Ddpca::GEMV(rota, tempVect, r_h_oh);
	//
	std::array<Ddpca::Real,3> omeg_h2_oh, omeg_2_oh, o_h2_oh, v_h2_oh;
	omeg_h2_oh = {0.0, i_2h, 1.0};
	omeg_2_oh = {0.0, - i_2h, 0.0};
	o_h2_oh = {- a_h2, 0.0, 0.0};
	tempVect = Ddpca::Cross(omeg_h2_oh, r_h_oh);
	v_h2_oh = Ddpca::Cross(omeg_2_oh, o_h2_oh);
	Ddpca::SCAL(-1.0, v_h2_oh);
	Ddpca::XPEY(v_h2_oh, tempVect);
	Ddpca::Real N_2xd = kapp_hxd * Ddpca::DOT(v_h2_oh, i_d_oh) + tau_hxd * Ddpca::DOT(v_h2_oh, j_d_oh)
		+ Ddpca::DOT(omeg_h2_oh, j_d_oh);
	Ddpca::Real N_2yd = tau_hxd * Ddpca::DOT(v_h2_oh, i_d_oh) + kapp_hyd * Ddpca::DOT(v_h2_oh, j_d_oh)
		- Ddpca::DOT(omeg_h2_oh, i_d_oh);
	std::array<Ddpca::Real,3> N_2_oh = i_d_oh;
	Ddpca::SCAL(N_2xd, N_2_oh);
	Ddpca::AXPY(N_2yd, j_d_oh, N_2_oh);
	Ddpca::Real B_11 = i_2h * x_d * std::cos(beta_c) - i_2h * a_1c * std::cos(beta_c) * std::cos(thet_c);
	Ddpca::Real B_12 = - i_2h * x_d * std::sin(beta_c) * std::sin(thet_c) 
		- i_2h * y_d * std::cos(thet_c) + i_2h * r_b2 * std::sin(beta_c) * std::cos(thet_c);
	Ddpca::Real PPhi_2Pthet_h = - B_11 * std::sin(thet_h - thet_1) + B_12 * std::cos(thet_h - thet_1);
	Ddpca::Real Psi_2 = Ddpca::DOT(N_2_oh, v_h2_oh) + PPhi_2Pthet_h;
	kapp_h2N = (N_2xd * N_2xd + N_2yd * N_2yd) / Psi_2;
	return Psi_2;
}

void DehwSurf::WORM_CURV_2_CART(Ddpca::Real xi_11, Ddpca::Real xi_12, 
	std::array<Ddpca::Real,3> &r_1_1, Ddpca::Real &thet_c){
// std::cout << "Here0\n";
	std::array<Ddpca::Real,2> x;//0 - thet_c, 1 - x_d
	x = {i_c1 * xi_11, d[1] / 2.0};//
	while(true){
// std::cout << "Here1\n";
		std::array<Ddpca::Real,2> func, deltX;
		Ddpca::DenseMatrix Dfunc(2,2);
		//function
		thet_c = x[0];
		Ddpca::Real thet_1 = i_1c * thet_c;
		Ddpca::Real x_d = x[1];
		Ddpca::Real y_d = - ((- std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c)) * x_d 
			+ (- r_b2 * std::sin(beta_c) * std::sin(thet_c) + a_1c * std::sin(beta_c))) / std::sin(thet_c);
		std::array<Ddpca::Real,3> r_c_c, T_o1_oc, r_1_o1;
		Ddpca::DenseMatrix R_oc_c(3,3), R_o1_oc(3,3), R_1_o1(3,3), tempMatr(3,3);
		r_c_c = {- x_d, r_b2 - y_d * std::sin(beta_c), y_d * std::cos(beta_c)};
		R_oc_c.Fill(
			std::cos(thet_c),std::sin(thet_c),0.0,
			- std::sin(thet_c),std::cos(thet_c),0.0,
			0.0,0.0,1.0
		);
		R_o1_oc.Fill(
			1.0,0.0,0.0,
			0.0,0.0,1.0,
			0.0,-1.0,0.0
		);
		R_1_o1.Fill(
			std::cos(thet_1),- std::sin(thet_1),0.0,
			std::sin(thet_1),std::cos(thet_1),0.0,
			0.0,0.0,1.0
		);
// std::cout << "Here2\n";
		T_o1_oc = {a_1c, 0.0, 0.0};
		Ddpca::GEMM(R_o1_oc, R_oc_c, tempMatr);
		Ddpca::GEMV(tempMatr, r_c_c, r_1_o1);
		Ddpca::XPEY(r_1_o1, T_o1_oc);
		Ddpca::GEMV(R_1_o1, r_1_o1, r_1_1);
		func = {thet_1 - std::atan2(r_1_o1[1], r_1_o1[0]) - xi_11,
			r_1_1[2] * r_1_1[2] 
			+ (a_h2 - std::sqrt(r_1_1[0] * r_1_1[0] + r_1_1[1] * r_1_1[1])) 
			* (a_h2 - std::sqrt(r_1_1[0] * r_1_1[0] + r_1_1[1] * r_1_1[1]))
			- xi_12 * xi_12};
// std::cout << "func = " << func[0] << " " << func[1] << "\n";
		//derivative 1
		Ddpca::Real Dy_d = - (((+ std::sin(beta_c) * std::sin(thet_c) - 0.0) * x_d 
			+ (- r_b2 * std::sin(beta_c) * std::cos(thet_c) + 0.0)) * std::sin(thet_c) 
			- ((- std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c)) * x_d 
			+ (- r_b2 * std::sin(beta_c) * std::sin(thet_c) + a_1c * std::sin(beta_c))) * std::cos(thet_c))
			/ (std::sin(thet_c) * std::sin(thet_c));
		std::array<Ddpca::Real,3> Dr_c_c, Dr_1_o1, Dr_1_1;
		Ddpca::DenseMatrix DR_oc_c(3,3), DR_1_o1(3,3);
		Dr_c_c = {0.0, 0.0 - Dy_d * std::sin(beta_c), Dy_d * std::cos(beta_c)};
		DR_oc_c.Fill(
			- std::sin(thet_c),std::cos(thet_c),0.0,
			- std::cos(thet_c),- std::sin(thet_c),0.0,
			0.0,0.0,0.0
		);
		DR_1_o1.Fill(
			- i_1c * std::sin(thet_1),- i_1c * std::cos(thet_1),0.0,
			i_1c * std::cos(thet_1),- i_1c * std::sin(thet_1),0.0,
			0.0,0.0,0.0
		);
// std::cout << "Here3\n";
		Ddpca::GEMM(R_o1_oc, DR_oc_c, tempMatr);
		Ddpca::GEMV(tempMatr, r_c_c, Dr_1_o1);
		Ddpca::GEMM(R_o1_oc, R_oc_c, tempMatr);
		std::array<Ddpca::Real,3> tempVect;
		Ddpca::GEMV(tempMatr, Dr_c_c, tempVect);
		Ddpca::XPEY(Dr_1_o1, tempVect);
		Ddpca::GEMV(DR_1_o1, r_1_o1, tempVect);
		Ddpca::GEMV(R_1_o1, Dr_1_o1, Dr_1_1);
		Ddpca::XPEY(Dr_1_1, tempVect);
		Dfunc(0,0) = i_1c - (Dr_1_o1[1] * r_1_o1[0] - r_1_o1[1] * Dr_1_o1[0])
			/ (r_1_o1[1] * r_1_o1[1] + r_1_o1[0] * r_1_o1[0]);
		Dfunc(1,0) = 2.0 * r_1_1[2] * Dr_1_1[2] 
			+ 2.0 * (a_h2 - std::sqrt(r_1_1[0] * r_1_1[0] + r_1_1[1] * r_1_1[1])) 
			* (0.0 - (r_1_1[0] * Dr_1_1[0] + r_1_1[1] * Dr_1_1[1]) 
			/ std::sqrt(r_1_1[0] * r_1_1[0] + r_1_1[1] * r_1_1[1]));
		//derivative 2
		Dy_d = - ((- std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c)) 
			+ 0.0) / std::sin(thet_c);
		Dr_c_c = {- 1.0, 0.0 - Dy_d * std::sin(beta_c), Dy_d * std::cos(beta_c)};
		Ddpca::GEMM(R_o1_oc, R_oc_c, tempMatr);
		Ddpca::GEMV(tempMatr, Dr_c_c, Dr_1_o1);
		Ddpca::GEMV(R_1_o1, Dr_1_o1, Dr_1_1);
		Dfunc(0,1) = - (Dr_1_o1[1] * r_1_o1[0] - r_1_o1[1] * Dr_1_o1[0])
			/ (r_1_o1[1] * r_1_o1[1] + r_1_o1[0] * r_1_o1[0]);
		Dfunc(1,1) = 2.0 * r_1_1[2] * Dr_1_1[2] 
			+ 2.0 * (a_h2 - std::sqrt(r_1_1[0] * r_1_1[0] + r_1_1[1] * r_1_1[1])) 
			* (0.0 - (r_1_1[0] * Dr_1_1[0] + r_1_1[1] * Dr_1_1[1]) 
			/ std::sqrt(r_1_1[0] * r_1_1[0] + r_1_1[1] * r_1_1[1]));
// std::cout << "Dfunc = " << Dfunc(0,0) << " " << Dfunc(0,1) << " " << Dfunc(1,0) << " " << Dfunc(1,1) << "\n";
		//iteration
		Ddpca::SCAL(-1.0, func);
		Ddpca::Solve(Dfunc,func,deltX);
// std::cout << "Here4\n";
// std::cout << "deltX = " << deltX[0] << " " << deltX[1] << "\n";
		if(Ddpca::NRM2(deltX) < 1.0E-12){
			if(a_h2 - std::sqrt(r_1_1[0] * r_1_1[0] + r_1_1[1] * r_1_1[1]) < 0.0){
				std::cout << "ERROR in DEHWSURF::WORM_CURV_2_CART!" << std::endl;
			}
			break;
		}
		Ddpca::XPEY(x, deltX);
// std::cout << "Here5\n"; Ddpca::I64 test; std::cin >> test;
	}
}

void DehwSurf::WHEE_G2L(std::array<Ddpca::Real,3> r_2_2, Ddpca::Real &angl_f, Ddpca::Real &radi_f, 
	Ddpca::Real &R_fmini, Ddpca::Real &R_fmaxi){
	Ddpca::Real radi_xi = a_h2 - std::sqrt(r_2_2[0] * r_2_2[0] + r_2_2[1] * r_2_2[1]);
	Ddpca::Real coor_zi = r_2_2[2];
	angl_f = std::atan2(coor_zi, radi_xi);
	radi_f = std::sqrt(radi_xi * radi_xi + coor_zi * coor_zi);
	Ddpca::Real angl_ai = angl_f - std::asin(offsR_a * std::sin(angl_f) / R_a[1]);
	R_fmini = (R_a[1] * std::cos(angl_ai) - offsR_a) / std::cos(angl_f);
	R_fmaxi = R_t[1];
}

void DehwSurf::WHEE_CURV_2_CART_1(Ddpca::Real xi_21, Ddpca::Real xi_22, std::array<Ddpca::Real,3> &r_2_2, 
	Ddpca::Real &thet_c, Ddpca::Real &thet_h, Ddpca::I64 f_lr, Ddpca::Real &x_d, Ddpca::Real &y_d){
	std::array<Ddpca::Real,2> x;
	x = {thet_c, thet_h};
	for(Ddpca::I64 ti = 0; ti < 1000; ti ++){
		std::array<Ddpca::Real,2> func, deltX;
		Ddpca::DenseMatrix Dfunc(2,2);
		thet_c = x[0];
		thet_h = x[1];
		//function
		Ddpca::Real thet_1 = i_1c * thet_c;
		FSME(thet_1, thet_h, x_d, y_d);
		Ddpca::DenseMatrix Dr_2_2(3,2);
		PD_WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2, Dr_2_2);
		Ddpca::Real coorX = a_h2 - std::sqrt(r_2_2[0] * r_2_2[0] + r_2_2[1] * r_2_2[1]);
		func = {std::atan2(r_2_2[2], coorX) - xi_21,
			r_2_2[2] * r_2_2[2] + coorX * coorX - xi_22 * xi_22};
		//derivative 1
		Ddpca::Real DcoorX = - (r_2_2[0] * Dr_2_2(0,0) + r_2_2[1] * Dr_2_2(1,0)) 
			/ std::sqrt(r_2_2[0] * r_2_2[0] + r_2_2[1] * r_2_2[1]);
		Dfunc(0,0) = (Dr_2_2(2,0) * coorX - r_2_2[2] * DcoorX) 
			/ (coorX * coorX + r_2_2[2] * r_2_2[2]);
		Dfunc(1,0) = 2.0 * r_2_2[2] * Dr_2_2(2,0) + 2.0 * coorX * DcoorX;
		Dfunc(0,0) *= i_1c;
		Dfunc(1,0) *= i_1c;
		//derivative 2
		DcoorX = - (r_2_2[0] * Dr_2_2(0,1) + r_2_2[1] * Dr_2_2(1,1)) 
			/ std::sqrt(r_2_2[0] * r_2_2[0] + r_2_2[1] * r_2_2[1]);
		Dfunc(0,1) = (Dr_2_2(2,1) * coorX - r_2_2[2] * DcoorX) 
			/ (coorX * coorX + r_2_2[2] * r_2_2[2]);
		Dfunc(1,1) = 2.0 * r_2_2[2] * Dr_2_2(2,1) + 2.0 * coorX * DcoorX;
		//iteration
		Ddpca::SCAL(-1.0, func);
		Ddpca::Solve(Dfunc,func,deltX);
		if(Ddpca::NRM2(deltX) < 1.0E-12){
			if(coorX < 0.0){
				std::cout << "ERROR in DEHWSURF::WHEE_CURV_2_CART_1!" << std::endl;
			}
			break;
		}
		Ddpca::Real rfac = 2.0;
		Ddpca::I64 rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			std::array<Ddpca::Real,2> x_test = x;
			Ddpca::AXPY(rfac, deltX, x_test);
			if(x_test[0] < 0.01 * Ddpca::PI || x_test[0] > 0.49 * Ddpca::PI){
				continue;
			}
			Ddpca::Real thet_hs, thet_hm;
			SINGULAR_C2H(x[0], thet_hs, thet_hm);
			if((f_lr == 1 && (x_test[1] <= thet_hs + 1.0E-14 || thet_hm - 1.0E-14 <= x_test[1]))
				|| (f_lr == 2 && (x_test[1] <= thet_hm + 1.0E-14 
				|| thet_hs + 2.0 * Ddpca::PI - 1.0E-14 <= x_test[1]))){
				continue;
			}
			Ddpca::Real thet_ct = x_test[0];
			Ddpca::Real thet_ht = x_test[1];
			//function
			Ddpca::Real thet_1t = i_1c * thet_ct;
			Ddpca::Real x_dt, y_dt;
			FSME(thet_1t, thet_ht, x_dt, y_dt);
			std::array<Ddpca::Real,3> r_2_2t;
			Ddpca::DenseMatrix Dr_2_2t(3,2);
			PD_WHEE_1H2R(x_dt, y_dt, thet_1t, thet_ht, r_2_2t, Dr_2_2t);
			Ddpca::Real cooX_t = a_h2 - std::sqrt(r_2_2t[0] * r_2_2t[0] + r_2_2t[1] * r_2_2t[1]);
			std::array<Ddpca::Real,2> func_t;
			func_t = {std::atan2(r_2_2t[2], cooX_t) - xi_21,
				r_2_2t[2] * r_2_2t[2] + cooX_t * cooX_t - xi_22 * xi_22};
			if(Ddpca::NRM2(func_t) < Ddpca::NRM2(func)){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
}

void DehwSurf::WHEE_CURV_2_CART_2(Ddpca::Real xi_21, Ddpca::Real xi_22, std::array<Ddpca::Real,3> &r_c_c, 
		Ddpca::Real &thet_c, Ddpca::Real &x_d, Ddpca::Real &y_d){
	std::array<Ddpca::Real,2> x;
	x = {thet_c, x_d};
	for(Ddpca::I64 ti = 0; ti < 1000; ti ++){
		std::array<Ddpca::Real,2> func, deltX;
		Ddpca::DenseMatrix Dfunc(2,2);
		thet_c = x[0];
		x_d = x[1];
		//function
		y_d = - ((- std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c)) * x_d 
			- r_b2 * std::sin(beta_c) * std::sin(thet_c) + a_1c * std::sin(beta_c)) 
			/ std::sin(thet_c);
		r_c_c = {- x_d, r_b2 - y_d * std::sin(beta_c), y_d * std::cos(beta_c)};
		Ddpca::Real coorX = a_h2 - std::sqrt(r_c_c[0] * r_c_c[0] + r_c_c[1] * r_c_c[1]);
		func = {std::atan2(r_c_c[2], coorX) - xi_21,
			r_c_c[2] * r_c_c[2] + coorX * coorX - xi_22 * xi_22};
		//derivative 1
		Ddpca::Real Dy_d = - (((+ std::sin(beta_c) * std::sin(thet_c) - 0.0) * x_d 
			- r_b2 * std::sin(beta_c) * std::cos(thet_c) + 0.0) * std::sin(thet_c) 
			- ((- std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c)) * x_d 
			- r_b2 * std::sin(beta_c) * std::sin(thet_c) + a_1c * std::sin(beta_c)) * std::cos(thet_c))
			/ (std::sin(thet_c) * std::sin(thet_c));
		std::array<Ddpca::Real,3> Dr_c_c;
		Dr_c_c = {- 0.0, 0.0 - Dy_d * std::sin(beta_c), Dy_d * std::cos(beta_c)};
		Ddpca::Real DcoorX = - (r_c_c[0] * Dr_c_c[0] + r_c_c[1] * Dr_c_c[1]) 
			/ std::sqrt(r_c_c[0] * r_c_c[0] + r_c_c[1] * r_c_c[1]);
		Dfunc(0,0) = (Dr_c_c[2] * coorX - r_c_c[2] * DcoorX) 
			/ (coorX * coorX + r_c_c[2] * r_c_c[2]);
		Dfunc(1,0) = 2.0 * r_c_c[2] * Dr_c_c[2] + 2.0 * coorX * DcoorX;
		//derivative 2
		Dy_d = - ((- std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c)) * 1.0 
			- 0.0 + 0.0) / std::sin(thet_c);
		Dr_c_c = {- 1.0, - Dy_d * std::sin(beta_c), Dy_d * std::cos(beta_c)};
		DcoorX = - (r_c_c[0] * Dr_c_c[0] + r_c_c[1] * Dr_c_c[1]) 
			/ std::sqrt(r_c_c[0] * r_c_c[0] + r_c_c[1] * r_c_c[1]);
		Dfunc(0,1) = (Dr_c_c[2] * coorX - r_c_c[2] * DcoorX) 
			/ (coorX * coorX + r_c_c[2] * r_c_c[2]);
		Dfunc(1,1) = 2.0 * r_c_c[2] * Dr_c_c[2] + 2.0 * coorX * DcoorX;
		//iteration
		Ddpca::SCAL(-1.0, func);
		Ddpca::Solve(Dfunc,func,deltX);
		if(Ddpca::NRM2(deltX) < 1.0E-12){
			if(coorX < 0.0){
				std::cout << "ERROR in DEHWSURF::WHEE_CURV_2_CART_2!" << std::endl;
			}
			break;
		}
		Ddpca::Real rfac = 2.0;
		Ddpca::I64 rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			std::array<Ddpca::Real,2> x_test = x;
			Ddpca::AXPY(rfac, deltX, x_test);
			if(x_test[0] < 0.01 * Ddpca::PI || x_test[0] > 0.49 * Ddpca::PI){
				continue;
			}
			Ddpca::Real thet_ct = x_test[0];
			Ddpca::Real x_dt = x_test[1];
			//function
			Ddpca::Real y_dt = - ((- std::sin(beta_c) * std::cos(thet_ct) - i_c1 * std::cos(beta_c)) * x_dt 
				- r_b2 * std::sin(beta_c) * std::sin(thet_ct) + a_1c * std::sin(beta_c)) 
				/ std::sin(thet_ct);
			std::array<Ddpca::Real,3> r_c_ct;
			r_c_ct = {- x_dt, r_b2 - y_dt * std::sin(beta_c), y_dt * std::cos(beta_c)};
			Ddpca::Real cooX_t = a_h2 - std::sqrt(r_c_ct[0] * r_c_ct[0] + r_c_ct[1] * r_c_ct[1]);
			std::array<Ddpca::Real,2> func_t;
			func_t = {std::atan2(r_c_ct[2], cooX_t) - xi_21,
				r_c_ct[2] * r_c_ct[2] + cooX_t * cooX_t - xi_22 * xi_22};
			if(Ddpca::NRM2(func_t) < Ddpca::NRM2(func)){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
}

void DehwSurf::WHEE_CURV_2_CART_3(Ddpca::Real xi_21, Ddpca::Real xi_22, std::array<Ddpca::Real,3> &r_2_2, 
	Ddpca::Real &thet_c, Ddpca::Real &thet_h, Ddpca::Real xi_11){
	std::array<Ddpca::Real,2> x;
	x = {thet_c, thet_h};
	for(Ddpca::I64 ti = 0; ti < 1000; ti ++){
		std::array<Ddpca::Real,2> func, deltX;
		Ddpca::DenseMatrix Dfunc(2,2);
		thet_c = x[0];
		thet_h = x[1];
		//function
		Ddpca::DenseMatrix Dr_2_2(3,2);
		WHEE_TRAN(thet_c, thet_h, xi_11, r_2_2, Dr_2_2);
		Ddpca::Real coorX = a_h2 - std::sqrt(r_2_2[0] * r_2_2[0] + r_2_2[1] * r_2_2[1]);
		func = {std::atan2(r_2_2[2], coorX) - xi_21,
			r_2_2[2] * r_2_2[2] + coorX * coorX - xi_22 * xi_22};
		//derivative 1
		Ddpca::Real DcoorX = - (r_2_2[0] * Dr_2_2(0,0) + r_2_2[1] * Dr_2_2(1,0)) 
			/ std::sqrt(r_2_2[0] * r_2_2[0] + r_2_2[1] * r_2_2[1]);
		Dfunc(0,0) = (Dr_2_2(2,0) * coorX - r_2_2[2] * DcoorX) 
			/ (coorX * coorX + r_2_2[2] * r_2_2[2]);
		Dfunc(1,0) = 2.0 * r_2_2[2] * Dr_2_2(2,0) + 2.0 * coorX * DcoorX;
		//derivative 2
		DcoorX = - (r_2_2[0] * Dr_2_2(0,1) + r_2_2[1] * Dr_2_2(1,1)) 
			/ std::sqrt(r_2_2[0] * r_2_2[0] + r_2_2[1] * r_2_2[1]);
		Dfunc(0,1) = (Dr_2_2(2,1) * coorX - r_2_2[2] * DcoorX) 
			/ (coorX * coorX + r_2_2[2] * r_2_2[2]);
		Dfunc(1,1) = 2.0 * r_2_2[2] * Dr_2_2(2,1) + 2.0 * coorX * DcoorX;
		//iteration
		Ddpca::SCAL(-1.0, func);
		Ddpca::Solve(Dfunc,func,deltX);
		if(Ddpca::NRM2(deltX) < 1.0E-12){
			if(coorX < 0.0){
				std::cout << "ERROR in DEHWSURF::WHEE_CURV_2_CART_3!" << std::endl;
			}
			break;
		}
		Ddpca::Real rfac = 2.0;
		Ddpca::I64 rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			std::array<Ddpca::Real,2> x_test = x;
			Ddpca::AXPY(rfac, deltX, x_test);
			if(x_test[0] < 0.01 * Ddpca::PI || x_test[0] > 0.49 * Ddpca::PI){
				continue;
			}
			Ddpca::Real thet_ct = x_test[0];
			Ddpca::Real thet_ht = x_test[1];
			//function
			std::array<Ddpca::Real,3> r_2_2t;
			WHEE_TRAN(thet_ct, thet_ht, xi_11, r_2_2t, Dr_2_2);
			Ddpca::Real cooX_t = a_h2 - std::sqrt(r_2_2t[0] * r_2_2t[0] + r_2_2t[1] * r_2_2t[1]);
			std::array<Ddpca::Real,2> func_t;
			func_t = {std::atan2(r_2_2t[2], cooX_t) - xi_21,
				r_2_2t[2] * r_2_2t[2] + cooX_t * cooX_t - xi_22 * xi_22};
			if(Ddpca::NRM2(func_t) < Ddpca::NRM2(func)){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
}

void DehwSurf::WHEE_TRAN(Ddpca::Real thet_c, Ddpca::Real thet_h, Ddpca::Real xi_11, 
	std::array<Ddpca::Real,3> &r_2_2, Ddpca::DenseMatrix &Dr_2_2){
	Ddpca::Real thet_1 = i_1c * thet_c;
	Ddpca::Real thet_2 = i_2h * thet_h;
	//
	Ddpca::Real C_1 = (std::tan(beta_c) * std::cos(thet_c) + i_c1) * std::cos(thet_1 - xi_11) 
		+ i_c1 * std::tan(beta_c) * std::sin(thet_c) * std::sin(thet_1 - xi_11) 
		- std::cos(thet_c) * std::sin(thet_c) * std::sin(thet_1 - xi_11);
	Ddpca::Real C_2 = i_c1 * r_b2 * std::sin(thet_c) - i_c1 * a_1c;
	Ddpca::Real x_a = - C_2 / C_1;
	Ddpca::Real z_a = ((std::tan(beta_c) * std::sin(thet_1 - xi_11) + std::sin(thet_c) * std::cos(thet_1 - xi_11)) 
		* x_a + r_b2 - a_1c * std::sin(thet_c)) / std::cos(thet_c);
	std::array<Ddpca::Real,3> r_1_1;
	r_1_1 = {x_a * std::cos(xi_11), - x_a * std::sin(xi_11), z_a};
	Ddpca::DenseMatrix R_oh_h(3,3,{
		std::cos(thet_h),std::sin(thet_h),0.0,
		- std::sin(thet_h),std::cos(thet_h),0.0,
		0.0,0.0,1.0
	});
	std::array<Ddpca::Real,3> r_h_oh;
	Ddpca::GEMV(R_oh_h, r_1_1, r_h_oh);
	Ddpca::DenseMatrix R_o2_oh(3,3,{
		1.0,0.0,0.0,
		0.0,0.0,-1.0,
		0.0,1.0,0.0
	});
	std::array<Ddpca::Real,3> T_o2_oh;
	T_o2_oh = {- a_h2, 0.0, 0.0}; 
	std::array<Ddpca::Real,3> r_h_o2;
	Ddpca::GEMV(R_o2_oh, r_h_oh, r_h_o2);
	Ddpca::XPEY(r_h_o2, T_o2_oh);
	Ddpca::DenseMatrix R_2_o2(3,3,{
		std::cos(thet_2),- std::sin(thet_2),0.0,
		std::sin(thet_2),std::cos(thet_2),0.0,
		0.0,0.0,1.0
	});
	Ddpca::GEMV(R_2_o2, r_h_o2, r_2_2);
	//derivative relative to thet_c
	Ddpca::Real DC_1 = (- std::tan(beta_c) * std::sin(thet_c) + 0.0) * std::cos(thet_1 - xi_11) 
		- (std::tan(beta_c) * std::cos(thet_c) + i_c1) * i_1c * std::sin(thet_1 - xi_11) 
		+ i_c1 * std::tan(beta_c) * std::cos(thet_c) * std::sin(thet_1 - xi_11) 
		+ i_c1 * std::tan(beta_c) * std::sin(thet_c) * i_1c * std::cos(thet_1 - xi_11) 
		+ std::sin(thet_c) * std::sin(thet_c) * std::sin(thet_1 - xi_11) 
		- std::cos(thet_c) * std::cos(thet_c) * std::sin(thet_1 - xi_11) 
		- std::cos(thet_c) * std::sin(thet_c) * i_1c * std::cos(thet_1 - xi_11);
	Ddpca::Real DC_2 = i_c1 * r_b2 * std::cos(thet_c) - 0.0;
	Ddpca::Real Dx_a = - (DC_2 * C_1 - C_2 * DC_1) / (C_1 * C_1);
	Ddpca::Real Dz_a = (((std::tan(beta_c) * i_1c * std::cos(thet_1 - xi_11) 
		+ std::cos(thet_c) * std::cos(thet_1 - xi_11) 
		- std::sin(thet_c) * i_1c * std::sin(thet_1 - xi_11)) * x_a
		+ (std::tan(beta_c) * std::sin(thet_1 - xi_11) + std::sin(thet_c) * std::cos(thet_1 - xi_11)) * Dx_a
		+ 0.0 - a_1c * std::cos(thet_c)) * std::cos(thet_c) 
		+ ((std::tan(beta_c) * std::sin(thet_1 - xi_11) + std::sin(thet_c) * std::cos(thet_1 - xi_11)) 
		* x_a + r_b2 - a_1c * std::sin(thet_c)) * std::sin(thet_c)) 
		/ (std::cos(thet_c) * std::cos(thet_c));
	std::array<Ddpca::Real,3> Dr_1_1;
	Dr_1_1 = {Dx_a * std::cos(xi_11), - Dx_a * std::sin(xi_11), Dz_a};
	std::array<Ddpca::Real,3> Dr_h_oh;
	std::array<Ddpca::Real,3> Dr_h_o2;
	Ddpca::GEMV(R_oh_h, Dr_1_1, Dr_h_oh);
	Ddpca::GEMV(R_o2_oh, Dr_h_oh, Dr_h_o2);
	std::array<Ddpca::Real,3> tempVect, tempVect_1;
	Ddpca::GEMV(R_2_o2, Dr_h_o2, tempVect);
	Dr_2_2(0,0) = tempVect[0];
	Dr_2_2(1,0) = tempVect[1];
	Dr_2_2(2,0) = tempVect[2];
	//derivative relative to thet_h
	Ddpca::DenseMatrix DR_oh_h(3,3,{
		- std::sin(thet_h),std::cos(thet_h),0.0,
		- std::cos(thet_h),- std::sin(thet_h),0.0,
		0.0,0.0,0.0
	});
	Ddpca::GEMV(DR_oh_h, r_1_1, Dr_h_oh);
	Ddpca::GEMV(R_o2_oh, Dr_h_oh, Dr_h_o2);
	Ddpca::DenseMatrix DR_2_o2(3,3,{
		- i_2h * std::sin(thet_2),- i_2h * std::cos(thet_2),0.0,
		i_2h * std::cos(thet_2),- i_2h * std::sin(thet_2),0.0,
		0.0,0.0,0.0
	});
	Ddpca::GEMV(DR_2_o2, r_h_o2, tempVect);
	Ddpca::GEMV(R_2_o2, Dr_h_o2, tempVect_1);
	Dr_2_2(0,1) = tempVect[0] + tempVect_1[0];
	Dr_2_2(1,1) = tempVect[1] + tempVect_1[1];
	Dr_2_2(2,1) = tempVect[2] + tempVect_1[2];
}

void DehwSurf::WHEE_PHAS(Ddpca::I64 ti, Ddpca::I64 tj, Ddpca::I64 f_ij, std::array<Ddpca::Real,3> r_2_2){
	if(fpha[ti][tj] == 0){
		fpha[ti][tj] = f_ij;
		cartCoor[1][ti][tj] = r_2_2;
	}
	else{
		Ddpca::Real phas_1 = std::atan2(cartCoor[1][ti][tj][1], cartCoor[1][ti][tj][0]);
		if(phas_1 < 0.0){
			phas_1 += 2.0 * Ddpca::PI;
		}
		Ddpca::Real phas_2 = std::atan2(r_2_2[1], r_2_2[0]);
		if(phas_2 < 0.0){
			phas_2 += 2.0 * Ddpca::PI;
		}
		if(phas_2 > phas_1){
			fpha[ti][tj] = f_ij;
			cartCoor[1][ti][tj] = r_2_2;
		}
	}
}

void DehwSurf::WORM_RELI(std::array<Ddpca::Real,3> &tempCoor, Ddpca::I64 ti, Ddpca::I64 tj){
	Ddpca::I64 reliLeng = 40;
	const Ddpca::I64 cc00Size = curvCoor[0][0].size();
	const Ddpca::I64 cc0Size = curvCoor[0].size();
	if(tj > cc00Size - reliLeng || ti < reliLeng || ti > cc0Size - reliLeng){
		std::array<Ddpca::Real,3> tempXYZ = tempCoor;
		std::vector<Ddpca::Real> reliAmou = {14.0E-6, 18.0E-6};//0 - tip, 1 - end
		Ddpca::Real reliExpo = 3.0;
		Ddpca::Real tempReli = 0.0;
		if(tj > cc00Size - reliLeng 
			&& ti >= reliLeng && ti <= cc0Size - reliLeng){
			tempReli = 
				std::pow((tj - (cc00Size - reliLeng)) / (Ddpca::Real)reliLeng, reliExpo) 
				* reliAmou[0];
		}
		else if(tj > cc00Size - reliLeng && ti < reliLeng){
			Ddpca::Real tempAngl = 
				std::atan((tj - (cc00Size - reliLeng)) / (Ddpca::Real)(reliLeng - ti));
			Ddpca::Real tempRati = tempAngl / (Ddpca::PI / 2.0);
			Ddpca::Real maxiAmou = reliAmou[1] 
				+ (- 1.0 + std::cos(tempRati * Ddpca::PI)) * (reliAmou[1] - reliAmou[0]) / 2.0;
			Ddpca::Real tempRadi = std::sqrt(std::pow(tj - (cc00Size - reliLeng), 2.0) 
				+ std::pow(reliLeng - ti, 2.0));
			tempReli = std::pow(tempRadi / (Ddpca::Real)reliLeng, reliExpo) * maxiAmou;
		}
		else if(tj > cc00Size - reliLeng && ti > cc0Size - reliLeng){
			Ddpca::Real tempAngl = std::atan((tj - (cc00Size - reliLeng)) 
				/ (Ddpca::Real)(ti - (cc0Size - reliLeng)));
			Ddpca::Real tempRati = tempAngl / (Ddpca::PI / 2.0);
			Ddpca::Real maxiAmou = reliAmou[1] 
				+ (- 1.0 + std::cos(tempRati * Ddpca::PI)) * (reliAmou[1] - reliAmou[0]) / 2.0;
			Ddpca::Real tempRadi = std::sqrt(std::pow(tj - (cc00Size - reliLeng), 2.0) 
				+ std::pow(ti - (cc0Size - reliLeng), 2.0));
			tempReli = std::pow(tempRadi / (Ddpca::Real)reliLeng, reliExpo) * maxiAmou;
		}
		else if(tj <= cc00Size - reliLeng && ti < reliLeng){
			tempReli = std::pow((reliLeng - ti) / (Ddpca::Real)reliLeng, reliExpo) * reliAmou[1];
		}
		else if(tj <= cc00Size - reliLeng && ti > cc0Size - reliLeng){
			tempReli = 
				std::pow((ti - (cc0Size - reliLeng)) / (Ddpca::Real)reliLeng, reliExpo) 
				* reliAmou[1];
		}
		if(std::abs(tempReli) > 1.0E-12){
			Ddpca::Real tempRadi_0 = std::sqrt(std::pow(tempXYZ[0], 2.0) + std::pow(tempXYZ[1], 2.0));
			Ddpca::Real tempRadi = a_h2 - tempRadi_0;
			tempRadi = std::sqrt(tempRadi * tempRadi + tempXYZ[2] * tempXYZ[2]);
			Ddpca::Real tempAngl = tempReli / tempRadi;
			Ddpca::Real tempThet_0 = std::asin(tempXYZ[2] / tempRadi);
			Ddpca::Real tempThet_1 = tempThet_0 + tempAngl;
			Ddpca::Real tempRadi_1 = a_h2 - tempRadi * std::cos(tempThet_1);
			Ddpca::Real tempFact = tempRadi_1 / tempRadi_0;
			tempXYZ[0] = tempFact * tempXYZ[0];
			tempXYZ[1] = tempFact * tempXYZ[1];
			tempXYZ[2] = tempXYZ[2] + tempRadi * (std::sin(tempThet_1) - std::sin(tempThet_0));
			tempCoor = tempXYZ;
		}
	}
}

void DehwSurf::WHEE_RELI(std::array<Ddpca::Real,3> &tempCoor, Ddpca::I64 ti, Ddpca::I64 tj){
	Ddpca::Real reliLeng = 40;
	const Ddpca::I64 cc1Size = curvCoor[1].size();
	if(tj < reliLeng || ti < reliLeng || ti > cc1Size - reliLeng){
		std::array<Ddpca::Real,3> tempXYZ = tempCoor;
		Ddpca::Real reliExpo = 3.0;
		std::vector<Ddpca::Real> reliAmou = {12.0E-6, 16.0E-6};//0 - tip, 1 - end
		//
		Ddpca::Real tempReli = 0.0;
		if(ti < reliLeng){
			if(tj >= reliLeng){
				tempReli = std::pow((reliLeng - ti) / (Ddpca::Real)reliLeng, reliExpo) * reliAmou[1];
			}
			else{
				Ddpca::Real tempAngl = std::atan((reliLeng - tj) / (Ddpca::Real)(reliLeng - ti));
				Ddpca::Real tempRati = tempAngl / (Ddpca::PI / 2.0);
				Ddpca::Real maxiAmou = reliAmou[1] 
					+ (- 1.0 + std::cos(tempRati * Ddpca::PI)) * (reliAmou[1] - reliAmou[0]) / 2.0;
				Ddpca::Real tempRadi = std::sqrt(std::pow((reliLeng - tj), 2.0) + std::pow((reliLeng - ti), 2.0));
				tempReli = std::pow((tempRadi / reliLeng), reliExpo) * maxiAmou;
			}
		}
		else if(ti > cc1Size - reliLeng){
			if(tj >= reliLeng){
				tempReli = 
					std::pow((ti - (cc1Size - reliLeng)) / (Ddpca::Real)reliLeng, reliExpo) 
					* reliAmou[1];
			}
			else{
				Ddpca::Real tempAngl = std::atan((reliLeng - tj) 
					/ (Ddpca::Real)(ti - (cc1Size - reliLeng)));
				Ddpca::Real tempRati = tempAngl / (Ddpca::PI / 2.0);
				Ddpca::Real maxiAmou = reliAmou[1] 
					+ (- 1.0 + std::cos(tempRati * Ddpca::PI)) * (reliAmou[1] - reliAmou[0]) / 2.0;
				Ddpca::Real tempRadi = std::sqrt(std::pow((reliLeng - tj), 2.0) 
					+ std::pow((ti - (cc1Size - reliLeng)), 2.0));
				tempReli = std::pow((tempRadi / reliLeng), reliExpo) * maxiAmou;
			}
		}
		else if(tj < reliLeng){
			tempReli = std::pow((reliLeng - tj) / (Ddpca::Real)reliLeng, reliExpo) * reliAmou[0];
		}
		if(std::abs(tempReli) > 1.0E-12){
			Ddpca::Real tempAngl = tempReli 
				/ std::sqrt(tempXYZ[0] * tempXYZ[0] + tempXYZ[1] * tempXYZ[1]);
			Ddpca::DenseMatrix tempMatr(3,3,{
				std::cos(tempAngl),std::sin(tempAngl),0.0,
				-std::sin(tempAngl),std::cos(tempAngl),0.0,
				0.0,0.0,1.0
			});
			Ddpca::GEMV(tempMatr,tempXYZ,tempCoor);
		}
	}
}

void DehwSurf::NEW_CONT_ZONE(Ddpca::I64 f_lr){
	//initial value
	Ddpca::Real thet_c, thet_h = 0.0;
	std::array<Ddpca::Real,3> r_2_2;
	Ddpca::Real angl_fi, radi_fi, R_fmini, R_fmaxi;
	Ddpca::I64 numb_c = 1000;
	Ddpca::I64 numb_h = 10000;
	Ddpca::Real thet_cL = 0.01 * Ddpca::PI;
	Ddpca::Real thet_cH = 0.49 * Ddpca::PI;
	Ddpca::Real epsl_t = 1.0E-8;
	Ddpca::I64 F_init = 0;
	for(Ddpca::I64 ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (Ddpca::Real)numb_c * ti;
		Ddpca::Real thet_hs, thet_hm;
		SINGULAR_C2H(thet_c, thet_hs, thet_hm);
		Ddpca::Real thet_hL, thet_hH;
		if(f_lr == 1){
			thet_hL = thet_hs + epsl_t;
			thet_hH = thet_hm - epsl_t;
		}
		else{
			thet_hL = thet_hm + epsl_t;
			thet_hH = thet_hs + 2.0 * Ddpca::PI - epsl_t;
		}
		if(thet_hL >= thet_hH){
			continue;
		}
		for(Ddpca::I64 tj = 0; tj <= numb_h && F_init == 0; tj ++){
			if(f_lr == 1){
				thet_h = thet_hH - (thet_hH - thet_hL) / (Ddpca::Real)numb_h * tj;
			}
			else{
				thet_h = thet_hL + (thet_hH - thet_hL) / (Ddpca::Real)numb_h * tj;
			}
			//????????????????????cillofe, !!!!!the choosing order of thet_h!!!!!
			Ddpca::Real thet_1 = i_1c * thet_c;
			Ddpca::Real x_d, y_d;
			FSME(thet_1, thet_h, x_d, y_d);
			WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2);
			Ddpca::Real kapp_h2N;
			if(CILFOSE_NI(thet_1, thet_h, kapp_h2N) > 0.0){
				WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
				if(-widtAngl <= angl_fi && angl_fi <= widtAngl 
					&& R_fmini <= radi_fi && radi_fi <= R_fmaxi){
					F_init = 1;
				}
			}
		}
	}
	if(F_init == 0){
		std::cout << "WARNING in DEHWSURF::NEW_CONT_ZONE" << f_lr << "!" << std::endl;
		return;
	}
	//closest point
	Ddpca::Real miniDist = 1.0E20;
	std::list<INDE_INIT> breaList;
	std::vector<std::vector<Ddpca::I64>> F_sear(curvCoor[1].size(), 
		std::vector<Ddpca::I64>(curvCoor[1][0].size(), 0));
	const Ddpca::I64 cc1Size = curvCoor[1].size();
	for(Ddpca::I64 ti = 0; ti < cc1Size - 1; ti ++){
		const Ddpca::I64 cc1tSize = curvCoor[1][ti].size();
		for(Ddpca::I64 tj = 0; tj < cc1tSize - 1; tj ++){
			Ddpca::Real epsl_x = (curvCoor[1][ti][tj][0] - curvCoor[1][ti + 1][tj][0]) / 4.0;
			Ddpca::Real epsl_y = (curvCoor[1][ti][tj + 1][1] - curvCoor[1][ti][tj][1]) / 4.0;
			if(curvCoor[1][ti + 1][tj][0] - epsl_x <= angl_fi && 
				angl_fi <= curvCoor[1][ti][tj][0] + epsl_x && 
				curvCoor[1][ti][tj][1] - epsl_y <= radi_fi && 
				radi_fi <= curvCoor[1][ti][tj + 1][1] + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, thet_h));
				F_sear[ti][tj] = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, thet_h));
				F_sear[ti + 1][tj] = 1;
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, thet_h));
				F_sear[ti][tj + 1] = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, thet_h));
				F_sear[ti + 1][tj + 1] = 1;
			}
			Ddpca::Real dist_ij = radi_fi * std::abs(angl_fi - curvCoor[1][ti][tj][0]) 
				+ std::abs(radi_fi - curvCoor[1][ti][tj][1]);
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	std::cout << "DEHWSURF::NEW_CONT_ZONE" << f_lr << ", miniDist = " << miniDist 
		<< ", initial number = " << breaList.size() << std::endl;
	//breadth-first search
	Ddpca::Real epsl_d = 1.0E-9;
	Ddpca::I64 coun_w = 0;
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		F_sear[temp_ij.ti][temp_ij.tj] = 2;
		//
		Ddpca::Real x_d, y_d;
		WHEE_CURV_2_CART_1(curvCoor[1][temp_ij.ti][temp_ij.tj][0], 
			curvCoor[1][temp_ij.ti][temp_ij.tj][1], 
			r_2_2, temp_ij.init_1, temp_ij.init_2, f_lr, x_d, y_d
		);
		//
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		Ddpca::Real dist_ij = 
			radi_fi * std::abs(angl_fi - curvCoor[1][temp_ij.ti][temp_ij.tj][0]) 
			+ std::abs(radi_fi - curvCoor[1][temp_ij.ti][temp_ij.tj][1]);
		if(dist_ij < epsl_d){
			//
			thet_c = temp_ij.init_1;
			Ddpca::Real thet_1 = i_1c * thet_c;
			std::array<Ddpca::Real,3> r_c_c, T_o1_oc, r_1_o1;
			Ddpca::DenseMatrix R_oc_c(3,3), R_o1_oc(3,3), R_1_o1(3,3);
			r_c_c = {- x_d, r_b2 - y_d * std::sin(beta_c), y_d * std::cos(beta_c)};
			R_oc_c.Fill(
				std::cos(thet_c),std::sin(thet_c),0.0,
				- std::sin(thet_c),std::cos(thet_c),0.0,
				0.0,0.0,1.0
			);
			R_o1_oc.Fill(
				1.0,0.0,0.0,
				0.0,0.0,1.0,
				0.0,-1.0,0.0
			);
			R_1_o1.Fill(
				std::cos(thet_1),- std::sin(thet_1),0.0,
				std::sin(thet_1),std::cos(thet_1),0.0,
				0.0,0.0,1.0
			);
			T_o1_oc = {a_1c, 0.0, 0.0};
			Ddpca::DenseMatrix tempMatr(3,3);
			Ddpca::GEMM(R_o1_oc, R_oc_c, tempMatr);
			Ddpca::GEMV(tempMatr, r_c_c, r_1_o1);
			Ddpca::XPEY(r_1_o1, T_o1_oc);
			Ddpca::Real woxi_11 =  thet_1 - std::atan2(r_1_o1[1], r_1_o1[0]);
			if(curvCoor[0][0][0][0] - 1.0E-12 <= woxi_11 
				&& woxi_11 <= curvCoor[0][curvCoor[0].size() - 1][0][0] + 1.0E-12){
				WHEE_PHAS(temp_ij.ti, temp_ij.tj, f_lr, r_2_2);
			}
			//
			std::array<std::array<Ddpca::I64,2>,8> inde = {{
				{{temp_ij.ti - 1, temp_ij.tj - 1}}, 
				{{temp_ij.ti, temp_ij.tj - 1}}, 
				{{temp_ij.ti + 1, temp_ij.tj - 1}}, 
				{{temp_ij.ti - 1, temp_ij.tj}}, 
				{{temp_ij.ti + 1, temp_ij.tj}}, 
				{{temp_ij.ti - 1, temp_ij.tj + 1}}, 
				{{temp_ij.ti, temp_ij.tj + 1}}, 
				{{temp_ij.ti + 1, temp_ij.tj + 1}}
			}};
			const Ddpca::I64 indeSize = inde.size();
			const Ddpca::I64 cc1Size = curvCoor[1].size();
			const Ddpca::I64 cc10Size = curvCoor[1][0].size();
			for(Ddpca::I64 ti = 0; ti < indeSize; ti ++){
				if(inde[ti][0] >= 0 && inde[ti][1] >= 0 
					&& inde[ti][0] < cc1Size && inde[ti][1] < cc10Size 
					&& F_sear[inde[ti][0]][inde[ti][1]] == 0){
					F_sear[inde[ti][0]][inde[ti][1]] = 1;
					breaList.push_back(INDE_INIT(inde[ti][0], inde[ti][1], 
						temp_ij.init_1, temp_ij.init_2)
					);
				}
			}
		}
		coun_w ++;
		if(coun_w % 4000 == 0){
			std::cout << coun_w << std::endl;
		}
	}
}

void DehwSurf::FORMER_CONT_ZONE(){
	//initial value
	Ddpca::Real thet_c, x_d;
	std::array<Ddpca::Real,3> r_c_c;
	Ddpca::Real angl_fi, radi_fi, R_fmini, R_fmaxi;
	Ddpca::I64 numb_c = 1000;
	Ddpca::Real thet_cL = 0.01 * Ddpca::PI;
	Ddpca::Real thet_cH = 0.49 * Ddpca::PI;
	Ddpca::I64 numb_d = 10000;
	Ddpca::Real x_dL = - 10.0 * a_1c;
	Ddpca::Real x_dH = 10.0 * a_1c;
	Ddpca::I64 F_init = 0;
	for(Ddpca::I64 ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (Ddpca::Real)numb_c * ti;
		for(Ddpca::I64 tj = 1; tj < numb_d && F_init == 0; tj ++){
			x_d = x_dL + (x_dH - x_dL) / (Ddpca::Real)numb_d * tj;
			Ddpca::Real y_d = - ((- std::sin(beta_c) * std::cos(thet_c) - i_c1 * std::cos(beta_c)) * x_d 
				- r_b2 * std::sin(beta_c) * std::sin(thet_c) + a_1c * std::sin(beta_c)) / std::sin(thet_c);
			r_c_c = {- x_d, r_b2 - y_d * std::sin(beta_c), y_d * std::cos(beta_c)};
			WHEE_G2L(r_c_c, angl_fi, radi_fi, R_fmini, R_fmaxi);
			if(-widtAngl <= angl_fi && angl_fi <= widtAngl 
				&& R_fmini <= radi_fi && radi_fi <= R_fmaxi){
				F_init = 1;
			}
		}
	}
	if(F_init == 0){
		std::cout << "WARNING in DEHWSURF::FORMER_CONT_ZONE!" << std::endl;
		return;
	}
	//closest point
	Ddpca::Real miniDist = 1.0E20;
	std::list<INDE_INIT> breaList;
	std::vector<std::vector<Ddpca::I64>> F_sear(curvCoor[1].size(),
		std::vector<Ddpca::I64>(curvCoor[1][0].size(), 0));
	const Ddpca::I64 cc1Size = curvCoor[1].size();
	for(Ddpca::I64 ti = 0; ti < cc1Size - 1; ti ++){
		const Ddpca::I64 cc10Size = curvCoor[1][0].size();
		for(Ddpca::I64 tj = 0; tj < cc10Size - 1; tj ++){
			Ddpca::Real epsl_x = (curvCoor[1][ti][tj][0] - curvCoor[1][ti + 1][tj][0]) / 4.0;
			Ddpca::Real epsl_y = (curvCoor[1][ti][tj + 1][1] - curvCoor[1][ti][tj][1]) / 4.0;
			if(curvCoor[1][ti + 1][tj][0] - epsl_x <= angl_fi && 
				angl_fi <= curvCoor[1][ti][tj][0] + epsl_x && 
				curvCoor[1][ti][tj][1] - epsl_y <= radi_fi && 
				radi_fi <= curvCoor[1][ti][tj + 1][1] + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, x_d));
				F_sear[ti][tj] = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, x_d));
				F_sear[ti + 1][tj] = 1;
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, x_d));
				F_sear[ti][tj + 1] = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, x_d));
				F_sear[ti + 1][tj + 1] = 1;
			}
			Ddpca::Real dist_ij = radi_fi * std::abs(angl_fi - curvCoor[1][ti][tj][0]) 
				+ std::abs(radi_fi - curvCoor[1][ti][tj][1]);
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	std::cout << "DEHWSURF::FORMER_CONT_ZONE, miniDist = " << miniDist 
		<< ", initial number = " << breaList.size() << std::endl;
	//breadth-first search
	Ddpca::Real epsl_d = 1.0E-9;
	Ddpca::I64 coun_w = 0;
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		F_sear[temp_ij.ti][temp_ij.tj] = 2;
		//
		Ddpca::Real y_d;
		WHEE_CURV_2_CART_2(curvCoor[1][temp_ij.ti][temp_ij.tj][0], 
			curvCoor[1][temp_ij.ti][temp_ij.tj][1], 
			r_c_c, temp_ij.init_1, temp_ij.init_2, y_d
		);
		//
		WHEE_G2L(r_c_c, angl_fi, radi_fi, R_fmini, R_fmaxi);
		Ddpca::Real dist_ij = 
			radi_fi * std::abs(angl_fi - curvCoor[1][temp_ij.ti][temp_ij.tj][0]) 
			+ std::abs(radi_fi - curvCoor[1][temp_ij.ti][temp_ij.tj][1]);
		if(dist_ij < epsl_d){
			//
			thet_c = temp_ij.init_1;
			Ddpca::Real thet_1 = i_1c * thet_c;
			std::array<Ddpca::Real,3> T_o1_oc, r_1_o1;
			Ddpca::DenseMatrix R_oc_c(3,3), R_o1_oc(3,3), R_1_o1(3,3), tempMatr(3,3);
			R_oc_c.Fill(
				std::cos(thet_c),std::sin(thet_c),0.0,
				- std::sin(thet_c),std::cos(thet_c),0.0,
				0.0,0.0,1.0
			);
			R_o1_oc.Fill(
				1.0,0.0,0.0,
				0.0,0.0,1.0,
				0.0,-1.0,0.0
			);
			R_1_o1.Fill(
				std::cos(thet_1),- std::sin(thet_1),0.0,
				std::sin(thet_1),std::cos(thet_1),0.0,
				0.0,0.0,1.0
			);
			T_o1_oc = {a_1c, 0.0, 0.0};
			Ddpca::GEMM(R_o1_oc, R_oc_c, tempMatr);
			Ddpca::GEMV(tempMatr, r_c_c, r_1_o1);
			Ddpca::XPEY(r_1_o1, T_o1_oc);
			Ddpca::Real woxi_11 =  thet_1 - std::atan2(r_1_o1[1], r_1_o1[0]);
			if(curvCoor[0][0][0][0] - 1.0E-12 <= woxi_11 
				&& woxi_11 <= curvCoor[0][curvCoor[0].size() - 1][0][0] + 1.0E-12){
				WHEE_PHAS(temp_ij.ti, temp_ij.tj, 3, r_c_c);
			}
			//
			std::array<std::array<Ddpca::I64,2>,8> inde = {{
				{{temp_ij.ti - 1, temp_ij.tj - 1}},
				{{temp_ij.ti, temp_ij.tj - 1}},
				{{temp_ij.ti + 1, temp_ij.tj - 1}},
				{{temp_ij.ti - 1, temp_ij.tj}},
				{{temp_ij.ti + 1, temp_ij.tj}},
				{{temp_ij.ti - 1, temp_ij.tj + 1}},
				{{temp_ij.ti, temp_ij.tj + 1}},
				{{temp_ij.ti + 1, temp_ij.tj + 1}}
			}};
			const Ddpca::I64 indeSize = inde.size();
			const Ddpca::I64 cc1Size = curvCoor[1].size();
			const Ddpca::I64 cc10Size = curvCoor[1][0].size();
			for(Ddpca::I64 ti = 0; ti < indeSize; ti ++){
				if(inde[ti][0] >= 0 && inde[ti][1] >= 0 
					&& inde[ti][0] < cc1Size && inde[ti][1] < cc10Size 
					&& F_sear[inde[ti][0]][inde[ti][1]] == 0){
					F_sear[inde[ti][0]][inde[ti][1]] = 1;
					breaList.push_back(INDE_INIT(inde[ti][0], inde[ti][1], 
						temp_ij.init_1, temp_ij.init_2)
					);
				}
			}
		}
		coun_w ++;
		if(coun_w % 4000 == 0){
			std::cout << coun_w << std::endl;
		}
	}
}

void DehwSurf::TRANSITION_ZONE(Ddpca::I64 f_hr){
	Ddpca::Real xi_11;
	if(f_hr == 1){
		xi_11 = curvCoor[0][0][0][0];//head
	}
	else{
		xi_11 = curvCoor[0][curvCoor[0].size() - 1][0][0];//rear
	}
	//initial value
	Ddpca::Real thet_c, thet_h;
	std::array<Ddpca::Real,3> r_2_2, r_1_1;
	Ddpca::Real angl_fi, radi_fi, R_fmini, R_fmaxi;
	Ddpca::I64 numb_c = 1000;
	Ddpca::Real thet_cL, thet_cH;
	WORM_CURV_2_CART(xi_11, a_h2 - d_f[0] / 2.0, r_1_1, thet_cL);
	WORM_CURV_2_CART(xi_11, d_f[1] / 2.0, r_1_1, thet_cH);
	thet_h = xi_11;
	Ddpca::I64 F_init = 0;
	for(Ddpca::I64 ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (Ddpca::Real)numb_c * ti;
		Ddpca::DenseMatrix Dr_2_2(3,2);
		WHEE_TRAN(thet_c, thet_h, xi_11, r_2_2, Dr_2_2);
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		if(-widtAngl <= angl_fi && angl_fi <= widtAngl 
			&& R_fmini <= radi_fi && radi_fi <= R_fmaxi){
			F_init = 1;
		}
	}
	if(F_init == 0){
		std::cout << "WARNING in DEHWSURF::TRANSITION_ZONE" << f_hr << "!" << std::endl;
		return;
	}
	//closest point
	Ddpca::Real miniDist = 1.0E20;
	std::list<INDE_INIT> breaList;
	std::vector<std::vector<Ddpca::I64>> F_sear(curvCoor[1].size(), 
		std::vector<Ddpca::I64>(curvCoor[1][0].size(), 0));
	const Ddpca::I64 cc1Size = curvCoor[1].size();
	for(Ddpca::I64 ti = 0; ti < cc1Size - 1; ti ++){
		const Ddpca::I64 cc1tSize = curvCoor[1][ti].size();
		for(Ddpca::I64 tj = 0; tj < cc1tSize - 1; tj ++){
			Ddpca::Real epsl_x = (curvCoor[1][ti][tj][0] - curvCoor[1][ti + 1][tj][0]) / 4.0;
			Ddpca::Real epsl_y = (curvCoor[1][ti][tj + 1][1] - curvCoor[1][ti][tj][1]) / 4.0;
			if(curvCoor[1][ti + 1][tj][0] - epsl_x <= angl_fi && 
				angl_fi <= curvCoor[1][ti][tj][0] + epsl_x && 
				curvCoor[1][ti][tj][1] - epsl_y <= radi_fi && 
				radi_fi <= curvCoor[1][ti][tj + 1][1] + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, thet_h));
				F_sear[ti][tj] = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, thet_h));
				F_sear[ti + 1][tj] = 1;
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, thet_h));
				F_sear[ti][tj + 1] = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, thet_h));
				F_sear[ti + 1][tj + 1] = 1;
			}
			Ddpca::Real dist_ij = radi_fi * std::abs(angl_fi - curvCoor[1][ti][tj][0]) 
				+ std::abs(radi_fi - curvCoor[1][ti][tj][1]);
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	std::cout << "DEHWSURF::TRANSITION_ZONE" << f_hr << ", miniDist = " << miniDist 
		<< ", initial number = " << breaList.size() << std::endl;
	//breadth-first search
	Ddpca::Real epsl_d = 1.0E-9;
	Ddpca::I64 coun_w = 0;
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		F_sear[temp_ij.ti][temp_ij.tj] = 2;
		//
		WHEE_CURV_2_CART_3(curvCoor[1][temp_ij.ti][temp_ij.tj][0], 
			curvCoor[1][temp_ij.ti][temp_ij.tj][1], 
			r_2_2, temp_ij.init_1, temp_ij.init_2, xi_11
		);
		//
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		Ddpca::Real dist_ij = 
			radi_fi * std::abs(angl_fi - curvCoor[1][temp_ij.ti][temp_ij.tj][0]) 
			+ std::abs(radi_fi - curvCoor[1][temp_ij.ti][temp_ij.tj][1]);
		if(dist_ij < epsl_d){
			WHEE_PHAS(temp_ij.ti, temp_ij.tj, 3 + f_hr, r_2_2);
			std::array<std::array<Ddpca::I64,2>,8> inde = {{
				{{temp_ij.ti - 1, temp_ij.tj - 1}},
				{{temp_ij.ti, temp_ij.tj - 1}},
				{{temp_ij.ti + 1, temp_ij.tj - 1}},
				{{temp_ij.ti - 1, temp_ij.tj}},
				{{temp_ij.ti + 1, temp_ij.tj}},
				{{temp_ij.ti - 1, temp_ij.tj + 1}},
				{{temp_ij.ti, temp_ij.tj + 1}},
				{{temp_ij.ti + 1, temp_ij.tj + 1}}
			}};
			const Ddpca::I64 indeSize = inde.size();
			const Ddpca::I64 cc1Size = curvCoor[1].size();
			const Ddpca::I64 cc10Size = curvCoor[1][0].size();
			for(Ddpca::I64 ti = 0; ti < indeSize; ti ++){
				if(inde[ti][0] >= 0 && inde[ti][1] >= 0 
					&& inde[ti][0] < cc1Size && inde[ti][1] < cc10Size 
					&& F_sear[inde[ti][0]][inde[ti][1]] == 0){
					F_sear[inde[ti][0]][inde[ti][1]] = 1;
					breaList.push_back(INDE_INIT(inde[ti][0], inde[ti][1], 
						temp_ij.init_1, temp_ij.init_2)
					);
				}
			}
		}
		coun_w ++;
		if(coun_w % 4000 == 0){
			std::cout << coun_w << std::endl;
		}
	}
}

void DehwSurf::WORM_ROOT_RADIUS(Ddpca::I64 flag, Ddpca::DenseMatrix tempPoin, 
	std::array<Ddpca::Real,2> &tempCent, Ddpca::Real &tempRadi, std::array<Ddpca::Real,2> &tempAngl){
	//
	std::array<Ddpca::Real,2> vect_1;
	vect_1[0] = tempPoin(0,1) - tempPoin(0,0);
	vect_1[1] = tempPoin(1,1) - tempPoin(1,0);
	Ddpca::NORMALIZE(vect_1);
	std::array<Ddpca::Real,2> vect_2;
	vect_2[0] = tempPoin(0,2) - tempPoin(0,0);
	vect_2[1] = tempPoin(1,2) - tempPoin(1,0);
	Ddpca::Real tempLeng_1 = Ddpca::DOT(vect_2, vect_1);
	Ddpca::Real tempLeng_2 = std::sqrt(vect_2[0] * vect_2[0] + vect_2[1] * vect_2[1] 
		- tempLeng_1 * tempLeng_1);
	Ddpca::Real targVari = tempLeng_1 / (R_f[0] - tempLeng_2);
	Ddpca::Real middAngl = std::asin(targVari / std::sqrt(1.0 + targVari * targVari)) - std::atan(1.0 / targVari);
	tempRadi = R_f[0] - tempLeng_1 / std::cos(middAngl);
	std::array<Ddpca::Real,2> tempVect;
	tempVect = {flag * vect_1[1], - flag * vect_1[0]};
	tempCent[0] = tempPoin(0,0) + tempRadi * tempVect[0];
	tempCent[1] = tempPoin(1,0) + tempRadi * tempVect[1];
	tempAngl[0] = std::atan2(- tempVect[1], - tempVect[0]);
	tempAngl[1] = tempAngl[0] + flag * (Ddpca::PI / 2.0 - middAngl);
}

void DehwSurf::WORM_ROOT(Ddpca::I64 indi, Ddpca::I64 flag, Ddpca::DenseMatrix &rootProf){
	//
	Ddpca::DenseMatrix tempPoin(3,3);
	if(flag == 1){
		tempPoin(0,0) = wormTosu.indexPoint[indi][0][0];
		tempPoin(1,0) = wormTosu.indexPoint[indi][0][1];
		tempPoin(2,0) = wormTosu.indexPoint[indi][0][2];
		tempPoin(0,1) = wormTosu.indexPoint[indi][1][0];
		tempPoin(1,1) = wormTosu.indexPoint[indi][1][1];
		tempPoin(2,1) = wormTosu.indexPoint[indi][1][2];
	}
	else{
		tempPoin(0,0) = wormToba.indexPoint[indi][0][0];
		tempPoin(1,0) = wormToba.indexPoint[indi][0][1];
		tempPoin(2,0) = wormToba.indexPoint[indi][0][2];
		tempPoin(0,1) = wormToba.indexPoint[indi][1][0];
		tempPoin(1,1) = wormToba.indexPoint[indi][1][1];
		tempPoin(2,1) = wormToba.indexPoint[indi][1][2];
	}
	Ddpca::Real tempXi11 = - curvCoor[0][indi][0][0];//(wormCurv[0] 
		// + (wormCurv[2] - wormCurv[0]) / (wormTosu.indexPoint.size() - 1) * indi);
	tempPoin(0,2) = a_h2 * std::cos(tempXi11);
	tempPoin(1,2) = a_h2 * std::sin(tempXi11);
	tempPoin(2,2) = 0.0;
	Ddpca::DenseMatrix tempPoin_(2,3);
	Ddpca::DenseMatrix tempMatr(3,3,{
		0.0,-std::cos(tempXi11),std::sin(tempXi11),
		0.0,-std::sin(tempXi11),-std::cos(tempXi11),
		1.0,0.0,0.0
	});
	for(Ddpca::I64 ti = 0; ti < 3; ti ++){
		std::array<Ddpca::Real,3> tempPoin_i, tempVect;
		tempVect[0] = tempPoin(0,ti);
		tempVect[1] = tempPoin(1,ti);
		tempVect[2] = tempPoin(2,ti);
		Ddpca::GEMV(tempMatr, tempVect, tempPoin_i);
		tempPoin_(0,ti) = tempPoin_i[0];
		tempPoin_(1,ti) = tempPoin_i[1] + a_h2;
	}
	std::array<Ddpca::Real,2> tempCent, tempAngl, tempPoinArce;
	Ddpca::Real tempRadi;
	WORM_ROOT_RADIUS(flag, tempPoin_, tempCent, tempRadi, tempAngl);
	tempPoinArce = {tempCent[0] + tempRadi * std::cos(tempAngl[1]),
		tempCent[1] + tempRadi * std::sin(tempAngl[1])};
	Ddpca::Real tempAnglArce = std::atan2(tempPoinArce[1], tempPoinArce[0]);
	//
	Ddpca::Real tempAnglRoot;
	Ddpca::Real tempAnglStar = std::acos(r_b2 / (d[1] / 2.0))
		- i_2h * tempXi11 - tootThicAngl[0] / 2.0;
	if(flag == 1){
		tempAnglRoot = tempAnglStar + pitcAngl / 2.0;
	}
	else{
		tempAnglRoot = tempAnglStar - pitcAngl / 2.0;
	}
	//
	Ddpca::Real sumLeng = flag * R_f[0] * (tempAnglRoot - tempAnglArce)
		+ flag * tempRadi * (tempAngl[1] - tempAngl[0]);
	std::array<Ddpca::Real,3> tempTran;
	tempTran = {a_h2 * std::cos(tempXi11), a_h2 * std::sin(tempXi11), 0.0};
	const Ddpca::I64 roprCols = rootProf.cols;
	for(Ddpca::I64 ti = 0; ti < roprCols; ti ++){
		Ddpca::Real leng_i = sumLeng / (roprCols - 1) * (Ddpca::Real)ti;
		std::array<Ddpca::Real,3> poin_i, tempVect;
		if(leng_i <= flag * R_f[0] * (tempAnglRoot - tempAnglArce)){
			Ddpca::Real angl_i = tempAnglRoot - flag * leng_i / R_f[0];
			poin_i = {R_f[0] * std::cos(angl_i), R_f[0] * std::sin(angl_i), 0.0};
		}
		else{
			leng_i = leng_i - flag * R_f[0] * (tempAnglRoot - tempAnglArce);
			Ddpca::Real angl_i = tempAngl[1] - flag * leng_i / tempRadi;
			poin_i = {tempCent[0] + tempRadi * std::cos(angl_i), 
				tempCent[1] + tempRadi * std::sin(angl_i), 0.0};
		}
		Ddpca::GEMTV(tempMatr, poin_i, tempVect);
		Ddpca::XPEY(tempVect, tempTran);
		rootProf(0,ti) = tempVect[0];
		rootProf(1,ti) = tempVect[1];
		rootProf(2,ti) = tempVect[2];
	}
}

std::array<Ddpca::Real,2> DehwSurf::WHEE_UNCONE(std::array<Ddpca::Real,3> tempXYZ, Ddpca::Real tempAlph_3){
	Ddpca::Real r_2 = std::sqrt(std::pow(tempXYZ[0], 2.0) + std::pow(tempXYZ[1], 2.0));
	Ddpca::Real r_1 = r_2 / std::cos(tempAlph_3);
	Ddpca::Real alph_2 = std::atan2(tempXYZ[1], tempXYZ[0]);
	Ddpca::Real alph_1 = r_2 * alph_2 / r_1;
	std::array<Ddpca::Real,2> tempXY;
	tempXY = {r_1 * std::cos(alph_1), r_1 * std::sin(alph_1)};
	return tempXY;
}

std::array<Ddpca::Real,3> DehwSurf::WHEE_CONE(std::array<Ddpca::Real,2> tempXY, Ddpca::Real tempAlph_3){
	Ddpca::Real r_1 = std::sqrt(std::pow(tempXY[0], 2.0) + std::pow(tempXY[1], 2.0));
	Ddpca::Real alph_1 = std::atan2(tempXY[1], tempXY[0]);
	Ddpca::Real r_2 = r_1 * std::cos(tempAlph_3);
	Ddpca::Real alph_2 = r_1 * alph_1 / r_2;
	Ddpca::Real r_3 = a_h2 / std::cos(tempAlph_3) - r_1;
	std::array<Ddpca::Real,3> tempXYZ;
	tempXYZ = {r_2 * std::cos(alph_2), r_2 * std::sin(alph_2), r_3 * std::sin(tempAlph_3)};
	return tempXYZ;
}

void DehwSurf::WHEE_ROOT(Ddpca::I64 indi, Ddpca::I64 flag, Ddpca::DenseMatrix &rootProf){
	//
	std::array<std::array<std::array<Ddpca::Real,2>,2>,2> tempProf;
	Ddpca::Real alph_3 = - curvCoor[1][indi][0][0];
		// (widtAngl - 2.0 * widtAngl / (wheeTosu.indexPoint.size() - 1) * (Ddpca::Real)indi);
	Ddpca::Real angl_ai = alph_3 - std::asin(offsR_a * std::sin(alph_3) / R_a[1]);
	Ddpca::Real R_fmini = R_t[1];
	Ddpca::Real R_fmaxi = (R_a[1] * std::cos(angl_ai) - offsR_a) / std::cos(alph_3);
	for(Ddpca::I64 ti = 0; ti <= 1; ti ++){//tooth surface or tooth back
		for(Ddpca::I64 tj = 0; tj <= 1; tj ++){
			std::array<Ddpca::Real,3> tempXYZ;
			if(ti == 0){
				tempXYZ = {wheeTosu.indexPoint[indi][tj][0], 
					wheeTosu.indexPoint[indi][tj][1], wheeTosu.indexPoint[indi][tj][2]};
			}
			else{
				tempXYZ = {wheeToba.indexPoint[indi][tj][0], 
					wheeToba.indexPoint[indi][tj][1], wheeToba.indexPoint[indi][tj][2]};
			}
			Ddpca::Real r_2 = std::sqrt(std::pow(tempXYZ[0], 2.0) + std::pow(tempXYZ[1], 2.0));
			Ddpca::Real alph_2 = std::atan2(tempXYZ[1], tempXYZ[0]);
			Ddpca::Real r_3 = R_fmini 
				+ (R_fmaxi - R_fmini) / (wheeTosu.indexPoint[0].size() - 1) * (Ddpca::Real)tj;
			Ddpca::Real r_1 = a_h2 / std::cos(alph_3) - r_3;
			Ddpca::Real alph_1 = r_2 * alph_2 / r_1;
			tempProf[tj][ti] = {r_1 * std::cos(alph_1), r_1 * std::sin(alph_1)};
		}
	}
	Ddpca::Real r_f = a_h2 / std::cos(alph_3) - (a_h2 - d_f[1] / 2.0);
	Ddpca::Real tempPitc = pitcAngl * std::cos(alph_3);
	//
	Ddpca::DenseMatrix tempPoin(2,3);
	tempPoin(0,0) = tempProf[0][flag][0];
	tempPoin(1,0) = tempProf[0][flag][1];
	tempPoin(0,1) = tempProf[1][flag][0];
	tempPoin(1,1) = tempProf[1][flag][1];
	tempPoin(0,2) = 0.0;
	tempPoin(1,2) = 0.0;
	std::array<Ddpca::Real,2> vect_1;
	vect_1 = {tempPoin(0,0) - tempPoin(0,1), tempPoin(1,0) - tempPoin(1,1)};
	Ddpca::NORMALIZE(vect_1);
	std::array<Ddpca::Real,2> vect_2;
	vect_2 = {tempPoin(0,2) - tempPoin(0,0), tempPoin(1,2) - tempPoin(1,0)};
	Ddpca::Real tempLeng_1 = Ddpca::DOT(vect_2, vect_1);
	Ddpca::Real tempLeng_2 = std::sqrt(vect_2[0] * vect_2[0] + vect_2[1] * vect_2[1] 
		- tempLeng_1 * tempLeng_1);
	Ddpca::Real targVari = tempLeng_1 / (r_f - tempLeng_2);
	Ddpca::Real middAngl = std::asin(targVari / std::sqrt(1.0 + targVari * targVari)) - std::atan(1.0 / targVari);
	Ddpca::Real tempRadi = tempLeng_1 / std::cos(middAngl) - r_f;
	Ddpca::Real tempSign = (flag == 0) ? 1.0 : -1.0;
	std::array<Ddpca::Real,2> tempVect, tempCent, tempAngl, tempPoinArce;
	tempVect = {-tempSign * vect_1[1], tempSign * vect_1[0]};
	tempCent[0] = tempPoin(0,0) + tempRadi * tempVect[0];
	tempCent[1] = tempPoin(1,0) + tempRadi * tempVect[1];
	tempAngl[0] = std::atan2(-tempVect[1], -tempVect[0]);
	tempAngl[1] = tempAngl[0] + tempSign * (Ddpca::PI / 2.0 - middAngl);
	tempPoinArce = {tempCent[0] + tempRadi * std::cos(tempAngl[1]),
		tempCent[1] + tempRadi * std::sin(tempAngl[1])};
	Ddpca::Real tempAnglArce = std::atan2(tempPoinArce[1], tempPoinArce[0]);
	//
	Ddpca::Real tempAnglRoot;
	std::array<Ddpca::Real,2> tempPoinOppo = tempProf[0][1 - flag];
	tempAnglRoot = (std::atan2(tempPoinOppo[1], tempPoinOppo[0]) 
		+ std::atan2(tempPoin(1,0), tempPoin(0,0))) / 2.0;
	tempAnglRoot -= tempSign * tempPitc / 2.0;
	//
	Ddpca::Real sumLeng = r_f * tempSign * (tempAnglArce - tempAnglRoot)
		+ tempRadi * tempSign * (tempAngl[1] - tempAngl[0]);
	const Ddpca::I64 roprCols = rootProf.cols;
	for(Ddpca::I64 ti = 0; ti < roprCols; ti ++){
		Ddpca::Real leng_i = sumLeng / (roprCols - 1) * ti;
		if(leng_i <= r_f * tempSign * (tempAnglArce - tempAnglRoot)){
			Ddpca::Real angl_i = tempAnglRoot + tempSign * leng_i / r_f;
			std::array<Ddpca::Real,2> tempResu;
			tempResu = {r_f * std::cos(angl_i), r_f * std::sin(angl_i)};
			std::array<Ddpca::Real,3> tempRopr = WHEE_CONE(tempResu, alph_3);
			rootProf(0,ti) = tempRopr[0];
			rootProf(1,ti) = tempRopr[1];
			rootProf(2,ti) = tempRopr[2];
		}
		else{
			leng_i = leng_i - r_f * tempSign * (tempAnglArce - tempAnglRoot);
			Ddpca::Real angl_i = tempAngl[1] - tempSign * leng_i / tempRadi;
			std::array<Ddpca::Real,2> tempResu;
			tempResu = {tempCent[0] + tempRadi * std::cos(angl_i), 
				tempCent[1] + tempRadi * std::sin(angl_i)};
			std::array<Ddpca::Real,3> tempRopr = WHEE_CONE(tempResu, alph_3);
			rootProf(0,ti) = tempRopr[0];
			rootProf(1,ti) = tempRopr[1];
			rootProf(2,ti) = tempRopr[2];
		}
	}
}

void DehwSurf::WORM_TS_GRID(){
	std::cout << "DehwSurf::WORM_TS_GRID\n";
	//
	Ddpca::Real domaCirc = 2.0 * Ddpca::PI / (Ddpca::Real)circNumb;
	Ddpca::Real deltTang = domaCirc / gridNumb[0][4];
	Ddpca::Real inteStar = wormCurv[1];
	while(inteStar - domaCirc >= wormCurv[0]){
		inteStar = inteStar - domaCirc;
	}
	gridNumb[0][5] = std::ceil((inteStar - wormCurv[0]) / deltTang);
	Ddpca::Real realStar = inteStar - gridNumb[0][5] * deltTang;
	Ddpca::Real inteEndi = wormCurv[1];
	while(inteEndi + domaCirc <= wormCurv[2]){
		inteEndi = inteEndi + domaCirc;
	}
	gridNumb[0][6] = std::floor((inteEndi - inteStar) / domaCirc + 1.0E-10) + 2;
	//curvilinear coordinate
	const Ddpca::I64 cc0Size_0 = (gridNumb[0][4] * (gridNumb[0][6] - 2) + gridNumb[0][5] * 2) 
		* (1 << (globInho + globHomo + locaLeve)) + 1;
	const Ddpca::I64 cc0Size_1 = gridNumb[0][3] * (1 << (globHomo + locaLeve)) + 1;
	curvCoor[0].resize(cc0Size_0);
	for(Ddpca::I64 ti = 0; ti < cc0Size_0; ++ ti){
		curvCoor[0][ti].resize(cc0Size_1);
	}
	deltTang /= (1 << (globInho + globHomo + locaLeve));
	for(Ddpca::I64 ti = 0; ti < cc0Size_0; ti ++){
		for(Ddpca::I64 tj = 0; tj < cc0Size_1; tj ++){
			curvCoor[0][ti][tj] = {realStar + ti * deltTang,
				R_t[0] + (R_a[0] - R_t[0]) / (cc0Size_1 - 1) * (Ddpca::Real)tj};
		}
	}
	//Cartesian coordinate
	cartCoor[0].resize(cc0Size_0);
	for(Ddpca::I64 ti = 0; ti < cc0Size_0; ++ ti){
		cartCoor[0][ti].resize(cc0Size_1);
	}
	for(Ddpca::I64 ti = 0; ti < cc0Size_0; ti ++){
		if(ti % 1000 == 0){
			std::cout << ti << "/" << cc0Size_0 << std::endl;
		}
		for(Ddpca::I64 tj = 0; tj < cc0Size_1; tj ++){
			Ddpca::Real thet_c;
			WORM_CURV_2_CART(curvCoor[0][ti][tj][0], 
				curvCoor[0][ti][tj][1], cartCoor[0][ti][tj], thet_c
			);
			if(reliSwit == 1){
				WORM_RELI(cartCoor[0][ti][tj], ti, tj);
			}
		}
	}
}

void DehwSurf::WHEE_TS_GRID(){
	std::cout << "DehwSurf::WHEE_TS_GRID\n";
	//curvilinear coordinate
	const Ddpca::I64 cc1Size_0 = gridNumb[1][4] * (1 << (globInho + globHomo + locaLeve)) + 1;
	const Ddpca::I64 cc1Size_1 = gridNumb[1][3] * (1 << (globHomo + locaLeve)) + 1;
	curvCoor[1].resize(cc1Size_0);
	for(Ddpca::I64 ti = 0; ti < cc1Size_0; ++ ti){
		curvCoor[1][ti].resize(cc1Size_1);
	}
	for(Ddpca::I64 ti = 0; ti < cc1Size_0; ti ++){
		Ddpca::Real angl_fi = widtAngl - 2.0 * widtAngl / (cc1Size_0 - 1) * (Ddpca::Real)ti;
		Ddpca::Real angl_ai = angl_fi - std::asin(offsR_a * std::sin(angl_fi) / R_a[1]);
		Ddpca::Real R_fmini = (R_a[1] * std::cos(angl_ai) - offsR_a) / std::cos(angl_fi);
		Ddpca::Real R_fmaxi = R_t[1];
		for(Ddpca::I64 tj = 0; tj < cc1Size_1; tj ++){
			curvCoor[1][ti][tj] = {angl_fi, 
				R_fmini + (R_fmaxi - R_fmini) / (cc1Size_1 - 1) * (Ddpca::Real)tj};
		}
	}
	//Cartesian coordinate	
	cartCoor[1].resize(cc1Size_0);
	fpha.resize(cc1Size_0);
	for(Ddpca::I64 ti = 0; ti < cc1Size_0; ++ ti){
		cartCoor[1][ti].resize(cc1Size_1);
		fpha[ti].assign(cc1Size_1,0);
	}
	NEW_CONT_ZONE(1);//left new contact zone
	NEW_CONT_ZONE(2);//right new contact zone
	if(modiTran == 0.0 && modiCent == 0.0){
		FORMER_CONT_ZONE();//former contact zone
	}
	TRANSITION_ZONE(1);//head transition zone
	TRANSITION_ZONE(2);//rear transition zone
	//tooth flank relief
	if(reliSwit == 1){
		for(Ddpca::I64 ti = 0; ti < cc1Size_0; ti ++){
			for(Ddpca::I64 tj = 0; tj < cc1Size_1; tj ++){
				WHEE_RELI(cartCoor[1][ti][tj], ti, tj);
			}
		}
	}
}

void DehwSurf::TOOT_SURF_GRID(){
	curvCoor.resize(2);
	cartCoor.resize(2);
	//worm tooth surface
	WORM_TS_GRID();
	const Ddpca::I64 cc0Size = curvCoor[0].size();
	const Ddpca::I64 cc00Size = curvCoor[0][0].size();
	wormTosu.Resize(cc0Size, cc00Size);
	for(Ddpca::I64 ti = 0; ti < cc0Size; ti ++){
		for(Ddpca::I64 tj = 0; tj < cc00Size; tj ++){
			std::array<Ddpca::Real,3> tempCoor = cartCoor[0][ti][tj];
			wormTosu.Insert(ti, tj, Ddpca::Coordinate(tempCoor[0], tempCoor[1], tempCoor[2]));
		}
	}
	//wheel tooth surface
	WHEE_TS_GRID();
	const Ddpca::I64 cc1Size = curvCoor[1].size();
	const Ddpca::I64 cc10Size = curvCoor[1][0].size();
	wheeTosu.Resize(cc1Size, cc10Size);
	for(Ddpca::I64 ti = 0; ti < cc1Size; ti ++){
		for(Ddpca::I64 tj = 0; tj < cc10Size; tj ++){
			std::array<Ddpca::Real,3> tempCoor = cartCoor[1][ti][tj];
			wheeTosu.Insert(cc1Size - 1 - ti, cc10Size - 1 - tj, 
				Ddpca::Coordinate(tempCoor[0], tempCoor[1], tempCoor[2])
			);
		}
	}
	//worm tooth back
	wormToba.Resize(cc0Size, cc00Size);
	Ddpca::DenseMatrix tempMatr(3,3,{
		std::cos(wormCurv[1]),std::sin(wormCurv[1]),0.0,
		-std::sin(wormCurv[1]),std::cos(wormCurv[1]),0.0,
		0.0,0.0,1.0
	});
	for(Ddpca::I64 ti = 0; ti < cc0Size; ti ++){
		for(Ddpca::I64 tj = 0; tj < cc00Size; tj ++){
			std::array<Ddpca::Real,3> tempXYZ = cartCoor[0][ti][tj];
			std::array<Ddpca::Real,3> tempXYZ_1;
			Ddpca::GEMV(tempMatr, tempXYZ, tempXYZ_1);
			tempXYZ = {tempXYZ_1[0], -tempXYZ_1[1], -tempXYZ[2]};
			Ddpca::GEMTV(tempMatr, tempXYZ, tempXYZ_1);
			wormToba.Insert(cc0Size - 1 - ti, tj, Ddpca::Coordinate(
				tempXYZ_1[0], tempXYZ_1[1], tempXYZ_1[2]
			));
		}
	}
	//wheel tooth back
	wheeToba.Resize(cc1Size, cc10Size);
	tempMatr.Fill(
		std::cos(backAngl[1]),-std::sin(backAngl[1]),0.0,
		-std::sin(backAngl[1]),-std::cos(backAngl[1]),0.0,
		0.0,0.0,-1.0
	);
	for(Ddpca::I64 ti = 0; ti < cc1Size; ti ++){
		for(Ddpca::I64 tj = 0; tj < cc10Size; tj ++){
			std::array<Ddpca::Real,3> tempXYZ;
			Ddpca::GEMV(tempMatr, cartCoor[1][cc1Size - 1 - ti][cc10Size - 1 - tj], tempXYZ);
			wheeToba.Insert(cc1Size - 1 - ti, tj, 
				Ddpca::Coordinate(tempXYZ[0], tempXYZ[1], tempXYZ[2]
			));
		}
	}
	// curvCoor[0].clear();
	cartCoor[0].clear();
	// curvCoor[1].clear();
	cartCoor[1].clear();
	fpha.clear();
}

void DehwSurf::ROOT_TRAN_GRID(){
	//worm
	Ddpca::I64 numb_0 = (curvCoor[0].size() - 1) / (1 << (locaLeve)) + 1;
	Ddpca::I64 numb_1 = (gridNumb[0][0] / 2) * (1 << (globHomo)) + 1;
	wormRtsu.Resize(numb_0, numb_1);
	wormRtba.Resize(numb_0, numb_1);
	for(Ddpca::I64 ti = 0; ti < numb_0; ti ++){
		Ddpca::DenseMatrix rootProf(3, numb_1);
		WORM_ROOT(ti * (1 << locaLeve), 1, rootProf);
		for(Ddpca::I64 tj = 0; tj < numb_1; tj ++){
			wormRtsu.Insert(ti, tj, Ddpca::Coordinate(rootProf(0,tj), rootProf(1,tj), rootProf(2,tj)));
		}
		WORM_ROOT(ti * (1 << locaLeve), -1, rootProf);
		for(Ddpca::I64 tj = 0; tj < numb_1; tj ++){
			wormRtba.Insert(ti, tj, Ddpca::Coordinate(rootProf(0,tj), rootProf(1,tj), rootProf(2,tj)));
		}
	}
	//wheel
	numb_0 = gridNumb[1][4] * (1 << (globInho + globHomo)) + 1;
	numb_1 = (gridNumb[1][0] / 2) * (1 << (globHomo)) + 1;
	wheeRtsu.Resize(numb_0, numb_1);
	wheeRtba.Resize(numb_0, numb_1);
	for(Ddpca::I64 ti = 0; ti < numb_0; ti ++){
		Ddpca::DenseMatrix rootProf(3, numb_1);
		WHEE_ROOT(ti * (1 << locaLeve), 0, rootProf);
		for(Ddpca::I64 tj = 0; tj < numb_1; tj ++){
			wheeRtsu.Insert(ti, tj, Ddpca::Coordinate(rootProf(0,tj), rootProf(1,tj), rootProf(2,tj)));
		}
		WHEE_ROOT(ti * (1 << locaLeve), 1, rootProf);
		for(Ddpca::I64 tj = 0; tj < numb_1; tj ++){
			wheeRtba.Insert(ti, tj, Ddpca::Coordinate(rootProf(0,tj), rootProf(1,tj), rootProf(2,tj)));
		}
	}
}

void DehwSurf::OUTPUT(std::string directoryPath){
	//
	std::vector<Ddpca::CurvedSurface*> wowhSurf = {
		&wormTosu, &wormToba, &wormRtsu, &wormRtba, 
		&wheeTosu, &wheeToba, &wheeRtsu, &wheeRtba
	};
	std::vector<std::string> wowhFina = {
		"/resuWOTS.txt", "/resuWOTB.txt", "/resuWORT.txt", "/resuWORB.txt", 
		"/resuWHTS.txt", "/resuWHTB.txt", "/resuWHRT.txt", "/resuWHRB.txt"};
	//
	const Ddpca::I64 wosuSize = wowhSurf.size();
	for(Ddpca::I64 tw = 0; tw < wosuSize; tw ++){
		std::ofstream tempOfst;
		tempOfst.open(directoryPath + wowhFina[tw], std::ios::out);
		tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
		const Ddpca::I64 woinSize = (* wowhSurf[tw]).indexPoint.size();
		for(Ddpca::I64 ti = 0; ti < woinSize; ti ++){
			const Ddpca::I64 woiiSize = (* wowhSurf[tw]).indexPoint[ti].size();
			for(Ddpca::I64 tj = 0; tj < woiiSize; tj ++){
				tempOfst << std::setw(30) << (* wowhSurf[tw]).indexPoint[ti][tj][0] 
					<< std::setw(30) << (* wowhSurf[tw]).indexPoint[ti][tj][1] 
					<< std::setw(30) << (* wowhSurf[tw]).indexPoint[ti][tj][2] << std::endl;
			}
		}
		tempOfst.close();
	}
}

void DehwSurf::ESTABLISH(std::string directoryPath){
	BASIC_PARAMETER();
	TOOT_SURF_GRID();
	ROOT_TRAN_GRID();
	#ifdef NDEBUG
		OUTPUT(directoryPath);
	#endif
}

#endif // _DehwSurf_hpp