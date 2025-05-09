
#ifndef _DEHWSURF_H
#define _DEHWSURF_H

#include "../CURVEDS.h"

class DEHWSURF{
public:
	//****************************************Input-1****************************************
	DEHWSURF();
	//****************************************Input-2****************************************
	long z[2];//teeth number
	double a_h2;//working center distance
	double modiTran;//modification of transmission ratio
	double modiCent;//modification of center distance
	double r_b2;//radius of base circle of worm wheel
	double beta_c;//inclination angle of generating plane
	//****************************************Input-3****************************************
	double z_k;//encircled teeth
	double d[2];//reference circle diameter at throat
	double h_a_s[2];//addendum coefficient
	double h_f_s[2];//dedendum coefficient
	double R_a[2];//tip arc
	double offsR_a;//offset of tip arc center of worm wheel
	double wheeWidt;//width of worm wheel
	double inneRadi[2];//inner rdius of hub
	double inpuTorq;//input torque
	//****************************************Input-4****************************************
	long globInho;//global inhomogeneous level
	long globHomo;//global homogeneous level
	long locaLeve;//local refinement level
	long reliSwit;//0 - no tooth flank relif, 1 - tooth flank relief
	VECTOR2L gridNumb;//number of grid division
	//***************************************************************************************
	double a_1c;//center distance of primary enveloping
	double i_1c;//transmission ratio of primary enveloping
	double i_c1;
	double i_h2;//transmission ratio of second enveloping
	double i_2h;
	double m_t;//transverse module
	double h_a[2];//addendum
	double h_f[2];//dedendum
	double d_f[2];//root circle
	double d_a[2];//tip circle
	double R_f[2];//root arc
	double R_t[2];//transition arc
	double alph;//nominal pressure angle
	double leadAngl;//nominal lead angle
	double pitcAngl;//pitch angle
	double tootThicCoef[2];//tooth thickness coefficient
	double halfAngl;//half working angle
	double starAngl;//starting angle
	double termAngl;//terminating angle
	double wormCurv[3];//curvilinear coordinate of worm
	double widtAngl;//face width angle of worm wheel
	double backlash;//backlash
	//axial tooth thickness of worm, transverse tooth thickness of worm wheel
	double tootThic[2];
	double tootThicAngl[2];//tooth thickness angle
	//simultaneous envelope of tooth surface and tooth back
	//the angle between the two rigidly connected coordinates
	double backAngl[2];
	long BASIC_PARAMETER();//calculate basic parameter
	//***************************************************************************************
	//singular thet_c to thet_h
	long SINGULAR_C2H(double thet_c, double &thet_hs, double &thet_hm);
	//first and second meshing equations
	long FSME(double thet_1, double thet_h, double &x_d, double &y_d);
	// //partial derivative of meshing equations
	long PD_FSME(double thet_1, double thet_h, 
		double &x_d, double &y_d, Eigen::Matrix<double,2,2> &Pxy_d);
	//x_d, y_d, thet_c to r_1
	long WORM_DC2R(double x_d, double y_d, double thet_c, Eigen::Vector3d &r_1_1);
	//x_d, y_d, thet_1, thet_h to r_2
	long WHEE_1H2R(double x_d, double y_d, double thet_1, double thet_h, 
		Eigen::Vector3d &r_2_2);
	//partial derivative of r_2_2 relative to thet_1, thet_h
	long PD_WHEE_1H2R(double x_d, double y_d, double thet_1, double thet_h, 
		Eigen::Vector3d &r_2_2, Eigen::Matrix<double,3,2> &Dr_2_2);
	//curvature interference limit function of first envelope
	long CILFOFE(double thet_1, double x_d, double y_d,
		double &Psi_1, double &kapp_cxd, double &kapp_cyd, double &tau_cxd);
	//curvature interference limit function of second envelope, non-interval version
	double CILFOSE_NI(double thet_1, double thet_h, double &kapp_h2N);
	//***************************************************************************************
	//worm tooth surface: curvilinear coordinate to Cartesian coordinate
	long WORM_CURV_2_CART(double xi_11, double xi_12, Eigen::Vector3d &r_1_1, double &thet_c);
	//global coordinate r_2_2 to local coordinate of worm wheel
	long WHEE_G2L(Eigen::Vector3d r_2_2, double &angl_f, double &radi_f, 
		double &R_fmini, double &R_fmaxi);
	//worm wheel tooth surface: curvilinear coordinate to Cartesian coordinate
	//thet_c, thet_h are initial values, new contact zone
	long WHEE_CURV_2_CART_1(double xi_21, double xi_22, Eigen::Vector3d &r_2_2, 
		double &thet_c, double &thet_h, long f_lr, double &x_d, double &y_d);
	//x_d, y_d are initial values, former contact zone
	long WHEE_CURV_2_CART_2(double xi_21, double xi_22, Eigen::Vector3d &r_c_c, 
		double &thet_c, double &x_d, double &y_d);
	//transition zone
	long WHEE_CURV_2_CART_3(double xi_21, double xi_22, Eigen::Vector3d &r_2_2, 
		double &thet_c, double &thet_h, double xi_11);
	//transition zone of worm wheel, xi_11 - head transition zone, rear transition zone
	long WHEE_TRAN(double thet_c, double thet_h, double xi_11, 
		Eigen::Vector3d &r_2_2, Eigen::Matrix<double,3,2> &Dr_2_2);
	//phase analysis of worm wheel tooth surface
	long WHEE_PHAS(long ti, long tj, long f_ij, Eigen::Vector3d r_2_2);
	//new contact zone, 1 - left, 2 - right
	long NEW_CONT_ZONE(long f_lr);
	long FORMER_CONT_ZONE();//former contact zone
	long TRANSITION_ZONE(long f_hr);//transition zone, 1 - head, 2 - rear
	//tooth flank relief of worm
	long WORM_RELI(Eigen::Vector3d &tempXYZ, long ti, long tj);
	//tooth flank relief of worm wheel
	long WHEE_RELI(Eigen::Vector3d &tempXYZ, long ti, long tj);
	//***************************************************************************************
	//curvilinear coordinate of tooth profile
	std::vector<std::vector<std::vector<Eigen::Vector2d>>> curvCoor;
	//Cartesian coordinate of tooth profile
	std::vector<std::vector<std::vector<Eigen::Vector3d>>> cartCoor;
	//flag of worm wheel tooth surface
	//1 - left new contact zone, 2 - right new contact zone, 3 - former contact zone
	//4 - head transition zone,  5 - rear transition zone,   0 - fail to solve
	VECTOR2L fpha;
	CURVEDS wormTosu;//worm tooth surface
	CURVEDS wormToba;//worm tooth back surface
	CURVEDS wormRtsu;//worm root transition surface
	CURVEDS wormRtba;//worm root transition back
	CURVEDS wheeTosu;//worm wheel tooth surface
	CURVEDS wheeToba;//worm wheel tooth back surface
	CURVEDS wheeRtsu;//worm wheel root transition surface
	CURVEDS wheeRtba;//worm wheel root transition back
	long circNumb;//number of worm domains in 2 * PI
	long WORM_TS_GRID();//grid discretization of worm tooth surface
	long WHEE_TS_GRID();//grid discretization of worm wheel tooth surface
	long TOOT_SURF_GRID();//grid discretization of tooth surface
	long ROOT_TRAN_GRID();//grid discretization of root transition surface
	//***************************************************************************************
	//calculate the radius of root transition arc
	long WORM_ROOT_RADIUS(long flag, Eigen::Matrix<double,2,3> tempPoin, 
		Eigen::Vector2d &tempCent, double &tempRadi, Eigen::Vector2d &tempAngl
	);
	//calculate nodes on root transition arc
	long WORM_ROOT(long indi, long flag, Eigen::MatrixXd &rootProf);
	Eigen::Vector2d WHEE_UNCONE(Eigen::Vector3d tempXYZ, double tempAlph_3);
	//transfer from "in unfolded cone surface" to r_2_2
	Eigen::Vector3d WHEE_CONE(Eigen::Vector2d tempXY, double tempAlph_3);
	//calculate nodes on root transition arc "in unfolded cone surface"
	long WHEE_ROOT(long indi, long flag, Eigen::MatrixXd &rootProf);
	long OUTPUT();//output wormTosu~wheeRtba to files
	long ESTABLISH();
};

class INDE_INIT{
public:
	long ti;
	long tj;
	double init_1;
	double init_2;
	INDE_INIT(long inde_1, long inde_2, double valu_1, double valu_2)
		: ti(inde_1), tj(inde_2), init_1(valu_1), init_2(valu_2){}
};

DEHWSURF::DEHWSURF(){
	//
	z[0] = 1;//must be 1
	z[1] = 40;
	a_h2 = 0.25;//m
	modiTran = 0.0;
	modiCent = 0.0;
	r_b2 = 0.158/2.0;//m
	beta_c = 11.0 * PI / 180.0;//rad
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
	gridNumb = {{4, 2, 2, 4, 4, 0, 0}, {4, 4, 2, 4, 8, 8 + z[0], 2}};
	globInho = 1;//direction of xi_11/facewidth
	globHomo = 2;//must >= 1
	locaLeve = 3;
	reliSwit = 1;
	circNumb = 8;
}

long DEHWSURF::BASIC_PARAMETER(){
	//
	a_1c = a_h2 + modiCent;
	i_h2 = (double)(z[1]) / z[0];
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
	alph = asin(2.0 * r_b2 / d[1]);
	leadAngl = atan(d[1] / i_h2 / d[0]);
	pitcAngl = 2.0 * PI / z[1];
	tootThicCoef[0] = 0.45;
	tootThicCoef[1] = 0.55;
	halfAngl = 0.5 * (z_k - tootThicCoef[0]) * pitcAngl;
	starAngl = alph - halfAngl;
	termAngl = starAngl + z_k * pitcAngl;
	wormCurv[0] = i_h2 * starAngl;
	wormCurv[2] = i_h2 * termAngl;
	wormCurv[1] = (wormCurv[0] + wormCurv[2]) / 2.0;
	while(wormCurv[1] - 2.0 * PI >= wormCurv[0]){
		wormCurv[1] = wormCurv[1] - 2.0 * PI;
	}
	//
	widtAngl = asin(wheeWidt / 2.0 / R_f[1]);
	backlash = 0.0;
	tootThic[0] = tootThicCoef[0] * PI * m_t - backlash;
	tootThic[1] = tootThicCoef[1] * PI * m_t;
	tootThicAngl[0] = tootThic[0] / (d[1] / 2.0);
	tootThicAngl[1] = tootThic[1] / (d[1] / 2.0);
	//figure 3.11 in [Zhou, L. Modification principle and manufacturing technology for hourglass 
	//worm drives (National University of Defense Technology Press, 2005)].
	backAngl[0] = 2.0 * alph + tootThicAngl[0];
	backAngl[1] = 2.0 * alph - tootThicAngl[1];
	return 1;
}

long DEHWSURF::SINGULAR_C2H(double thet_c, double &thet_hs, double &thet_hm){
	//0 < thet_c < PI / 2.0
	double thet_1 = i_1c * thet_c;
	double C_m11 = - i_2h * cos(beta_c) * sin(thet_c);
	double C_m12 = i_c1 * i_2h * cos(beta_c) * cos(thet_c)
		+ i_2h * sin(beta_c);
	double C_m13 = i_c1 * cos(beta_c) * sin(thet_c);
	double a2CC = atan2(C_m11, C_m12);
	if(C_m13 > sqrt(C_m11 * C_m11 + C_m12 * C_m12)){
		thet_hs = thet_1 - a2CC - PI / 2.0;
		thet_hm = thet_hs;
	}
	else{
		thet_hs = thet_1 - PI - a2CC + asin(C_m13 / sqrt(C_m11 * C_m11 + C_m12 * C_m12));
		thet_hm = thet_1 - a2CC - asin(C_m13 / sqrt(C_m11 * C_m11 + C_m12 * C_m12));
	}
	return 1;
}

long DEHWSURF::FSME(double thet_1, double thet_h, double &x_d, double &y_d){
	//
	double thet_c = i_c1 * thet_1;
	//
	Eigen::Matrix<double,2,2> coefA;
	coefA << - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c), sin(thet_c),
		sin(beta_c) * cos(thet_c)
		+ i_2h * cos(beta_c) * cos(thet_h - thet_1) 
		- i_2h * sin(beta_c) * sin(thet_c) * sin(thet_h - thet_1),
		-sin(thet_c) - i_2h * cos(thet_c) * sin(thet_h -thet_1);
	Eigen::Vector2d coefB;
	coefB << - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c),
		+ r_b2 * sin(beta_c) * sin(thet_c) 
		+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * sin(thet_h - thet_1)
		- a_1c * sin(beta_c)
		- i_2h * a_1c * cos(beta_c) * cos(thet_c) * cos(thet_h - thet_1)
		+ i_2h * a_h2 * cos(beta_c) * cos(thet_c);
	Eigen::Vector2d xy_d = coefA.fullPivLu().solve(-coefB);
	x_d = xy_d(0);
	y_d = xy_d(1);
	return 1;
}

long DEHWSURF::PD_FSME(double thet_1, double thet_h, double &x_d, double &y_d, 
	Eigen::Matrix<double,2,2> &Pxy_d){
	//
	double thet_c = i_c1 * thet_1;
	//
	Eigen::Matrix<double,2,2> coefA;
	coefA << - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c), sin(thet_c),
		sin(beta_c) * cos(thet_c)
		+ i_2h * cos(beta_c) * cos(thet_h - thet_1) 
		- i_2h * sin(beta_c) * sin(thet_c) * sin(thet_h - thet_1),
		-sin(thet_c) - i_2h * cos(thet_c) * sin(thet_h -thet_1);
	Eigen::Vector2d coefB;
	coefB << - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c),
		+ r_b2 * sin(beta_c) * sin(thet_c) 
		+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * sin(thet_h - thet_1)
		- a_1c * sin(beta_c)
		- i_2h * a_1c * cos(beta_c) * cos(thet_c) * cos(thet_h - thet_1)
		+ i_2h * a_h2 * cos(beta_c) * cos(thet_c);
	Eigen::Vector2d xy_d = coefA.fullPivLu().solve(-coefB);
	x_d = xy_d(0);
	y_d = xy_d(1);
	//derivative relative to thet_1
	Eigen::Matrix<double,2,2> PcoefA;
	PcoefA << sin(beta_c) * sin(thet_c) * i_c1, cos(thet_c) * i_c1,
		- sin(beta_c) * sin(thet_c) * i_c1
		- i_2h * cos(beta_c) * sin(thet_h - thet_1) * -1.0
		- i_2h * sin(beta_c) * cos(thet_c) * i_c1 * sin(thet_h - thet_1)
		- i_2h * sin(beta_c) * sin(thet_c) * cos(thet_h - thet_1) * -1.0,
		- cos(thet_c) * i_c1 
		+ i_2h * sin(thet_c) * i_c1 * sin(thet_h -thet_1) 
		- i_2h * cos(thet_c) * cos(thet_h -thet_1) * -1.0;
	Eigen::Vector2d PcoefB;
	PcoefB << - r_b2 * sin(beta_c) * cos(thet_c) * i_c1,
		r_b2 * sin(beta_c) * cos(thet_c) * i_c1
		- i_2h * r_b2 * sin(beta_c) * sin(thet_c) * i_c1 * sin(thet_h - thet_1)
		+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * cos(thet_h - thet_1) * -1.0
		+ i_2h * a_1c * cos(beta_c) * sin(thet_c) * i_c1 * cos(thet_h - thet_1)
		+ i_2h * a_1c * cos(beta_c) * cos(thet_c) * sin(thet_h - thet_1) * -1.0
		- i_2h * a_h2 * cos(beta_c) * sin(thet_c) * i_c1;
	Pxy_d.block(0,0,2,1) = coefA.fullPivLu().solve(- PcoefB - PcoefA * xy_d);
	//derivative relative to thet_h
	PcoefA << 0.0, 0.0,
		0.0 
		- i_2h * cos(beta_c) * sin(thet_h - thet_1) 
		- i_2h * sin(beta_c) * sin(thet_c) * cos(thet_h - thet_1),
		- 0.0 - i_2h * cos(thet_c) * cos(thet_h -thet_1);
	PcoefB << 0.0 + 0.0,
		+ 0.0 
		+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * cos(thet_h - thet_1)
		- 0.0
		+ i_2h * a_1c * cos(beta_c) * cos(thet_c) * sin(thet_h - thet_1)
		+ 0.0;
	Pxy_d.block(0,1,2,1) = coefA.fullPivLu().solve(- PcoefB - PcoefA * xy_d);
	return 1;
}

long DEHWSURF::WORM_DC2R(double x_d, double y_d, double thet_c, Eigen::Vector3d &r_1_1){
	double thet_1 = i_1c * thet_c;
	Eigen::Vector3d xyz;
	xyz << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
	Eigen::Matrix<double,3,3> rota;
	//R_oc,c
	rota << cos(thet_c), -sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	xyz = rota * xyz;
	//R_o1,oc
	rota << 1.0, 0.0, 0.0,
		0.0, 0.0, -1.0,
		0.0, 1.0, 0.0;
	xyz = rota * xyz;
	//T_o1,oc
	xyz(0) = xyz(0) + a_1c;
	//R_1,o1
	rota << cos(thet_1), sin(thet_1), 0.0,
		-sin(thet_1), cos(thet_1), 0.0,
		0.0, 0.0, 1.0;
	r_1_1 = rota * xyz;
	return 1;
}

long DEHWSURF::WHEE_1H2R(double x_d, double y_d, double thet_1, double thet_h, 
	Eigen::Vector3d &r_2_2){
	//
	double thet_c = i_c1 * thet_1;
	double thet_2 = i_2h * thet_h;
	//
	Eigen::Vector3d xyz;
	Eigen::Matrix<double,3,3> rota;
	WORM_DC2R(x_d, y_d, thet_c, xyz);
	//R_oh,h
	rota << cos(thet_h), -sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	xyz = rota * xyz;
	//R_o2,oh
	rota << 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, -1.0, 0.0;
	xyz = rota * xyz;
	//T_o2,oh
	xyz(0) = xyz(0) - a_h2;
	//R_2,o2
	rota << cos(thet_2), sin(thet_2), 0.0,
		-sin(thet_2), cos(thet_2), 0.0,
		0.0, 0.0, 1.0;
	r_2_2 = rota * xyz;
	return 1;
}

long DEHWSURF::PD_WHEE_1H2R(double x_d, double y_d, double thet_1, double thet_h, 
	Eigen::Vector3d &r_2_2, Eigen::Matrix<double,3,2> &Dr_2_2){
	//
	double thet_c = i_c1 * thet_1;
	double thet_2 = i_2h * thet_h;
	//relative to thet_1
	Eigen::Matrix<double,2,2> Dxy_d;
	PD_FSME(thet_1, thet_h, x_d, y_d, Dxy_d);
	double Dx_d = Dxy_d(0,0);
	double Dy_d = Dxy_d(1,0);
	Eigen::Vector3d r_c_c, Dr_c_c;
	r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
	Dr_c_c << - Dx_d, - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
	Eigen::Matrix<double,3,3> R_oc_c, DR_oc_c;
	R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	DR_oc_c << - i_c1 * sin(thet_c), - i_c1 * cos(thet_c), 0.0,
		i_c1 * cos(thet_c), - i_c1 * sin(thet_c), 0.0,
		0.0, 0.0, 0.0;
	Eigen::Vector3d r_c_oc = R_oc_c * r_c_c;
	Eigen::Vector3d Dr_c_oc = DR_oc_c * r_c_c + R_oc_c * Dr_c_c;
	Eigen::Matrix<double,3,3> R_o1_oc;
	Eigen::Vector3d T_o1_oc;
	R_o1_oc << 1.0, 0.0, 0.0,
		0.0, 0.0, -1.0,
		0.0, 1.0, 0.0;
	T_o1_oc << a_1c, 0.0, 0.0;
	Eigen::Vector3d r_c_o1 = R_o1_oc * r_c_oc + T_o1_oc;
	Eigen::Vector3d Dr_c_o1 = R_o1_oc * Dr_c_oc;
	Eigen::Matrix<double,3,3> R_1_o1, DR_1_o1;
	R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
		-sin(thet_1), cos(thet_1), 0.0,
		0.0, 0.0, 1.0;
	DR_1_o1 << -sin(thet_1), cos(thet_1), 0.0,
		-cos(thet_1), -sin(thet_1), 0.0,
		0.0, 0.0, 0.0;
	Eigen::Vector3d r_1_1 = R_1_o1 * r_c_o1;
	Eigen::Vector3d Dr_1_1 = DR_1_o1 * r_c_o1 + R_1_o1 * Dr_c_o1;
	Eigen::Matrix<double,3,3> R_oh_h;
	R_oh_h << cos(thet_h), - sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	Eigen::Vector3d r_h_oh = R_oh_h * r_1_1;
	Eigen::Vector3d Dr_h_oh = R_oh_h * Dr_1_1;
	Eigen::Matrix<double,3,3> R_o2_oh;
	Eigen::Vector3d T_o2_oh;
	R_o2_oh << 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, -1.0, 0.0;
	T_o2_oh << - a_h2, 0.0, 0.0;
	Eigen::Vector3d r_h_o2 = R_o2_oh * r_h_oh + T_o2_oh;
	Eigen::Vector3d Dr_h_o2 = R_o2_oh * Dr_h_oh;
	Eigen::Matrix<double,3,3> R_2_o2;
	R_2_o2 << cos(thet_2), sin(thet_2), 0.0,
		-sin(thet_2), cos(thet_2), 0.0,
		0.0, 0.0, 1.0;
	r_2_2 = R_2_o2 * r_h_o2;
	Dr_2_2.block(0,0,3,1) = R_2_o2 * Dr_h_o2;
	//relative to thet_h
	Dx_d = Dxy_d(0,1);
	Dy_d = Dxy_d(1,1);
	Dr_c_c << - Dx_d, - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
	Dr_c_oc = R_oc_c * Dr_c_c;
	Dr_c_o1 = R_o1_oc * Dr_c_oc;
	Dr_1_1 = R_1_o1 * Dr_c_o1;
	Eigen::Matrix<double,3,3> DR_oh_h;
	DR_oh_h << -sin(thet_h), -cos(thet_h), 0.0,
		cos(thet_h), -sin(thet_h), 0.0,
		0.0, 0.0, 0.0;
	Dr_h_oh = DR_oh_h * r_1_1 + R_oh_h * Dr_1_1;
	Dr_h_o2 = R_o2_oh * Dr_h_oh;
	Eigen::Matrix<double,3,3> DR_2_o2;
	DR_2_o2 << -i_2h * sin(thet_2), i_2h * cos(thet_2), 0.0,
		-i_2h * cos(thet_2), -i_2h * sin(thet_2), 0.0,
		0.0, 0.0, 0.0;
	Dr_2_2.block(0,1,3,1) = DR_2_o2 * r_h_o2 + R_2_o2 * Dr_h_o2;
	return 1;
}

long DEHWSURF::CILFOFE(double thet_1, double x_d, double y_d, 
	double &Psi_1, double &kapp_1xd, double &kapp_1yd, double &tau_1xd){
	//
	double thet_c = thet_1 / i_1c;
	//
	double kapp_cxd = 0.0;
	double kapp_cyd = 0.0;
	double tau_cxd = 0.0;
	Eigen::Vector3d i_d_c, j_d_c, i_d_oc, j_d_oc, omeg_c1_oc, v_c1_oc;
	i_d_c << -1.0, 0.0, 0.0;
	j_d_c << 0.0, - sin(beta_c), cos(beta_c);
	//R_oc,c
	Eigen::Matrix<double,3,3> rota;
	rota << cos(thet_c), -sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	i_d_oc = rota * i_d_c;
	j_d_oc = rota * j_d_c;
	omeg_c1_oc << 0.0, -1.0, i_c1;
	v_c1_oc << - y_d * cos(beta_c) 
		- i_c1 * (- x_d * sin(thet_c) + cos(thet_c) * (r_b2 - y_d * sin(beta_c))),
		i_c1 * (- x_d * cos(thet_c) - sin(thet_c) * (r_b2 - y_d * sin(beta_c))),
		- x_d * cos(thet_c) - sin(thet_c) * (r_b2 - y_d * sin(beta_c)) + a_1c;
	double N_1xd = kapp_cxd * v_c1_oc.dot(i_d_oc) + tau_cxd * v_c1_oc.dot(j_d_oc)
		+ omeg_c1_oc.dot(j_d_oc);
	double N_1yd = tau_cxd * v_c1_oc.dot(i_d_oc) + kapp_cyd * v_c1_oc.dot(j_d_oc)
		- omeg_c1_oc.dot(i_d_oc);
	Eigen::Vector3d N_1_oc = N_1xd * i_d_oc + N_1yd * j_d_oc;
	double PPhi_1Pthet_1 = x_d * sin(beta_c) * sin(thet_c) / i_1c
		+ y_d * cos(thet_c) / i_1c
		- r_b2 * sin(beta_c) * cos(thet_c) / i_1c;
	Psi_1 = N_1_oc.dot(v_c1_oc) + PPhi_1Pthet_1;
	double kapp_c1xd = N_1xd * N_1xd / Psi_1;
	double kapp_c1yd = N_1yd * N_1yd / Psi_1;	
	double tau_c1xd = N_1xd * N_1yd / Psi_1;
	kapp_1xd = kapp_cxd - kapp_c1xd;
	kapp_1yd = kapp_cyd - kapp_c1yd;
	tau_1xd = tau_cxd - tau_c1xd;
	return 1;
}

double DEHWSURF::CILFOSE_NI(double thet_1, double thet_h, double &kapp_h2N){
	//
	double thet_c = thet_1 / i_1c;
	double thet_2 = thet_h / i_h2;
	//
	double x_d, y_d;
	FSME(thet_1, thet_h, x_d, y_d);
	double Psi_1, kapp_1xd, kapp_1yd, tau_1xd;
	CILFOFE(thet_1, x_d, y_d, Psi_1, kapp_1xd, kapp_1yd, tau_1xd);
	double kapp_hxd = kapp_1xd;
	double kapp_hyd = kapp_1yd;
	double tau_hxd = tau_1xd;
	//
	Eigen::Vector3d i_d_c, j_d_c, i_d_oh, j_d_oh;
	i_d_c << -1.0, 0.0, 0.0;
	j_d_c << 0.0, - sin(beta_c), cos(beta_c);
	Eigen::Vector3d r_c_c, r_h_oh;
	r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
	//R_oc,c
	Eigen::Matrix<double,3,3> rota;
	rota << cos(thet_c), -sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	i_d_oh = rota * i_d_c;
	j_d_oh = rota * j_d_c;
	r_h_oh = rota * r_c_c;
	//R_o1,oc
	rota << 1.0, 0.0, 0.0,
		0.0, 0.0, -1.0,
		0.0, 1.0, 0.0;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//T_o1,oc
	r_h_oh(0) = r_h_oh(0) + a_1c;
	//R_1,o1
	rota << cos(thet_1), sin(thet_1), 0.0,
		-sin(thet_1), cos(thet_1), 0.0,
		0.0, 0.0, 1.0;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//R_oh,h
	rota << cos(thet_h), -sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//
	Eigen::Vector3d omeg_h2_oh, omeg_2_oh, o_h2_oh, v_h2_oh;
	omeg_h2_oh << 0.0, i_2h, 1.0;
	omeg_2_oh << 0.0, - i_2h, 0.0;
	o_h2_oh << - a_h2, 0.0, 0.0;
	v_h2_oh = omeg_h2_oh.cross(r_h_oh) - omeg_2_oh.cross(o_h2_oh);
	double N_2xd = kapp_hxd * v_h2_oh.dot(i_d_oh) + tau_hxd * v_h2_oh.dot(j_d_oh)
		+omeg_h2_oh.dot(j_d_oh);
	double N_2yd = tau_hxd * v_h2_oh.dot(i_d_oh) + kapp_hyd * v_h2_oh.dot(j_d_oh)
		-omeg_h2_oh.dot(i_d_oh);
	Eigen::Vector3d N_2_oh = N_2xd * i_d_oh + N_2yd * j_d_oh;
	double B_11 = i_2h * x_d * cos(beta_c) - i_2h * a_1c * cos(beta_c) * cos(thet_c);
	double B_12 = - i_2h * x_d * sin(beta_c) * sin(thet_c) 
		- i_2h * y_d * cos(thet_c) + i_2h * r_b2 * sin(beta_c) * cos(thet_c);
	double PPhi_2Pthet_h = - B_11 * sin(thet_h - thet_1) + B_12 * cos(thet_h - thet_1);
	double Psi_2 = N_2_oh.dot(v_h2_oh) + PPhi_2Pthet_h;
	kapp_h2N = (N_2xd * N_2xd + N_2yd * N_2yd) / Psi_2;
	return Psi_2;
}

long DEHWSURF::WORM_CURV_2_CART(double xi_11, double xi_12, 
	Eigen::Vector3d &r_1_1, double &thet_c){
	Eigen::Vector2d x;//0 - thet_c, 1 - x_d
	x << i_c1 * xi_11, d[1] / 2.0;//
	while(true){
		Eigen::Vector2d func, deltX;
		Eigen::Matrix<double,2,2> Dfunc;
		//function
		thet_c = x(0);
		double thet_1 = i_1c * thet_c;
		double x_d = x(1);
		double y_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
			+ (- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c))) / sin(thet_c);
		Eigen::Vector3d r_c_c, T_o1_oc, r_1_o1;
		Eigen::Matrix<double,3,3> R_oc_c, R_o1_oc, R_1_o1;
		r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
		R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
			sin(thet_c), cos(thet_c), 0.0,
			0.0, 0.0, 1.0;
		R_o1_oc << 1.0, 0.0, 0.0,
			0.0, 0.0, -1.0,
			0.0, 1.0, 0.0;
		R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
			- sin(thet_1), cos(thet_1), 0.0,
			0.0, 0.0, 1.0;
		T_o1_oc << a_1c, 0.0, 0.0;
		r_1_o1 = R_o1_oc * R_oc_c * r_c_c + T_o1_oc;
		r_1_1 = R_1_o1 * r_1_o1;
		func << thet_1 - atan2(r_1_o1(1), r_1_o1(0)) - xi_11,
			r_1_1(2) * r_1_1(2) 
			+ (a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1))) 
			* (a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)))
			- xi_12 * xi_12;
		//derivative 1
		double Dy_d = - (((+ sin(beta_c) * sin(thet_c) - 0.0) * x_d 
			+ (- r_b2 * sin(beta_c) * cos(thet_c) + 0.0)) * sin(thet_c) 
			- ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
			+ (- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c))) * cos(thet_c))
			/ (sin(thet_c) * sin(thet_c));
		Eigen::Vector3d Dr_c_c, Dr_1_o1, Dr_1_1;
		Eigen::Matrix<double,3,3> DR_oc_c, DR_1_o1;
		Dr_c_c << 0.0, 0.0 - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
		DR_oc_c << - sin(thet_c), - cos(thet_c), 0.0,
			cos(thet_c), - sin(thet_c), 0.0,
			0.0, 0.0, 0.0;
		DR_1_o1 << - i_1c * sin(thet_1), i_1c * cos(thet_1), 0.0,
			- i_1c * cos(thet_1), - i_1c * sin(thet_1), 0.0,
			0.0, 0.0, 0.0;
		Dr_1_o1 = R_o1_oc * DR_oc_c * r_c_c + R_o1_oc * R_oc_c * Dr_c_c;
		Dr_1_1 = DR_1_o1 * r_1_o1 + R_1_o1 * Dr_1_o1;
		Dfunc(0,0) = i_1c - (Dr_1_o1(1) * r_1_o1(0) - r_1_o1(1) * Dr_1_o1(0))
			/ (r_1_o1(1) * r_1_o1(1) + r_1_o1(0) * r_1_o1(0));
		atan2(r_1_o1(1), r_1_o1(0));
		Dfunc(1,0) = 2.0 * r_1_1(2) * Dr_1_1(2) 
			+ 2.0 * (a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1))) 
			* (0.0 - (r_1_1(0) * Dr_1_1(0) + r_1_1(1) * Dr_1_1(1)) 
			/ sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)));
		//derivative 2
		Dy_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) 
			+ 0.0) / sin(thet_c);
		Dr_c_c << - 1.0, 0.0 - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
		Dr_1_o1 = R_o1_oc * R_oc_c * Dr_c_c;
		Dr_1_1 = R_1_o1 * Dr_1_o1;
		Dfunc(0,1) = - (Dr_1_o1(1) * r_1_o1(0) - r_1_o1(1) * Dr_1_o1(0))
			/ (r_1_o1(1) * r_1_o1(1) + r_1_o1(0) * r_1_o1(0));
		Dfunc(1,1) = 2.0 * r_1_1(2) * Dr_1_1(2) 
			+ 2.0 * (a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1))) 
			* (0.0 - (r_1_1(0) * Dr_1_1(0) + r_1_1(1) * Dr_1_1(1)) 
			/ sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)));
		//iteration
		deltX = Dfunc.fullPivLu().solve(-func);
		if(deltX.norm() < 1.0E-12){
			if(a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)) < 0.0){
				std::cout << "ERROR in DEHWSURF::WORM_CURV_2_CART!" << std::endl;
			}
			break;
		}
		x += deltX;
	}
	return 1;
}

long DEHWSURF::WHEE_G2L(Eigen::Vector3d r_2_2, double &angl_f, double &radi_f, 
	double &R_fmini, double &R_fmaxi){
	double radi_xi = a_h2 - sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
	double coor_zi = r_2_2(2);
	angl_f = atan2(coor_zi, radi_xi);
	radi_f = sqrt(radi_xi * radi_xi + coor_zi * coor_zi);
	double angl_ai = angl_f - asin(offsR_a * sin(angl_f) / R_a[1]);
	R_fmini = (R_a[1] * cos(angl_ai) - offsR_a) / cos(angl_f);
	R_fmaxi = R_t[1];
	return 1;
}

long DEHWSURF::WHEE_CURV_2_CART_1(double xi_21, double xi_22, Eigen::Vector3d &r_2_2, 
	double &thet_c, double &thet_h, long f_lr, double &x_d, double &y_d){
	Eigen::Vector2d x;
	x << thet_c, thet_h;
	for(long ti = 0; ti < 1000; ti ++){
		Eigen::Vector2d func, deltX;
		Eigen::Matrix<double,2,2> Dfunc;
		thet_c = x(0);
		thet_h = x(1);
		//function
		double thet_1 = i_1c * thet_c;
		FSME(thet_1, thet_h, x_d, y_d);
		Eigen::Matrix<double,3,2> Dr_2_2;
		PD_WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2, Dr_2_2);
		double coorX = a_h2 - sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		func << atan2(r_2_2(2), coorX) - xi_21,
			r_2_2(2) * r_2_2(2) + coorX * coorX - xi_22 * xi_22;
		//derivative 1
		double DcoorX = - (r_2_2(0) * Dr_2_2(0,0) + r_2_2(1) * Dr_2_2(1,0)) 
			/ sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		Dfunc(0,0) = (Dr_2_2(2,0) * coorX - r_2_2(2) * DcoorX) 
			/ (coorX * coorX + r_2_2(2) * r_2_2(2));
		Dfunc(1,0) = 2.0 * r_2_2(2) * Dr_2_2(2,0) + 2.0 * coorX * DcoorX;
		Dfunc.block(0,0,2,1) = i_1c * Dfunc.block(0,0,2,1);
		//derivative 2
		DcoorX = - (r_2_2(0) * Dr_2_2(0,1) + r_2_2(1) * Dr_2_2(1,1)) 
			/ sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		Dfunc(0,1) = (Dr_2_2(2,1) * coorX - r_2_2(2) * DcoorX) 
			/ (coorX * coorX + r_2_2(2) * r_2_2(2));
		Dfunc(1,1) = 2.0 * r_2_2(2) * Dr_2_2(2,1) + 2.0 * coorX * DcoorX;
		//iteration
		deltX = Dfunc.fullPivLu().solve(-func);
		if(deltX.norm() < 1.0E-12){
			if(coorX < 0.0){
				std::cout << "ERROR in DEHWSURF::WHEE_CURV_2_CART_1!" << std::endl;
			}
			break;
		}
		double rfac = 2.0;
		long rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			Eigen::Vector2d x_test = x + rfac * deltX;
			if(x_test(0) < 0.01 * PI || x_test(0) > 0.49 * PI){
				continue;
			}
			double thet_hs, thet_hm;
			SINGULAR_C2H(x(0), thet_hs, thet_hm);
			if((f_lr == 1 && (x_test(1) <= thet_hs + 1.0E-14 || thet_hm - 1.0E-14 <= x_test(1)))
				|| (f_lr == 2 && (x_test(1) <= thet_hm + 1.0E-14 
				|| thet_hs + 2.0 * PI - 1.0E-14 <= x_test(1)))){
				continue;
			}
			double thet_ct = x_test(0);
			double thet_ht = x_test(1);
			//function
			double thet_1t = i_1c * thet_ct;
			double x_dt, y_dt;
			FSME(thet_1t, thet_ht, x_dt, y_dt);
			Eigen::Vector3d r_2_2t;
			Eigen::Matrix<double,3,2> Dr_2_2t;
			PD_WHEE_1H2R(x_dt, y_dt, thet_1t, thet_ht, r_2_2t, Dr_2_2t);
			double cooX_t = a_h2 - sqrt(r_2_2t(0) * r_2_2t(0) + r_2_2t(1) * r_2_2t(1));
			Eigen::Vector2d func_t;
			func_t << atan2(r_2_2t(2), cooX_t) - xi_21,
				r_2_2t(2) * r_2_2t(2) + cooX_t * cooX_t - xi_22 * xi_22;
			if(func_t.norm() < func.norm()){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
	return 1;
}

long DEHWSURF::WHEE_CURV_2_CART_2(double xi_21, double xi_22, Eigen::Vector3d &r_c_c, 
		double &thet_c, double &x_d, double &y_d){
	Eigen::Vector2d x;
	x << thet_c, x_d;
	for(long ti = 0; ti < 1000; ti ++){
		Eigen::Vector2d func, deltX;
		Eigen::Matrix<double,2,2> Dfunc;
		thet_c = x(0);
		x_d = x(1);
		//function
		y_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
			- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c)) 
			/ sin(thet_c);
		r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
		double coorX = a_h2 - sqrt(r_c_c(0) * r_c_c(0) + r_c_c(1) * r_c_c(1));
		func << atan2(r_c_c(2), coorX) - xi_21,
			r_c_c(2) * r_c_c(2) + coorX * coorX - xi_22 * xi_22;
		//derivative 1
		double Dy_d = - (((+ sin(beta_c) * sin(thet_c) - 0.0) * x_d 
			- r_b2 * sin(beta_c) * cos(thet_c) + 0.0) * sin(thet_c) 
			- ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
			- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c)) * cos(thet_c))
			/ (sin(thet_c) * sin(thet_c));
		Eigen::Vector3d Dr_c_c;
		Dr_c_c << - 0.0, 0.0 - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
		double DcoorX = - (r_c_c(0) * Dr_c_c(0) + r_c_c(1) * Dr_c_c(1)) 
			/ sqrt(r_c_c(0) * r_c_c(0) + r_c_c(1) * r_c_c(1));
		Dfunc(0,0) = (Dr_c_c(2) * coorX - r_c_c(2) * DcoorX) 
			/ (coorX * coorX + r_c_c(2) * r_c_c(2));
		Dfunc(1,0) = 2.0 * r_c_c(2) * Dr_c_c(2) + 2.0 * coorX * DcoorX;
		//derivative 2
		Dy_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * 1.0 
			- 0.0 + 0.0) / sin(thet_c);
		Dr_c_c << - 1.0, - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
		DcoorX = - (r_c_c(0) * Dr_c_c(0) + r_c_c(1) * Dr_c_c(1)) 
			/ sqrt(r_c_c(0) * r_c_c(0) + r_c_c(1) * r_c_c(1));
		Dfunc(0,1) = (Dr_c_c(2) * coorX - r_c_c(2) * DcoorX) 
			/ (coorX * coorX + r_c_c(2) * r_c_c(2));
		Dfunc(1,1) = 2.0 * r_c_c(2) * Dr_c_c(2) + 2.0 * coorX * DcoorX;
		//iteration
		deltX = Dfunc.fullPivLu().solve(-func);
		if(deltX.norm() < 1.0E-12){
			if(coorX < 0.0){
				std::cout << "ERROR in DEHWSURF::WHEE_CURV_2_CART_2!" << std::endl;
			}
			break;
		}
		double rfac = 2.0;
		long rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			Eigen::Vector2d x_test = x + rfac * deltX;
			if(x_test(0) < 0.01 * PI || x_test(0) > 0.49 * PI){
				continue;
			}
			double thet_ct = x_test(0);
			double x_dt = x_test(1);
			//function
			double y_dt = - ((- sin(beta_c) * cos(thet_ct) - i_c1 * cos(beta_c)) * x_dt 
				- r_b2 * sin(beta_c) * sin(thet_ct) + a_1c * sin(beta_c)) 
				/ sin(thet_ct);
			Eigen::Vector3d r_c_ct;
			r_c_ct << - x_dt, r_b2 - y_dt * sin(beta_c), y_dt * cos(beta_c);
			double cooX_t = a_h2 - sqrt(r_c_ct(0) * r_c_ct(0) + r_c_ct(1) * r_c_ct(1));
			Eigen::Vector2d func_t;
			func_t << atan2(r_c_ct(2), cooX_t) - xi_21,
				r_c_ct(2) * r_c_ct(2) + cooX_t * cooX_t - xi_22 * xi_22;
			if(func_t.norm() < func.norm()){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
	return 1;
}

long DEHWSURF::WHEE_CURV_2_CART_3(double xi_21, double xi_22, Eigen::Vector3d &r_2_2, 
	double &thet_c, double &thet_h, double xi_11){
	Eigen::Vector2d x;
	x << thet_c, thet_h;
	for(long ti = 0; ti < 1000; ti ++){
		Eigen::Vector2d func, deltX;
		Eigen::Matrix<double,2,2> Dfunc;
		thet_c = x(0);
		thet_h = x(1);
		//function
		Eigen::Matrix<double,3,2> Dr_2_2;
		WHEE_TRAN(thet_c, thet_h, xi_11, r_2_2, Dr_2_2);
		double coorX = a_h2 - sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		func << atan2(r_2_2(2), coorX) - xi_21,
			r_2_2(2) * r_2_2(2) + coorX * coorX - xi_22 * xi_22;
		//derivative 1
		double DcoorX = - (r_2_2(0) * Dr_2_2(0,0) + r_2_2(1) * Dr_2_2(1,0)) 
			/ sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		Dfunc(0,0) = (Dr_2_2(2,0) * coorX - r_2_2(2) * DcoorX) 
			/ (coorX * coorX + r_2_2(2) * r_2_2(2));
		Dfunc(1,0) = 2.0 * r_2_2(2) * Dr_2_2(2,0) + 2.0 * coorX * DcoorX;
		//derivative 2
		DcoorX = - (r_2_2(0) * Dr_2_2(0,1) + r_2_2(1) * Dr_2_2(1,1)) 
			/ sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		Dfunc(0,1) = (Dr_2_2(2,1) * coorX - r_2_2(2) * DcoorX) 
			/ (coorX * coorX + r_2_2(2) * r_2_2(2));
		Dfunc(1,1) = 2.0 * r_2_2(2) * Dr_2_2(2,1) + 2.0 * coorX * DcoorX;
		//iteration
		deltX = Dfunc.fullPivLu().solve(-func);
		if(deltX.norm() < 1.0E-11){
			if(coorX < 0.0){
				std::cout << "ERROR in DEHWSURF::WHEE_CURV_2_CART_3!" << std::endl;
			}
			break;
		}
		double rfac = 2.0;
		long rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			Eigen::Vector2d x_test = x + rfac * deltX;
			if(x_test(0) < 0.01 * PI || x_test(0) > 0.49 * PI){
				continue;
			}
			double thet_ct = x_test(0);
			double thet_ht = x_test(1);
			//function
			Eigen::Vector3d r_2_2t;
			WHEE_TRAN(thet_ct, thet_ht, xi_11, r_2_2t, Dr_2_2);
			double cooX_t = a_h2 - sqrt(r_2_2t(0) * r_2_2t(0) + r_2_2t(1) * r_2_2t(1));
			Eigen::Vector2d func_t;
			func_t << atan2(r_2_2t(2), cooX_t) - xi_21,
				r_2_2t(2) * r_2_2t(2) + cooX_t * cooX_t - xi_22 * xi_22;
			if(func_t.norm() < func.norm()){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
	return 1;
}

long DEHWSURF::WHEE_TRAN(double thet_c, double thet_h, double xi_11, 
	Eigen::Vector3d &r_2_2, Eigen::Matrix<double,3,2> &Dr_2_2){
	double thet_1 = i_1c * thet_c;
	double thet_2 = i_2h * thet_h;
	//
	double C_1 = (tan(beta_c) * cos(thet_c) + i_c1) * cos(thet_1 - xi_11) 
		+ i_c1 * tan(beta_c) * sin(thet_c) * sin(thet_1 - xi_11) 
		- cos(thet_c) * sin(thet_c) * sin(thet_1 - xi_11);
	double C_2 = i_c1 * r_b2 * sin(thet_c) - i_c1 * a_1c;
	double x_a = - C_2 / C_1;
	double z_a = ((tan(beta_c) * sin(thet_1 - xi_11) + sin(thet_c) * cos(thet_1 - xi_11)) 
		* x_a + r_b2 - a_1c * sin(thet_c)) / cos(thet_c);
	Eigen::Vector3d r_1_1;
	r_1_1 << x_a * cos(xi_11), - x_a * sin(xi_11), z_a;
	Eigen::Matrix<double,3,3> R_oh_h;
	R_oh_h << cos(thet_h), - sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	Eigen::Vector3d r_h_oh = R_oh_h * r_1_1;
	Eigen::Matrix<double,3,3> R_o2_oh;
	R_o2_oh << 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, -1.0, 0.0;
	Eigen::Vector3d T_o2_oh;
	T_o2_oh << - a_h2, 0.0, 0.0;
	Eigen::Vector3d r_h_o2 = R_o2_oh * r_h_oh + T_o2_oh;
	Eigen::Matrix<double,3,3> R_2_o2;
	R_2_o2 << cos(thet_2), sin(thet_2), 0.0,
		- sin(thet_2), cos(thet_2), 0.0,
		0.0, 0.0, 1.0;
	r_2_2 = R_2_o2 * r_h_o2;
	//derivative relative to thet_c
	double DC_1 = (- tan(beta_c) * sin(thet_c) + 0.0) * cos(thet_1 - xi_11) 
		- (tan(beta_c) * cos(thet_c) + i_c1) * i_1c * sin(thet_1 - xi_11) 
		+ i_c1 * tan(beta_c) * cos(thet_c) * sin(thet_1 - xi_11) 
		+ i_c1 * tan(beta_c) * sin(thet_c) * i_1c * cos(thet_1 - xi_11) 
		+ sin(thet_c) * sin(thet_c) * sin(thet_1 - xi_11) 
		- cos(thet_c) * cos(thet_c) * sin(thet_1 - xi_11) 
		- cos(thet_c) * sin(thet_c) * i_1c * cos(thet_1 - xi_11);
	double DC_2 = i_c1 * r_b2 * cos(thet_c) - 0.0;
	double Dx_a = - (DC_2 * C_1 - C_2 * DC_1) / (C_1 * C_1);
	double Dz_a = (((tan(beta_c) * i_1c * cos(thet_1 - xi_11) 
		+ cos(thet_c) * cos(thet_1 - xi_11) 
		- sin(thet_c) * i_1c * sin(thet_1 - xi_11)) * x_a
		+ (tan(beta_c) * sin(thet_1 - xi_11) + sin(thet_c) * cos(thet_1 - xi_11)) * Dx_a
		+ 0.0 - a_1c * cos(thet_c)) * cos(thet_c) 
		+ ((tan(beta_c) * sin(thet_1 - xi_11) + sin(thet_c) * cos(thet_1 - xi_11)) 
		* x_a + r_b2 - a_1c * sin(thet_c)) * sin(thet_c)) 
		/ (cos(thet_c) * cos(thet_c));
	Eigen::Vector3d Dr_1_1;
	Dr_1_1 << Dx_a * cos(xi_11), - Dx_a * sin(xi_11), Dz_a;
	Eigen::Vector3d Dr_h_oh = R_oh_h * Dr_1_1;
	Eigen::Vector3d Dr_h_o2 = R_o2_oh * Dr_h_oh;
	Dr_2_2.block(0,0,3,1) = R_2_o2 * Dr_h_o2;
	//derivative relative to thet_h
	Eigen::Matrix<double,3,3> DR_oh_h;
	DR_oh_h << - sin(thet_h), - cos(thet_h), 0.0,
		cos(thet_h), - sin(thet_h), 0.0,
		0.0, 0.0, 0.0;
	Dr_h_oh = DR_oh_h * r_1_1;
	Dr_h_o2 = R_o2_oh * Dr_h_oh;
	Eigen::Matrix<double,3,3> DR_2_o2;
	DR_2_o2 << - i_2h * sin(thet_2), i_2h * cos(thet_2), 0.0,
		- i_2h * cos(thet_2), - i_2h * sin(thet_2), 0.0,
		0.0, 0.0, 0.0;
	Dr_2_2.block(0,1,3,1) = DR_2_o2 * r_h_o2 + R_2_o2 * Dr_h_o2;
	return 1;
}

long DEHWSURF::WHEE_PHAS(long ti, long tj, long f_ij, Eigen::Vector3d r_2_2){
	if(fpha[ti][tj] == 0){
		fpha[ti][tj] = f_ij;
		cartCoor[1][ti][tj] = r_2_2;
	}
	else{
		double phas_1 = atan2(cartCoor[1][ti][tj](1), cartCoor[1][ti][tj](0));
		if(phas_1 < 0.0){
			phas_1 += 2.0 * PI;
		}
		double phas_2 = atan2(r_2_2(1), r_2_2(0));
		if(phas_2 < 0.0){
			phas_2 += 2.0 *PI;
		}
		if(phas_2 > phas_1){
			fpha[ti][tj] = f_ij;
			cartCoor[1][ti][tj] = r_2_2;
		}
	}
	return 1;
}

long DEHWSURF::WORM_RELI(Eigen::Vector3d &tempCoor, long ti, long tj){
	long reliLeng = 40;
	if(tj > curvCoor[0][0].size() - reliLeng 
		|| ti < reliLeng || ti > curvCoor[0].size() - reliLeng){
		Eigen::Vector3d tempXYZ = tempCoor;
		std::vector<double> reliAmou = {14.0E-6, 18.0E-6};//0 - tip, 1 - end
		double reliExpo = 3.0;
		double tempReli = 0.0;
		if(tj > curvCoor[0][0].size() - reliLeng 
			&& ti >= reliLeng && ti <= curvCoor[0].size() - reliLeng){
			tempReli = 
				pow((tj - (curvCoor[0][0].size() - reliLeng)) / (double)reliLeng, reliExpo) 
				* reliAmou[0];
		}
		else if(tj > curvCoor[0][0].size() - reliLeng && ti < reliLeng){
			double tempAngl = 
				atan((tj - (curvCoor[0][0].size() - reliLeng)) / (double)(reliLeng - ti));
			double tempRati = tempAngl / (PI / 2.0);
			double maxiAmou = reliAmou[1] 
				+ (- 1.0 + cos(tempRati * PI)) * (reliAmou[1] - reliAmou[0]) / 2.0;
			double tempRadi = sqrt(pow(tj - (curvCoor[0][0].size() - reliLeng), 2.0) 
				+ pow(reliLeng - ti, 2.0));
			tempReli = pow(tempRadi / (double)reliLeng, reliExpo) * maxiAmou;
		}
		else if(tj > curvCoor[0][0].size() - reliLeng && ti > curvCoor[0].size() - reliLeng){
			double tempAngl = atan((tj - (curvCoor[0][0].size() - reliLeng)) 
				/ (double)(ti - (curvCoor[0].size() - reliLeng)));
			double tempRati = tempAngl / (PI / 2.0);
			double maxiAmou = reliAmou[1] 
				+ (- 1.0 + cos(tempRati * PI)) * (reliAmou[1] - reliAmou[0]) / 2.0;
			double tempRadi = sqrt(pow(tj - (curvCoor[0][0].size() - reliLeng), 2.0) 
				+ pow(ti - (curvCoor[0].size() - reliLeng), 2.0));
			tempReli = pow(tempRadi / (double)reliLeng, reliExpo) * maxiAmou;
		}
		else if(tj <= curvCoor[0][0].size() - reliLeng && ti < reliLeng){
			tempReli = pow((reliLeng - ti) / (double)reliLeng, reliExpo) * reliAmou[1];
		}
		else if(tj <= curvCoor[0][0].size() - reliLeng && ti > curvCoor[0].size() - reliLeng){
			tempReli = 
				pow((ti - (curvCoor[0].size() - reliLeng)) / (double)reliLeng, reliExpo) 
				* reliAmou[1];
		}
		if(abs(tempReli) > 1.0E-12){
			double tempRadi_0 = sqrt(pow(tempXYZ(0), 2.0) + pow(tempXYZ(1), 2.0));
			double tempRadi = a_h2 - tempRadi_0;
			tempRadi = sqrt(tempRadi * tempRadi + tempXYZ(2) * tempXYZ(2));
			double tempAngl = tempReli / tempRadi;
			double tempThet_0 = asin(tempXYZ(2) / tempRadi);
			double tempThet_1 = tempThet_0 + tempAngl;
			double tempRadi_1 = a_h2 - tempRadi * cos(tempThet_1);
			double tempFact = tempRadi_1 / tempRadi_0;
			tempXYZ(0) = tempFact * tempXYZ(0);
			tempXYZ(1) = tempFact * tempXYZ(1);
			tempXYZ(2) = tempXYZ(2) + tempRadi * (sin(tempThet_1) - sin(tempThet_0));
			tempCoor = tempXYZ;
		}
	}
	return 1;
}

long DEHWSURF::WHEE_RELI(Eigen::Vector3d &tempCoor, long ti, long tj){
	double reliLeng = 40;
	if(tj < reliLeng || ti < reliLeng || ti > curvCoor[1].size() - reliLeng){
		Eigen::Vector3d tempXYZ = tempCoor;
		double reliExpo = 3.0;
		std::vector<double> reliAmou = {12.0E-6, 16.0E-6};//0 - tip, 1 - end
		//
		double tempReli = 0.0;
		if(ti < reliLeng){
			if(tj >= reliLeng){
				tempReli = pow((reliLeng - ti) / (double)reliLeng, reliExpo) * reliAmou[1];
			}
			else{
				double tempAngl = atan((reliLeng - tj) / (double)(reliLeng - ti));
				double tempRati = tempAngl / (PI / 2.0);
				double maxiAmou = reliAmou[1] 
					+ (- 1.0 + cos(tempRati * PI)) * (reliAmou[1] - reliAmou[0]) / 2.0;
				double tempRadi = sqrt(pow((reliLeng - tj), 2.0) + pow((reliLeng - ti), 2.0));
				tempReli = pow((tempRadi / reliLeng), reliExpo) * maxiAmou;
			}
		}
		else if(ti > curvCoor[1].size() - reliLeng){
			if(tj >= reliLeng){
				tempReli = 
					pow((ti - (curvCoor[1].size() - reliLeng)) / (double)reliLeng, reliExpo) 
					* reliAmou[1];
			}
			else{
				double tempAngl = atan((reliLeng - tj) 
					/ (double)(ti - (curvCoor[1].size() - reliLeng)));
				double tempRati = tempAngl / (PI / 2.0);
				double maxiAmou = reliAmou[1] 
					+ (- 1.0 + cos(tempRati * PI)) * (reliAmou[1] - reliAmou[0]) / 2.0;
				double tempRadi = sqrt(pow((reliLeng - tj), 2.0) 
					+ pow((ti - (curvCoor[1].size() - reliLeng)), 2.0));
				tempReli = pow((tempRadi / reliLeng), reliExpo) * maxiAmou;
			}
		}
		else if(tj < reliLeng){
			tempReli = pow((reliLeng - tj) / (double)reliLeng, reliExpo) * reliAmou[0];
		}
		if(abs(tempReli) > 1.0E-12){
			double tempAngl = tempReli 
				/ sqrt(tempXYZ(0) * tempXYZ(0) + tempXYZ(1) * tempXYZ(1));
			Eigen::Matrix<double,3,3> tempMatr;
			tempMatr << cos(tempAngl), -sin(tempAngl), 0.0,
				sin(tempAngl), cos(tempAngl), 0.0,
				0.0, 0.0, 1.0;
			tempXYZ = tempMatr * tempXYZ;
			tempCoor = tempXYZ;
		}
	}
	return 1;
}

long DEHWSURF::NEW_CONT_ZONE(long f_lr){
	//initial value
	double thet_c, thet_h;
	Eigen::Vector3d r_2_2;
	double angl_fi, radi_fi, R_fmini, R_fmaxi;
	long numb_c = 1000;
	long numb_h = 10000;
	double thet_cL = 0.01 * PI;
	double thet_cH = 0.49 * PI;
	double epsl_t = 1.0E-8;
	long F_init = 0;
	for(long ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (double)numb_c * ti;
		double thet_hs, thet_hm;
		SINGULAR_C2H(thet_c, thet_hs, thet_hm);
		double thet_hL, thet_hH;
		if(f_lr == 1){
			thet_hL = thet_hs + epsl_t;
			thet_hH = thet_hm - epsl_t;
		}
		else{
			thet_hL = thet_hm + epsl_t;
			thet_hH = thet_hs + 2.0 * PI - epsl_t;
		}
		if(thet_hL >= thet_hH){
			continue;
		}
		for(long tj = 0; tj <= numb_h && F_init == 0; tj ++){
			if(f_lr == 1){
				thet_h = thet_hH - (thet_hH - thet_hL) / (double)numb_h * tj;
			}
			else{
				thet_h = thet_hL + (thet_hH - thet_hL) / (double)numb_h * tj;
			}
			//????????????????????cillofe, !!!!!the choosing order of thet_h!!!!!
			double thet_1 = i_1c * thet_c;
			double x_d, y_d;
			FSME(thet_1, thet_h, x_d, y_d);
			WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2);
			double kapp_h2N;
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
		return -1;
	}
	//closest point
	double miniDist = 1.0E20;
	std::list<INDE_INIT> breaList;
	Eigen::MatrixXi F_sear = Eigen::MatrixXi::Zero(curvCoor[1].size(), 
		curvCoor[1][0].size()
	);
	for(long ti = 0; ti < curvCoor[1].size() - 1; ti ++){
		for(long tj = 0; tj < curvCoor[1][ti].size() - 1; tj ++){
			double epsl_x = (curvCoor[1][ti][tj](0) - curvCoor[1][ti + 1][tj](0)) / 4.0;
			double epsl_y = (curvCoor[1][ti][tj + 1](1) - curvCoor[1][ti][tj](1)) / 4.0;
			if(curvCoor[1][ti + 1][tj](0) - epsl_x <= angl_fi && 
				angl_fi <= curvCoor[1][ti][tj](0) + epsl_x && 
				curvCoor[1][ti][tj](1) - epsl_y <= radi_fi && 
				radi_fi <= curvCoor[1][ti][tj + 1](1) + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, thet_h));
				F_sear(ti, tj) = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, thet_h));
				F_sear(ti + 1, tj) = 1;
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, thet_h));
				F_sear(ti, tj + 1) = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, thet_h));
				F_sear(ti + 1, tj + 1) = 1;
			}
			double dist_ij = radi_fi * abs(angl_fi - curvCoor[1][ti][tj](0)) 
				+ abs(radi_fi - curvCoor[1][ti][tj](1));
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	std::cout << "DEHWSURF::NEW_CONT_ZONE" << f_lr << ", miniDist = " << miniDist 
		<< ", initial number = " << breaList.size() << std::endl;
	//breadth-first search
	double epsl_d = 1.0E-9;
	long coun_w = 0;
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		F_sear(temp_ij.ti, temp_ij.tj) = 2;
		//
		double x_d, y_d;
		WHEE_CURV_2_CART_1(curvCoor[1][temp_ij.ti][temp_ij.tj](0), 
			curvCoor[1][temp_ij.ti][temp_ij.tj](1), 
			r_2_2, temp_ij.init_1, temp_ij.init_2, f_lr, x_d, y_d
		);
		//
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		double dist_ij = 
			radi_fi * abs(angl_fi - curvCoor[1][temp_ij.ti][temp_ij.tj](0)) 
			+ abs(radi_fi - curvCoor[1][temp_ij.ti][temp_ij.tj](1));
		if(dist_ij < epsl_d){
			//
			thet_c = temp_ij.init_1;
			double thet_1 = i_1c * thet_c;
			Eigen::Vector3d r_c_c, T_o1_oc, r_1_o1;
			Eigen::Matrix<double,3,3> R_oc_c, R_o1_oc, R_1_o1;
			r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
			R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
				sin(thet_c), cos(thet_c), 0.0,
				0.0, 0.0, 1.0;
			R_o1_oc << 1.0, 0.0, 0.0,
				0.0, 0.0, -1.0,
				0.0, 1.0, 0.0;
			R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
				- sin(thet_1), cos(thet_1), 0.0,
				0.0, 0.0, 1.0;
			T_o1_oc << a_1c, 0.0, 0.0;
			r_1_o1 = R_o1_oc * R_oc_c * r_c_c + T_o1_oc;
			double woxi_11 =  thet_1 - atan2(r_1_o1(1), r_1_o1(0));
			if(curvCoor[0][0][0](0) - 1.0E-12 <= woxi_11 
				&& woxi_11 <= curvCoor[0][curvCoor[0].size() - 1][0](0) + 1.0E-12){
				WHEE_PHAS(temp_ij.ti, temp_ij.tj, f_lr, r_2_2);
			}
			//
			Eigen::Matrix<long,8,2> inde;
			inde << temp_ij.ti - 1, temp_ij.tj - 1, 
				temp_ij.ti, temp_ij.tj - 1, 
				temp_ij.ti + 1, temp_ij.tj - 1, 
				temp_ij.ti - 1, temp_ij.tj, 
				temp_ij.ti + 1, temp_ij.tj, 
				temp_ij.ti - 1, temp_ij.tj + 1, 
				temp_ij.ti, temp_ij.tj + 1, 
				temp_ij.ti + 1, temp_ij.tj + 1;
			for(long ti = 0; ti < inde.rows(); ti ++){
				if(inde(ti,0) >= 0 && inde(ti,1) >= 0 
					&& inde(ti,0) < curvCoor[1].size() && inde(ti,1) < curvCoor[1][0].size() 
					&& F_sear(inde(ti,0), inde(ti,1)) == 0){
					F_sear(inde(ti,0), inde(ti,1)) = 1;
					breaList.push_back(INDE_INIT(inde(ti,0), inde(ti,1), 
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
	return 1;
}

long DEHWSURF::FORMER_CONT_ZONE(){
	//initial value
	double thet_c, x_d;
	Eigen::Vector3d r_c_c;
	double angl_fi, radi_fi, R_fmini, R_fmaxi;
	long numb_c = 1000;
	double thet_cL = 0.01 * PI;
	double thet_cH = 0.49 * PI;
	long numb_d = 10000;
	double x_dL = - 10.0 * a_1c;
	double x_dH = 10.0 * a_1c;
	long F_init = 0;
	for(long ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (double)numb_c * ti;
		for(long tj = 1; tj < numb_d && F_init == 0; tj ++){
			x_d = x_dL + (x_dH - x_dL) / (double)numb_d * tj;
			double y_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
				- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c)) / sin(thet_c);
			r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
			WHEE_G2L(r_c_c, angl_fi, radi_fi, R_fmini, R_fmaxi);
			if(-widtAngl <= angl_fi && angl_fi <= widtAngl 
				&& R_fmini <= radi_fi && radi_fi <= R_fmaxi){
				F_init = 1;
			}
		}
	}
	if(F_init == 0){
		std::cout << "WARNING in DEHWSURF::FORMER_CONT_ZONE!" << std::endl;
		return -1;
	}
	//closest point
	double miniDist = 1.0E20;
	std::list<INDE_INIT> breaList;
	Eigen::MatrixXi F_sear = Eigen::MatrixXi::Zero(curvCoor[1].size(), curvCoor[1][0].size());
	for(long ti = 0; ti < curvCoor[1].size() - 1; ti ++){
		for(long tj = 0; tj < curvCoor[1][0].size() - 1; tj ++){
			double epsl_x = (curvCoor[1][ti][tj](0) - curvCoor[1][ti + 1][tj](0)) / 4.0;
			double epsl_y = (curvCoor[1][ti][tj + 1](1) - curvCoor[1][ti][tj](1)) / 4.0;
			if(curvCoor[1][ti + 1][tj](0) - epsl_x <= angl_fi && 
				angl_fi <= curvCoor[1][ti][tj](0) + epsl_x && 
				curvCoor[1][ti][tj](1) - epsl_y <= radi_fi && 
				radi_fi <= curvCoor[1][ti][tj + 1](1) + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, x_d));
				F_sear(ti, tj) = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, x_d));
				F_sear(ti + 1, tj) = 1;
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, x_d));
				F_sear(ti, tj + 1) = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, x_d));
				F_sear(ti + 1, tj + 1) = 1;
			}
			double dist_ij = radi_fi * abs(angl_fi - curvCoor[1][ti][tj](0)) 
				+ abs(radi_fi - curvCoor[1][ti][tj](1));
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	std::cout << "DEHWSURF::FORMER_CONT_ZONE, miniDist = " << miniDist 
		<< ", initial number = " << breaList.size() << std::endl;
	//breadth-first search
	double epsl_d = 1.0E-9;
	long coun_w = 0;
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		F_sear(temp_ij.ti, temp_ij.tj) = 2;
		//
		double y_d;
		WHEE_CURV_2_CART_2(curvCoor[1][temp_ij.ti][temp_ij.tj](0), 
			curvCoor[1][temp_ij.ti][temp_ij.tj](1), 
			r_c_c, temp_ij.init_1, temp_ij.init_2, y_d
		);
		//
		WHEE_G2L(r_c_c, angl_fi, radi_fi, R_fmini, R_fmaxi);
		double dist_ij = 
			radi_fi * abs(angl_fi - curvCoor[1][temp_ij.ti][temp_ij.tj](0)) 
			+ abs(radi_fi - curvCoor[1][temp_ij.ti][temp_ij.tj](1));
		if(dist_ij < epsl_d){
			//
			thet_c = temp_ij.init_1;
			double thet_1 = i_1c * thet_c;
			Eigen::Vector3d T_o1_oc, r_1_o1;
			Eigen::Matrix<double,3,3> R_oc_c, R_o1_oc, R_1_o1;
			R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
				sin(thet_c), cos(thet_c), 0.0,
				0.0, 0.0, 1.0;
			R_o1_oc << 1.0, 0.0, 0.0,
				0.0, 0.0, -1.0,
				0.0, 1.0, 0.0;
			R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
				- sin(thet_1), cos(thet_1), 0.0,
				0.0, 0.0, 1.0;
			T_o1_oc << a_1c, 0.0, 0.0;
			r_1_o1 = R_o1_oc * R_oc_c * r_c_c + T_o1_oc;
			double woxi_11 =  thet_1 - atan2(r_1_o1(1), r_1_o1(0));
			if(curvCoor[0][0][0](0) - 1.0E-12 <= woxi_11 
				&& woxi_11 <= curvCoor[0][curvCoor[0].size() - 1][0](0) + 1.0E-12){
				WHEE_PHAS(temp_ij.ti, temp_ij.tj, 3, r_c_c);
			}
			//
			Eigen::Matrix<long,8,2> inde;
			inde << temp_ij.ti - 1, temp_ij.tj - 1, 
				temp_ij.ti, temp_ij.tj - 1, 
				temp_ij.ti + 1, temp_ij.tj - 1, 
				temp_ij.ti - 1, temp_ij.tj, 
				temp_ij.ti + 1, temp_ij.tj, 
				temp_ij.ti - 1, temp_ij.tj + 1, 
				temp_ij.ti, temp_ij.tj + 1, 
				temp_ij.ti + 1, temp_ij.tj + 1;
			for(long ti = 0; ti < inde.rows(); ti ++){
				if(inde(ti,0) >= 0 && inde(ti,1) >= 0 
					&& inde(ti,0) < curvCoor[1].size() && inde(ti,1) < curvCoor[1][0].size() 
					&& F_sear(inde(ti,0), inde(ti,1)) == 0){
					F_sear(inde(ti,0), inde(ti,1)) = 1;
					breaList.push_back(INDE_INIT(inde(ti,0), inde(ti,1), 
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
	return 1;
}

long DEHWSURF::TRANSITION_ZONE(long f_hr){
	double xi_11;
	if(f_hr == 1){
		xi_11 = curvCoor[0][0][0](0);//head
	}
	else{
		xi_11 = curvCoor[0][curvCoor[0].size() - 1][0](0);//rear
	}
	//initial value
	double thet_c, thet_h;
	Eigen::Vector3d r_2_2, r_1_1;
	double angl_fi, radi_fi, R_fmini, R_fmaxi;
	long numb_c = 1000;
	double thet_cL, thet_cH;
	WORM_CURV_2_CART(xi_11, a_h2 - d_f[0] / 2.0, r_1_1, thet_cL);
	WORM_CURV_2_CART(xi_11, d_f[1] / 2.0, r_1_1, thet_cH);
	thet_h = xi_11;
	long F_init = 0;
	for(long ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (double)numb_c * ti;
		Eigen::Matrix<double,3,2> Dr_2_2;
		WHEE_TRAN(thet_c, thet_h, xi_11, r_2_2, Dr_2_2);
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		if(-widtAngl <= angl_fi && angl_fi <= widtAngl 
			&& R_fmini <= radi_fi && radi_fi <= R_fmaxi){
			F_init = 1;
		}
	}
	if(F_init == 0){
		std::cout << "WARNING in DEHWSURF::TRANSITION_ZONE" << f_hr << "!" << std::endl;
		return -1;
	}
	//closest point
	double miniDist = 1.0E20;
	std::list<INDE_INIT> breaList;
	Eigen::MatrixXi F_sear = Eigen::MatrixXi::Zero(curvCoor[1].size(), curvCoor[1][0].size());
	for(long ti = 0; ti < curvCoor[1].size() - 1; ti ++){
		for(long tj = 0; tj < curvCoor[1][ti].size() - 1; tj ++){
			double epsl_x = (curvCoor[1][ti][tj](0) - curvCoor[1][ti + 1][tj](0)) / 4.0;
			double epsl_y = (curvCoor[1][ti][tj + 1](1) - curvCoor[1][ti][tj](1)) / 4.0;
			if(curvCoor[1][ti + 1][tj](0) - epsl_x <= angl_fi && 
				angl_fi <= curvCoor[1][ti][tj](0) + epsl_x && 
				curvCoor[1][ti][tj](1) - epsl_y <= radi_fi && 
				radi_fi <= curvCoor[1][ti][tj + 1](1) + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, thet_h));
				F_sear(ti, tj) = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, thet_h));
				F_sear(ti + 1, tj) = 1;
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, thet_h));
				F_sear(ti, tj + 1) = 1;
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, thet_h));
				F_sear(ti + 1, tj + 1) = 1;
			}
			double dist_ij = radi_fi * abs(angl_fi - curvCoor[1][ti][tj](0)) 
				+ abs(radi_fi - curvCoor[1][ti][tj](1));
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	std::cout << "DEHWSURF::TRANSITION_ZONE" << f_hr << ", miniDist = " << miniDist 
		<< ", initial number = " << breaList.size() << std::endl;
	//breadth-first search
	double epsl_d = 1.0E-9;
	long coun_w = 0;
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		F_sear(temp_ij.ti, temp_ij.tj) = 2;
		//
		WHEE_CURV_2_CART_3(curvCoor[1][temp_ij.ti][temp_ij.tj](0), 
			curvCoor[1][temp_ij.ti][temp_ij.tj](1), 
			r_2_2, temp_ij.init_1, temp_ij.init_2, xi_11
		);
		//
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		double dist_ij = 
			radi_fi * abs(angl_fi - curvCoor[1][temp_ij.ti][temp_ij.tj](0)) 
			+ abs(radi_fi - curvCoor[1][temp_ij.ti][temp_ij.tj](1));
		if(dist_ij < epsl_d){
			WHEE_PHAS(temp_ij.ti, temp_ij.tj, 3 + f_hr, r_2_2);
			Eigen::Matrix<long,8,2> inde;
			inde << temp_ij.ti - 1, temp_ij.tj - 1, 
				temp_ij.ti, temp_ij.tj - 1, 
				temp_ij.ti + 1, temp_ij.tj - 1, 
				temp_ij.ti - 1, temp_ij.tj, 
				temp_ij.ti + 1, temp_ij.tj, 
				temp_ij.ti - 1, temp_ij.tj + 1, 
				temp_ij.ti, temp_ij.tj + 1, 
				temp_ij.ti + 1, temp_ij.tj + 1;
			for(long ti = 0; ti < inde.rows(); ti ++){
				if(inde(ti,0) >= 0 && inde(ti,1) >= 0 
					&& inde(ti,0) < curvCoor[1].size() && inde(ti,1) < curvCoor[1][0].size() 
					&& F_sear(inde(ti,0), inde(ti,1)) == 0){
					F_sear(inde(ti,0), inde(ti,1)) = 1;
					breaList.push_back(INDE_INIT(inde(ti,0), inde(ti,1), 
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
	return 1;
}

long DEHWSURF::WORM_ROOT_RADIUS(long flag, Eigen::Matrix<double,2,3> tempPoin, 
	Eigen::Vector2d &tempCent, double &tempRadi, Eigen::Vector2d &tempAngl){
	Eigen::Vector2d vect_1 = tempPoin.block(0,1,2,1) - tempPoin.block(0,0,2,1);
	vect_1 = vect_1.eval() / vect_1.norm();
	Eigen::Vector2d vect_2 = tempPoin.block(0,2,2,1) - tempPoin.block(0,0,2,1);
	double tempLeng_1 = vect_2.dot(vect_1);
	double tempLeng_2 = sqrt(vect_2(0) * vect_2(0) + vect_2(1) * vect_2(1) 
		- tempLeng_1 * tempLeng_1);
	double targVari = tempLeng_1 / (R_f[0] - tempLeng_2);
	double middAngl = asin(targVari / sqrt(1.0 + targVari * targVari)) - atan(1.0 / targVari);
	tempRadi = R_f[0] - tempLeng_1 / cos(middAngl);
	Eigen::Vector2d tempVect;
	tempVect << flag * vect_1(1), - flag * vect_1(0);
	tempCent = tempPoin.block(0,0,2,1) + tempRadi * tempVect;
	tempAngl(0) = atan2(- tempVect(1), - tempVect(0));
	tempAngl(1) = tempAngl(0) + flag * (PI / 2.0 - middAngl);
	return 1;
}

long DEHWSURF::WORM_ROOT(long indi, long flag, Eigen::MatrixXd &rootProf){
	//
	Eigen::Matrix<double,3,3> tempPoin;
	if(flag == 1){
		tempPoin.block(0,0,3,1) << wormTosu.indiPoin[indi][0][0], 
			wormTosu.indiPoin[indi][0][1], wormTosu.indiPoin[indi][0][2];
		tempPoin.block(0,1,3,1) << wormTosu.indiPoin[indi][1][0], 
			wormTosu.indiPoin[indi][1][1], wormTosu.indiPoin[indi][1][2];
	}
	else{
		tempPoin.block(0,0,3,1) << wormToba.indiPoin[indi][0][0], 
			wormToba.indiPoin[indi][0][1], wormToba.indiPoin[indi][0][2];
		tempPoin.block(0,1,3,1) << wormToba.indiPoin[indi][1][0], 
			wormToba.indiPoin[indi][1][1], wormToba.indiPoin[indi][1][2];
	}
	double tempXi11 = - curvCoor[0][indi][0](0);//(wormCurv[0] 
		// + (wormCurv[2] - wormCurv[0]) / (wormTosu.indiPoin.size() - 1) * indi);
	tempPoin.block(0,2,3,1) << a_h2 * cos(tempXi11), a_h2 * sin(tempXi11), 0.0;
	Eigen::Matrix<double,2,3> tempPoin_;
	Eigen::Matrix<double,3,3> tempMatr;
	tempMatr << 0.0, 0.0, 1.0,
		-cos(tempXi11), -sin(tempXi11), 0.0,
		sin(tempXi11), -cos(tempXi11), 0.0;
	for(long ti = 0; ti < 3; ti ++){
		Eigen::Matrix<double,3,1> tempPoin_i = tempMatr * tempPoin.block(0,ti,3,1).eval();
		tempPoin_(0,ti) = tempPoin_i(0);
		tempPoin_(1,ti) = tempPoin_i(1) + a_h2;
	}
	Eigen::Vector2d tempCent;
	double tempRadi;
	Eigen::Vector2d tempAngl;
	WORM_ROOT_RADIUS(flag, tempPoin_, tempCent, tempRadi, tempAngl);
	Eigen::Vector2d tempPoinArce;
	tempPoinArce << tempCent(0) + tempRadi * cos(tempAngl(1)),
		tempCent(1) + tempRadi * sin(tempAngl(1));
	double tempAnglArce = atan2(tempPoinArce(1), tempPoinArce(0));
	//
	double tempAnglRoot;
	double tempAnglStar = acos(r_b2 / (d[1] / 2.0))
		- i_2h * tempXi11 - tootThicAngl[0] / 2.0;
	if(flag == 1){
		tempAnglRoot = tempAnglStar + pitcAngl / 2.0;
	}
	else{
		tempAnglRoot = tempAnglStar - pitcAngl / 2.0;
	}
	//
	double sumLeng = flag * R_f[0] * (tempAnglRoot - tempAnglArce)
		+ flag * tempRadi * (tempAngl(1) - tempAngl(0));
	Eigen::Vector3d tempTran;
	tempTran << a_h2 * cos(tempXi11), a_h2 * sin(tempXi11), 0.0;
	for(long ti = 0; ti < rootProf.cols(); ti ++){
		double leng_i = sumLeng / (rootProf.cols() - 1) * (double)ti;
		Eigen::Vector3d poin_i;
		if(leng_i <= flag * R_f[0] * (tempAnglRoot - tempAnglArce)){
			double angl_i = tempAnglRoot - flag * leng_i / R_f[0];
			poin_i << R_f[0] * cos(angl_i), R_f[0] * sin(angl_i), 0.0;
		}
		else{
			leng_i = leng_i - flag * R_f[0] * (tempAnglRoot - tempAnglArce);
			double angl_i = tempAngl(1) - flag * leng_i / tempRadi;
			poin_i << tempCent(0) + tempRadi * cos(angl_i), 
				tempCent(1) + tempRadi * sin(angl_i), 0.0;
		}
		rootProf.block(0,ti,3,1) = tempMatr.transpose() * poin_i.eval() + tempTran;
	}
	return 1;
}

Eigen::Vector2d DEHWSURF::WHEE_UNCONE(Eigen::Vector3d tempXYZ, double tempAlph_3){
	Eigen::Vector2d tempXY;
	double r_2 = sqrt(pow(tempXYZ(0), 2.0) + pow(tempXYZ(1), 2.0));
	double r_1 = r_2 / cos(tempAlph_3);
	double alph_2 = atan2(tempXYZ(1), tempXYZ(0));
	double alph_1 = r_2 * alph_2 / r_1;
	tempXY << r_1 * cos(alph_1), r_1 * sin(alph_1);
	return tempXY;
}

Eigen::Vector3d DEHWSURF::WHEE_CONE(Eigen::Vector2d tempXY, double tempAlph_3){
	Eigen::Vector3d tempXYZ;
	double r_1 = sqrt(pow(tempXY(0), 2.0) + pow(tempXY(1), 2.0));
	double alph_1 = atan2(tempXY(1), tempXY(0));
	double r_2 = r_1 * cos(tempAlph_3);
	double alph_2 = r_1 * alph_1 / r_2;
	double r_3 = a_h2 / cos(tempAlph_3) - r_1;
	tempXYZ << r_2 * cos(alph_2), r_2 * sin(alph_2), r_3 * sin(tempAlph_3);
	return tempXYZ;
}

long DEHWSURF::WHEE_ROOT(long indi, long flag, Eigen::MatrixXd &rootProf){
	//
	Eigen::Matrix<Eigen::Vector2d,2,2> tempProf;
	double alph_3 = - curvCoor[1][indi][0](0);
		// (widtAngl - 2.0 * widtAngl / (wheeTosu.indiPoin.size() - 1) * (double)indi);
	double angl_ai = alph_3 - asin(offsR_a * sin(alph_3) / R_a[1]);
	double R_fmini = R_t[1];
	double R_fmaxi = (R_a[1] * cos(angl_ai) - offsR_a) / cos(alph_3);
	for(long ti = 0; ti <= 1; ti ++){//tooth surface or tooth back
		for(long tj = 0; tj <= 1; tj ++){
			Eigen::Vector3d tempXYZ;
			if(ti == 0){
				tempXYZ << wheeTosu.indiPoin[indi][tj][0], 
					wheeTosu.indiPoin[indi][tj][1], wheeTosu.indiPoin[indi][tj][2];
			}
			else{
				tempXYZ << wheeToba.indiPoin[indi][tj][0], 
					wheeToba.indiPoin[indi][tj][1], wheeToba.indiPoin[indi][tj][2];
			}
			double r_2 = sqrt(pow(tempXYZ(0), 2.0) + pow(tempXYZ(1), 2.0));
			double alph_2 = atan2(tempXYZ(1), tempXYZ(0));
			double r_3 = R_fmini 
				+ (R_fmaxi - R_fmini) / (wheeTosu.indiPoin[0].size() - 1) * (double)tj;
			double r_1 = a_h2 / cos(alph_3) - r_3;
			double alph_1 = r_2 * alph_2 / r_1;
			tempProf(tj,ti) << r_1 * cos(alph_1), r_1 * sin(alph_1);
		}
	}
	double r_f = a_h2 / cos(alph_3) - (a_h2 - d_f[1] / 2.0);
	double tempPitc = pitcAngl * cos(alph_3);
	//
	Eigen::Matrix<double,2,3> tempPoin;
	tempPoin.block(0,0,2,1) = tempProf(0,flag);
	tempPoin.block(0,1,2,1) = tempProf(1,flag);
	tempPoin.block(0,2,2,1) << 0.0, 0.0;
	Eigen::Vector2d vect_1 = tempPoin.block(0,0,2,1) - tempPoin.block(0,1,2,1);
	vect_1 = vect_1.eval() / vect_1.norm();
	Eigen::Vector2d vect_2 = tempPoin.block(0,2,2,1) - tempPoin.block(0,0,2,1);
	double tempLeng_1 = vect_2.dot(vect_1);
	double tempLeng_2 = sqrt(vect_2(0) * vect_2(0) + vect_2(1) * vect_2(1) 
		- tempLeng_1 * tempLeng_1);
	double targVari = tempLeng_1 / (r_f - tempLeng_2);
	double middAngl = asin(targVari / sqrt(1.0 + targVari * targVari)) - atan(1.0 / targVari);
	double tempRadi = tempLeng_1 / cos(middAngl) - r_f;
	double tempSign = (flag == 0) ? 1.0 : -1.0;
	Eigen::Vector2d tempVect;
	tempVect << -tempSign * vect_1(1), tempSign * vect_1(0);
	Eigen::Vector2d tempCent = tempPoin.block(0,0,2,1) + tempRadi * tempVect;
	Eigen::Vector2d tempAngl;
	tempAngl(0) = atan2(-tempVect(1), -tempVect(0));
	tempAngl(1) = tempAngl(0) + tempSign * (PI / 2.0 - middAngl);
	Eigen::Vector2d tempPoinArce;
	tempPoinArce << tempCent(0) + tempRadi * cos(tempAngl(1)),
		tempCent(1) + tempRadi * sin(tempAngl(1));
	double tempAnglArce = atan2(tempPoinArce(1), tempPoinArce(0));
	//
	double tempAnglRoot;
	Eigen::Vector2d tempPoinOppo = tempProf(0, 1 - flag);
	tempAnglRoot = (atan2(tempPoinOppo(1), tempPoinOppo(0)) 
		+ atan2(tempPoin(1,0), tempPoin(0,0))) / 2.0;
	tempAnglRoot -= tempSign * tempPitc / 2.0;
	//
	double sumLeng = r_f * tempSign * (tempAnglArce - tempAnglRoot)
		+ tempRadi * tempSign * (tempAngl(1) - tempAngl(0));
	for(long ti = 0; ti < rootProf.cols(); ti ++){
		double leng_i = sumLeng / (rootProf.cols() - 1) * ti;
		if(leng_i <= r_f * tempSign * (tempAnglArce - tempAnglRoot)){
			double angl_i = tempAnglRoot + tempSign * leng_i / r_f;
			Eigen::Vector2d tempResu;
			tempResu << r_f * cos(angl_i), r_f * sin(angl_i);
			rootProf.block(0,ti,3,1) = WHEE_CONE(tempResu, alph_3);
		}
		else{
			leng_i = leng_i - r_f * tempSign * (tempAnglArce - tempAnglRoot);
			double angl_i = tempAngl(1) - tempSign * leng_i / tempRadi;
			Eigen::Vector2d tempResu;
			tempResu << tempCent(0) + tempRadi * cos(angl_i), 
				tempCent(1) + tempRadi * sin(angl_i);
			rootProf.block(0,ti,3,1) = WHEE_CONE(tempResu, alph_3);
		}
	}
	return 1;
}

long DEHWSURF::WORM_TS_GRID(){
	OUTPUT_TIME("DEHWSURF::WORM_TS_GRID");
	//
	double domaCirc = 2.0 * PI / (double)circNumb;
	double deltTang = domaCirc / gridNumb[0][4];
	double inteStar = wormCurv[1];
	while(inteStar - domaCirc >= wormCurv[0]){
		inteStar = inteStar - domaCirc;
	}
	gridNumb[0][5] = ceil((inteStar - wormCurv[0]) / deltTang);
	double realStar = inteStar - gridNumb[0][5] * deltTang;
	double inteEndi = wormCurv[1];
	while(inteEndi + domaCirc <= wormCurv[2]){
		inteEndi = inteEndi + domaCirc;
	}
	gridNumb[0][6] = floor((inteEndi - inteStar) / domaCirc + 1.0E-10) + 2;
	//curvilinear coordinate
	VECT_RESI(curvCoor[0], 
		(gridNumb[0][4] * (gridNumb[0][6] - 2) + gridNumb[0][5] * 2) 
		* (1 << (globInho + globHomo + locaLeve)) + 1, 
		gridNumb[0][3] * (1 << (globHomo + locaLeve)) + 1
	);
	deltTang /= (1 << (globInho + globHomo + locaLeve));
	for(long ti = 0; ti < curvCoor[0].size(); ti ++){
		for(long tj = 0; tj < curvCoor[0][ti].size(); tj ++){
			curvCoor[0][ti][tj] << realStar + ti * deltTang,
				R_t[0] + (R_a[0] - R_t[0]) / (curvCoor[0][ti].size() - 1) * (double)tj;
		}
	}
	//Cartesian coordinate
	VECT_RESI(cartCoor[0], curvCoor[0].size(), curvCoor[0][0].size());
	for(long ti = 0; ti < curvCoor[0].size(); ti ++){
		if(ti % 1000 == 0){
			std::cout << ti << "/" << curvCoor[0].size() << std::endl;
		}
		for(long tj = 0; tj < curvCoor[0][ti].size(); tj ++){
			double thet_c;
			WORM_CURV_2_CART(curvCoor[0][ti][tj](0), 
				curvCoor[0][ti][tj](1), cartCoor[0][ti][tj], thet_c
			);
			if(reliSwit == 1){
				WORM_RELI(cartCoor[0][ti][tj], ti, tj);
			}
		}
	}
	return 1;
}

long DEHWSURF::WHEE_TS_GRID(){
	OUTPUT_TIME("DEHWSURF::WHEE_TS_GRID");
	//curvilinear coordinate
	VECT_RESI(curvCoor[1], gridNumb[1][4] * (1 << (globInho + globHomo + locaLeve)) + 1, 
		gridNumb[1][3] * (1 << (globHomo + locaLeve)) + 1
	);
	for(long ti = 0; ti < curvCoor[1].size(); ti ++){
		double angl_fi = widtAngl - 2.0 * widtAngl / (curvCoor[1].size() - 1) * (double)ti;
		double angl_ai = angl_fi - asin(offsR_a * sin(angl_fi) / R_a[1]);
		double R_fmini = (R_a[1] * cos(angl_ai) - offsR_a) / cos(angl_fi);
		double R_fmaxi = R_t[1];
		for(long tj = 0; tj < curvCoor[1][ti].size(); tj ++){
			curvCoor[1][ti][tj] << angl_fi, 
				R_fmini + (R_fmaxi - R_fmini) / (curvCoor[1][ti].size() - 1) * (double)tj;
		}
	}
	//Cartesian coordinate	
	VECT_RESI(cartCoor[1], curvCoor[1].size(), curvCoor[1][0].size());
	VECT_ASSI(fpha, curvCoor[1].size(), curvCoor[1][0].size(), 0);
	NEW_CONT_ZONE(1);//left new contact zone
	NEW_CONT_ZONE(2);//right new contact zone
	if(modiTran == 0.0 && modiCent == 0.0){
		FORMER_CONT_ZONE();//former contact zone
	}
	TRANSITION_ZONE(1);//head transition zone
	TRANSITION_ZONE(2);//rear transition zone
	//tooth flank relief
	if(reliSwit == 1){
		for(long ti = 0; ti < curvCoor[1].size(); ti ++){
			for(long tj = 0; tj < curvCoor[1][0].size(); tj ++){
				WHEE_RELI(cartCoor[1][ti][tj], ti, tj);
			}
		}
	}
	return 1;
}

long DEHWSURF::TOOT_SURF_GRID(){
	curvCoor.resize(2);
	cartCoor.resize(2);
	//worm tooth surface
	WORM_TS_GRID();
	VECT_RESI(wormTosu.indiPoin, curvCoor[0].size(), curvCoor[0][0].size());
	for(long ti = 0; ti < curvCoor[0].size(); ti ++){
		for(long tj = 0; tj < curvCoor[0][ti].size(); tj ++){
			Eigen::Vector3d tempCoor = cartCoor[0][ti][tj];
			wormTosu.INSERT(ti, tj, COOR(tempCoor(0), tempCoor(1), tempCoor(2)));
		}
	}
	//wheel tooth surface
	WHEE_TS_GRID();
	VECT_RESI(wheeTosu.indiPoin, curvCoor[1].size(), curvCoor[1][0].size());
	for(long ti = 0; ti < curvCoor[1].size(); ti ++){
		for(long tj = 0; tj < curvCoor[1][ti].size(); tj ++){
			Eigen::Vector3d tempCoor = cartCoor[1][ti][tj];
			wheeTosu.INSERT(curvCoor[1].size() - 1 - ti, curvCoor[1][ti].size() - 1 - tj, 
				COOR(tempCoor(0), tempCoor(1), tempCoor(2))
			);
		}
	}
	//worm tooth back
	VECT_RESI(wormToba.indiPoin, curvCoor[0].size(), curvCoor[0][0].size());
	Eigen::Matrix<double,3,3> tempMatr;
	tempMatr << cos(wormCurv[1]), -sin(wormCurv[1]), 0.0,
		sin(wormCurv[1]), cos(wormCurv[1]), 0.0,
		0.0, 0.0, 1.0;
	for(long ti = 0; ti < curvCoor[0].size(); ti ++){
		for(long tj = 0; tj < curvCoor[0][ti].size(); tj ++){
			Eigen::Vector3d tempXYZ = cartCoor[0][ti][tj];
			Eigen::Vector3d tempXYZ_1 = tempMatr * tempXYZ;
			tempXYZ << tempXYZ_1(0), -tempXYZ_1(1), -tempXYZ(2);
			tempXYZ_1 = tempMatr.transpose() * tempXYZ;
			wormToba.INSERT(curvCoor[0].size() - 1 - ti, tj, COOR(
				tempXYZ_1(0), tempXYZ_1(1), tempXYZ_1(2)
			));
		}
	}
	//wheel tooth back
	VECT_RESI(wheeToba.indiPoin, curvCoor[1].size(), curvCoor[1][0].size());
	tempMatr << cos(backAngl[1]), -sin(backAngl[1]), 0.0,
		-sin(backAngl[1]), -cos(backAngl[1]), 0.0,
		0.0, 0.0, -1.0;
	for(long ti = 0; ti < curvCoor[1].size(); ti ++){
		for(long tj = 0; tj < curvCoor[1][ti].size(); tj ++){
			Eigen::Vector3d tempXYZ = tempMatr 
				* cartCoor[1][curvCoor[1].size() - 1 - ti][curvCoor[1][ti].size() - 1 - tj];
			wheeToba.INSERT(curvCoor[1].size() - 1 - ti, tj, 
				COOR(tempXYZ(0), tempXYZ(1), tempXYZ(2)
			));
		}
	}
	// curvCoor[0].clear();
	cartCoor[0].clear();
	// curvCoor[1].clear();
	cartCoor[1].clear();
	fpha.clear();
	return 1;
}

long DEHWSURF::ROOT_TRAN_GRID(){
	//worm
	long numb_0 = (curvCoor[0].size() - 1) / (1 << (locaLeve)) + 1;
	long numb_1 = (gridNumb[0][0] / 2) * (1 << (globHomo)) + 1;
	VECT_RESI(wormRtsu.indiPoin, numb_0, numb_1);
	VECT_RESI(wormRtba.indiPoin, numb_0, numb_1);
	for(long ti = 0; ti < numb_0; ti ++){
		Eigen::MatrixXd rootProf;
		rootProf.resize(3, numb_1);
		WORM_ROOT(ti * (1 << locaLeve), 1, rootProf);
		for(long tj = 0; tj < numb_1; tj ++){
			wormRtsu.INSERT(ti, tj, COOR(rootProf(0,tj), rootProf(1,tj), rootProf(2,tj)));
		}
		WORM_ROOT(ti * (1 << locaLeve), -1, rootProf);
		for(long tj = 0; tj < numb_1; tj ++){
			wormRtba.INSERT(ti, tj, COOR(rootProf(0,tj), rootProf(1,tj), rootProf(2,tj)));
		}
	}
	//wheel
	numb_0 = gridNumb[1][4] * (1 << (globInho + globHomo)) + 1;
	numb_1 = (gridNumb[1][0] / 2) * (1 << (globHomo)) + 1;
	VECT_RESI(wheeRtsu.indiPoin, numb_0, numb_1);
	VECT_RESI(wheeRtba.indiPoin, numb_0, numb_1);
	for(long ti = 0; ti < numb_0; ti ++){
		Eigen::MatrixXd rootProf;
		rootProf.resize(3, numb_1);
		WHEE_ROOT(ti * (1 << locaLeve), 0, rootProf);
		for(long tj = 0; tj < numb_1; tj ++){
			wheeRtsu.INSERT(ti, tj, COOR(rootProf(0,tj), rootProf(1,tj), rootProf(2,tj)));
		}
		WHEE_ROOT(ti * (1 << locaLeve), 1, rootProf);
		for(long tj = 0; tj < numb_1; tj ++){
			wheeRtba.INSERT(ti, tj, COOR(rootProf(0,tj), rootProf(1,tj), rootProf(2,tj)));
		}
	}
	return 1;
}

long DEHWSURF::OUTPUT(){
	//
	std::vector<CURVEDS*> wowhSurf = {&wormTosu, &wormToba, &wormRtsu, &wormRtba, 
		&wheeTosu, &wheeToba, &wheeRtsu, &wheeRtba
	};
	std::vector<std::string> wowhFina = {
		"resuWOTS.txt", "resuWOTB.txt", "resuWORT.txt", "resuWORB.txt", 
		"resuWHTS.txt", "resuWHTB.txt", "resuWHRT.txt", "resuWHRB.txt"};
	//
	for(long tw = 0; tw < wowhSurf.size(); tw ++){
		std::ofstream tempOfst;
		tempOfst.open(DIRECTORY(wowhFina[tw]), std::ios::out);
		tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
		for(long ti = 0; ti < (* wowhSurf[tw]).indiPoin.size(); ti ++){
			for(long tj = 0; tj < (* wowhSurf[tw]).indiPoin[ti].size(); tj ++){
				tempOfst << std::setw(30) << (* wowhSurf[tw]).indiPoin[ti][tj][0] 
					<< std::setw(30) << (* wowhSurf[tw]).indiPoin[ti][tj][1] 
					<< std::setw(30) << (* wowhSurf[tw]).indiPoin[ti][tj][2] << std::endl;
			}
		}
		tempOfst.close();
	}
	return 1;
}

long DEHWSURF::ESTABLISH(){
	BASIC_PARAMETER();
	TOOT_SURF_GRID();
	ROOT_TRAN_GRID();
	if(debuMode >= 1){
		OUTPUT();
	}
	return 1;
}

#endif
