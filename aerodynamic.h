#pragma once
#include<string>
#include"math_missile.h"

const double ATMO_RHO0 = 1.225; //海平面大气密度(kg/m3)
const double ATMO_R0 = 6356.766; //大气常数
const double ATMO_P0 = 101325.0; //标准大气压
const double ATMO_SEASURF_TEMP = 25;//海平面气温(摄氏度)
const double ATMO_AIRGAMMA = 1.4; //空气比热比
const double  ATMO_RGAS = 287.14;//普适气体常数J/(kg*K)
const double  ATMO_K = 273.15;//摄氏开氏温度互换
//const double  ATMO_STBO = 5.6697e-12;//斯蒂芬-玻尔兹曼常数


/****************************************************************************************/
//说明：大气密度计算模块
//输入：高度
//输出：当前大气密度、当前温度、当前声速
//注意：输入高度单位是米
 /**************************************************************************************/
class Atmosphere_Model
{
public:  
	Atmosphere_Model();
	~Atmosphere_Model();

public:
	double phi;//当前大气密度(kg/m3)
	double Tempreture;//当前温度(摄氏度)
	double a_s; //当前声速(m/s)
	double  Kp_air;//大气压强
	double m_CalAtmosphere(double IN_dHeight);//输入的单位是km ,计算声速，返回胜诉，且更新温度
	double g_CalAtomspherePressure(double IN_dH);//根据高度和温度计算大气压强
};



/*********************************************************************************************************/
//说明：气动插值模块
//输入：飞行状态、高度、舵偏、参考面积、参考长度、插值表
//输出：气动力、气动力矩、气动力矩系数、气动力系数、马赫数、动压
//使用说明：直接调用CalAero函数即可，如果想拉偏气动系数，请先init_lapian，再调用CalAero函数
//注意：插值模块用的是角度，需要转化一下
/********************************************************************************************************/

class Aero_F
{
public:
	Aero_F();
	Aero_F(string filepath);
	~Aero_F();

public:
	double Ma; //马赫数
	double q_dyn; //动压
	Vector3d Cx; //插值输出的气动力系数
	Vector3d Cm; //插值输出的气动力矩系数

	Vector3d F_Aero; //气动力
	Vector3d M_Aero; //气动力矩


	Vector4d Rudder_physics; //物理舵偏角  X字舵 有耦合,此处的气动舵，不能跟 燃气舵混合
	Vector3d Rudder_math = Vector3d(0, 0, 0); //数学舵偏角
	bool iflapian; //拉偏开关
	Vector3d F_lapian; //气动系数拉偏，后续做模拟仿真用
	Vector3d M_lapian; //气动力矩系数拉偏，后续做模拟仿真用
	FlyStatus status; //飞行状态

	Atmosphere_Model atmos_m;
private:
	//这部分再修改一下
	
	double phi; //大气密度
	double a_s; //声速
	double S_ref; //参考面积
	double L_ref; //参考长度

	//以下是插值表
	bool iftable;
	//vector<double> C_alpha;
	//vector<double> C_beta;
	//vector<double> C_Ma;
	//vector<double> C_delta;
	vector<double> C_alpha = { -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30 };
	vector<double> C_beta = { -30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30 };
	vector<double> C_Ma = { 0.4,0.8,1.2,2.0,3.0,4.0,5.0,6.0,8.0,10.0 };
	vector<double> C_delta = { -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30 };
	vector<double> CN;
	vector<double> CA;
	vector<double> CM;
	vector<double> C_Mx_delta;
	vector<double> C_Mx;

public:
	void test(string filepath);
	void update_Aero(const FlyStatus& in_status, const Vector3d& in_Rudder);
	
	void init_lapian(Vector3d in_F_lapian, Vector3d in_M_lapian);

	void interpolation();
	Vector3d inverse_interpolation(Vector3d Cx, vector<double> Xs, vector<double> Ys, vector<double> Zs, Vector3d Rudder_math, Vector3d Ma);
	Vector3d F_current();
	Vector3d M_current();
	void output_FM();
	//void interpolation_F();

	void test_F();

	void readtable(vector <double> in_C_alpha, vector <double> in_C_beta, vector <double> in_C_Ma, vector <double> in_C_delta, vector <double> in_C_Mx_delta, vector <double> in_CN, vector <double> in_CA, vector <double> in_C_Mx, vector <double> in_C_Mz);

	void table(string filepath);

	void CalAero(const FlyStatus& in_status, const Vector3d& in_Rudder);

	double interp1(size_t n_Xs, vector<double> Xs, vector<double> Data, double X);

	double interp2(size_t n_Xs, vector<double> Xs, size_t n_Ys, vector<double> Ys, vector<double> Data, double X, double Y);

	double interp3(size_t n_Xs, vector<double> Xs, size_t n_Ys, vector<double> Ys, size_t n_Zs, vector<double> Zs, vector<double> Data, double X, double Y, double Z);

	double interp3_inverse(vector<double> Xs, vector<double> Ys, vector<double> Zs, vector<double> Data, double X, double Y, double Z);

};

