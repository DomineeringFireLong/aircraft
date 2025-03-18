#pragma once
#include"missile.h"
#include"math_missile.h"
#include"aerodynamic.h"
#include"controller.h"
#include"guide.h"
using namespace Eigen;
using namespace std;

class powered_phase_stand_ballastic
{
public:
	//标准弹道
	Vector3d  F,K1, K2, k_rel;//，空气力 重 推 科里奥利力
	double G, P;//注意主动段标准弹道与导弹有同名变量
	Vector3d P_t;
	Vector3d PM,AM,Expected_M,Expected_Acc;
	Vector3d program_angle = Vector3d((Pi-0.001) / 2, 0, 0);
	Vector3d program_omega = Vector3d(0, 0, 0);
	missile mymissile; 
	Vector3d expected_Aero_F{ 0,0,0 };
	controller mycontroller;
	Guidance_law myguidar;

	double engine_m = 750, gas_m = 507, mass_flow = 56.3333, specific_impulse=3089.0,first_engine_worktime=9, second_engine_worktime = 44;//此处给值进入不到构造函数，先执行构造函数
	double S=0.01, Pa=3.0e3;//发动机喷口的面积和压强,估计 todo
	
	Vector4d gas_vane_physical = Vector4d(0, 0, 0,0);
	Vector3d gas_vane_math = Vector3d(0, 0, 0);
	
	Aero_F aero_model;
	Vector3d x_c, x_g;//尾部到质心距离 和 弹体半径
	bool Ifdisplay = 0;
	//改模块融入了 发动机特性
	powered_phase_stand_ballastic(missile m1,string areo_data= "F:/SM_3/", double P = 174e3, Vector3d K1 = Vector3d(0, 0, 0), \
		Vector3d K2 = Vector3d(0, 0, 0), Vector3d PM = Vector3d(0, 0, 0), Vector3d AM = Vector3d(0, 0, 0), Vector3d point = Vector3d(3000,10000, 0),double second_engine_worktime1=44);//前面有默认初始值 后面必须有初始值
	void updata_parameter();

	Vector3d F_sum();
	Vector3d M_sum(Vector3d gas_vane_deviation, Vector3d rudder_deviation);
	//每段的控制力也就是推力公式不同,其他的力写到math库中

	//助推器产生的推力,n=0，1，2代表不同 数学舵
	Vector3d thrust(Vector3d gas_vane_deviation);//先直接考虑数学舵
	//高度和质量
	double Gravity(double h, double m);

	// 发射非惯性系、变质量产生的
	Vector3d K_1();
	Vector3d K_2();
	//力矩
	Vector3d thrust_Moment(Vector3d gas_vane_deviation);
	Vector3d Aero_Moment(Vector3d rudder_deviation);
	Vector3d Aero(Vector3d rudder_deviation);
	
	Vector3d Program_angle_turning(double time_current);
	//第一级 助推段
	void boost_turning_phase(bool Ifdisplay);
	//执行机构，燃气舵建模 期望M->舵偏
	Vector3d Cal_gas_vane_deviation(Vector3d Expected_torque);
	
	//第二级 主动段
	void powered_stage(bool Ifdisplay);

	//制导模块与控制模块的连接：a->F空气->攻角
	Vector3d Calculate_Expected_control_angle(Vector3d F_sum);

};