#pragma once
#include"missile.h"
#include"math_missile.h"
#include<fstream>
#include<random>
//分模块写成类：动力学（导弹）、制导、控制、整个模拟、
//各模块之间的信息传递，以及导弹对象的传递
//不同段的动力学方程的多态 融合


using namespace Eigen;
using namespace std;
class Guidance_law
{
	double t_final = 44;//预测制导的末端时间
	double t_start = 0;//开始时间
	double  Q = 10;//二次型性能指标的末端项的加权值
	double time_work=0;
public:
	Vector3d R, R_dot;//LOS坐标系
	double r, r_dot, theta, theta_dot,sigma, sigma_dot;//LOS坐标系
	Vector3d expect_point,expect_accelerate;
	missile SM_3,target;//可以扩展期望点为missile移动
	Vector3d ZCMD = Vector3d(10000, 80000, 0);
	Guidance_law(Vector3d point, double t_start1, double t_final1);

	//Vector3d expect_point;
	//powered_phase_guidance
	//主动段助推段段采用方案弹道，不用导引
	
	//大气层内中段采用预测制导
	Vector3d powered_stage(double Q1, double  now_time);//输入导弹参数和 制导参数，输出期望过载
	//末端
	void terminal_phase_guidance(missile m1, missile target);
	void terminal_phase_guidance(missile m1, Eigen::Vector3d expect_point1);
	Vector3d Integration_G(double a, double b, std::function<Vector3d(double)> f, double interval = 0.01);
	double get_theta_dot();
	Vector3d VL_Integrataion_2D(double x0, double x1, std::function<Vector3d(double)> f, double interval);
	virtual ~Guidance_law() {};//此处为纯虚；通过基类的指针来删除派生类的对象时，基类的析构函数应该是虚的。否则其删除效果将无法实现。
	//missile slidemode(int k1 = 6, int k2 = 0.1, int k3 = 6, int k4 = 0.1, double theta_e = PI / 6, double eta_e = PI / 6, string root = "./3D_slidemode_guide.txt");
};



//制导器继承弹道，包含其信息,末制导律是基于los系的
//class terminal_Guider:terminal_guidance {
//public:
//	Data_Buffer dandao_data;
//	missile SM_3;
//	terminal_Guider(missile m1, missile target);
//	terminal_Guider(missile m1, Vector3d expect_point);
//	//missile proportional_guide(int Ky = 4, int Kz = 4);
//	//missile slidemode(int k1 = 6, int k2 = 0.1, int k3 = 6, int k4 = 0.1, double theta_e = PI / 6, double eta_e = PI / 6, string root = "./Vertical_rising_section_guide.txt");
//
//};






