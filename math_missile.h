//#pragma once
#ifndef MATH_MISSILE
#define MATH_MISSILE


#include<iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <random>
#include <assert.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <direct.h>
#include <numeric>
#include"matplotlibcpp.h"
namespace plt = matplotlibcpp;

//#include "matplotlibcpp.h"
//namespace plt = matplotlibcpp;
//#include<matplotlibcpp.h>

#define Pi 3.14159265358979
#define PI 3.14159265358979
#define NATURE_E 2.718281828454904
//声明变量和函数作为外部变量或者外部函数供其它文件使用

#define LARGENUM  1e10
#define SMALLNUM  1e-10
#define a57 (180.0 / PI)
#define DEG2RAD (PI/180.0)
#define RAD2DEG (180.0/PI)
#define g_nature 9.8
#define Deg2Rad(a)        (a)*PI/180.0
#define Rad2Deg(a)        (a)*180.0/PI
//角度计算一律用弧度，输出时转化成度

#define iszero(x)         abs(x)<SMALLNUM
//地球参数
#define earth_r0  6371.0087714e3
#define earth_g0  9.80665
#define EPS		 1e-10
#define EARTH_RAVG  6371000	//地球平均半径(m)
#define EARTH_AE  6378140.0	//地球长半轴(m)
#define EARTH_BE  6356752.0	//地球短半轴(m)
#define earth_omega  7.292115e-5	//地球自转角速度(rad/s)
#define EARTH_FM  3.986e14	//地球引力常量(m3/s2)
#define EARTH_J2  1.08263e-3	//地球引力修正J2项

#define ATMO_TA 1	//平均大气透过率

using namespace Eigen;
using namespace std;
extern Vector3d default_Angle_expect;
Vector3d G2T(Vector3d gxoy, double theta1, double eta1);//地面->弹道
Vector3d T2G(Vector3d txoy, double theta1, double eta1);//弹道->地面

Vector3d B2T(Vector3d bxoy, double alpha1, double beta1);//弹体->弹道
Vector3d T2B(Vector3d txoy, double alpha1, double beta1);
Vector3d rao_z(Vector3d xoy, double theta1);
Vector3d rao_y(Vector3d xoy, double theta1);
//标量取符号
template<class T1>
T1 sign_scalar(T1 a)
{
	if (a > 0)
		return 1;
	else
	{
		if (a < 0)
			return -1;
		else
			return 0;
	}
}
//标量饱和
template<class T>
T sat(T x, T limit)
{
	if (x > limit)
		return limit;
	else {
		if (x < -limit)
			return -1 * limit;
		else
			return x;
	}
}
//向量
Vector3d Sign(Vector3d x);
double modulo(VectorXd x);
VectorXd abs(VectorXd x); 
VectorXd sat(VectorXd x, double limit);
VectorXd sign_vector(VectorXd x);
Vector3d calAlphaBetanu(Vector3d Eularangle, double theta, double psi);

Vector3d sature(Vector3d x, double limit);
Matrix3d Cross_product(Vector3d x);

Vector3d sign_sat(Vector3d x, double limit);

struct FlyStatus //飞行状态
{
	double m; //导弹质量
	double v_dandao; //弹道坐标系下的速度
	Vector3d v; //弹体坐标系下的速度
	Vector3d EularAngle; //导弹的姿态角：俯仰偏航滚转
	Vector3d attitude_omega; //导弹的转动角速度在弹体坐标轴的分量
	Vector3d Pos; //导弹位置，地面坐标系
	double alpha; //攻角
	double beta;  //侧滑角
	double theta; //弹道倾角
	double segma; //弹道偏角
	double nu; //倾侧角
	double S_ref;
	double L_ref;
	Vector3d J;
	double t;
	//FlyStatus();//构造函数
};

Vector3d Integration(double a, double b, Vector3d(*f)(double x), double interval = 0.01);
Vector3d SIN_V3(double x);
Vector3d  DefiniteIntegration_2D(double x0, double x1, double y0, double y1, Vector3d(*f)(double x, double y), double interval = 0.01);
Vector3d Variablelimit_Integrataion_2D_(double x0, double x1, Vector3d(*f)(double x), double interval = 0.005);
Vector3d Eye_V3(double x, double y);
void Integration_test();
//Vector3d Launch2LaunchIv(double B0, double A0, double t, const Vector3d v);


//数据处理的类,功能：绘制图像

class Data_Processor {
	ifstream infile;
	ofstream outfile;// 未知重写说明符	可能是需要 std::
	vector<double> Data;
public:
	Data_Processor();
	vector<double> read_file(string root);//读取文件
	void plot_curve_2();//绘制曲线
	void plot_curve_3(vector<double> x, vector<double> y, vector<double> z);
	void test();
};

#endif // !MATH_MISSILE