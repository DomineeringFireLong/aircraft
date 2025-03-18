#pragma once
#include<stdexcept>
#include"missile.h"
#include"math_missile.h"

using namespace Eigen;
using namespace std;
//#define Pi  3.14159265358979323846
//#define default_Angle_expect Vector3d(30 / Pi, 30 / Pi, 0)
//h文件只声明，不定义，编译器容易重复定义报错
//------------------------------------------------------//
//包含 滑模控制器，和执行机构模块(控制分配)，解算实际控制量

//extern double Pi;


//double sign(double x);
//
//double mo(Vector3d x);
//
//VectorXd sat(VectorXd x, double limit);
//
//double sat(double x, double limit);
//
//Vector3d abs(Vector3d x);
//
//VectorXd sign(VectorXd x);



//控制器，输入导弹的角度状态，输出期望的力矩
class controller
{
public:
	Vector3d F, M;	//VectorXd f;//此处初始化会报错：应输入类型说明符  在类中定义 如 vector<int> a(10);//会出现错误
	Vector3d V, g, f_, f2_;
	Matrix3d E, eta, k1,k2;//MatrixXd E=MatrixXd::Zero(3, 3);//Matrix<double, 3, 3> E;或Matrix3d E
	Vector3d Angle_expectation, omega_expectation;
	FlyStatus state;

	controller();
	Vector3d control_instruct_abr(FlyStatus state, Vector3d Angle_expectation,Vector3d Omg_expectation);
	Vector3d control_instruct_eular(FlyStatus state, Vector3d Angle_expectation, Vector3d Omg_expectation);
	//Vector3d control_instruct(double a, double b, double r, double wx, double wy, double wz, double ix, double iy, double iz, Vector3d Omg_expectation);
};




//（期望力矩和实际力矩不同）

