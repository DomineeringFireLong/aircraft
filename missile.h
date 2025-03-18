#pragma once
#ifndef _MISSILE_
#define _MISSILE_
#include <vector>
#include<iostream>
#include <cmath>
#include<string>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <random>
#include"math_missile.h"
//否则会循环依赖
//#include"controller.h"
//#include"guide.h"

//为了让函数可以在各个.cpp中共享，正确的做法就是在.h中只声明函数，并在另一个.cpp中实现这个函数。这样就不会冲突了。
using namespace Eigen;
using namespace std;
////把运动模块嵌入到导弹类中,运动学方程（质心，绕质心转动）

class missile
{
public:
	//质心运动学、绕质心运动学
	Eigen::Vector3d pos,vel,acc; //位置、速度、加速度 弹道坐标系  vel是速度坐标系 stt 弹道和速度相同
	Eigen::Vector3d Eular_angle, Eular_omega;//我们所谓的“姿态角”其实是欧拉角(是个过程量，三个角度的定义都不在同一个坐标系下)、加速度，角加速度
	Eigen::Vector3d attitude_omega, attitude_omega_acc;//真正的姿态角速度（绕弹体的角速度） //力矩直接作用
	//质心动力学、绕质心动力学
	Eigen::Vector3d F,P,G,A= Vector3d::Zero();//合力、推力、重力、空气动力(合力 基于弹道系)  4个轨控发动机（算入推力中）
	Eigen::Vector3d M, M_p, M_a;//合力矩、推力矩(发动机)、空气动力产生的力矩 (合力 基于弹体系)  6个姿控发动机（算入推力矩中）
	
	//角度自由度
	double theta, sigma, theta_dot, sigma_dot;//弹道俯仰、偏角 地面->弹道  theta sigma
	double alpha, beta, nu;//速度->弹体  速度->弹体 弹道->速度
	double alpha_dot, beta_dot, nu_dot;
	
	//三维取出来一维
	double velocity = this->vel[0]; // 弹道系  后者为地面系sqrt(pow(this->vel[0], 2) + pow(this->vel[1], 2) + pow(this->vel[2], 2));
	double x = pos[0], y = pos[1], z = pos[2];
	//针对发射惯性系的姿态角 todo
	double pitch = Eular_angle[0], yaw = Eular_angle[1], roll = Eular_angle[2];//俯仰角(pitch),滚动角(roll),偏航角(yaw)	 fai psi gama

	double wx = attitude_omega[0], wy = attitude_omega[1], wz = attitude_omega[2];
	double ax = acc[0], ay = acc[1], az = acc[2];//加速度 弹道坐标系 下
	double Fx = F[0],Fy = F[1],Fz = F[2];//动力、升力、侧力 弹道坐标系
	//-----弹体总体参数-----------
	double m;
	Vector3d inertia;
	double ix = inertia[0], iy = inertia[1], iz = inertia[2];
	double L_body = 6.58;
	double L_ref = 0.3261;
	double S_ref = 0.223961;
	double time=0.0;//计时器(程序时间，不是实际实际)
	double X_G=3.4;	//质心位置 todo
	//带默认值的默认构造函数,构造函数的参数是按顺序传递的，无法跳过前面的参数直接给后面的参数赋值

	missile(double m = 1501, double xg=3.4, Vector3d inertia= Vector3d(23.64544063, 5427.480754, 5427.480754), Vector3d pos = Vector3d::Zero(), \
		Vector3d vel = Vector3d::Zero(), Vector3d acc = Vector3d::Zero(), \
		Vector3d Eular_angle = Vector3d(Pi/2,0,0), double theta = Pi / 2, double sigma = 0, Vector3d Eular_omega = Vector3d::Zero(), \
		Vector3d attitude_omega = Vector3d::Zero(),Vector3d attitude_omega_acc = Vector3d::Zero(), double theta_dot = 0,double sigma_dot = 0,\
		double alpha=0, double  alpha_dot = 0, double  beta = 0, double  beta_dot = 0, double  nu = 0, double  nu_dot = 0);
	/*
	missile(Vector3d pos, Vector3d vel, Vector3d acc,Vector3d Eular_angle, Vector3d Eular_omega, \
		Vector3d attitude_omega, Vector3d attitude_omega_acc, double theta, double sigma, double theta_dot, double sigma_dot, double m);
	*/
	//拷贝
	missile(const missile& m);
	//赋值运算符重载,形参和返回的都是 对象的引用
	//missile operator=(missile& m);
	virtual ~missile() {};
	//void update_data();//todo

	double dt = 0.01;
	void dynamic(Vector3d F);//弹道坐标系
	void rotate(Vector3d M);//弹体坐标系  输入M 更新 w w_dot 返回 attitude_omega
	Vector3d get_Eular_angle();
	Vector3d get_w();
	Vector3d get_w_dot();
	Vector3d change_J(Vector3d I);
	Vector3d get_accelerate(double t);

	FlyStatus Extract_status();
};






#endif
