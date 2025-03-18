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
//�����ѭ������
//#include"controller.h"
//#include"guide.h"

//Ϊ���ú��������ڸ���.cpp�й�����ȷ������������.h��ֻ����������������һ��.cpp��ʵ����������������Ͳ����ͻ�ˡ�
using namespace Eigen;
using namespace std;
////���˶�ģ��Ƕ�뵽��������,�˶�ѧ���̣����ģ�������ת����

class missile
{
public:
	//�����˶�ѧ���������˶�ѧ
	Eigen::Vector3d pos,vel,acc; //λ�á��ٶȡ����ٶ� ��������ϵ  vel���ٶ�����ϵ stt �������ٶ���ͬ
	Eigen::Vector3d Eular_angle, Eular_omega;//������ν�ġ���̬�ǡ���ʵ��ŷ����(�Ǹ��������������ǶȵĶ��嶼����ͬһ������ϵ��)�����ٶȣ��Ǽ��ٶ�
	Eigen::Vector3d attitude_omega, attitude_omega_acc;//��������̬���ٶȣ��Ƶ���Ľ��ٶȣ� //����ֱ������
	//���Ķ���ѧ�������Ķ���ѧ
	Eigen::Vector3d F,P,G,A= Vector3d::Zero();//��������������������������(���� ���ڵ���ϵ)  4����ط����������������У�
	Eigen::Vector3d M, M_p, M_a;//�����ء�������(������)�������������������� (���� ���ڵ���ϵ)  6���˿ط������������������У�
	
	//�Ƕ����ɶ�
	double theta, sigma, theta_dot, sigma_dot;//����������ƫ�� ����->����  theta sigma
	double alpha, beta, nu;//�ٶ�->����  �ٶ�->���� ����->�ٶ�
	double alpha_dot, beta_dot, nu_dot;
	
	//��άȡ����һά
	double velocity = this->vel[0]; // ����ϵ  ����Ϊ����ϵsqrt(pow(this->vel[0], 2) + pow(this->vel[1], 2) + pow(this->vel[2], 2));
	double x = pos[0], y = pos[1], z = pos[2];
	//��Է������ϵ����̬�� todo
	double pitch = Eular_angle[0], yaw = Eular_angle[1], roll = Eular_angle[2];//������(pitch),������(roll),ƫ����(yaw)	 fai psi gama

	double wx = attitude_omega[0], wy = attitude_omega[1], wz = attitude_omega[2];
	double ax = acc[0], ay = acc[1], az = acc[2];//���ٶ� ��������ϵ ��
	double Fx = F[0],Fy = F[1],Fz = F[2];//���������������� ��������ϵ
	//-----�����������-----------
	double m;
	Vector3d inertia;
	double ix = inertia[0], iy = inertia[1], iz = inertia[2];
	double L_body = 6.58;
	double L_ref = 0.3261;
	double S_ref = 0.223961;
	double time=0.0;//��ʱ��(����ʱ�䣬����ʵ��ʵ��)
	double X_G=3.4;	//����λ�� todo
	//��Ĭ��ֵ��Ĭ�Ϲ��캯��,���캯���Ĳ����ǰ�˳�򴫵ݵģ��޷�����ǰ��Ĳ���ֱ�Ӹ�����Ĳ�����ֵ

	missile(double m = 1501, double xg=3.4, Vector3d inertia= Vector3d(23.64544063, 5427.480754, 5427.480754), Vector3d pos = Vector3d::Zero(), \
		Vector3d vel = Vector3d::Zero(), Vector3d acc = Vector3d::Zero(), \
		Vector3d Eular_angle = Vector3d(Pi/2,0,0), double theta = Pi / 2, double sigma = 0, Vector3d Eular_omega = Vector3d::Zero(), \
		Vector3d attitude_omega = Vector3d::Zero(),Vector3d attitude_omega_acc = Vector3d::Zero(), double theta_dot = 0,double sigma_dot = 0,\
		double alpha=0, double  alpha_dot = 0, double  beta = 0, double  beta_dot = 0, double  nu = 0, double  nu_dot = 0);
	/*
	missile(Vector3d pos, Vector3d vel, Vector3d acc,Vector3d Eular_angle, Vector3d Eular_omega, \
		Vector3d attitude_omega, Vector3d attitude_omega_acc, double theta, double sigma, double theta_dot, double sigma_dot, double m);
	*/
	//����
	missile(const missile& m);
	//��ֵ���������,�βκͷ��صĶ��� ���������
	//missile operator=(missile& m);
	virtual ~missile() {};
	//void update_data();//todo

	double dt = 0.01;
	void dynamic(Vector3d F);//��������ϵ
	void rotate(Vector3d M);//��������ϵ  ����M ���� w w_dot ���� attitude_omega
	Vector3d get_Eular_angle();
	Vector3d get_w();
	Vector3d get_w_dot();
	Vector3d change_J(Vector3d I);
	Vector3d get_accelerate(double t);

	FlyStatus Extract_status();
};






#endif
