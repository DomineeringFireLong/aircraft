#include<iostream>
#include<cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include<fstream>
#include<random>
#include"missile.h"//����ŵ�ͬһ��·���� �����Ѳ���



using namespace Eigen;
using namespace std;
//missile::missile(double m , Vector3d pos , Vector3d vel , Vector3d acc ,Vector3d Eular_angle, double theta, double sigma,\
//	Vector3d Eular_omega = Vector3d::Zero(), Vector3d attitude_omega = Vector3d::Zero(), \
//	Vector3d attitude_omega_acc = Vector3d::Zero(), double theta_dot = 0, double sigma_dot = 0)

missile::missile(double m , double xg, Vector3d inertia, Vector3d pos, Vector3d vel , Vector3d acc , \
	Vector3d Eular_angle , double theta , double sigma , Vector3d Eular_omega , \
	Vector3d attitude_omega , Vector3d attitude_omega_acc, double theta_dot, double sigma_dot,\
	double alpha, double  alpha_dot, double  beta, double  beta_dot, double  nu, double  nu_dot)

{
	this->m = m;
	this->X_G = xg;
	this->inertia = inertia;
	this->pos = pos;
	this->vel = vel;
	this->acc = acc;

	this->Eular_angle = Eular_angle;//��̬��
	this->Eular_omega = Eular_omega;//��̬���ٶ�
	this->attitude_omega = attitude_omega;//�Ƶ���Ľ��ٶ�
	this->attitude_omega_acc = attitude_omega_acc;

	this->x = this->pos[0];
	this->y = this->pos[1];
	this->z = this->pos[2];

	this->ix = inertia[0];
	this->iy = inertia[1];
	this->iz = inertia[2];

	this->velocity = this->vel[0];
	this->ax = this->acc[0];
	this->ay = this->acc[1];
	this->az = this->acc[2];

	this->theta = theta;
	this->sigma= sigma;
	this->sigma_dot = sigma_dot;
	this->theta_dot = theta_dot;

	this->pitch = Eular_angle[0];//���� wz
	this->yaw = Eular_angle[1];//ƫ��  wy
	this->roll = Eular_angle[2];// wx

	this->alpha = alpha;
	this->alpha_dot = alpha_dot;
	this->beta = beta;
	this->beta_dot = beta_dot;
	this->nu = nu;
	this->nu_dot = nu_dot;
}
//�������캯�����ֳƸ��ƹ��캯��,�βα��������ã�����������Ϊconst
missile::missile(const missile& other) {
	this->pos = other.pos;
	this->vel = other.vel;
	this->acc = other.acc;

	this->Eular_angle = other.Eular_angle;
	this->Eular_omega = other.Eular_omega;
	this->attitude_omega = other.attitude_omega;
	this->attitude_omega_acc = other.attitude_omega_acc;

	this->x = this->pos[0];
	this->y = this->pos[1];
	this->z = this->pos[2];
	
	this->m = other.m;
	this->inertia = other.inertia;
	this->ix = other.inertia[0];
	this->iy = other.inertia[1];
	this->iz = other.inertia[2];
	this->X_G= other.X_G;

	this->velocity = this->vel[0];

	this->ax = this->acc[0];
	this->ay = this->acc[1];
	this->az = this->acc[2];

	this->theta = other.theta;
	this->sigma = other.sigma;
	this->sigma_dot = other.sigma_dot;
	this->theta_dot = other.theta_dot;

	this->pitch = Eular_angle[0];//���� wz
	this->yaw = Eular_angle[1];//ƫ��  wy
	this->roll = Eular_angle[2];// wx

	this->alpha = other.alpha;
	this->alpha_dot = other.alpha_dot;
	this->beta= other.beta;
	this->beta_dot = other.beta_dot;
	this->nu = other.nu;
	this->nu_dot = other.nu_dot;;
	cout << "�������캯��" << endl;//ֻҪ��ʼ��ʱ����ܸ�ֵ���������������޸� ֵ
}

//��ֵ���������,ע�ⲻ�ǹ��캯��������ͨ���������
//missile missile::operator=(missile& const other) {
//
//	missile m1(other);
//
//	return m1;//������ʱ���� 
//}


FlyStatus missile::Extract_status() {

	FlyStatus current_status = { m,velocity,vel,Eular_angle,attitude_omega,pos,alpha,beta,theta,sigma,nu,S_ref,L_ref,inertia,time};
	return current_status;

}
//------------------------------------------���Ķ���ѧ----------------------------------------//
void missile::dynamic(Vector3d F) {
	velocity = vel[0];//��������ϵ
	acc = F/ m;
	theta_dot = F[1] / (m * velocity);////�������ʱΪȡ ���� ��double�;���С�������� %ȡ��
	sigma_dot = -1 * F[2] / (m * velocity * cos(theta));
	
	//ע�������任,���ⲿ��ı�
	//ʵ���ǿ��Ʒ�����أ��ı��ٶȷ����ƫ��	//���ĳ�������任�ǳ��󣬿����ǻ���ʱ������û��dt
	vel += Vector3d(F[0] / m,0, 0) * dt;
	theta += theta_dot * dt;
	sigma += sigma_dot * dt;
	velocity = vel[0];
	//cout << vel << endl;
	//vel���ٶ�����ϵ��ת������������ϵ�����ڵ�������ϵ��xyz���㣬������Ŀ��Ե�λ�ø���
	x += velocity * cos(theta) * cos(sigma) * dt;
	y += velocity * sin(theta) * dt;
	z += -1 * velocity * cos(theta) * sin(sigma) * dt;
	//cout << "theta" << theta << endl;
	//cout <<"���Ƹ߶�" << y << endl;
	pos = Vector3d(x, y, z);


}

//---------------------------------------������ת������ѧ---------------------------------------//
//������ת�����˶�ѧģ��    ����M�����W
void missile::rotate(Vector3d M)//Vector3d  ����M ���� w w_dot  ֻ�ܿ��ƹ��ǡ��໬�ǡ��ٶȹ�ת��
{
	wx = this->attitude_omega[0];
	wy = this->attitude_omega[1];
	wz = this->attitude_omega[2];
	attitude_omega_acc[0] = (M[0] - (this->iz - this->iy) * this->wz * this->wy) / this->ix;
	attitude_omega_acc[1] = (M[1] - (this->ix - this->iz) * this->wx * this->wz) / this->iy;
	attitude_omega_acc[2] = (M[2] - (this->iy - this->ix) * this->wy * this->wx) / this->iz;
	attitude_omega += this->attitude_omega_acc * dt;//���ٶ�
	//���Ƕ�Ӧ��ϵ,�Ҵ˴�����Ϊ��̬�ǣ������ǹ��ǣ��໬�ǣ��ٶ����ǣ�Ӧ������̬��
	//this->r += this->attitude_omega[0] * dt;//�Ƕ�
	//this->b += this->attitude_omega[1] * dt;
	//this->a += this->attitude_omega[2] * dt;
	wx = this->attitude_omega[0];
	wy = this->attitude_omega[1];
	wz = this->attitude_omega[2];
	

	//�ٶȵ�����
	alpha_dot = wz - tan(beta) * (wx * cos(alpha) - wy * sin(alpha));
	beta_dot = wx * sin(alpha) + wy * cos(alpha);
	//�������ٶ�
	nu = (wx * cos(alpha) - wy * sin(alpha)) / cos(beta);
	//�˴��Ƕ�Ӧ��ϵ
	alpha += alpha_dot * dt;//�Ƕ�
	beta += beta_dot * dt;
	nu += nu_dot * dt;

	//�״�
	//this->roll += this->attitude_omega[0] * dt;//�Ƕ�
	//this->pothi += this->attitude_omega[1] * dt;
	//this->sigma += this->attitude_omega[2] * dt;
	//���� ���ٶ��� ��̬�ǵĵ��� ��ϵ ����ֱ�ӱ任 
	//���ܳ��ָ���90��ʱ ����:�ı�˳��  todo
	//��ֱ��theta=90
	//����תƫ�����ٸ����͹�ת
	//Eular_angle[0] += (sin(roll) * wy + cos(roll) * wz) * dt;
	//Eular_angle[1] += (cos(roll) * wy  - sin(roll) * wz) /cos(theta) * dt;
	//Eular_angle[2] += (wx+tan(theta)*( -cos(roll) * wy + sin(roll) * wz)) * dt;
	
	//����ת��������ƫ���͹�ת
	Eular_angle[0] += (sin(roll) * wy + cos(roll) * wz) / cos(sigma) * dt;
	Eular_angle[1] += (cos(roll) * wy - sin(roll) * wz) * dt;
	Eular_angle[2] += (wx + tan(sigma) * (sin(roll) * wy + cos(roll) * wz)) * dt;

	attitude_omega<< wx, wy, wz; //+= this->attitude_omega * dt;
	//�������Ƕ�ͨ��΢�ַ��̺ͼ��ι�ϵ���������ں�
	Vector3d abg_geometry =calAlphaBetanu(Eular_angle, theta, sigma);
	alpha = (alpha + abg_geometry[0]) / 2;
	beta= (beta + abg_geometry[1]) / 2;
	nu= (nu + abg_geometry[2]) / 2;
	//alpha =abg_geometry[0];
	//beta =abg_geometry[1];
	//nu =abg_geometry[2];


	//return Vector3d(alpha,beta,nu);// attitude_omega;
}

Vector3d missile::get_Eular_angle()
{
	return Eular_angle;
}


Vector3d missile::get_w()
{
	return this->attitude_omega;
}

Vector3d missile::get_w_dot()
{
	return this->attitude_omega_acc;
}
Vector3d missile::change_J(Vector3d I)
{
	ix = I[0];
	iy = I[1];
	iz = I[2];
	inertia = I;
	Vector3d J(ix, iy, iz);
	return J;
}
Vector3d missile::get_accelerate(double t)
{
	return acc;

}


//std::variant �����庯���ķ������� ��������һ�������ڷ��ز�ͬ���͵�ֵ


/*
void control_test()
{

	Vector3d z;
	z << 0, 0, 0;
	Vector3d Omega__expect;
	Omega__expect << 0.5, 0.5, 0;
	missile_rotate MR(0, 0, 0, 0.01, 0.01, 0.01, z, z, 1, 1, 1, Omega__expect);//default_Omega_expect
	Vector3d angle = MR.get_angle();
	Vector3d w = MR.get_w();
	Vector3d J = MR.get_J();
	Vector3d ae = MR.qd_control.Angle_expectation;
	int count = 0;
	double temp;
	while ((abs(angle[0] - ae[0]) > 0.1) || (abs(angle[1] - ae[1]) > 0.1) || (abs(angle[2] - ae[2]) < 0.1))//mo(w - MR.qd_control.Omg_expectation) > 0.5
	{
		count += 1;
		std::cout << "--����error =  " << (angle - ae) << std::endl;
		//std::cout << mo(MR.get_angle() - MR.qd_control.Omg_expectation) << std::endl;
		//�Ӷ���ѧģ���л�ò���		//�Ѳ������������ģ�飬������������
		Vector3d M = MR.qd_control.control_instruct(angle[0], angle[1], angle[2], w[0], w[1], w[2], J[0], J[1], J[2], ae);
		//���ظ�������ѧģ�飬����״̬�����ĸ���

		temp = angle[0];
		MR.YD(M);
		w = MR.get_w();
		angle = MR.get_angle();//���²���
		if (abs(temp - ae[0]) < (abs(angle[0] - ae[0])))
			break;

		if (count > 500)
		{
			std::cout << "-------���Ƴ�ʱ��5s����----------" << std::endl;
			break;
		}
	}
	std::cout << "\n----------����ʱ��----------:" << count * 0.01 << "------------\n" << std::endl;

	angle = MR.get_angle();

	std::cout << "----------���ƽǶ�error---------:" << (angle - ae).transpose() * 180 / Pi << "------------\n" << std::endl;
	std::cout << "---------��̬��Ϊ--------:" << MR.sigma * 180 / Pi << " " << MR.pothi * 180 / Pi << " " << MR.roll * 180 / Pi << "------------\n" << std::endl;

}


*/