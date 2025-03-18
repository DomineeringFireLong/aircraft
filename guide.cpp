#include<iostream>
#include<cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include<graphics.h>
#include<fstream>
#include<random>
#include"guide.h"

using namespace std;
using namespace Eigen;

Guidance_law::Guidance_law(Vector3d point,double t_start1, double t_final1):SM_3()
{
	
	expect_point = point;
	t_start = t_start1;
	t_final = t_final1;
	time_work = 0.0;
	r = sqrt(pow(expect_point[0] - SM_3.pos[0], 2) + pow(expect_point[1] - SM_3.pos[1], 2) + pow(expect_point[2] - SM_3.pos[2], 2));

	theta = atan((expect_point[1] - SM_3.pos[1]) / sqrt(pow(expect_point[0] - SM_3.pos[0], 2) + pow(expect_point[2] - SM_3.pos[2], 2)));
	sigma = atan(-(expect_point[2] - SM_3.pos[2]) / expect_point[0] - SM_3.pos[0]);
	theta_dot = 0;
	sigma_dot = 0;
	ZCMD = Vector3d(5000, 10000, 0);
	R = expect_point - SM_3.pos;
	//R[0] = expect_point[0] - SM_3.pos[0];
	//R[1] = expect_point[1] - SM_3.pos[1];
	//R[2] = expect_point[2] - SM_3.pos[2];

	r_dot = G2T(T2G(SM_3.vel, SM_3.theta, SM_3.sigma), theta, sigma)[0];
	R_dot = G2T(T2G(SM_3.vel, SM_3.theta, SM_3.sigma), theta, sigma);

}

Vector3d Guidance_law::powered_stage(double Q1, double  now_time) {
	//����LQR���Զ��������ţ�Ԥ���Ƶ�?
	Q = Q1;
	//SM_3 = m1;
	time_work = now_time;
	//����LOS
	theta = atan((expect_point[1] - SM_3.pos[1]) / sqrt(pow(expect_point[0] - SM_3.pos[0], 2) + pow(expect_point[2] - SM_3.pos[2], 2)));
	sigma = atan((expect_point[2] - SM_3.pos[2]) / expect_point[0] - SM_3.pos[0]);//atan �����-pi/2  0;// 
	if (sigma + PI / 2 < 0.01)
		sigma = 0;
	//�����Ƶ���������λ�������ǻ���los������������Ҫ���� los���ٶ�

	r_dot = -1*G2T(T2G(SM_3.vel, SM_3.theta, SM_3.sigma), theta, sigma)[0];
	R = expect_point - SM_3.pos;
	r = modulo(R);
	R_dot = -1*G2T(T2G(SM_3.vel, SM_3.theta, SM_3.sigma), theta, sigma);//�����ٶ�ת��losϵ
	theta_dot = - (SM_3.velocity * cos(theta) * sin(SM_3.theta) - SM_3.velocity * sin(theta) * cos(SM_3.theta) * cos((sigma- SM_3.sigma))) / r;
	sigma_dot = (SM_3.velocity * cos(SM_3.theta) * sin(sigma - SM_3.sigma)) / (r * cos(theta));

	////����Ѱ���	ZEM

	////�˴�����time_work����0��������һ��������ֻ����ǰһ���๤��0.01s
	std::function<Vector3d(double)> ptr = std::bind(&missile::get_accelerate,&SM_3,std::placeholders::_1);//Vector3d (Guidance_law::*ptr)(double t);//��ͨ����ָ��:��������(*ָ�������)(������);��ָ�����ָ������������ͬ�ĺ���
	//ZCMD = R+R_dot*(t_final - time_work)- VL_Integrataion_2D(time_work, t_final, ptr, 0.01);//��Ա����ָ�� ����� �β�ָ������ͬһ����
	//expect_accelerate =((t_final - time_work) /(1+ pow((t_final - time_work), 3) / 3)) * ZCMD;//���Ƶ��ɣ�����Ч����һ��


	//����Ч������ ����qС�ǶȲų���
	//expect_accelerate = Vector3d(0, 12 * SM_3.velocity * theta_dot + 6 * SM_3.velocity * (SM_3.theta - 60/180*PI) / (r / r_dot),0);
	//	expect_accelerate = Vector3d(0, -3 * r_dot* theta_dot , 3 * r_dot * sigma_dot);//���ǵ�ת��ƽ�����z����̫��
	
	//�������Ƶ��ɣ������������ȫ�ȼ�,����Ч����
	//ZCMD = R + R_dot * (t_final - time_work) - VL_Integrataion_2D(time_work, t_final, ptr, 0.01);
	//double rtf = modulo(ZCMD);
	//expect_accelerate = Vector3d(0, 3 * pow(r, 3) / (pow(r, 3) - pow(rtf, 3)) * r_dot * theta_dot, 0);
	
	//��������
	//expect_accelerate = Vector3d(0,  - 4 * r_dot * theta_dot, 0);//3*r_dot*sigma_dot*cos(theta)
	

	//�����ά�����Ƶ���,������

	expect_accelerate[0] = 0;
	expect_accelerate[1] = -2 * r_dot * theta_dot - 1 * r * pow(sigma_dot, 2) * sin(theta) * cos(theta) + 6 * (theta - 60 / 180 * PI) / pow((44-now_time),2)+4*r* theta_dot/(44 - now_time);
	expect_accelerate[2] = 0;

	//-2*r_dot*sigma_dot*cos(theta)+2*r*theta_dot*sigma_dot*sin(theta)+r*cos(theta)*6*(sigma-0)/ pow((44 - now_time), 2)+ 4*r*cos(theta) * sigma_dot / (44 - now_time);
	//expect_accelerate[0] = 0;
	//expect_accelerate[1] = -1 * r * pow(sigma_dot, 2) * sin(theta) * cos(theta) + 6 * pow(r_dot, 2) * (theta - 60 / 180 * PI) / r - 6 * theta_dot * r_dot;
	//expect_accelerate[2] = 0;

	//2*r*theta_dot*sigma_dot*sin(theta)+6*pow(r_dot,2)*(sigma)*cos(theta)/r-6*sigma_dot*r_dot*cos(theta);


//	double(Guidance_law:: *p)() = get_theta_dot;
	//std::function<double()> p = get_theta_dot;


	//expect_accelerate[0] = 0;
	//expect_accelerate[1] = -5* theta_dot/ (44 - now_time) - 10 *theta / pow((44 - now_time), 2) -10 * Integration(0, now_time,p,0.01) / pow((44-now_time),3);
	//expect_accelerate[2] = 0;


	expect_accelerate = sature(expect_accelerate, 300);
	
	
	return expect_accelerate;
}

//���λ���

Vector3d Guidance_law::VL_Integrataion_2D(double x0, double x1, std::function<Vector3d(double)> f, double interval)// Vector3d(Guidance_law::*f)(double t)
{
	double lenx = x1 - x0;
	unsigned long Nx = lenx / interval;
	Vector3d result = Vector3d(0, 0, 0);
	for (unsigned long i = 0; i < Nx; ++i)
	{
		result += Integration_G(x0, x0 + i * interval, f, interval) * interval;//ֱ�Ӹ���ַ����
	}
	return result;
}

Vector3d Guidance_law::Integration_G(double a, double b, std::function<Vector3d(double)> f, double interval)
{
	double len = b - a;
	int N = len / interval;
	Vector3d result{ 0.0,0.0,0.0 };
	//for (unsigned long i = 0; i < N; ++i) {
	//	result += f(a + i * interval) * interval;
	//}
	result = f(a) * len;
	return result;
}


double Guidance_law::get_theta_dot(){
	return theta_dot;
}

void Guidance_law::terminal_phase_guidance(missile m1, missile target1)//ֱ���к��в����ĵ�����Ŀ�����й������
{
	//ֻ�ǳ�ʼ�������ܸ��£�vel���ڵ�������ϵ���׸��£������ٶ�����ϵ
	this->SM_3 = m1;
	this->target = target1;
	this->r_dot = 0;
	this->theta_dot = 0;
	this->sigma_dot = 0;
	this->r = sqrt(pow(target.pos[0] - SM_3 .pos[0], 2) + pow(target.pos[1] - SM_3 .pos[1], 2) + pow(target.pos[2] - SM_3 .pos[2], 2));
	//1.ע�⣺��������ϵ���ڵ����ƫ�ǡ���� ���㹫ʽ,���ݳ�ʼ�ٶ�ʸ�������ʼ�ǣ���Ŀ����ϵ�Ƕ������λ�Ʒ��򣬵����ĽǶ����ٶȷ���Ķ���
	this->theta = atan((this->target.y - this->SM_3 .y) / sqrt(pow(this->target.x - this->SM_3 .x, 2) + pow(this->target.z - this->SM_3 .z, 2)));
	this->sigma = atan(-(this->target.z) / this->target.x);
	
}

void Guidance_law::terminal_phase_guidance(missile m1, Eigen::Vector3d target)
{
	//ֻ�ǳ�ʼ�������ܸ��£�vel���ڵ�������ϵ���׸��£������ٶ�����ϵ
	this->SM_3 = m1;
	this->expect_point = target;
	this->r_dot = 0;
	this->theta_dot = 0;
	this->sigma_dot = 0;
	this->r = sqrt(pow(expect_point[0] - expect_point[0], 2) + pow(expect_point[1] - SM_3.pos[1], 2) + pow(expect_point[2] - SM_3.pos[2], 2));
	//1.ע�⣺��������ϵ���ڵ����ƫ�ǡ���� ���㹫ʽ,���ݳ�ʼ�ٶ�ʸ�������ʼ�ǣ���Ŀ����ϵ�Ƕ������λ�Ʒ��򣬵����ĽǶ����ٶȷ���Ķ���
	this->theta = atan((this->expect_point[1] - this->SM_3.y) / sqrt (pow(this->expect_point[0] - this->SM_3.x, 2) + pow(this->expect_point[2] - this->SM_3.z, 2)));
	this->sigma = atan(-(this->expect_point[2]) / this->expect_point[0]);

}


//֮ǰд���Ƶ��͵�������ѧ/�˶�ѧģ���ں���һ���ˣ���Ҫ�ֿ�

/*
missile Guider::proportional_guide(int Ky = 4, int Kz = 4, Data_Buffer dandao_data)
{
		//*******************�������****************
	//ofstream outfile;
	//outfile.open("./3D_proportional_guide.txt");
	string r1 = "./3D_proportional_guide.txt";
	//outfile << "����λ��X Y Z �ٶ�v  �������theta ����ƫ��sigma  " << "Ŀ��λ��X Y Z �ٶ�v �������theta ����ƫ��sigma  " << endl; //�Ƕ�theta eta
	string r2 = "����λ��X Y Z �ٶ�v  �������theta ����ƫ��sigma Ŀ��λ��X Y Z �ٶ�v �������theta ����ƫ��sigma ";
	//*****************�Ƶ�����ѭ��******************
	//		//�������ʴ�С����
	//this->SM_3.vel = Vector3d(this->SM_3.velocity, 0.0, 0.0);
	//this->expect_point.vel = Vector3d(this->expect_point.velocity, 0.0, 0.0);
	double dt = 0.01;
	double count = 0;
	while (this->r > 10)
	{
		count += dt;
		outfile << this->SM_3.x << " " << this->SM_3.y << " " << this->SM_3.z << " " << this->SM_3.vel[0] << " " << this->SM_3.theta * 180 / 3.1415926 << " " << this->SM_3.sigma * 180 / 3.1415926
			<< " " << this->expect_point.x << " " << this->expect_point.y << " " << this->expect_point.z << " " << this->expect_point.vel[0] << " " << this->expect_point.theta * 180 / 3.1415926 << " " << this->expect_point.sigma * 180 / 3.1415926 << endl;

		//����Ƕ�ʱ������ת�Ƕȣ�this->SM_3.theta * 180 / 3.1415926 << this->SM_3.sigma * 180 / 3.1415926
		//*******������Ŀ��λ��״̬��ƫ�ǡ���ǡ�����*************
		
		this->r_dot = (G2T(T2G(this->expect_point.vel, this->expect_point.theta, this->expect_point.sigma), this->theta, this->sigma) - G2T(T2G(this->SM_3.vel, this->SM_3.theta, this->SM_3.sigma), this->theta, this->sigma))[0];
		this->theta_dot = (G2T(T2G(this->expect_point.vel, this->expect_point.theta, this->expect_point.sigma), this->theta, this->sigma) - G2T(T2G(this->SM_3.vel, this->SM_3.theta, this->SM_3.sigma), this->theta, this->sigma))[1] / this->r;
		this->sigma_dot = (-1 * (G2T(T2G(this->expect_point.vel, this->expect_point.theta, this->expect_point.sigma), this->theta, this->sigma)[2]) + (G2T(T2G(this->SM_3.vel, this->SM_3.theta, this->SM_3.sigma), this->theta, this->sigma))[2]) / (this->r * cos(this->theta));

		//������������
		this->r = sqrt(pow(this->expect_point.x - this->SM_3.x, 2) + pow(this->expect_point.y - this->SM_3.y, 2) + pow(this->expect_point.z - this->SM_3.z, 2));// 
		//std::cout << this->r << endl;
		//this->r += this->r_dot*dt;
		this->sigma += this->sigma_dot * dt;
		this->theta += this->theta_dot * dt;
		//ֻ������������"�����Ƕȱ仯��"��ֻ�ʺ��㷨��������û��"��������������"��ʵ�ʵĸı��ٶȷ���
		//(�����Ϣ-(�Ƶ���)->�����Ƕ�->�������->���õ�������-(������)->��ƫ-(����ѧģ�͡��˶�ѧ��ʽ)->ʵ������ʵ�ʹ���->λ�ñ任->�����Ϣ
		this->SM_3.ay = Ky * -1 * (this->r_dot) * this->theta_dot;
		this->SM_3.az = -Kz * -1 * (this->r_dot) * this->sigma_dot * cos(this->theta);

		//this->SM_3.ay = Ky * (this->r_dot * this->theta_dot+this->r+theta_dot2);
		//this->SM_3.az = -Kz * abs(this->r_dot) * this->sigma_dot * cos(this->theta);
		//***������***
		this->SM_3.Fx = 0;
		this->SM_3.Fy = this->SM_3.m * this->SM_3.ay;
		this->SM_3.Fz = this->SM_3.m * this->SM_3.az;
		//cout << this->SM_3.ay << "    " << this->SM_3.az << endl;//ay���ٶȴ���� 40   az -40
		this->SM_3.acc = Vector3d(this->SM_3.Fx / this->SM_3.m, 0, 0);
		this->SM_3.theta_dot = this->SM_3.Fy / (this->SM_3.m * this->SM_3.velocity);
		this->SM_3.sigma_dot = -1 * this->SM_3.Fz / (this->SM_3.m * this->SM_3.velocity * cos(this->SM_3.theta));

		//ʵ���ǿ��Ʒ�����أ��ı��ٶȷ����ƫ��	//���ĳ�������任�ǳ��󣬿����ǻ���ʱ������û��dt
		this->SM_3.vel += this->SM_3.acc * dt;
		this->SM_3.theta += this->SM_3.theta_dot * dt;
		this->SM_3.sigma += this->SM_3.sigma_dot * dt;
		this->SM_3.velocity = this->SM_3.vel[0];
		//cout << this->SM_3.vel << endl;
		//vel���ٶ�����ϵ��ת������������ϵ�����ڵ�������ϵ��xyz���㣬������Ŀ��Ե�λ�ø���
		this->SM_3.x += this->SM_3.velocity * cos(this->SM_3.theta) * cos(this->SM_3.sigma) * dt;
		this->SM_3.y += this->SM_3.velocity * sin(this->SM_3.theta) * dt;
		this->SM_3.z += -1 * this->SM_3.velocity * cos(this->SM_3.theta) * sin(this->SM_3.sigma) * dt;
		this->expect_point.x += this->expect_point.velocity * cos(this->expect_point.theta) * cos(this->expect_point.sigma) * dt;
		this->expect_point.y += this->expect_point.velocity * sin(this->expect_point.theta) * dt;
		this->expect_point.z += -1 * this->expect_point.velocity * cos(this->expect_point.theta) * sin(this->expect_point.sigma) * dt;
		this->SM_3.pos = Vector3d(this->SM_3.x, this->SM_3.y, this->SM_3.z);
		this->expect_point.pos = Vector3d(this->expect_point.x, this->expect_point.y, this->expect_point.z);
		//xyzһֱ�ڸ��£�����û�и�pos��ֵ������pos���ǳ�ֵ
		//c++��ifֻ��һ����䣬û{}��ֱ���˳�
		if (this->r_dot > 0 || count >= 60)//|| count >= 60
		{
			std::cout << "δ�ܸ���" << endl;
			break;
		}
	}
	outfile.close();

	std::cout << this->expect_point.pos[0] - this->SM_3.pos[0] << std::endl;
	std::cout << this->expect_point.pos[1] - this->SM_3.pos[1] << std::endl;
	std::cout << this->expect_point.pos[2] - this->SM_3.pos[2] << std::endl;
	std::cout << "�Ѱ���:" << sqrt(pow(this->expect_point.pos[0] - this->SM_3.pos[0], 2) + pow(this->expect_point.pos[1] - this->SM_3.pos[1], 2) + pow(this->expect_point.pos[2] - this->SM_3.pos[2], 2));
	return SM_3;
}

*/


/*
missile terminal_phase_guidance::slidemode(int k1 = 6, int k2 = 0.1, int k3 = 6, int k4 = 0.1, double theta_e = PI / 6, double eta_e = PI / 6, string root = "./3D_slidemode_guide.txt" )//��ģ���ܲ�����Ӱ�죬���·�ɢ,һ��������һ�������� �Ϳ��ܵ��������ͷ�ɢ
{//���ɷ��ʱ�������Ҫ��public
	
	string root = "./3D_slidemode_guide.txt";
	//------------------------�������--------------
	ofstream outfile;  outfile.open("./3D_slidemode_guide.txt");
	
	outfile << "����λ��X Y Z �ٶ�v  �������theta ����ƫ��sigma  " << "Ŀ��λ��X Y Z �ٶ�v �������theta ����ƫ��sigma  " << endl; //�Ƕ�theta eta
	ofstream accfile;
	//accfile.open("./3D_acc_guide.txt");

	//*****************�Ƶ�����ѭ��******************
	//		//�������ʴ�С����
	//this->SM_3.vel = Vector3d(this->SM_3.velocity, 0.0, 0.0);
	//this->expect_point.vel = Vector3d(this->expect_point.velocity, 0.0, 0.0);
	double dt = 0.01;
	double count = 0;
	double s1 = 0;
	double s2 = 0;
	double s1_dot = 0;
	double s2_dot = 0;
	double ay_dm = 0;
	double az_dm = 0;
	double ep1 = 0.001;
	double ep2 = 0.001;
	bool flag = 1;
	double temp = this->r;
	while (this->r > 4)
	{
		count += dt;	//����Ƕ�ʱ������ת�Ƕȣ�this->SM_3.theta * 180 / 3.1415926 << this->SM_3.sigma * 180 / 3.1415926
		outfile << this->SM_3.x << " " << this->SM_3.y << " " << this->SM_3.z << " " << this->SM_3.vel[0] << " " << this->SM_3.theta * 180 / 3.1415926 << " " << this->SM_3.sigma * 180 / 3.1415926
			<< " " << this->expect_point.x << " " << this->expect_point.y << " " << this->expect_point.z << " " << this->expect_point.vel[0] << " " << this->expect_point.theta * 180 / 3.1415926 << " " << this->expect_point.sigma * 180 / 3.1415926 << endl;
		//outfile << this->SM_3.x << " " << this->SM_3.y << " " << this->SM_3.z << "   " << this->expect_point.x << " " << this->expect_point.y << " " << this->expect_point.z << " " << endl;
		//����˶�ѧģ�ͣ�r��theta��eta������ M��T��V��theta��eta����
		this->r_dot = (this->expect_point.velocity * cos(this->expect_point.theta) * cos(this->theta) * cos(this->expect_point.sigma - this->sigma) + this->expect_point.velocity * sin(this->expect_point.theta) * sin(this->theta)) - (this->SM_3.velocity * cos(this->SM_3.theta) * cos(this->theta) * cos(this->SM_3.sigma - this->sigma) + this->SM_3.velocity * sin(this->SM_3.theta) * sin(this->theta));
		this->theta_dot = ((this->expect_point.velocity * cos(this->theta) * sin(this->expect_point.theta) - this->expect_point.velocity * sin(this->theta) * cos(this->expect_point.theta) * cos((this->expect_point.sigma - this->sigma))) - (this->SM_3.velocity * cos(this->theta) * sin(this->SM_3.theta) - this->SM_3.velocity * sin(this->theta) * cos(this->SM_3.theta) * cos((this->SM_3.sigma - this->sigma)))) / this->r;
		this->sigma_dot = (this->SM_3.velocity * cos(this->SM_3.theta) * sin(this->sigma - this->SM_3.sigma) - this->expect_point.velocity * cos(this->expect_point.theta) * sin(this->sigma - this->expect_point.sigma)) / (this->r * cos(this->theta));
		//������������
		this->r = sqrt(pow(this->expect_point.x - this->SM_3.x, 2) + pow(this->expect_point.y - this->SM_3.y, 2) + pow(this->expect_point.z - this->SM_3.z, 2));
		this->sigma += this->sigma_dot * dt;
		this->theta += this->theta_dot * dt;
		//���������ɣ�����ay��az  
		//epsilion = y0-(this->expect_point.y - this->SM_3.y);
		//t_go = this->r / this->r_dot + 0.0000001;
		s1 = this->theta_dot + k1 * (this->theta - theta_e);
		s2 = this->sigma_dot + k3 * (this->sigma - eta_e);
		ay_dm = sat(-2 * this->r_dot * this->theta_dot - this->r * cos(this->theta) * sin(this->theta) * this->sigma_dot * this->sigma_dot + this->r * k1 * this->theta_dot + this->r * ep1 * sign(s1) + this->r * k2 * s1, 250.0);//ģ�庯�����Զ�������ע������������ģ���������
		az_dm = sat(2 * this->r_dot * this->sigma_dot * cos(this->theta) - 2 * this->r * this->sigma_dot * this->theta_dot * sin(this->theta) - this->r * cos(this->theta) * k3 * this->sigma_dot - this->r * cos(this->theta) * (ep2 * sign(s2) + k4 * s2), 250.0);
		//ay_dm = -2 * this->r_dot * this->theta_dot - this->r * cos(this->theta) * sin(this->theta) * this->sigma_dot * this->sigma_dot + this->r * k1 * this->theta_dot + this->r * ep1 * sign(s1) + this->r * k2 * s1;
		//az_dm = 2 * this->r_dot * this->sigma_dot * cos(this->theta) - 2 * this->r * this->sigma_dot * this->theta_dot * sin(this->theta) - this->r * cos(this->theta) * k3 * this->sigma_dot - this->r * cos(this->theta) * (ep2 * sign(s2) + k4 * s2);

		//ay_dm = -this->r * this->sigma_dot * this->sigma_dot * sin(this->theta) * cos(this->theta) + 6 * (this->r_dot * this->r_dot * this->theta) / this->r - 6 * this->r_dot * this->theta_dot;
		//az_dm = 2 * this->r * this->theta_dot * this->sigma_dot * sin(this->theta) + 6 * (this->r_dot * this->r_dot * this->sigma * cos(this->theta)) / this->r - 6 * this->sigma_dot * this->r_dot * cos(this->theta_dot);
		this->SM_3.acc = G2T(T2G(Vector3d(0, ay_dm, az_dm), this->theta, this->sigma), this->SM_3.theta, this->SM_3.sigma);
		accfile << this->SM_3.acc[0] << " " << this->SM_3.acc[1] << " " << this->SM_3.acc[2] << endl;
		//ay az����Թ�������ϵ�ļ��ٶ�   ���ǵ�Ŀ��������ϵ �ǵ���ϵ
		//this->SM_3.acc = G2T(Vector3d(0, ay_dm, az_dm),this->SM_3.theta,this->SM_3.sigma);//ay az����ڵ��� �ǶԵ�
		//cout << this->SM_3.acc[0]<< "  "<<this->SM_3.acc[1] << "  " << this->SM_3.acc[2] << endl;
		//this->SM_3.acc =Vector3d(0, ay_dm, az_dm);
		this->SM_3.F = this->SM_3.m * this->SM_3.acc;
		this->SM_3.Fx = 0;// this->SM_3.F[0];//ʵ�ʵ���û�취���ٶ�
		this->SM_3.Fy = this->SM_3.F[1];
		this->SM_3.Fz = this->SM_3.F[2];
		//cout << this->SM_3.acc << "    " << endl;
		//ay��az->Fy Fz
		//this->SM_3.Fy = this->SM_3.m * this->SM_3.ay;// ;ay_dm
		//this->SM_3.Fz = this->SM_3.m * this->SM_3.az;// ;az_dm
		//�������Ķ���ѧģ��
		this->SM_3.theta_dot = this->SM_3.Fy / (this->SM_3.m * this->SM_3.velocity);
		this->SM_3.sigma_dot = -1 * this->SM_3.Fz / (this->SM_3.m * this->SM_3.velocity * cos(this->SM_3.theta));
		//this->SM_3.vel+=Vector3d(this->SM_3.Fx / this->SM_3.m, 0, 0);
		//�ٰѼ��ٶ�ת�����ٶ�����ϵ��
		this->SM_3.acc = Vector3d(this->SM_3.Fx / this->SM_3.m, 0, 0);

		//�����˶�ѧģ��
		//if (this->SM_3.velocity >=800)
		//	this->SM_3.acc = Vector3d(0, 0, 0);
		this->SM_3.vel += this->SM_3.acc * dt;
		//cout << this->SM_3.vel << endl;
		this->SM_3.theta += this->SM_3.theta_dot * dt;
		this->SM_3.sigma += this->SM_3.sigma_dot * dt;
		this->SM_3.velocity = this->SM_3.vel[0];
		//�ٶ�����ϵת������������ϵ
		this->SM_3.x += this->SM_3.velocity * cos(this->SM_3.theta) * cos(this->SM_3.sigma) * dt;
		this->SM_3.y += this->SM_3.velocity * sin(this->SM_3.theta) * dt;
		this->SM_3.z += -1 * this->SM_3.velocity * cos(this->SM_3.theta) * sin(this->SM_3.sigma) * dt;
		this->expect_point.x += this->expect_point.velocity * cos(this->expect_point.theta) * cos(this->expect_point.sigma) * dt;
		this->expect_point.y += this->expect_point.velocity * sin(this->expect_point.theta) * dt;
		this->expect_point.z += -1 * this->expect_point.velocity * cos(this->expect_point.theta) * sin(this->expect_point.sigma) * dt;
		this->SM_3.pos = Vector3d(this->SM_3.x, this->SM_3.y, this->SM_3.z);
		this->expect_point.pos = Vector3d(this->expect_point.x, this->expect_point.y, this->expect_point.z);
		//���������˳����
		if (this->r > temp)//r��С�ˣ�����r_dot����0��   todo
		{
			flag = 1;
			if (temp > 15)
				flag = 0;
			break;
		}
		temp = this->r;
	}
	outfile.close();
	accfile.close();
	if (!flag)
		std::cout << "�Ѱ���>15��δ�ܸ���" << endl;
	std::cout << this->expect_point.pos[0] - this->SM_3.pos[0] << std::endl;
	std::cout << this->expect_point.pos[1] - this->SM_3.pos[1] << std::endl;
	std::cout << this->expect_point.pos[2] - this->SM_3.pos[2] << std::endl;
	std::cout << "�Ѱ���:" << (temp + this->r) / 2;

	return SM_3;
}
*/






//��Ա��������ͨ������������������������е�һ����Ա������ϵͳҪ������3�����涼Ҫƥ�䣺
//(1) �������������ͺͲ����ĸ�����	(2) ��������ֵ���͡�(3) �����������ࡣ
//����ָ���ó�Ա������ָ�������һ����ʽΪ������������������::*ָ����������������б���
//void (student::* pf)(); //���庯��ָ��pf
//void student::* (pf()); //���Ƿ���ֵΪvoid ����ָ��ĺ���
//ָ�����ָ��һ�����ó�Ա������һ����ʽΪ��ָ������� = &����������Ա��������
//����ж��ͬ��Ķ������ǹ���ͬһ����������Σ��ඨ��ĳ�Ա��������������ĳһ�����󣬶������ڸ��࣬�����е��������