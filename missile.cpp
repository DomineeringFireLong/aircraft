#include<iostream>
#include<cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include<fstream>
#include<random>
#include"missile.h"//必须放到同一个路径下 否则搜不到



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

	this->Eular_angle = Eular_angle;//姿态角
	this->Eular_omega = Eular_omega;//姿态角速度
	this->attitude_omega = attitude_omega;//绕弹体的角速度
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

	this->pitch = Eular_angle[0];//俯仰 wz
	this->yaw = Eular_angle[1];//偏航  wy
	this->roll = Eular_angle[2];// wx

	this->alpha = alpha;
	this->alpha_dot = alpha_dot;
	this->beta = beta;
	this->beta_dot = beta_dot;
	this->nu = nu;
	this->nu_dot = nu_dot;
}
//拷贝构造函数，又称复制构造函数,形参必须是引用，但并不限制为const
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

	this->pitch = Eular_angle[0];//俯仰 wz
	this->yaw = Eular_angle[1];//偏航  wy
	this->roll = Eular_angle[2];// wx

	this->alpha = other.alpha;
	this->alpha_dot = other.alpha_dot;
	this->beta= other.beta;
	this->beta_dot = other.beta_dot;
	this->nu = other.nu;
	this->nu_dot = other.nu_dot;;
	cout << "拷贝构造函数" << endl;//只要初始化时候才能赋值，否则调用完后不能修改 值
}

//赋值运算符重载,注意不是构造函数，是普通运算符函数
//missile missile::operator=(missile& const other) {
//
//	missile m1(other);
//
//	return m1;//返回临时变量 
//}


FlyStatus missile::Extract_status() {

	FlyStatus current_status = { m,velocity,vel,Eular_angle,attitude_omega,pos,alpha,beta,theta,sigma,nu,S_ref,L_ref,inertia,time};
	return current_status;

}
//------------------------------------------质心动力学----------------------------------------//
void missile::dynamic(Vector3d F) {
	velocity = vel[0];//弹道坐标系
	acc = F/ m;
	theta_dot = F[1] / (m * velocity);////整数相除时为取 整数 ，double型就是小数除法； %取余
	sigma_dot = -1 * F[2] / (m * velocity * cos(theta));
	
	//注意质量变换,在外部类改变
	//实际是控制法向过载，改变速度方向的偏角	//如果某个变量变换非常大，可能是积分时，导数没乘dt
	vel += Vector3d(F[0] / m,0, 0) * dt;
	theta += theta_dot * dt;
	sigma += sigma_dot * dt;
	velocity = vel[0];
	//cout << vel << endl;
	//vel是速度坐标系，转换到地面坐标系，对于地面坐标系对xyz运算，导弹、目标对地位置更新
	x += velocity * cos(theta) * cos(sigma) * dt;
	y += velocity * sin(theta) * dt;
	z += -1 * velocity * cos(theta) * sin(sigma) * dt;
	//cout << "theta" << theta << endl;
	//cout <<"近似高度" << y << endl;
	pos = Vector3d(x, y, z);


}

//---------------------------------------绕质心转动动力学---------------------------------------//
//绕质心转动的运动学模型    输入M，输出W
void missile::rotate(Vector3d M)//Vector3d  输入M 更新 w w_dot  只能控制攻角、侧滑角、速度滚转角
{
	wx = this->attitude_omega[0];
	wy = this->attitude_omega[1];
	wz = this->attitude_omega[2];
	attitude_omega_acc[0] = (M[0] - (this->iz - this->iy) * this->wz * this->wy) / this->ix;
	attitude_omega_acc[1] = (M[1] - (this->ix - this->iz) * this->wx * this->wz) / this->iy;
	attitude_omega_acc[2] = (M[2] - (this->iy - this->ix) * this->wy * this->wx) / this->iz;
	attitude_omega += this->attitude_omega_acc * dt;//角速度
	//不是对应关系,且此处积分为姿态角，并不是攻角，侧滑角，速度倾侧角，应该是姿态角
	//this->r += this->attitude_omega[0] * dt;//角度
	//this->b += this->attitude_omega[1] * dt;
	//this->a += this->attitude_omega[2] * dt;
	wx = this->attitude_omega[0];
	wy = this->attitude_omega[1];
	wz = this->attitude_omega[2];
	

	//速度到弹体
	alpha_dot = wz - tan(beta) * (wx * cos(alpha) - wy * sin(alpha));
	beta_dot = wx * sin(alpha) + wy * cos(alpha);
	//弹道到速度
	nu = (wx * cos(alpha) - wy * sin(alpha)) / cos(beta);
	//此处是对应关系
	alpha += alpha_dot * dt;//角度
	beta += beta_dot * dt;
	nu += nu_dot * dt;

	//易错处
	//this->roll += this->attitude_omega[0] * dt;//角度
	//this->pothi += this->attitude_omega[1] * dt;
	//this->sigma += this->attitude_omega[2] * dt;
	//弹体 角速度与 姿态角的导数 关系 不能直接变换 
	//可能出现俯仰90度时 奇异:改变顺序  todo
	//垂直段theta=90
	//先旋转偏航，再俯仰和滚转
	//Eular_angle[0] += (sin(roll) * wy + cos(roll) * wz) * dt;
	//Eular_angle[1] += (cos(roll) * wy  - sin(roll) * wz) /cos(theta) * dt;
	//Eular_angle[2] += (wx+tan(theta)*( -cos(roll) * wy + sin(roll) * wz)) * dt;
	
	//先旋转俯仰，再偏航和滚转
	Eular_angle[0] += (sin(roll) * wy + cos(roll) * wz) / cos(sigma) * dt;
	Eular_angle[1] += (cos(roll) * wy - sin(roll) * wz) * dt;
	Eular_angle[2] += (wx + tan(sigma) * (sin(roll) * wy + cos(roll) * wz)) * dt;

	attitude_omega<< wx, wy, wz; //+= this->attitude_omega * dt;
	//将三个角度通过微分方程和几何关系进行数据融合
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


//std::variant 来定义函数的返回类型 允许你在一个函数内返回不同类型的值


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
		std::cout << "--过程error =  " << (angle - ae) << std::endl;
		//std::cout << mo(MR.get_angle() - MR.qd_control.Omg_expectation) << std::endl;
		//从动力学模型中获得参数		//把参数传入控制器模块，产生期望力矩
		Vector3d M = MR.qd_control.control_instruct(angle[0], angle[1], angle[2], w[0], w[1], w[2], J[0], J[1], J[2], ae);
		//力矩给到动力学模块，产生状态变量的更新

		temp = angle[0];
		MR.YD(M);
		w = MR.get_w();
		angle = MR.get_angle();//更新参数
		if (abs(temp - ae[0]) < (abs(angle[0] - ae[0])))
			break;

		if (count > 500)
		{
			std::cout << "-------控制超时：5s以上----------" << std::endl;
			break;
		}
	}
	std::cout << "\n----------控制时间----------:" << count * 0.01 << "------------\n" << std::endl;

	angle = MR.get_angle();

	std::cout << "----------控制角度error---------:" << (angle - ae).transpose() * 180 / Pi << "------------\n" << std::endl;
	std::cout << "---------姿态角为--------:" << MR.sigma * 180 / Pi << " " << MR.pothi * 180 / Pi << " " << MR.roll * 180 / Pi << "------------\n" << std::endl;

}


*/