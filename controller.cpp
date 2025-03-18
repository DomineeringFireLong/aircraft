#include"controller.h"
using namespace Eigen;
using namespace std;
//内联函数的定义应该放在头文件.h中,内联函数类外定义的时候，定义应该和声明放在一起，都放在头文件.h中，不然会出现LNK2019错误，无法解析外部命令！！！

controller::controller() {
	//控制角度
	double a = 0;
	double b = 0;
	double r = 0;

	double ix = 0;
	double iy = 0;// iy1;
	double iz = 0;//iz1;

	this->V = Vector3d::Ones(3);//注意要赋值3 否则出现维数报错
	this->eta = Matrix3d::Identity(3, 3);	
	this->k1 = Matrix3d::Identity(3, 3);
	this->k2 = Matrix3d::Identity(3, 3);

	this->f_[0] = 0;// wz - tan(b) * (wx * cos(a) - wy * sin(a));
	this->f_[1] = 0;//wx * sin(a) + wy * cos(a);
	this->f_[2] = 0;// (wx * cos(a) - wy * sin(a)) / cos(b);
	this->f2_[0] = 0;
	this->f2_[1] = 0;
	this->f2_[2] = 0;

	g[0] = 0;//1/ix;
	g[1] = 0;//1/iy;
	g[2] = 0;//1/iz;
	F[0] = 0;//(wx * sin(a) * tan(b) + wy * cos(a) * tan(b)) * f[0] - (wx * cos(a) / (pow(cos(b), 2)) - wy * sin(a) / (pow(cos(b), 2))) * f[1];
	F[1] = 0;//(-wy * sin(a) + wx * cos(a)) * f[0];
	F[2] = 0;//(-wx * sin(a) / cos(b) - wy * cos(a) / cos(b)) * f[0] + (wx * cos(a) * tan(b) / cos(b) - wy * sin(a) * tan(b) / cos(b)) * f[1];	

	//不执行的程序有错误，也影响main函数
	E(0, 0) = -tan(b) * cos(a) / ix;
	E(0, 1) = tan(b) * sin(a) / iy;
	E(0, 2) = 1 / iz;
	E(1, 0) = sin(a) / ix;
	E(1, 1) = cos(a) / iy;
	E(1, 2) = 0;
	E(2, 0) = cos(a) / (ix * cos(b));
	E(2, 1) = -sin(a) / (iy * cos(b));
	E(2, 2) = 0;
	Angle_expectation << 0, 0, 0;
	omega_expectation << 0, 0, 0;
	M << 0, 0, 0;
}
//有bug lnk2019
Vector3d controller::control_instruct_abr(FlyStatus state, Vector3d Angle_expectation,Vector3d Omg_expectation)//todo
{
	this->eta << 0.2, 0, 0, 0, 5, 0, 0, 0, 23;// 第i的通道 增大 使得该通道参数 更收敛 ，但其他发散 eta 第三个通道大了，使得 r收敛
	this->k1 << 0.2, 0, 0, 0, 10, 0, 0, 0, 15;//第i的通道 增大 使得该通道参数 更收敛 ，但其他发散
	this->Angle_expectation = Angle_expectation;//期望的角度
	this->omega_expectation = Omg_expectation;//期望的角速度
	this->state = state;
	double a = state.alpha;
	double b = state.beta;
	double r = state.nu;
	double wx = state.attitude_omega[0];
	double wy = state.attitude_omega[1];
	double wz = state.attitude_omega[2];
	Vector3d Angle_current(a, b, r);//能直接赋值		//Omg_current << a, b, r;
	Vector3d omega = Vector3d::Zero();
	omega<< wx, wy, wz;
	//Vector3d omega(wx, wy, wz);
	double ix = state.J[0];
	double iy = state.J[1];// iy1;
	double iz = state.J[2];//iz1;
	//x=a,b,r,wx,wy,wz
	//x_dot=fx+gM
	//y=hx
	//---------------------------------------------耿洁 p wx r -wy q wz  u rv-----------------------------------------------
	f_[0] = wz - tan(b) * (wx * cos(a) - wy * sin(a));//alpha_dot
	f_[1] = wx * sin(a) + wy * cos(a);
	f_[2] = (wx * cos(a) - wy * sin(a)) / cos(b);
	f2_[0] = -(iz - iy) / ix * wz * wy;
	f2_[1] = -(ix - iz) / iy * wx * wz;
	f2_[2] = -(iy - ix) / iz * wy * wx;
	g[0] = 1 / ix;
	g[1] = 1 / iy;
	g[2] = 1 / iz;
	//倒数关系	tanα ·cotα＝1		sinα ·cscα＝1		cosα ·secα＝1
	F[0] = (wx * sin(a) * tan(b) + wy * cos(a) * tan(b)) * f_[0] - (wx * cos(a) - wy * sin(a)) / (pow(cos(b), 2)) * f_[1] - (cos(a) * tan(b)) * f2_[0] + f2_[1] - (sin(a) * tan(b)) * f2_[2];
	F[1] = (-wy * sin(a) + wx * cos(a)) * f_[0] + sin(a) * f2_[0] - cos(a) * f2_[2];
	F[2] = (-wx * sin(a) / cos(b) - wy * cos(a) / cos(b)) * f_[0] + (wx * cos(a) * tan(b) / cos(b) - wy * sin(a) * tan(b) / cos(b)) * f_[1] + (cos(a) / cos(b)) * f2_[0] + (sin(a) / cos(b) * f2_[2]);
	//此处“线性化”只是把 角度二阶导的控制量拆成M表示，并实现解耦
	E(0, 0) = -tan(b) * cos(a) / ix;
	E(0, 1) = tan(b) * sin(a) / iy;
	E(0, 2) = 1 / iz;
	E(1, 0) = sin(a) / ix;
	E(1, 1) = cos(a) / iy;
	E(1, 2) = 0;
	E(2, 0) = cos(a) / (ix * cos(b));
	E(2, 1) = -sin(a) / (iy * cos(b));
	E(2, 2) = 0;

	Vector3d s = (omega - omega_expectation) + k1 * (Angle_current - Angle_expectation);

	//todo
	// V= -k * (omega - omega_expectation) - eta * sign(s);//todo
	// M = (E.inverse() * (V - F));//-sat(,15);//气动力矩是有限的，必须限幅  E.inverse()可能出错
	return M;
}
/*
Vector3d controller::control_instruct(double a, double b, double r, double wx, double wy, double wz, double ix, double iy, double iz, Vector3d Omg_expectation)
{

	Vector3d Angle_current(a, b, r);//能直接赋值		//Omg_current << a, b, r;
	Vector3d omega(wx, wy, wz);

	f_[0] = wz - tan(b) * (wx * cos(a) - wy * sin(a));
	f_[1] = wx * sin(a) + wy * cos(a);
	f_[2] = (wx * cos(a) - wy * sin(a)) / cos(b);
	f2_[0] = -(iz - iy) / ix * wz * wy;
	f2_[1] = -(ix - iz) / iy * wx * wz;
	f2_[2] = -(iy - ix) / iz * wy * wx;
	g[0] = 1 / ix;
	g[1] = 1 / iy;
	g[2] = 1 / iz;
	//F[0] = (wx * sin(a) * tan(b) + wy * cos(a) * tan(b)) * f_[0] - (wx * cos(a) - wy * sin(a) )/ (pow(cos(b), 2)) * f_[1];
	//F[1] = (-wy * sin(a) + wx * cos(a)) * f_[0];
	//F[2] = (-wx * sin(a) / cos(b) - wy * cos(a) / cos(b)) * f_[0] + (wx * cos(a) * tan(b) / cos(b) - wy * sin(a) * tan(b) / cos(b)) * f_[1];

	F[0] = (wx * sin(a) * tan(b) + wy * cos(a) * tan(b)) * f_[0] - (wx * cos(a) - wy * sin(a)) / (pow(cos(b), 2)) * f_[1] - (cos(a) * tan(b)) * f2_[0] + f2_[1] - (sin(a) * tan(b)) * f2_[2];
	F[1] = (-wy * sin(a) + wx * cos(a)) * f_[0] + sin(a) * f2_[0] - cos(a) * f2_[2];
	F[2] = (-wx * sin(a) / cos(b) - wy * cos(a) / cos(b)) * f_[0] + (wx * cos(a) * tan(b) / cos(b) - wy * sin(a) * tan(b) / cos(b)) * f_[1] + (cos(a) / cos(b)) * f2_[0] + (sin(a) / cos(b) * f2_[2]);
	E(0, 0) = -tan(b) * cos(a) / ix;
	E(0, 1) = tan(b) * sin(a) / iy;
	E(0, 2) = 1 / iz;
	E(1, 0) = sin(a) / ix;
	E(1, 1) = cos(a) / iy;
	E(1, 2) = 0;
	E(2, 0) = cos(a) / (ix * cos(b));
	E(2, 1) = -sin(a) / (iy * cos(b));
	E(2, 2) = 0;

	Vector3d s = (omega - omega_expectation) + k * (Angle_current - Angle_expectation);

	V = -k * (omega - omega_expectation) - eta * sign(s);//todo
	M = (E.inverse() * (V - F));//-sat(,15);//气动力矩是有限的，必须限幅  E.inverse()可能出错
	//std::cout << "逆运算正确性 " << E.inverse()*E << std::endl;
	return M;

}

*/

/*
missile_rotate::missile_rotate(double fai1 = 0, double pothi1 = 0, double gamma1 = 0, double a1 = 0, double b1 = 0, double r1 = 0, Vector3d Omega1 = Vector3d::Zero(), Vector3d Omega_dot1 = Vector3d::Zero(), double ix1 = 1, double iy1 = 1, \
	double iz1 = 1, Vector3d Angle_expect1 = default_Angle_expect) :qd_control(Angle_expect1)//c++的类/函数调用 不能加形式参数  Vector3d::Zero(), Omega_dot1[0], Omega_dot1[1], Omega_dot1[2], ix1, iy1, iz1,
{

	this->a = a1;
	this->b = b1;
	this->r = r1;
	this->Contorl_Angle << a1, b1, r1;
	this->Omega << Omega1;
	this->wx = this->Omega[0];
	this->wy = this->Omega[1];
	this->wz = this->Omega[2];
	this->Omega_dot << Omega_dot1;
	this->ix = ix1;
	this->iy = iy1;
	this->iz = iz1;
	this->fai = fai1;
	this->pothi = pothi1;
	this->gamma = gamma1;
	//this->qd_control.Angle_expectation = Angle_expect1;
	//att << 0, 0, 0;
}

//绕质心转动的运动学模型，输入M，输出W
void missile_rotate::YD(Vector3d M)//Vector3d  输入M 更新 w w_dot
{
	this->Omega_dot[0] = M[0] / ix;
	this->Omega_dot[1] = M[1] / iy;
	this->Omega_dot[2] = M[2] / iz;
	this->Omega += this->Omega_dot * dt;//角速度
	//不是对应关系,且此处积分为姿态角，并不是攻角，侧滑角，速度倾侧角，应该是姿态角
	//this->r += this->Omega[0] * dt;//角度
	//this->b += this->Omega[1] * dt;
	//this->a += this->Omega[2] * dt;
	wx = this->Omega[0];
	wy = this->Omega[1];
	wz = this->Omega[2];
	qd_control.f_[0] = wz - tan(b) * (wx * cos(a) - wy * sin(a));
	qd_control.f_[1] = wx * sin(a) + wy * cos(a);
	qd_control.f_[2] = (wx * cos(a) - wy * sin(a)) / cos(b);
	//此处是对应关系
	this->a += qd_control.f_[0] * dt;//角度
	this->b += qd_control.f_[1] * dt;
	this->r += qd_control.f_[2] * dt;
	//易错处
	//this->gamma += this->Omega[0] * dt;//角度
	//this->pothi += this->Omega[1] * dt;
	//this->fai += this->Omega[2] * dt;
	//弹体 角速度与 姿态角的导数 关系 不能直接变换
	//有问题？？？
	this->gamma += (wx + (-tan(fai) * cos(gamma) * wy + (tan(fai) * sin(gamma) * wz))) * dt;
	this->pothi += (cos(gamma) * wy / cos(fai) - sin(gamma) * wz / cos(fai)) * dt;
	this->fai += (sin(gamma) * wx + cos(r) * wz) * dt;
	this->Contorl_Angle << a, b, r; //+= this->Omega * dt;

	//return Omega;
}

Vector3d missile_rotate::get_angle()
{
	return this->Contorl_Angle;
}


Vector3d missile_rotate::get_w()
{
	return this->Omega;
}

Vector3d missile_rotate::get_w_dot()
{
	return this->Omega_dot;
}
Vector3d missile_rotate::get_J()
{
	Vector3d J(ix, iy, iz);
	return J;
}

//std::variant 来定义函数的返回类型 允许你在一个函数内返回不同类型的值



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
	std::cout << "---------姿态角为--------:" << MR.fai * 180 / Pi << " " << MR.pothi * 180 / Pi << " " << MR.gamma * 180 / Pi << "------------\n" << std::endl;

}




*/
//为了降低代码耦合，提高各模块的独立性。控制器不设置导弹类，只传递信息和计算力矩，实际导弹参数更新由导弹类完成
Vector3d controller::control_instruct_eular(FlyStatus state1, Vector3d Angle_expectation1, Vector3d Omg_expectation1) {

	Angle_expectation = Angle_expectation1;//期望的角度
	omega_expectation = Omg_expectation1;//期望的角速度
	state = state1;

	//
	double theta = state.EularAngle[0];
	double sigma = state.EularAngle[1];
	double gm = state.EularAngle[2];
	double wx = state.attitude_omega[0];
	double wy = state.attitude_omega[1];
	double wz = state.attitude_omega[2];
	Vector3d Angle_current(theta, sigma, gm);//能直接赋值		//Omg_current << a, b, r;
	Vector3d omega = Vector3d::Zero();

	//Vector3d omega(wx, wy, wz);
	double ix = state.J[0];
	double iy = state.J[1];// iy1;
	double iz = state.J[2];//iz1;
	//
	double theta_dot, sigma_dot, gm_dot;
	//double theta_dot2, sigma_dot2, gm_dot2;//控制器这里的姿态角微分方程只是为了求力矩用的，而参数更新在动力学模型实现


	//90->60->60
	//eta << 0.1, 0, 0,//参数增大，可以加快收敛速度，减少位置偏移,但是参数大了，导致震荡，角度收敛不准  1/20/30但姿态角不准
	//	0, 5, 0,
	//	0, 0, 10;// 0.1 1 3
	//k1 << 2, 0, 0, //2 5 20
	//	0, 5, 0,
	//	0, 0, 10;
	//k2 << 2, 0, 0,//第i的通道 增大 使得该通道参数 更收敛 ，但其他发散  1 2 30
	//	0, 10, 0,
	//	0, 0, 50;

	//90->60->45->45
	eta <<0.1, 0, 0,//参数增大，可以加快收敛速度，减少位置偏移,但是参数大了，导致震荡，角度收敛不准  1/20/30但姿态角不准
		0, 5, 0, 
		0, 0, 20;// 0.1 1 3
	k1 << 2, 0, 0, //2 5 20
		0, 5, 0, 
		0, 0, 20;
	k2 << 2, 0, 0,//第i的通道 增大 使得该通道参数 更收敛 ，但其他发散  1 2 30
		0, 10, 0, 
		0, 0,100;

	Vector3d slide_surface = (state.attitude_omega - omega_expectation) + k1 * (state.EularAngle - Angle_expectation);
	V = -k1 * (state.attitude_omega - omega_expectation) - eta * sign_sat(slide_surface,0.01) - k2 * slide_surface;//lnk-2019报错语句
	//V = -k1 * (state.attitude_omega - omega_expectation) - eta * Sign(slide_surface) - k2 * slide_surface;//
	//V <<0,0,0;
	//旋转的运动学模型
	//一阶导,由于theta=90，交换了旋转循序
	theta_dot = wy * sin(gm)/cos(sigma) + wz * cos(gm)/cos(sigma);
	sigma_dot = wy * cos(gm) - wz * sin(gm);
	gm_dot = wx + wy * tan(sigma) * sin(gm) + wz * tan(sigma) * cos(gm);
	//二阶导
	//theta_dot2 = wy * cos(gm) * gm_dot - wz * sin(gm) * gm_dot;//SECX=1/COSX
	//sigma_dot2 = wy * cos(gm) * tan(theta) / cos(theta) * theta_dot - wz * sin(gm) * tan(theta) / cos(theta) * theta_dot - wy * sin(gm) / cos(theta) * gm_dot - wz * cos(gm) / cos(theta) * gm_dot;
	//gm_dot2 = wy * sin(gm) * tan(theta) * gm_dot + wz * cos(gm) * tan(theta) * gm_dot - wy * cos(gm) / cos(theta) / cos(theta) * theta_dot + wz * sin(gm) / cos(theta) / cos(theta) * theta_dot;

	//theta_dot2 = wy * cos(gm) * gm_dot - wz * sin(gm) * gm_dot +sin(gm)* (iz - ix) * wx * wz / iy+cos(gm)*(ix-iy)*wx*wy/iz;//SECX=1/COSX
	//sigma_dot2 = wy * cos(gm) * tan(theta) / cos(theta) * theta_dot - wz * sin(gm) * tan(theta) / cos(theta) * theta_dot - wy * sin(gm) / cos(theta) * gm_dot - wz * cos(gm) / cos(theta) * gm_dot + cos(gm)/cos(theta) * (iz - ix) * wx * wz / iy -sin(gm)/cos(theta) * (ix - iy) * wx * wy / iz;
	//gm_dot2 = wy * sin(gm) * tan(theta) * gm_dot + wz * cos(gm) * tan(theta) * gm_dot - wy * cos(gm) / cos(theta) / cos(theta) * theta_dot + wz * sin(gm) / cos(theta) / cos(theta) * theta_dot + (iy - iz) * wz * wy / ix - cos(gm) *tan(theta) * (iz - ix) * wx * wz / iy + sin(gm) *tan(theta) * (ix - iy) * wx * wy / iz;

	F[0] = wy * cos(gm) * gm_dot - wz * sin(gm) * gm_dot + sin(gm) * (iz - ix) * wx * wz / iy + cos(gm) * (ix - iy) * wx * wy / iz;//SECX=1/COSX
	F[1] = wy * cos(gm) * tan(theta) / cos(theta) * theta_dot - wz * sin(gm) * tan(theta) / cos(theta) * theta_dot - wy * sin(gm) / cos(theta) * gm_dot - wz * cos(gm) / cos(theta) * gm_dot + cos(gm) / cos(theta) * (iz - ix) * wx * wz / iy - sin(gm) / cos(theta) * (ix - iy) * wx * wy / iz;
	F[2] = wy * sin(gm) * tan(theta) * gm_dot + wz * cos(gm) * tan(theta) * gm_dot - wy * cos(gm) / cos(theta) / cos(theta) * theta_dot + wz * sin(gm) / cos(theta) / cos(theta) * theta_dot + (iy - iz) * wz * wy / ix - cos(gm) * tan(theta) * (iz - ix) * wx * wz / iy + sin(gm) * tan(theta) * (ix - iy) * wx * wy / iz;
	
	E(0, 0) = 0;
	E(0, 1) = sin(gm) / iy;
	E(0, 2) = cos(gm) / iz;
	
	E(1, 0) =	0 ;
	E(1, 1) = cos(gm)/cos(theta) / iy;
	E(1, 2) = -1 * sin(gm) / cos(theta) / iz;
	
	E(2, 0) = 1 / ix;
	E(2, 1) = -tan(theta)*cos(gm) / iy;
	E(2, 2) = tan(theta) * sin(gm) / iz;
	//cout << E.determinant();
	try
	{
		if (E.determinant() == 0)
			throw 0;
	}
	catch (int)
	{
		cout << "出现奇异值 " << endl;
		return Vector3d{ 0,0,0 };
	}
	M = E.inverse() * (V - F);
	M = sature(M, 5000);//对控制力矩世家饱和
	M[0] = 0;
	M[1] = 0;
	return M;
}