#include"controller.h"
using namespace Eigen;
using namespace std;
//���������Ķ���Ӧ�÷���ͷ�ļ�.h��,�����������ⶨ���ʱ�򣬶���Ӧ�ú���������һ�𣬶�����ͷ�ļ�.h�У���Ȼ�����LNK2019�����޷������ⲿ�������

controller::controller() {
	//���ƽǶ�
	double a = 0;
	double b = 0;
	double r = 0;

	double ix = 0;
	double iy = 0;// iy1;
	double iz = 0;//iz1;

	this->V = Vector3d::Ones(3);//ע��Ҫ��ֵ3 �������ά������
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

	//��ִ�еĳ����д���ҲӰ��main����
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
//��bug lnk2019
Vector3d controller::control_instruct_abr(FlyStatus state, Vector3d Angle_expectation,Vector3d Omg_expectation)//todo
{
	this->eta << 0.2, 0, 0, 0, 5, 0, 0, 0, 23;// ��i��ͨ�� ���� ʹ�ø�ͨ������ ������ ����������ɢ eta ������ͨ�����ˣ�ʹ�� r����
	this->k1 << 0.2, 0, 0, 0, 10, 0, 0, 0, 15;//��i��ͨ�� ���� ʹ�ø�ͨ������ ������ ����������ɢ
	this->Angle_expectation = Angle_expectation;//�����ĽǶ�
	this->omega_expectation = Omg_expectation;//�����Ľ��ٶ�
	this->state = state;
	double a = state.alpha;
	double b = state.beta;
	double r = state.nu;
	double wx = state.attitude_omega[0];
	double wy = state.attitude_omega[1];
	double wz = state.attitude_omega[2];
	Vector3d Angle_current(a, b, r);//��ֱ�Ӹ�ֵ		//Omg_current << a, b, r;
	Vector3d omega = Vector3d::Zero();
	omega<< wx, wy, wz;
	//Vector3d omega(wx, wy, wz);
	double ix = state.J[0];
	double iy = state.J[1];// iy1;
	double iz = state.J[2];//iz1;
	//x=a,b,r,wx,wy,wz
	//x_dot=fx+gM
	//y=hx
	//---------------------------------------------���� p wx r -wy q wz  u rv-----------------------------------------------
	f_[0] = wz - tan(b) * (wx * cos(a) - wy * sin(a));//alpha_dot
	f_[1] = wx * sin(a) + wy * cos(a);
	f_[2] = (wx * cos(a) - wy * sin(a)) / cos(b);
	f2_[0] = -(iz - iy) / ix * wz * wy;
	f2_[1] = -(ix - iz) / iy * wx * wz;
	f2_[2] = -(iy - ix) / iz * wy * wx;
	g[0] = 1 / ix;
	g[1] = 1 / iy;
	g[2] = 1 / iz;
	//������ϵ	tan�� ��cot����1		sin�� ��csc����1		cos�� ��sec����1
	F[0] = (wx * sin(a) * tan(b) + wy * cos(a) * tan(b)) * f_[0] - (wx * cos(a) - wy * sin(a)) / (pow(cos(b), 2)) * f_[1] - (cos(a) * tan(b)) * f2_[0] + f2_[1] - (sin(a) * tan(b)) * f2_[2];
	F[1] = (-wy * sin(a) + wx * cos(a)) * f_[0] + sin(a) * f2_[0] - cos(a) * f2_[2];
	F[2] = (-wx * sin(a) / cos(b) - wy * cos(a) / cos(b)) * f_[0] + (wx * cos(a) * tan(b) / cos(b) - wy * sin(a) * tan(b) / cos(b)) * f_[1] + (cos(a) / cos(b)) * f2_[0] + (sin(a) / cos(b) * f2_[2]);
	//�˴������Ի���ֻ�ǰ� �Ƕȶ��׵��Ŀ��������M��ʾ����ʵ�ֽ���
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
	// M = (E.inverse() * (V - F));//-sat(,15);//�������������޵ģ������޷�  E.inverse()���ܳ���
	return M;
}
/*
Vector3d controller::control_instruct(double a, double b, double r, double wx, double wy, double wz, double ix, double iy, double iz, Vector3d Omg_expectation)
{

	Vector3d Angle_current(a, b, r);//��ֱ�Ӹ�ֵ		//Omg_current << a, b, r;
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
	M = (E.inverse() * (V - F));//-sat(,15);//�������������޵ģ������޷�  E.inverse()���ܳ���
	//std::cout << "��������ȷ�� " << E.inverse()*E << std::endl;
	return M;

}

*/

/*
missile_rotate::missile_rotate(double fai1 = 0, double pothi1 = 0, double gamma1 = 0, double a1 = 0, double b1 = 0, double r1 = 0, Vector3d Omega1 = Vector3d::Zero(), Vector3d Omega_dot1 = Vector3d::Zero(), double ix1 = 1, double iy1 = 1, \
	double iz1 = 1, Vector3d Angle_expect1 = default_Angle_expect) :qd_control(Angle_expect1)//c++����/�������� ���ܼ���ʽ����  Vector3d::Zero(), Omega_dot1[0], Omega_dot1[1], Omega_dot1[2], ix1, iy1, iz1,
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

//������ת�����˶�ѧģ�ͣ�����M�����W
void missile_rotate::YD(Vector3d M)//Vector3d  ����M ���� w w_dot
{
	this->Omega_dot[0] = M[0] / ix;
	this->Omega_dot[1] = M[1] / iy;
	this->Omega_dot[2] = M[2] / iz;
	this->Omega += this->Omega_dot * dt;//���ٶ�
	//���Ƕ�Ӧ��ϵ,�Ҵ˴�����Ϊ��̬�ǣ������ǹ��ǣ��໬�ǣ��ٶ����ǣ�Ӧ������̬��
	//this->r += this->Omega[0] * dt;//�Ƕ�
	//this->b += this->Omega[1] * dt;
	//this->a += this->Omega[2] * dt;
	wx = this->Omega[0];
	wy = this->Omega[1];
	wz = this->Omega[2];
	qd_control.f_[0] = wz - tan(b) * (wx * cos(a) - wy * sin(a));
	qd_control.f_[1] = wx * sin(a) + wy * cos(a);
	qd_control.f_[2] = (wx * cos(a) - wy * sin(a)) / cos(b);
	//�˴��Ƕ�Ӧ��ϵ
	this->a += qd_control.f_[0] * dt;//�Ƕ�
	this->b += qd_control.f_[1] * dt;
	this->r += qd_control.f_[2] * dt;
	//�״�
	//this->gamma += this->Omega[0] * dt;//�Ƕ�
	//this->pothi += this->Omega[1] * dt;
	//this->fai += this->Omega[2] * dt;
	//���� ���ٶ��� ��̬�ǵĵ��� ��ϵ ����ֱ�ӱ任
	//�����⣿����
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

//std::variant �����庯���ķ������� ��������һ�������ڷ��ز�ͬ���͵�ֵ



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
	std::cout << "---------��̬��Ϊ--------:" << MR.fai * 180 / Pi << " " << MR.pothi * 180 / Pi << " " << MR.gamma * 180 / Pi << "------------\n" << std::endl;

}




*/
//Ϊ�˽��ʹ�����ϣ���߸�ģ��Ķ����ԡ������������õ����ֻ࣬������Ϣ�ͼ������أ�ʵ�ʵ������������ɵ��������
Vector3d controller::control_instruct_eular(FlyStatus state1, Vector3d Angle_expectation1, Vector3d Omg_expectation1) {

	Angle_expectation = Angle_expectation1;//�����ĽǶ�
	omega_expectation = Omg_expectation1;//�����Ľ��ٶ�
	state = state1;

	//
	double theta = state.EularAngle[0];
	double sigma = state.EularAngle[1];
	double gm = state.EularAngle[2];
	double wx = state.attitude_omega[0];
	double wy = state.attitude_omega[1];
	double wz = state.attitude_omega[2];
	Vector3d Angle_current(theta, sigma, gm);//��ֱ�Ӹ�ֵ		//Omg_current << a, b, r;
	Vector3d omega = Vector3d::Zero();

	//Vector3d omega(wx, wy, wz);
	double ix = state.J[0];
	double iy = state.J[1];// iy1;
	double iz = state.J[2];//iz1;
	//
	double theta_dot, sigma_dot, gm_dot;
	//double theta_dot2, sigma_dot2, gm_dot2;//�������������̬��΢�ַ���ֻ��Ϊ���������õģ������������ڶ���ѧģ��ʵ��


	//90->60->60
	//eta << 0.1, 0, 0,//�������󣬿��Լӿ������ٶȣ�����λ��ƫ��,���ǲ������ˣ������𵴣��Ƕ�������׼  1/20/30����̬�ǲ�׼
	//	0, 5, 0,
	//	0, 0, 10;// 0.1 1 3
	//k1 << 2, 0, 0, //2 5 20
	//	0, 5, 0,
	//	0, 0, 10;
	//k2 << 2, 0, 0,//��i��ͨ�� ���� ʹ�ø�ͨ������ ������ ����������ɢ  1 2 30
	//	0, 10, 0,
	//	0, 0, 50;

	//90->60->45->45
	eta <<0.1, 0, 0,//�������󣬿��Լӿ������ٶȣ�����λ��ƫ��,���ǲ������ˣ������𵴣��Ƕ�������׼  1/20/30����̬�ǲ�׼
		0, 5, 0, 
		0, 0, 20;// 0.1 1 3
	k1 << 2, 0, 0, //2 5 20
		0, 5, 0, 
		0, 0, 20;
	k2 << 2, 0, 0,//��i��ͨ�� ���� ʹ�ø�ͨ������ ������ ����������ɢ  1 2 30
		0, 10, 0, 
		0, 0,100;

	Vector3d slide_surface = (state.attitude_omega - omega_expectation) + k1 * (state.EularAngle - Angle_expectation);
	V = -k1 * (state.attitude_omega - omega_expectation) - eta * sign_sat(slide_surface,0.01) - k2 * slide_surface;//lnk-2019�������
	//V = -k1 * (state.attitude_omega - omega_expectation) - eta * Sign(slide_surface) - k2 * slide_surface;//
	//V <<0,0,0;
	//��ת���˶�ѧģ��
	//һ�׵�,����theta=90����������תѭ��
	theta_dot = wy * sin(gm)/cos(sigma) + wz * cos(gm)/cos(sigma);
	sigma_dot = wy * cos(gm) - wz * sin(gm);
	gm_dot = wx + wy * tan(sigma) * sin(gm) + wz * tan(sigma) * cos(gm);
	//���׵�
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
		cout << "��������ֵ " << endl;
		return Vector3d{ 0,0,0 };
	}
	M = E.inverse() * (V - F);
	M = sature(M, 5000);//�Կ����������ұ���
	M[0] = 0;
	M[1] = 0;
	return M;
}