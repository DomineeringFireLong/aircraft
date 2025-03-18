#include<iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <random>
#include"math_missile.h"

//extern const double Pi = 3.14159265358979323846;//应用于一个全局变量，函数或模板声明，说明该符号具有外部链接(external linkage)属性,一般而言，C++全局变量的作用范围仅限于当前的文件
using namespace Eigen;
using namespace std;
Vector3d default_Angle_expect(30 / Pi, 30 / Pi, 0);
//地面坐标系转换到弹道坐标系
Vector3d G2T(Vector3d gxoy, double theta1, double eta1)
{
	Matrix3d transform;
	transform << cos(theta1) * cos(eta1), sin(theta1), -cos(theta1) * sin(eta1),
		-sin(theta1) * cos(eta1), cos(theta1), sin(theta1)* sin(eta1),
		sin(eta1), 0, cos(eta1);

	return transform * gxoy;
}
//弹道坐标系转换到地面坐标系
Vector3d T2G(Vector3d txoy, double theta1, double eta1)
{
	Matrix3d transform1;
	transform1 << cos(theta1) * cos(eta1), -sin(theta1) * cos(eta1), sin(eta1),
		sin(theta1), cos(theta1), 0,
		-cos(theta1) * sin(eta1), sin(theta1)* sin(eta1), cos(eta1);
	//transform1 = transform1.transpose();
	return  transform1 * txoy;
}
//绕转
Vector3d rao_z(Vector3d xoy, double theta1)
{
	Matrix3d transform2;
	transform2 << cos(theta1), sin(theta1), 0.0,
		-sin(theta1), cos(theta1), 0,
		0.0, 0.0, 1.0;
	return transform2 * xoy;
}
//绕y转
Vector3d rao_y(Vector3d xoy, double theta1)
{
	Matrix3d transform3;
	transform3 << cos(theta1), 0.0, -sin(theta1),
		0.0, 1.0, 0.0,
		sin(theta1), 0, cos(theta1);
	return transform3 * xoy;
}

//
Vector3d T2B(Vector3d txoy, double alpha1, double beta1)
{
	Matrix3d transform1;
	transform1 << cos(alpha1) * cos(beta1), -sin(alpha1) * cos(beta1), sin(beta1),
		sin(alpha1), cos(alpha1), 0,
		-cos(alpha1) * sin(beta1), sin(alpha1)* sin(beta1), cos(beta1);
	transform1 = transform1.transpose();
	return  transform1 * txoy;
}
//
Vector3d B2T(Vector3d bxoy, double alpha1, double beta1)
{
	Matrix3d transform1;
	transform1 << cos(alpha1) * cos(beta1), -sin(alpha1) * cos(beta1), sin(beta1),
		sin(alpha1), cos(alpha1), 0,
		-cos(alpha1) * sin(beta1), sin(alpha1)* sin(beta1), cos(beta1);
	return  transform1 * bxoy;

}
Vector3d Sign(Vector3d x)
{
	for (int i = 0; i < 3; i++)
	{
		if (x[i] > 0)
			x[i] = 1;
		else
		{
			if (x[i] < 0)
				x[i] = -1;
			else
				x[i] = 0;
		}
	}
	return x;
}



//向量求模
double modulo(VectorXd x)
{
	return sqrt(x.transpose() * x);
}

Vector3d sature(Vector3d x, double limit)
{
	for (int i = 0; i < 3; i++)
	{
		if (x[i] > limit)
			x[i] = limit;
		else {
			if (x[i] < -limit)
				x[i] = -limit;
		}
	}
	return x;
}

//向量饱和
VectorXd sat(VectorXd x, double limit)
{
	for (int i = 0; i < x.size(); i++)
	{
		if (x[i] > limit)
			x[i] = limit;
		else {
			if (x[i] < -limit)
				x[i] = -limit;
		}
	}
	return x;
}

//向量取绝对值
VectorXd abs(VectorXd x)
{
	VectorXd absx;
	for(int i=0;i<absx.size();i++)
		absx << abs(x[i]);//abs(x[0]), abs(x[1]), abs(x[2])
	return absx;
}

//向量取符号
VectorXd sign_vector(VectorXd x)
{
	for (int i = 0; i < x.size(); i++)
	{
		if (x[i] > 0)
			x[i] = 1;
		else
		{
			if (x[i] < 0)
				x[i] = -1;
			else
				x[i] = 0;
		}
	}
	return x;
}

Matrix3d Cross_product(Vector3d x)
{
	Matrix3d X;
	X << 0, -x[2], x[1], 
		x[2], 0, -x[0], 
		-x[1], x[0], 0;
	return X;
}





Vector3d calAlphaBetanu(Vector3d Eularangle,double theta,double psi)//Vector3d ballistic_angle
{
	// 知五求三
	double pitch = Eularangle[0];
	double yaw = Eularangle[1];
	double roll = Eularangle[2];

	//double nu = ballistic_angle[2];

	double sin_beta = sin(yaw - psi) * cos(theta) * cos(roll)
		+ cos(yaw - psi) * cos(theta) * sin(pitch) * sin(roll)
		- sin(theta) * sin(roll) * cos(pitch);
	double beta = asin(sin_beta);

	double sin_alpha = (cos(theta) * sin(pitch) * cos(roll) * cos(yaw - psi) - cos(theta) * sin(roll) * sin(yaw - psi) - sin(theta) * cos(pitch) * cos(roll)) / cos(beta);
	double alpha = asin(sin_alpha);

	double sin_nu = (cos(alpha) *sin(beta) * sin(pitch) - sin(alpha) * sin(beta) *cos(roll) * cos(pitch)+cos(beta)*sin(roll)*cos(pitch)) / cos(theta);
	double nu = asin(sin_nu);

	return {alpha, beta,nu};
}





////--------------------------欧拉角类-----------------------------//
/*
using namespace std;
Angle::Angle() {
}

JAngle::~Angle()
{
}
  
double Angle::calEta(double alpha, double beta)
{
	// 《远程火箭弹道学》, Eq.(6-1-11)
	return acos(cos(alpha) * cos(beta));
}

Vector3d Angle::calThetaSigmaGammav(const Vector3d& v)
{
	// 《弹道导弹弹道学》, P25
	double vx = v.x, vy = v.y, vz = v.z;

	//计算发射系弹道倾角
	double theta;
	if (vx < -SMALLNUM)
	{
		if (vy > 0)
			theta = PI - atan(fabs(vy / vx));
		else
			theta = -PI + atan(fabs(vy / vx));
	}
	else if (vx < SMALLNUM)
	{
		if (vy >= 0)
			theta = PI / 2 - SMALLNUM;//垂直发射，初始值
		else
			theta = -PI / 2 + SMALLNUM;
	}
	else
	{
		theta = atan(vy / vx);
	}

	//计算发射系弹道偏角
	double vel = v.norm();
	double sigma;
	if (vel < SMALLNUM)
		sigma = 0;
	else
		sigma = -asin(vz / vel); // 弹道偏角，以左为正

	double gammav = 0;
	return { theta, sigma, gammav };
}

Vector3d Angle::calPhiPsiGamma(Vector3d theta_sigma_gammav, Vector3d alpha_beta)
{
	double theta = theta_sigma_gammav.x;
	double sigma = theta_sigma_gammav.y;
	double gammav = theta_sigma_gammav.z;
	double alpha = alpha_beta.x;
	double beta = alpha_beta.y;

	double psi = asin(cos(alpha) * cos(beta) * sin(sigma) - sin(alpha) * cos(sigma) * sin(gammav) + cos(alpha) * sin(beta) * cos(sigma) * sin(gammav));

	double fais = (cos(alpha) * cos(beta) * sin(theta) * cos(sigma) + sin(alpha) * sin(theta) * sin(sigma) * sin(gammav) + sin(alpha) * cos(theta) * cos(gammav) - cos(alpha) * sin(beta) * sin(theta) * sin(sigma) * cos(gammav) + cos(alpha) * sin(beta) * cos(theta) * sin(gammav)) / cos(psi);
	double faic = (cos(alpha) * cos(beta) * cos(theta) * cos(sigma) + sin(alpha) * cos(theta) * sin(sigma) * sin(gammav) - sin(alpha) * sin(theta) * cos(gammav) - cos(alpha) * sin(beta) * cos(theta) * sin(sigma) * cos(gammav) - cos(alpha) * sin(beta) * sin(theta) * sin(gammav)) / cos(psi);
	if ((faic * faic + fais * fais - 1 >= 0.1) || ((faic * faic + fais * fais - 1 <= -0.1)))
	{
		cout << "角度计算错误" << endl;
	}

	double fai;
	if (fais > 1)
		fais = 1 - SMALLNUM;
	if (fais < -1)
		fais = -1 + SMALLNUM;

	if (faic >= 0)
		fai = asin(fais);
	else
	{
		if (fais >= 0)
			fai = PI - asin(fais);
		else
			fai = -PI - asin(fais);
	}
	return { fai, psi, 0 };
}

Vector3d Angle::calPhiPsiGamma(double B0, double A0, double t, const Vector3d& phi_psi_gamma_a)
{
	double X = cos(B0) * cos(A0);
	double Y = sin(B0);
	double Z = -cos(B0) * sin(A0);
	double S = sin(Earth::OMEAGA * t);
	double C = cos(Earth::OMEAGA * t);
	Matrix m1(3, 3);           // 发射到发惯系
	m1[0][0] = X * X * (1 - C) + C;
	m1[0][1] = X * Y * (1 - C) - Z * S;
	m1[0][2] = X * Z * (1 - C) + Y * S;
	m1[1][0] = X * Y * (1 - C) + Z * S;
	m1[1][1] = Y * Y * (1 - C) + C;
	m1[1][2] = Y * Z * (1 - C) - X * S;
	m1[2][0] = X * Z * (1 - C) - Y * S;
	m1[2][1] = Y * Z * (1 - C) + X * S;
	m1[2][2] = Z * Z * (1 - C) + C;

	// 发惯系姿态角
	double phia = phi_psi_gamma_a[0];
	double psia = phi_psi_gamma_a[1];
	double gammaa = phi_psi_gamma_a[2];
	Matrix m2(3, 3);                 // 发惯到弹体系转换矩阵
	m2[0][0] = cos(phia) * cos(psia);
	m2[0][1] = sin(phia) * cos(psia);
	m2[0][2] = -sin(psia);
	m2[1][0] = cos(phia) * sin(psia) * sin(gammaa) - sin(phia) * cos(gammaa);
	m2[1][1] = sin(phia) * sin(psia) * sin(gammaa) + cos(phia) * cos(gammaa);
	m2[1][2] = cos(psia) * sin(gammaa);
	m2[2][0] = cos(phia) * sin(psia) * cos(gammaa) + sin(phia) * sin(gammaa);
	m2[2][1] = sin(phia) * sin(psia) * cos(gammaa) - cos(phia) * sin(gammaa);
	m2[2][2] = cos(psia) * cos(gammaa);

	auto m = m2 * m1; // 发射到弹体系转换矩阵

	// 计算发射系姿态角
	double psi = asin(-m[0][2]);
	double sin_phi = m[0][1] / cos(psi);
	double cos_phi = m[0][0] / cos(psi);
	double sin_gamma = m[1][2] / cos(psi);
	double cos_gamma = m[2][2] / cos(psi);
	double phi =Math::calAngle2(sin_phi, cos_phi);
	double gamma = Math::calAngle2(sin_gamma, cos_gamma);
	return Vector3d(phi, psi, gamma);
}

Vector3d Angle::calPhiPsiGammaa(double B0, double A0, double t, const Vector3d& phi_psi_gamma)
{
	double X = cos(B0) * cos(A0);
	double Y = sin(B0);
	double Z = -cos(B0) * sin(A0);
	double S = sin(ZJEarth::OMEAGA * t);
	double C = cos(ZJEarth::OMEAGA * t);
	ZJMatrix m1(3, 3);
	m1[0][0] = X * X * (1 - C) + C;
	m1[0][1] = X * Y * (1 - C) - Z * S;
	m1[0][2] = X * Z * (1 - C) + Y * S;
	m1[1][0] = X * Y * (1 - C) + Z * S;
	m1[1][1] = Y * Y * (1 - C) + C;
	m1[1][2] = Y * Z * (1 - C) - X * S;
	m1[2][0] = X * Z * (1 - C) - Y * S;
	m1[2][1] = Y * Z * (1 - C) + X * S;
	m1[2][2] = Z * Z * (1 - C) + C;
	m1 = m1.Transpose(); // 发惯系到发射系

	double phi = phi_psi_gamma[0];
	double psi = phi_psi_gamma[1];
	double gamma = phi_psi_gamma[2];
	ZJMatrix m2(3, 3);

	//        m2.SetElement(0, 0, cos(phi)*cos(psi));
	//        m2.SetElement(0, 1, cos(phi)*sin(psi)*sin(gamma) - sin(phi)*cos(gamma));
	//        m2.SetElement(0, 2, cos(phi)*sin(psi)*cos(gamma) + sin(phi)*sin(gamma));
	//        m2.SetElement(1, 0, sin(phi)*cos(psi));
	//        m2.SetElement(1, 1, sin(phi)*sin(psi)*sin(gamma) + cos(phi)*cos(gamma));
	//        m2.SetElement(1, 2, sin(phi)*sin(psi)*cos(gamma) - cos(phi)*sin(gamma));
	//        m2.SetElement(2, 0, -sin(psi));
	//        m2.SetElement(2, 1, cos(psi)*sin(gamma));
	//        m2.SetElement(2, 2, cos(psi)*cos(gamma));

	double c_phi = cos(phi);
	double s_phi = sin(phi);
	double c_psi = cos(psi);
	double s_psi = sin(psi);
	double s_gamma = sin(gamma);
	double c_gamma = cos(gamma);
	m2.SetElement(0, 0, c_phi * c_psi);
	m2.SetElement(0, 1, c_phi * s_psi * s_gamma - s_phi * c_gamma);
	m2.SetElement(0, 2, c_phi * s_psi * c_gamma + s_phi * s_gamma);
	m2.SetElement(1, 0, s_phi * c_psi);
	m2.SetElement(1, 1, s_phi * s_psi * s_gamma + c_phi * c_gamma);
	m2.SetElement(1, 2, s_phi * s_psi * c_gamma - c_phi * s_gamma);
	m2.SetElement(2, 0, -s_psi);
	m2.SetElement(2, 1, c_psi * s_gamma);
	m2.SetElement(2, 2, c_psi * c_gamma);

	m2 = m2.Transpose(); // 发射系到弹体系
	auto m = m2 * m1; // 发惯系到弹体系

	// 计算发惯系姿态角
	double psia = asin(-m[0][2]);
	double sin_phia = m[0][1] / cos(psia);
	double cos_phia = m[0][0] / cos(psia);
	double sin_gammaa = m[1][2] / cos(psia);
	double cos_gammaa = m[2][2] / cos(psia);
	double phia = ZJMath::calAngle2(sin_phia, cos_phia);
	double gammaa = ZJMath::calAngle2(sin_gammaa, cos_gammaa);
	return Vector3d(phia, psia, gammaa);
}

Vector3d ZJAngle::calAlphaBetaNu(Vector3d theta_sigma_gammav, Vector3d phi_psi_gamma)
{
	// 《远程火箭飞行动力学与制导》陈克俊.Eq.3-2-32
	double theta = theta_sigma_gammav[0];
	double sigma = theta_sigma_gammav[1];
	double phi = phi_psi_gamma[0];
	double psi = phi_psi_gamma[1];
	double gamma = phi_psi_gamma[2];

	double sin_beta = cos(phi - theta) * cos(sigma) * sin(psi) * cos(gamma)
		+ sin(phi - theta) * cos(sigma) * sin(gamma)
		- sin(sigma) * cos(psi) * cos(gamma);
	double beta = asin(sin_beta);

	double cos_alpha = (cos(sigma) * cos(theta - phi) * cos(psi) + sin(sigma) * sin(psi)) / cos(beta);
	double sin_alpha = -(-cos(psi) * sin(gamma) * sin(sigma) + cos(sigma) * (cos(gamma) * sin(theta - phi) + cos(theta - phi) * sin(gamma) * sin(psi))) / cos(beta);
	double alpha = ZJMath::calAngle2(sin_alpha, cos_alpha);

	double sin_nu = (cos(alpha) * cos(psi) * sin(gamma) - sin(alpha) * sin(psi)) / cos(sigma);
	double cos_nu = (cos(beta) * cos(psi) * cos(gamma) + sin(beta) * (sin(gamma) * cos(psi) * sin(alpha) + cos(alpha) * sin(psi))) / cos(sigma);
	double nu = ZJMath::calAngle2(sin_nu, cos_nu);

	return { alpha, beta, nu };
}
*/


Vector3d sign_sat(Vector3d x,double limit)
{
	for (int i = 0; i < 3; i++)
	{
		if (x[i] > limit)
			x[i] = 1;//sign_scalar(x[i])
		else if (x[i] < -limit)
		{
			x[i] = -1;
		}
		else{
			x[i] = x[i]/limit;
		}
	}
	return x;
}

Vector3d SIN_V3(double x)
{
	return Vector3d(1, 1, 1);//test
	//return Vector3d(sin(x[0]), sin(x[1]), sin(x[2]));
}
//定积分:一重积分是面积

//Vector3d Integration(double a, double b, Vector3d(*f)(Vector3d x), double interval)
//{
//	double len = b - a;
//	Vector3d t0 = Vector3d(a, a, a);
//	Vector3d dt = Vector3d(interval, interval, interval);
//	//int N = (int)(len / interval);//	int N = (int)len / interval; 这样会出错，先对len取整了，而不是除完之后的结果
//	int N = (len / interval);//	
//	Vector3d result{ 0.0,0.0,0.0 };
//	for (unsigned long i = 0; i < N; ++i) {
//		result += f(t0 + i * dt) * interval;
//	}
//	return result;
//}
Vector3d Integration(double a, double b, Vector3d (*f)(double x), double interval)
{
	double len = b - a;
	int N = len / interval;
	//int N = (int)(len / interval);//	int N = (int)len / interval; 这样会出错，先对len取整了，而不是除完之后的结果
	
	Vector3d result{0.0,0.0,0.0};
	for (unsigned long i = 0; i < N; ++i) {
		result += f(a + i * interval) * interval;
	}
	return result;
}

//计算二元函数的定积分，二重定积分是体积,由于上下限独立，可以拆成曲面积分

Vector3d DefiniteIntegration_2D(double x0, double x1, double y0, double y1,Vector3d(*f)(double x, double y), double interval)
{
	double lenx = x1 - x0, leny = y1 - y0;
	unsigned long Nx = lenx / interval, Ny = leny / interval;
	Vector3d result = Vector3d(0,0,0);
	for (unsigned long i = 0; i < Nx; ++i) {
		for (unsigned long j = 0; j < Ny; ++j) {
			result += f(x0 + i * interval, y0 + j * interval) * (interval * interval);
		}
	}
	return result;

}
//二次积分
Vector3d Variablelimit_Integrataion_2D_(double x0, double x1, Vector3d(*f)(double x), double interval)
{
	double lenx = x1 - x0;
	unsigned long Nx = lenx / interval;
	Vector3d result = Vector3d(0,0, 0);
	for (unsigned long i = 0; i < Nx; ++i) 
	{	
		result += Integration(x0, x0+ i * interval, *f, interval)*interval;
	}
	return result;
}


Vector3d Eye_V3(double x, double y)
{
	return Vector3d(1, 1, 1);
}

void Integration_test()
{  
//	cout << Integration(0, PI / 2, SIN_V3) << endl;
	//cout << DefiniteIntegration_2D(0, 2.5, 0, 1.4, Eye_V3) << endl;;
	cout << Variablelimit_Integrataion_2D_ (0,2, SIN_V3)<< endl;
} 
// 发射到发惯系转换矩阵
//Vector3d Launch2LaunchIv(double B0, double A0, double t, const Vector3d v)
//{
//	double X = cos(B0) * cos(A0);
//	double Y = sin(B0);
//	double Z = -cos(B0) * sin(A0);
//	double S = sin(earth_omega * t);
//	double C = cos(earth_omega * t);
//
//	Matrix3d m1;   
//	m1[0,0] = X * X * (1 - C) + C;
//	m1[0,1] = X * Y * (1 - C) - Z * S;
//	m1[0,2] = X * Z * (1 - C) + Y * S;
//
//	m1[1,0] = X * Y * (1 - C) + Z * S;
//	m1[1,1] = Y * Y * (1 - C) + C;
//	m1[1,2] = Y * Z * (1 - C) - X * S;
//
//	m1[2,0] = X * Z * (1 - C) - Y * S;
//	m1[2,1] = Y * Z * (1 - C) + X * S;
//	m1[2,2] = Z * Z * (1 - C) + C;
//
//	return m1 * v;
//}



Data_Processor::Data_Processor() {




}

vector<double>  Data_Processor::read_file(string root) {

	vector<double> indata;
	double temp;
	infile.open(root);
	//防止读文件出错
	if (!infile)
	{
		cout << "没有找到文件" << endl;
		return indata;
	}
	else
		cout << "成功找到文件" << endl;
	//
	outfile.open(".\out.txt");//需要输出时
	while (!infile.eof())
	{
		infile >> temp;//读入,输入流对象 必须通过>>运算符 才会迭代次数
		indata.push_back(temp);//保存到容器中
		//outfile<<temp//保存到输出文件中
	}
	infile.close();
	outfile.close();
	return indata;
}

void Data_Processor::plot_curve_2() {




}

void Data_Processor::plot_curve_3(vector<double> x, vector<double> y, vector<double> z) {

	vector<vector<double>> X, Y, Z;
	vector<double> x_row, y_row, z_row;
	for (size_t i = 0; i < x.size(); i++)
	{

		for (size_t i = 0; i < x.size(); i++)
		{

			x_row.push_back(x[i]);//如果显示"变量已被优化掉，因而不可用。"需要设置 编译器的优化选项
			y_row.push_back(y[i]);
			z_row.push_back(z[i]);
		
		}

		X.push_back(x_row);
		Y.push_back(x_row);
		Z.push_back(x_row);
	}

	//plt::plot_surface(x, y, z);//该函数的形参是 vector<vector<double>>
	plt::plot_surface(X, Y, Z); 
	plt::show();

}

void Data_Processor::test()
{
	vector<double> x(100), y(100), z(100);
	for (size_t i = 0; i < 100; i++)
	{
		x.at(i) = i * 0.1;
		y.at(i) = i * 0.1;
		z.at(i) = i * 0.1;
	}
	plot_curve_3(x, y, z);

}