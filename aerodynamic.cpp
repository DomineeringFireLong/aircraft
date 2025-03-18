//代码格式默认是GB2312 ,不能改成其他格式，再改回来中文出错
#include "aerodynamic.h"

//气动插值部分
Aero_F::Aero_F()
{
	iflapian = false;
	iftable = false;
}

Aero_F::Aero_F(string filepath)
{
	iflapian = false;
	iftable = false;
	table(filepath);
	cout << "气动数据表导入完成" << endl;
}


Aero_F::~Aero_F()
{

}

void Aero_F::init_lapian(Vector3d in_F_lapian, Vector3d in_M_lapian)
{
	iflapian = true;
	F_lapian = in_F_lapian;
	M_lapian = in_M_lapian;
	cout << "气动参数已拉偏" << endl;
	return;
}
void Aero_F::test(string filepath)
{
	string file1 = "气动数据/"; //
	vector <string> filename;
	filename = { "CX.txt" };
	size_t iDataNum = 0;
	double dTemp = 0;//

	vector<vector<double>> output;
	output.resize(filename.size());

	ifstream input;
	
	for (size_t i = 0; i < filename.size(); i++)
	{
		iDataNum = 0;
		string file = filepath + file1 + filename[i];
		input.open(file);

		if (!input)
		{
			cout << "没有读到文件" << filename[i] << endl;
			return;
		}
		else
			cout << "成功读到" << filename[i] << "文件" << endl;
		while (!input.eof())
		{
			input >> dTemp;
			iDataNum++;
			//cout << "iDataNum "<< iDataNum <<": " << dTemp << endl;
			//input >> output[i][iDataNum];//锟斤拷锟斤拷岜拷锟斤拷锟斤拷为锟斤拷一维锟斤拷 元锟斤拷 锟斤拷锟斤拷锟斤拷小锟斤拷确锟斤拷
		}


		input.close();

		output[i].resize(iDataNum);

		input.open(file);

		for (size_t j = 0; j < iDataNum; j++)
		{
			input >> output[i][j];
		}

		input.close();
	}		

	//for (size_t j = 0; j < output[0].size(); j++)
	//{
	//	ofstream outf("./out_test.txt");
	//	outf << output[0][j];
	//	outf.close();
	//}
	cout << output.size() << endl;// ouput 锟斤拷一维 output[0]锟斤拷338
	CA = output[0];
	cout << CA.size() << endl;

}
void Aero_F::update_Aero(const FlyStatus& in_status, const Vector3d& in_Rudder)  //导弹的信息传递给气动模块，产生空气动力信息
{
	//先直接给数学气动舵，转化之后再加todo
	status = in_status;
	//Rudder_physics
	Rudder_math = in_Rudder;
	//Rudder_math[0] = -1 * (Rudder_physics[0] + Rudder_physics[1]) + (Rudder_physics[2] + Rudder_physics[3]);//俯仰舵
	//Rudder_math[1] = (Rudder_physics[0] + Rudder_physics[3]) + (Rudder_physics[1] + Rudder_physics[2]);//偏航舵
	//Rudder_math[2] = -1 * (Rudder_physics[0] + Rudder_physics[1]) + (Rudder_physics[2] + Rudder_physics[3]);//滚转舵
	atmos_m.m_CalAtmosphere(status.Pos[1]);

	//C_alpha={ -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30 };
	//C_beta = { -30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30 };
	//C_Ma= {0.4,0.8,1.2,2.0,3.0,4.0,5.0,6.0,8.0,10.0};
	//C_delta= { -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30 };

	phi = atmos_m.phi;
	a_s = atmos_m.a_s;

	S_ref = status.S_ref;
	L_ref = status.L_ref;

	Ma = status.v_dandao / a_s;
	q_dyn = 0.5 * phi * status.v_dandao * status.v_dandao;

	return;
}
// 
void Aero_F::interpolation()
{
	double alpha = status.alpha * RAD2DEG;
	double beta = status.beta * RAD2DEG;

	//Rudder_math[0] = -1 * (Rudder_physics[0] + Rudder_physics[1]) + (Rudder_physics[2] + Rudder_physics[3]);
	//Rudder_math[1] = -1 * (Rudder_physics[0] + Rudder_physics[3]) + (Rudder_physics[1] + Rudder_physics[2]);
	//Rudder_math[2] = -1 * (Rudder_physics[0] + Rudder_physics[1] + Rudder_physics[2] + Rudder_physics[3]);
	
	Cx[0]= interp3(C_alpha.size(), C_alpha, C_delta.size(), C_delta, C_Ma.size(), C_Ma, CA, alpha, Rudder_math[0], Ma); 
	if (Cx[0] < 0)
		Cx[0] = 0;
	Cx[1] = interp3(C_alpha.size(), C_alpha, C_delta.size(), C_delta, C_Ma.size(), C_Ma, CN, alpha, Rudder_math[1], Ma); 
	Cx[2] = -interp3(C_beta.size(), C_beta, C_delta.size(), C_delta, C_Ma.size(), C_Ma, CN, beta, Rudder_math[2], Ma); 

	Cm[0] = 0;
	Cm[1] = interp3(C_beta.size(), C_beta, C_delta.size(), C_delta, C_Ma.size(), C_Ma, CM, beta, Rudder_math[1], Ma); 
	Cm[2] = interp3(C_alpha.size(), C_alpha, C_delta.size(), C_delta, C_Ma.size(), C_Ma, CM, alpha, Rudder_math[2], Ma); 

	//if (iflapian)
	//{
	//	Cx[0] += 0.1 * abs(Cx[0]) * F_lapian.x;
	//	Cx[1] += 0.1 * abs(Cx[1]) * F_lapian.y;
	//	Cx[2] += 0.1 * abs(Cx[2]) * F_lapian.z;
	//	Cm.x += 0.2 * abs(Cm.x) * M_lapian.x;
	//	Cm.y += 0.2 * abs(Cm.y) * M_lapian.y;
	//	Cm.z += 0.2 * abs(Cm.z) * M_lapian.z;
	//}
	//return;
}

/*
void Aero_F::interpolation_F()
{
	double alpha = status.alpha * RAD2DEG;
	double beta = status.beta * RAD2DEG;
	

	Rudder_math[0] = -1 * (Rudder_physics[0] + Rudder_physics[1]) + (Rudder_physics[2] + Rudder_physics[3]);
	Rudder_math[1] = (Rudder_physics[0] + Rudder_physics[3]) + (Rudder_physics[1] + Rudder_physics[2]);
	Rudder_math[2] = -1 * (Rudder_physics[0] + Rudder_physics[1]) + (Rudder_physics[2] + Rudder_physics[3]);

	///////////////////////////////////////////////////////////////////////////////
	//遍历法
	///////////////////////////////////////////////////////////////////////////////
	//Cx[0] = -interp3(C_alpha.size(), C_alpha, C_Ma.size(), C_Ma, C_delta.size(), C_delta, CA, alpha, Ma, delta); //阻力（负值）
	//Cx[0]= -interp3(C_alpha.size(), C_alpha, C_Ma.size(), C_Ma, C_delta.size(), C_delta, CA, alpha, Ma, rubber.z); //阻力（负值）
	//Cx[1]= interp3(C_alpha.size(), C_alpha, C_Ma.size(), C_Ma, C_delta.size(), C_delta, CN, alpha, Ma, rubber.z); //升力
	//Cx[2] = -interp3(C_beta.size(), C_beta, C_Ma.size(), C_Ma, C_delta.size(), C_delta, CN, beta, Ma, rubber.y); //侧向力

	////暂时没有用到滚转，将滚转系数锁定为0
	////Cm.x = BDInterpolation::interp2(n_5, C_Mx_delta, n_3, C_Ma, C_Mx, delta, Ma);
	//Cm.x = 0;//滚转
	//Cm.y = interp3(C_beta.size(), C_beta,  C_Ma.size(), C_Ma, C_delta.size(), C_delta, CM,  beta, Ma, rubber.y); //偏航
	//Cm.z = interp3(C_alpha.size(), C_alpha, C_Ma.size(), C_Ma, C_delta.size(), C_delta, CM,  alpha, Ma, rubber.z); //俯仰


	///////////////////////////////////////////////////////////////////////////////
	//二分法
	///////////////////////////////////////////////////////////////////////////////
	Cx.x = -interp3_breeze(C_alpha, C_Ma, C_delta, CA, alpha, Ma, rubber.z); //阻力（负值）
	Cx.y = interp3_breeze(C_alpha, C_Ma, C_delta, CN, alpha, Ma, rubber.z); //升力
	Cx.z = -interp3_breeze(C_beta, C_Ma, C_delta, CN, beta, Ma, rubber.y); //侧向力

	//暂时没有用到滚转，将滚转系数锁定为0
	//Cm.x = BDInterpolation::interp2(n_5, C_Mx_delta, n_3, C_Ma, C_Mx, delta, Ma);
	Cm.x = 0;//滚转
	Cm.y = interp3_breeze(C_beta, C_Ma, C_delta, CM, beta, Ma, rubber.y); //偏航
	Cm.z = interp3_breeze(C_alpha, C_Ma, C_delta, CM, alpha, Ma, rubber.z); //俯仰


	//后续做蒙塔卡洛仿真需要将气动系数拉偏


	if (iflapian)
	{
		Cx[0] += 0.1 * abs(Cx[0]) * F_lapian.x;
		Cx[1] += 0.1 * abs(Cx[1]) * F_lapian.y;
		Cx[2] += 0.1 * abs(Cx[2]) * F_lapian.z;
		Cm.x += 0.2 * abs(Cm.x) * M_lapian.x;
		Cm.y += 0.2 * abs(Cm.y) * M_lapian.y;
		Cm.z += 0.2 * abs(Cm.z) * M_lapian.z;
	}
	return;
}
*/

void Aero_F::test_F()
{
	double alpha = status.alpha * RAD2DEG;
	double beta = status.beta * RAD2DEG;
	//Rudder_math[0] = -1 * (Rudder_physics[0] + Rudder_physics[1]) + (Rudder_physics[2] + Rudder_physics[3]);
	//Rudder_math[1] = (Rudder_physics[0] + Rudder_physics[3]) + (Rudder_physics[1] + Rudder_physics[2]);
	//Rudder_math[2] = -1 * (Rudder_physics[0] + Rudder_physics[1]) + (Rudder_physics[2] + Rudder_physics[3]);

	if (Ma < 1)
	{
		Cx[0] = -(0.00061 * alpha * alpha + 0.33);
	}
	else
	{
		Cx[0] = -0.3;
	}
	Cx[1] = 0.00014 * alpha * alpha * alpha + 0.2 * alpha + 0.022 * Rudder_math[2];
	Cx[2] = -(0.00014 * beta * beta * beta + 0.2 * beta + 0.022 * Rudder_math[1]);
	Cm[0] = 0;
	Cm[1] = -1.6e-5 * beta * beta * beta - 0.018 * beta - 0.03 * Rudder_math[1];
	Cm[2] = -1.6e-5 * alpha * alpha * alpha - 0.018 * alpha - 0.03 * Rudder_math[2];

	if (iflapian)
	{
		Cx[0] += 0.1 * abs(Cx[0]) * F_lapian[0];
		Cx[1] += 0.1 * abs(Cx[1]) * F_lapian[1];
		Cx[2] += 0.1 * abs(Cx[2]) * F_lapian[2];
		Cm[0] += 0.2 * abs(Cm[0]) * M_lapian[0];
		Cm[1] += 0.2 * abs(Cm[1]) * M_lapian[1];
		Cm[2] += 0.2 * abs(Cm[2]) * M_lapian[2];
	}
}

//此函数暂未使用
void Aero_F::readtable(vector <double> in_C_alpha, vector <double> in_C_beta, vector <double> in_C_Ma, vector <double> in_C_delta, vector <double> in_C_Mx_delta, vector <double> in_CN, vector <double> in_CA, vector <double> in_C_Mx, vector <double> in_C_Mz)
{
	CN = in_CN;
	CA = in_CA;
	C_Mx_delta = in_C_Mx_delta;
	C_Mx = in_C_Mx;
	CM = in_C_Mz;
	C_alpha = in_C_alpha;
	C_beta = in_C_beta;
	C_Ma = in_C_Ma;
	C_delta = in_C_delta;
	cout << "气动参数表导入完成" << endl;

	iftable = true;

	return;
}

//读取气动表
void Aero_F::table(string filepath)
{
	string file1 = "气动数据/"; //子文件夹
	vector <string> filename;
	//filename = { "CN.txt","CA.txt","CM.txt","C_alpha.txt","C_beta.txt", "C_Ma.txt", "C_delta.txt" };
	filename = { "CX.txt","CY.txt","CM.txt"};
	int iDataNum = 0;//数据个数
	double dTemp = 0;//缓存

	vector<vector<double>> output;//
	output.resize(filename.size());//[3][130]

	ifstream input;
	//有几个文件夹循环几次，output是一个二维的向量容器，第一个变量连得容器存储该文件夹的数据
	for (size_t i = 0; i < filename.size(); i++)
	{
		iDataNum = 0;
		string file = filepath + file1 + filename[i];
		input.open(file);

		if (!input)
		{
			cout << "没有找到文件" << filename[i] << endl;
			return;
		}
		else
			cout << "成功找到" << filename[i] << "文件" << endl;
		//可能出现死循环现象
		while (!input.eof())
		{
			input >> dTemp;
			
			//if (iDataNum % 169 == 0)
			//{
			//	cout << "第" << iDataNum / 169 << "组数据读入" << dTemp << endl;
			//}
			iDataNum++;
		}
		cout<<filename[i]<<"成功读入 " << endl;
		input.close();

		output[i].resize(iDataNum);

		input.open(file);
		for (size_t j = 0; j < iDataNum; j++)
		{

			input >> output[i][j];
		}

		input.close();
	}

	CA = output[0];
	CN = output[1];
	CM = output[2];
	//C_alpha = output[3];
	//C_beta = output[4];
	//C_Ma = output[5];
	//C_delta = output[6];

	iftable = true;
	return;

}

Vector3d Aero_F::F_current() {

	return F_Aero;
}
Vector3d Aero_F::M_current() {

	return M_Aero;
}

void Aero_F::output_FM() {
	cout << F_Aero.transpose() << endl;
	cout << M_Aero.transpose() << endl;
}

void Aero_F::CalAero(const FlyStatus& in_status, const Vector3d& in_Rudder)
{
	update_Aero(in_status, in_Rudder);
	iftable = 1;

	if (!iftable)
	{
		//cout << "气动参数表没有导入，请导入气动参数表" << endl;
		test_F();
		F_Aero = q_dyn * S_ref * Cx;
		M_Aero = q_dyn * S_ref * L_ref * Cm;
		return;
	}

	interpolation();
	Cx[0] = -Cx[0];//阻力系数是正，加符号
	F_Aero = q_dyn * S_ref * Cx; //弹体系下的气动力
	M_Aero = q_dyn * S_ref * L_ref * Cm;
	//cout << q_dyn << endl;
	//cout << S_ref << endl;
	//cout << L_ref << endl;
	//cout << Cm << endl;
	return;
}

//double interp1(int n_Xs, double *Xs, double *Data, double X)
//{
//	// Step1：确定Xs向位置
//	int id;
//	for (int i = 0; i != n_Xs - 1; i++)
//	{
//		id = i;
//		if (X >= Xs[i] && X <= Xs[i + 1])
//			break;
//		else
//			continue;
//	}
//
//	// Step2：Xs向插值
//	double ratio;
//
//	if (X < Xs[0])
//	{
//		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
//		return Data[0] - ratio * (Data[1] - Data[0]);
//	}
//	else if (X >= Xs[n_Xs - 1])
//	{
//		ratio = (X - Xs[n_Xs - 1]) / (Xs[n_Xs - 1] - Xs[n_Xs - 2]);
//		return Data[n_Xs - 1] + ratio * (Data[n_Xs - 1] - Data[n_Xs - 2]);
//	}
//	else
//	{
//		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
//		return Data[id] + (Data[id + 1] - Data[id])*ratio;
//	}
//}


//插值函数部分
double Aero_F::interp1(size_t n_Xs, vector<double> Xs, vector<double> Data, double X)
{
	// Step1：确定Xs向位置
	size_t id = 0;
	for (size_t i = 0; i != n_Xs - 1; i++)
	{
		id = i;
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		return Data[0] - ratio * (Data[1] - Data[0]);
	}
	else if (X >= Xs[n_Xs - 1])
	{
		ratio = (X - Xs[n_Xs - 1]) / (Xs[n_Xs - 1] - Xs[n_Xs - 2]);
		return Data[n_Xs - 1] + ratio * (Data[n_Xs - 1] - Data[n_Xs - 2]);
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		return Data[id] + (Data[id + 1] - Data[id]) * ratio;
	}

}

double Aero_F::interp2(size_t n_Xs, vector<double> Xs, size_t n_Ys, vector<double> Ys, vector<double> Data, double X, double Y)
{
	// Step1：确定Xs向位置
	size_t id = 0;
	for (size_t i = 0; i != n_Xs - 1; i++)
	{
		id = i;
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	vector<double> Y_Data;
	Y_Data.resize(n_Ys);
	/*double *Y_Data = new double[n_Ys];*/
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		for (size_t i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs] - ratio * (Data[i * n_Xs + 1] - Data[i * n_Xs]);
	}
	else if (X >= Xs[n_Xs - 1])
	{
		ratio = (X - Xs[n_Xs - 1]) / (Xs[n_Xs - 1] - Xs[n_Xs - 2]);
		for (size_t i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + n_Xs - 1] + ratio * (Data[i * n_Xs + n_Xs - 1] - Data[i * n_Xs + n_Xs - 2]);
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		for (size_t i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + id] + (Data[i * n_Xs + id + 1] - Data[i * n_Xs + id]) * ratio;
	}

	// Step3：确定Ys向位置
	for (size_t i = 0; i != n_Ys - 1; i++)
	{
		id = i;
		if (Y >= Ys[i] && Y <= Ys[i + 1])
			break;
		else
			continue;
	}

	// Step4：Ys向插值
	double	ResultVal;
	if (Y < Ys[0])
	{
		ratio = (Ys[0] - Y) / (Ys[1] - Ys[0]);
		ResultVal = Y_Data[0] - ratio * (Y_Data[1] - Y_Data[0]);
	}

	else if (Y >= Ys[n_Ys - 1])
	{
		ratio = (Y - Ys[n_Ys - 1]) / (Ys[n_Ys - 1] - Ys[n_Ys - 2]);
		ResultVal = Y_Data[n_Ys - 1] + ratio * (Y_Data[n_Ys - 1] - Y_Data[n_Ys - 2]);
	}
	else
	{
		ratio = (Y - Ys[id]) / (Ys[id + 1] - Ys[id]);
		ResultVal = Y_Data[id] + (Y_Data[id + 1] - Y_Data[id]) * ratio;
	}
	//delete[] Y_Data;
	return ResultVal;
}

double Aero_F::interp3(size_t n_Xs, vector<double> Xs, size_t n_Ys, vector<double> Ys, size_t n_Zs, vector<double> Zs, vector<double> Data, double X, double Y, double Z)
{
	// Step1：Xs、Ys向二维插值
	vector<double> Z_Data;
	Z_Data.resize(n_Zs);

	vector<double> data1;
	data1.resize(n_Xs * n_Ys);
	
	for (size_t i = 0; i != n_Zs; i++)
	{
		for (size_t j = 0; j < data1.size(); j++)
		{
			data1[j] = Data[i * n_Xs * n_Ys + j];
		}
		Z_Data[i] = interp2(n_Xs, Xs, n_Ys, Ys, data1, X, Y);
	}

	// Step2：确定Zs向位置
	size_t id;
	for (size_t i = 0; i != n_Zs - 1; i++)
	{
		id = i;
		if (Z >= Zs[i] && Z <= Zs[i + 1])
			break;
		else
			continue;
	}

	// Step3：Zs向插值
	double ResultVal;
	double ratio;
	if (Z < Zs[0])
	{
		ratio = (Zs[0] - Z) / (Zs[1] - Zs[0]);
		ResultVal = Z_Data[0] - ratio * (Z_Data[1] - Z_Data[0]);
	}
	else if (Z >= Zs[n_Zs - 1])
	{
		ratio = (Z - Zs[n_Zs - 1]) / (Zs[n_Zs - 1] - Zs[n_Zs - 2]);
		ResultVal = Z_Data[n_Zs - 1] + ratio * (Z_Data[n_Zs - 1] - Z_Data[n_Zs - 2]);
	}

	else
	{
		ratio = (Z - Zs[id]) / (Zs[id + 1] - Zs[id]);
		ResultVal = Z_Data[id] + (Z_Data[id + 1] - Z_Data[id]) * ratio;
	}
	//delete[] Z_Data;
	return ResultVal;
}

//todo
double Aero_F::interp3_inverse(vector<double> Xs, vector<double> Ys, vector<double> Zs, vector<double> Data, double X, double Y, double Z)
{
	vector<double> Z_Data;
	Z_Data.resize(Zs.size());

	vector<double> data1;
	data1.resize(Xs.size() * Ys.size());

	for (size_t i = 0; i != Zs.size(); i++)
	{
		for (size_t j = 0; j < data1.size(); j++)
		{
			data1[j] = Data[i * Xs.size() * Ys.size() + j];
		}
		Z_Data[i] = interp2(Xs.size(), Xs, Ys.size(), Ys, data1, X, Y);
	}

	return 0.0;
}

Vector3d Aero_F::inverse_interpolation(Vector3d Cx, vector<double> Xs, vector<double> Ys, vector<double> Zs, Vector3d Rudder_math, Vector3d Ma)
{
	vector<double> Z_Data;//存储double数据的容器
	Z_Data.resize(Zs.size());

}

//大气密度计算部分
Atmosphere_Model::Atmosphere_Model()
{
	phi = ATMO_RHO0;
	a_s = 340.0;
	Tempreture = ATMO_SEASURF_TEMP + ATMO_K;
}

Atmosphere_Model::~Atmosphere_Model()
{

}

double Atmosphere_Model::m_CalAtmosphere(double IN_dHeight)//m
{
	
	double dW = 0;
	double dHeightKM = 0;
	double dH = 0;

	dHeightKM = IN_dHeight/1000 ;
	dH = dHeightKM / (1 + dHeightKM / ATMO_R0);

	//计算大气密度
	if (dHeightKM >= 0 && dHeightKM <= 11.0191)
	{
		dW = 1 - dH / 44.3308;
		phi = pow(dW, 4.2559) * ATMO_RHO0;
	}
	else if (dHeightKM > 11.0191 && dHeightKM <= 20.0631)
	{
		dW = exp((14.9647 - dH) / 6.3416);
		phi = 1.5898 * 0.1 * dW * ATMO_RHO0;
	}
	else if (dHeightKM > 20.0631 && dHeightKM <= 32.3619)
	{
		dW = 1 + (dH - 24.0921) / 221.552;
		phi = 3.2722 * 0.01 * pow(dW, -35.1629) * ATMO_RHO0;
	}
	else if (dHeightKM > 32.3619 && dHeightKM <= 47.3501)
	{
		dW = 1 + (dH - 39.7499) / 89.4107;
		phi = 3.2618 * 0.001 * pow(dW, -13.2011) * ATMO_RHO0;
	}
	else if (dHeightKM > 47.3501 && dHeightKM <= 51.4125)
	{
		dW = exp((48.6252 - dH) / 7.9223);
		phi = 9.4920 * 0.0001 * dW * ATMO_RHO0;
	}
	else if (dHeightKM > 51.4125 && dHeightKM <= 71.8020)
	{
		dW = 1 - (dH - 59.4390) / 88.2218;
		phi = 2.5280 * 0.0001 * pow(dW, 11.2011) * ATMO_RHO0;
	}
	else if (dHeightKM > 71.8020 && dHeightKM <= 86.0000)
	{
		dW = 1 - (dH - 78.0303) / 100.2950;
		phi = 1.7632 * 0.00001 * pow(dW, 16.0816) * ATMO_RHO0;
	}
	else if (dHeightKM > 86.0000 && dHeightKM <= 91.0000)
	{
		dW = exp((87.2848 - dH) / 5.4700);
		phi = 3.6411 * 0.000001 * dW * ATMO_RHO0;
	}
	else
	{
		phi = 0;
	}


	//计算气温
	if (dHeightKM < 0)
	{
		Tempreture = ATMO_SEASURF_TEMP + ATMO_K;
	}
	else if (dHeightKM <= 11.0191)
	{
		Tempreture = 288.15 * (1 - dH / 44.3308);
	}
	else if (dHeightKM <= 20.0631)
	{
		Tempreture = 216.650;
	}
	else if (dHeightKM <= 32.1619)
	{
		Tempreture = 221.552 * (1 + (dH - 24.9021) / 221.552);
	}
	else if (dHeightKM > 32.3619 && dHeightKM <= 47.3501)
	{
		Tempreture = 250.350 * (1 + (dH - 39.7499) / 89.4107);
	}
	else if (dHeightKM > 47.3501 && dHeightKM <= 51.4125)
	{
		Tempreture = 270.650;
	}
	else if (dHeightKM > 51.4125 && dHeightKM <= 71.8020)
	{
		Tempreture = 247.021 * (1 - (dH - 59.4390) / 88.2218);
	}
	else if (dHeightKM > 71.8020 && dHeightKM <= 86.0000)
	{
		Tempreture = 200.590 * (1 - (dH - 78.0303) / 100.2950);
	}
	else if (dHeightKM > 86.0000 && dHeightKM <= 91.0000)
	{
		Tempreture = 186.8700;
	}

	//计算声速
	a_s = sqrt(ATMO_AIRGAMMA * ATMO_RGAS * Tempreture);

	return a_s;
}

double Atmosphere_Model::g_CalAtomspherePressure(double IN_dH)//m
{
	m_CalAtmosphere(IN_dH);
	const double ATMO_P0 = 101325.0;                                 // 标准大气压(pa)
	const double ATMO_K = 273.15;                                   // 摄氏开氏温度互换
	Kp_air = ATMO_P0*0.001 / exp(IN_dH*1000 / (18400 * (1 + Tempreture / ATMO_K)));
	return  Kp_air;
}
