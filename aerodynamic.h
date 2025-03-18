#pragma once
#include<string>
#include"math_missile.h"

const double ATMO_RHO0 = 1.225; //��ƽ������ܶ�(kg/m3)
const double ATMO_R0 = 6356.766; //��������
const double ATMO_P0 = 101325.0; //��׼����ѹ
const double ATMO_SEASURF_TEMP = 25;//��ƽ������(���϶�)
const double ATMO_AIRGAMMA = 1.4; //�������ȱ�
const double  ATMO_RGAS = 287.14;//�������峣��J/(kg*K)
const double  ATMO_K = 273.15;//���Ͽ����¶Ȼ���
//const double  ATMO_STBO = 5.6697e-12;//˹�ٷ�-������������


/****************************************************************************************/
//˵���������ܶȼ���ģ��
//���룺�߶�
//�������ǰ�����ܶȡ���ǰ�¶ȡ���ǰ����
//ע�⣺����߶ȵ�λ����
 /**************************************************************************************/
class Atmosphere_Model
{
public:  
	Atmosphere_Model();
	~Atmosphere_Model();

public:
	double phi;//��ǰ�����ܶ�(kg/m3)
	double Tempreture;//��ǰ�¶�(���϶�)
	double a_s; //��ǰ����(m/s)
	double  Kp_air;//����ѹǿ
	double m_CalAtmosphere(double IN_dHeight);//����ĵ�λ��km ,�������٣�����ʤ�ߣ��Ҹ����¶�
	double g_CalAtomspherePressure(double IN_dH);//���ݸ߶Ⱥ��¶ȼ������ѹǿ
};



/*********************************************************************************************************/
//˵����������ֵģ��
//���룺����״̬���߶ȡ���ƫ���ο�������ο����ȡ���ֵ��
//��������������������ء���������ϵ����������ϵ�������������ѹ
//ʹ��˵����ֱ�ӵ���CalAero�������ɣ��������ƫ����ϵ��������init_lapian���ٵ���CalAero����
//ע�⣺��ֵģ���õ��ǽǶȣ���Ҫת��һ��
/********************************************************************************************************/

class Aero_F
{
public:
	Aero_F();
	Aero_F(string filepath);
	~Aero_F();

public:
	double Ma; //�����
	double q_dyn; //��ѹ
	Vector3d Cx; //��ֵ�����������ϵ��
	Vector3d Cm; //��ֵ�������������ϵ��

	Vector3d F_Aero; //������
	Vector3d M_Aero; //��������


	Vector4d Rudder_physics; //�����ƫ��  X�ֶ� �����,�˴��������棬���ܸ� ȼ������
	Vector3d Rudder_math = Vector3d(0, 0, 0); //��ѧ��ƫ��
	bool iflapian; //��ƫ����
	Vector3d F_lapian; //����ϵ����ƫ��������ģ�������
	Vector3d M_lapian; //��������ϵ����ƫ��������ģ�������
	FlyStatus status; //����״̬

	Atmosphere_Model atmos_m;
private:
	//�ⲿ�����޸�һ��
	
	double phi; //�����ܶ�
	double a_s; //����
	double S_ref; //�ο����
	double L_ref; //�ο�����

	//�����ǲ�ֵ��
	bool iftable;
	//vector<double> C_alpha;
	//vector<double> C_beta;
	//vector<double> C_Ma;
	//vector<double> C_delta;
	vector<double> C_alpha = { -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30 };
	vector<double> C_beta = { -30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30 };
	vector<double> C_Ma = { 0.4,0.8,1.2,2.0,3.0,4.0,5.0,6.0,8.0,10.0 };
	vector<double> C_delta = { -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30 };
	vector<double> CN;
	vector<double> CA;
	vector<double> CM;
	vector<double> C_Mx_delta;
	vector<double> C_Mx;

public:
	void test(string filepath);
	void update_Aero(const FlyStatus& in_status, const Vector3d& in_Rudder);
	
	void init_lapian(Vector3d in_F_lapian, Vector3d in_M_lapian);

	void interpolation();
	Vector3d inverse_interpolation(Vector3d Cx, vector<double> Xs, vector<double> Ys, vector<double> Zs, Vector3d Rudder_math, Vector3d Ma);
	Vector3d F_current();
	Vector3d M_current();
	void output_FM();
	//void interpolation_F();

	void test_F();

	void readtable(vector <double> in_C_alpha, vector <double> in_C_beta, vector <double> in_C_Ma, vector <double> in_C_delta, vector <double> in_C_Mx_delta, vector <double> in_CN, vector <double> in_CA, vector <double> in_C_Mx, vector <double> in_C_Mz);

	void table(string filepath);

	void CalAero(const FlyStatus& in_status, const Vector3d& in_Rudder);

	double interp1(size_t n_Xs, vector<double> Xs, vector<double> Data, double X);

	double interp2(size_t n_Xs, vector<double> Xs, size_t n_Ys, vector<double> Ys, vector<double> Data, double X, double Y);

	double interp3(size_t n_Xs, vector<double> Xs, size_t n_Ys, vector<double> Ys, size_t n_Zs, vector<double> Zs, vector<double> Data, double X, double Y, double Z);

	double interp3_inverse(vector<double> Xs, vector<double> Ys, vector<double> Zs, vector<double> Data, double X, double Y, double Z);

};

