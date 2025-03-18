#pragma once
#include"missile.h"
#include"math_missile.h"
#include<fstream>
#include<random>
//��ģ��д���ࣺ����ѧ�����������Ƶ������ơ�����ģ�⡢
//��ģ��֮�����Ϣ���ݣ��Լ���������Ĵ���
//��ͬ�εĶ���ѧ���̵Ķ�̬ �ں�


using namespace Eigen;
using namespace std;
class Guidance_law
{
	double t_final = 44;//Ԥ���Ƶ���ĩ��ʱ��
	double t_start = 0;//��ʼʱ��
	double  Q = 10;//����������ָ���ĩ����ļ�Ȩֵ
	double time_work=0;
public:
	Vector3d R, R_dot;//LOS����ϵ
	double r, r_dot, theta, theta_dot,sigma, sigma_dot;//LOS����ϵ
	Vector3d expect_point,expect_accelerate;
	missile SM_3,target;//������չ������Ϊmissile�ƶ�
	Vector3d ZCMD = Vector3d(10000, 80000, 0);
	Guidance_law(Vector3d point, double t_start1, double t_final1);

	//Vector3d expect_point;
	//powered_phase_guidance
	//���������ƶζβ��÷������������õ���
	
	//���������жβ���Ԥ���Ƶ�
	Vector3d powered_stage(double Q1, double  now_time);//���뵼�������� �Ƶ������������������
	//ĩ��
	void terminal_phase_guidance(missile m1, missile target);
	void terminal_phase_guidance(missile m1, Eigen::Vector3d expect_point1);
	Vector3d Integration_G(double a, double b, std::function<Vector3d(double)> f, double interval = 0.01);
	double get_theta_dot();
	Vector3d VL_Integrataion_2D(double x0, double x1, std::function<Vector3d(double)> f, double interval);
	virtual ~Guidance_law() {};//�˴�Ϊ���飻ͨ�������ָ����ɾ��������Ķ���ʱ���������������Ӧ������ġ�������ɾ��Ч�����޷�ʵ�֡�
	//missile slidemode(int k1 = 6, int k2 = 0.1, int k3 = 6, int k4 = 0.1, double theta_e = PI / 6, double eta_e = PI / 6, string root = "./3D_slidemode_guide.txt");
};



//�Ƶ����̳е�������������Ϣ,ĩ�Ƶ����ǻ���losϵ��
//class terminal_Guider:terminal_guidance {
//public:
//	Data_Buffer dandao_data;
//	missile SM_3;
//	terminal_Guider(missile m1, missile target);
//	terminal_Guider(missile m1, Vector3d expect_point);
//	//missile proportional_guide(int Ky = 4, int Kz = 4);
//	//missile slidemode(int k1 = 6, int k2 = 0.1, int k3 = 6, int k4 = 0.1, double theta_e = PI / 6, double eta_e = PI / 6, string root = "./Vertical_rising_section_guide.txt");
//
//};






