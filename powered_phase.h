#pragma once
#include"missile.h"
#include"math_missile.h"
#include"aerodynamic.h"
#include"controller.h"
#include"guide.h"
using namespace Eigen;
using namespace std;

class powered_phase_stand_ballastic
{
public:
	//��׼����
	Vector3d  F,K1, K2, k_rel;//�������� �� �� ���������
	double G, P;//ע�������α�׼�����뵼����ͬ������
	Vector3d P_t;
	Vector3d PM,AM,Expected_M,Expected_Acc;
	Vector3d program_angle = Vector3d((Pi-0.001) / 2, 0, 0);
	Vector3d program_omega = Vector3d(0, 0, 0);
	missile mymissile; 
	Vector3d expected_Aero_F{ 0,0,0 };
	controller mycontroller;
	Guidance_law myguidar;

	double engine_m = 750, gas_m = 507, mass_flow = 56.3333, specific_impulse=3089.0,first_engine_worktime=9, second_engine_worktime = 44;//�˴���ֵ���벻�����캯������ִ�й��캯��
	double S=0.01, Pa=3.0e3;//��������ڵ������ѹǿ,���� todo
	
	Vector4d gas_vane_physical = Vector4d(0, 0, 0,0);
	Vector3d gas_vane_math = Vector3d(0, 0, 0);
	
	Aero_F aero_model;
	Vector3d x_c, x_g;//β�������ľ��� �� ����뾶
	bool Ifdisplay = 0;
	//��ģ�������� ����������
	powered_phase_stand_ballastic(missile m1,string areo_data= "F:/SM_3/", double P = 174e3, Vector3d K1 = Vector3d(0, 0, 0), \
		Vector3d K2 = Vector3d(0, 0, 0), Vector3d PM = Vector3d(0, 0, 0), Vector3d AM = Vector3d(0, 0, 0), Vector3d point = Vector3d(3000,10000, 0),double second_engine_worktime1=44);//ǰ����Ĭ�ϳ�ʼֵ ��������г�ʼֵ
	void updata_parameter();

	Vector3d F_sum();
	Vector3d M_sum(Vector3d gas_vane_deviation, Vector3d rudder_deviation);
	//ÿ�εĿ�����Ҳ����������ʽ��ͬ,��������д��math����

	//����������������,n=0��1��2����ͬ ��ѧ��
	Vector3d thrust(Vector3d gas_vane_deviation);//��ֱ�ӿ�����ѧ��
	//�߶Ⱥ�����
	double Gravity(double h, double m);

	// ����ǹ���ϵ��������������
	Vector3d K_1();
	Vector3d K_2();
	//����
	Vector3d thrust_Moment(Vector3d gas_vane_deviation);
	Vector3d Aero_Moment(Vector3d rudder_deviation);
	Vector3d Aero(Vector3d rudder_deviation);
	
	Vector3d Program_angle_turning(double time_current);
	//��һ�� ���ƶ�
	void boost_turning_phase(bool Ifdisplay);
	//ִ�л�����ȼ���潨ģ ����M->��ƫ
	Vector3d Cal_gas_vane_deviation(Vector3d Expected_torque);
	
	//�ڶ��� ������
	void powered_stage(bool Ifdisplay);

	//�Ƶ�ģ�������ģ������ӣ�a->F����->����
	Vector3d Calculate_Expected_control_angle(Vector3d F_sum);

};