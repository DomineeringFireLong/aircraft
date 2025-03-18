//������_����ת�䣺��׼�������ɣ�����ʸ��+pid����
#include"powered_phase.h"
using namespace std;
powered_phase_stand_ballastic::powered_phase_stand_ballastic(missile m1,string aero_data, double P, Vector3d K1,\
	Vector3d K2, Vector3d PM, Vector3d AM, Vector3d point, double second_engine_worktime1):mymissile(m1),myguidar(point,0, second_engine_worktime1), mycontroller()
{
	//mymissile=m1;//������, ��Ķ����Ա �����ó�Ա��ʼ���б�����ֵ���ڹ��캯���︳ֵ���õ���Ȼ��Ĭ��!!!
	this->G = mymissile.m*9.8;
	this->P = P;
	this->K1 = K1;
	this->K2= K2;
	this->x_g << 3.18, 0, 0;
	this->x_c << 0.1775, 0, 0;
	this->aero_model = Aero_F(aero_data);
	this->aero_model.status = m1.Extract_status();//������mymissile���壬���׳���,�˴�������Ϣ����
	this->PM=PM;
	this->AM = AM;
	//mycontroller;//Ĭ�Ϲ��캯��,��Ĺ��캯������ڳ�ʼ���б�
	
}

double powered_phase_stand_ballastic::Gravity(double h,double m) {

	double gn = earth_g0;
	double r = earth_r0;
	//cout << "�߶�" << h << endl;
	//cout << "R" <<r << endl;
	double g = gn * pow(r / (r + h), 2);
	G = g * m;
	//cout << "�������ٶ�" << g << endl;
	return G;
}


Vector3d powered_phase_stand_ballastic::thrust(Vector3d gas_vane_deviation) {

	//����ʵ������ �ڶ�������ʱ�ٷ��� ,��ʱ��Ϊ��,����ΪX��
	//gas_vane_math[0] = ((gas_vane_physical[0] + gas_vane_physical[1]) - (gas_vane_physical[2] + gas_vane_physical[3]))/4;//������
	//gas_vane_math[1] = ((gas_vane_physical[0] + gas_vane_physical[3]) - (gas_vane_physical[1] + gas_vane_physical[2])) / 4;//ƫ����
	//gas_vane_math[2] = ((gas_vane_physical[0] + gas_vane_physical[1] + gas_vane_physical[2] + gas_vane_physical[3])) / 4;//��ת��
	//P_t = Vector3d(P, P * sqrt(2) / 2 * gas_vane_deviation[2], -1 * P * sqrt(2) / 2 * gas_vane_deviation[1]);
	
	
	//SM3�İڶ������+�ͣ�����ȼ���棩
	gas_vane_math= gas_vane_deviation;
	P = mass_flow * specific_impulse+S*(Pa- aero_model.atmos_m.g_CalAtomspherePressure(mymissile.pos[1]));
	P_t = Vector3d(P, P / 2 * gas_vane_deviation[2], -1 * P / 2 * gas_vane_deviation[1]);//���ط�����Ϊ�淽��Ϊ�˲�����ƫ���˶����ڶ���[1]=0
	return P_t;
}


//--------------------------------------------������ϵ/��ת�����ķǹ�����������----------------------------------------------------------
Vector3d powered_phase_stand_ballastic::K_1() {
	Matrix3d w_ = Cross_product(Vector3d(0,earth_omega,0));
	K1 = -2 * mymissile.m * w_ * mymissile.vel;
	return K1;
}
Vector3d powered_phase_stand_ballastic::K_2() {
	
	Matrix3d w_ = Cross_product(mymissile.attitude_omega);
	K2 = -2 * mass_flow * w_ * x_g;
	return K2;
}
//---------------------------------------------------------�����������ؼ���------------------------------------------------------------
Vector3d powered_phase_stand_ballastic::Aero(Vector3d rudder_deviation) {
	//��Ȼ�����Լ���status�ṹ�壬���ǲ���Ҫ���ݵ������и���
	aero_model.status= mymissile.Extract_status();//��ɾ����update����
	aero_model.CalAero(aero_model.status, rudder_deviation);//����ģ��Ƕ���˵�����״̬��Ϣ���������Բ�
	mymissile.A = aero_model.F_Aero;
	return mymissile.A;
}

Vector3d powered_phase_stand_ballastic::thrust_Moment(Vector3d gas_vane_deviation) {

	gas_vane_math = gas_vane_deviation;
	// X��
	//PM[0] = -1 * P* sqrt(2) * x_c[0] *  gas_vane_math[0];
	//PM[1] = -1 * P* sqrt(2) / 2 * x_g[0] *  gas_vane_math[1];
	//PM[2] = -1 * P* sqrt(2) / 2 * x_g[0] *  gas_vane_math[2];
	//ʮ����
	PM[0] = -1 * P  * x_c[0] * gas_vane_math[0];//btt���ƣ���������ת���Ҳ�����ƫ�������� [0] [1]��Ϊ0
	PM[1] = -1 * P / 2 * x_g[0] * gas_vane_math[1];
	PM[2] = -1 * P / 2 * x_g[0] * gas_vane_math[2];
	//��ֱ����Σ�ֻ���Ƹ���ͨ�����ɣ�[2]
	return PM;
}
Vector3d powered_phase_stand_ballastic::Aero_Moment(Vector3d rudder_deviation) {

	aero_model.status = mymissile.Extract_status();//��ɾ
	aero_model.CalAero(aero_model.status, rudder_deviation);
	AM = aero_model.M_Aero;
	return AM;
}

//����->����/�ٶ�����ϵ
Vector3d powered_phase_stand_ballastic::F_sum() {
	//���� ����->��������ϵ
	mymissile.G =Vector3d(0,-1*Gravity(mymissile.pos[1], mymissile.m),0);
	Vector3d G_trajectory = G2T(mymissile.G, mymissile.theta, mymissile.sigma);
	//cout << "G_trajectory " << G_trajectory << endl;
	//���� ����->��������ϵ
	mymissile.P = thrust(aero_model.Rudder_math);
	Vector3d P_trajectory = B2T(mymissile.P, mymissile.alpha, mymissile.beta);
	//cout << "P_trajectory" << P_trajectory << endl;
	mymissile.A = Aero(aero_model.Rudder_math);
	//cout << "mymissile.A" << mymissile.A << endl;
	K1 = K_1();// ������
	K2 = K_2();//���ӿ�����
	F=G_trajectory + P_trajectory + mymissile.A + K1 + K2;
	F[2]=0;//ǿ�Ʋ�����Ϊ0
	return F;
}
//������->��������ϵ
Vector3d powered_phase_stand_ballastic::M_sum(Vector3d gas_vane_deviation, Vector3d rudder_deviation) {
	//�������������أ��������������˿ط�������������
	PM = thrust_Moment(gas_vane_deviation);
	AM = Aero_Moment(rudder_deviation);
	return PM + AM;
}
//����״̬�������任
void powered_phase_stand_ballastic::updata_parameter() {

	mymissile.m -= mass_flow*mymissile.dt;
	aero_model.status = mymissile.Extract_status();

}

//90->60->60 
//
//Vector3d powered_phase_stand_ballastic::Program_angle_turning(double time_current) {
//	int turn_phase;
//	if (time_current <= 3)
//		turn_phase = 1;
//	else if (time_current <= 6)
//	{
//		turn_phase = 2;
//	}
//	else
//		turn_phase = 3;//�ֳ�����ת��
//	switch (turn_phase)
//	{
//	case 1:
//		program_angle = Vector3d(PI / 2, 0, 0);
//		program_omega = Vector3d(0, 0, 0);
//		break;
//	case 2:
//
//		program_angle = Vector3d((PI / 2) * exp(-1 * (time_current - 3) * log(3.0 / 2.0) / 3), 0, 0);//log���� ��eΪ�� log10 ��10Ϊ��
//		//cout << exp(-1.0 * (time_current - 1.0) * log(3.0 / 2.0) / 2.0) << endl;//ע��Ҫ��С����������int���� һֱ��1
//		program_omega = Vector3d((PI / 2) * (-1 * log(3 / 2) / 3) * exp(-1 * (time_current - 3) * log(3 / 2) /3), 0, 0);
//		break;
//	case 3:
//		program_angle = Vector3d(PI / 3, 0, 0);
//		program_omega = Vector3d(0, 0, 0);
//		break;
//	default:
//		break;
//	}
//	return program_angle;
//}


//90->60->45->45    1 3 0.03 km
Vector3d powered_phase_stand_ballastic::Program_angle_turning(double time_current) {
	int turn_phase;
	if (time_current <= 1)
		turn_phase = 1;
	else if (time_current <= 3)
	{
		turn_phase = 2;		
	}
	else if (time_current <= 6)
	{
		turn_phase = 3;
	}
	else
		turn_phase = 4;//�ֳ��Ķ�ת��
	switch (turn_phase)
	{
	case 1:
		program_angle = Vector3d(PI / 2, 0, 0);
		program_omega = Vector3d(0, 0, 0);
		break;
	case 2:
		
		program_angle = Vector3d((PI/2)*exp(-1*(time_current-1)*log(3.0/2.0)/2.0), 0, 0);//log���� ��eΪ�� log10 ��10Ϊ��
		//cout << exp(-1.0 * (time_current - 1.0) * log(3.0 / 2.0) / 2.0) << endl;//ע��Ҫ��С����������int���� һֱ��1
		program_omega = Vector3d((PI / 2)*(-1*log(3 / 2) / 2) * exp(-1 * (time_current - 1) * log(3 / 2) / 2), 0, 0);
		break;
	case 3:
		program_angle = Vector3d((PI / 2) / 54 * pow((time_current - 6), 2) + PI / 4, 0, 0);//todo
		program_omega = Vector3d((PI / 2) /27 * (time_current - 6) / 27, 0, 0);
		break;
	case 4:
		program_angle = Vector3d(PI / 4, 0, 0);
		program_omega = Vector3d(0, 0, 0);
		break;
	default:
		break;
	}
	return program_angle;
}


Vector3d powered_phase_stand_ballastic::Cal_gas_vane_deviation(Vector3d Expected_torque)
{
	//�����������ط�����
	gas_vane_math[0] = Expected_torque[0] / (-1 * P * x_c[0]);
	gas_vane_math[1] = Expected_torque[1] / (-1 * P / 2 * x_g[0]);
	gas_vane_math[2] = Expected_torque[2] / (-1 * P / 2 * x_g[0]);
	//�޷�
	gas_vane_math = sature(gas_vane_math, 7.5);//todo 
	return gas_vane_math;
}
//MK72����ʸ������
void powered_phase_stand_ballastic::boost_turning_phase(bool Ifdisplay) {
	int count = 0;
	while (mymissile.time < first_engine_worktime-0.01)
	{
		updata_parameter();
		Program_angle_turning(mymissile.time);
		Expected_M = mycontroller.control_instruct_eular(aero_model.status, program_angle, program_omega);
		//������������������
		gas_vane_math= Cal_gas_vane_deviation(mycontroller.M);//�����ȼ����ƫ
		aero_model.Rudder_math = Vector3d(0, 0, 0);//�������ȸ�0
		mymissile.M = M_sum(gas_vane_math, aero_model.Rudder_math);
		 //cout << "������ǵ��� " << mymissile.theta_dot << endl;
		//cout << "������� " << mymissile.theta << endl;
		//cout << "������� " << mymissile.theta << endl;
		//cout << "����ƫ�ǵ��� " << mymissile.sigma_dot << endl;
		//cout << "����ƫ�� " << mymissile.sigma << endl;

		//����3s���һ��,ȡ�����о������
		//cout << "������: " << mymissile.A << endl;
		mymissile.rotate(mycontroller.M);
		//todo��Ҫ������������Ķ���ܵ�ƫ�ƣ�Ȼ���ټ����ʵ�ʵ������˴�����ԭ���������ŵ��壬�������ص������ǣ�
		//����ִ�л���ģ�飬����->������ܽǶ�->ʵ�ʵ�����(ʵ�ʵ��Ƶ�ϵͳ�Ϳ���ϵͳ�����)
	
		
		mymissile.dynamic(F_sum());

		//cout <<"��̬��: " << mymissile.Eular_angle << endl;
	
		//if (int(mymissile.time * 100) % 100 == 0)//�����������е�ƫ��
		if (count%100 == 0&&(count>=100||count==0)) {
			
			//cout << mymissile.time << "ʱ�̣�������ǰλ��: " << mymissile.pos.transpose() << endl;
			cout << mymissile.time << "ʱ�̣������Ƕ�: " << program_angle.transpose() * RAD2DEG << endl;
			cout << "ʵ����̬��: " << mymissile.Eular_angle.transpose() * RAD2DEG << endl;

			cout << "������� " << mymissile.theta * RAD2DEG << endl;
			cout << "����ƫ�� " << mymissile.sigma * RAD2DEG << endl;

			//cout << "�������: " << mymissile.theta<< endl;
			//cout << mymissile.time << "ʱ�̣�����: " << mymissile.M.transpose() << endl;
			if (Ifdisplay == 1)
			{
				cout << "-----------------------------------------------" << endl;
				cout << "��������: " << mymissile.m << endl;
				cout << "������ǵ���: " << mymissile.theta_dot * RAD2DEG << endl;
				cout << "�������: " << mymissile.theta * RAD2DEG << endl;
				cout << "����ƫ�ǵ���: " << mymissile.sigma_dot * RAD2DEG << endl;
				cout << "����ƫ��: " << mymissile.sigma * RAD2DEG << endl;

				cout << "����: " << mymissile.G.transpose() << endl;
				cout << "��������: " << mymissile.A.transpose() << endl;
				cout << "����: " << mymissile.P.transpose() << endl;
				cout << "������1: " << K1.transpose() << endl;
				cout << "������2: " << K2.transpose() << endl;
				cout << "�������ܺ���" << F_sum().transpose() << endl;
				cout << "��̬��: " << mymissile.Eular_angle.transpose() * RAD2DEG << endl;
				cout << "������ٶ�: " << mymissile.attitude_omega.transpose() * RAD2DEG << endl;
				cout << "���ٶ�: " << mymissile.acc.transpose() << endl;
				cout << "�ٶ�: " << mymissile.vel.transpose() << endl;
				cout << "-----------------------------------------------" << endl;
			}

		}

		count+=1;
		mymissile.time += mymissile.dt;
	}
	cout <<"-----------"<< mymissile.time << "ʱ�̵�һ�����ƶν�����������ǰλ��" << "-----------" << mymissile.pos.transpose() << endl;
	cout <<"-------------------------------------------------------------------------------------" << endl;

}

//MK104˫��������������������������е���Ϳ��ƶ��棬����������

Vector3d powered_phase_stand_ballastic::Calculate_Expected_control_angle(Vector3d F_sum) {

	//��������ϵ
	mymissile.G = Vector3d(0, -1 * Gravity(mymissile.pos[1], mymissile.m), 0);
	Vector3d G_trajectory = G2T(mymissile.G, mymissile.theta, mymissile.sigma);
	
	mymissile.P = thrust(aero_model.Rudder_math);//���� ����->��������ϵ
	Vector3d P_trajectory = B2T(mymissile.P, mymissile.alpha, mymissile.beta);
	K1 = K_1();// ������
	K2 = K_2();//���ӿ�����

	expected_Aero_F = F_sum- G_trajectory- P_trajectory-K1-K2;
	aero_model.update_Aero(mymissile.Extract_status(), aero_model.Rudder_math);//������ƫ��ʱ��Ϊ��0
	Vector3d expected_Cx = expected_Aero_F / (aero_model.q_dyn * mymissile.S_ref);
	//���� ����ϵ�� �� ��ֵ�� ����

}


void powered_phase_stand_ballastic::powered_stage(bool Ifdisplay) {
	//��ȥ��һ��������MK72(750,507)���ڶ���������MK104(550 422)����
	//���µ�������
	mymissile.m = 1501 - engine_m;//750
	//���·���������
	engine_m = 550,
	gas_m = 422,
	mass_flow = 9.5925, 
	specific_impulse = 2294;
	second_engine_worktime = 44;//44
	myguidar.SM_3 = mymissile;//�ѵ����������뵽�Ƶ�ģ��
	//�����������º͵������ĵĸı�,֮���ٿ���todo


	//�Ƶ��Ϳ��Ʒ���Ҳ��ͬ
	int count = 0;
	while (myguidar.SM_3.time - first_engine_worktime < second_engine_worktime-0.01)
	{
	
		updata_parameter();
		myguidar.expect_point = Vector3d{11000,25000,0};
		Expected_Acc = myguidar.powered_stage(9999, myguidar.SM_3.time - first_engine_worktime);
		myguidar.SM_3.M = Vector3d(0, 0, 0);
		//	Expected_M = mycontroller.control_instruct_eular(aero_model.status, program_angle, program_omega);
		//������������������
		//gas_vane_math = Cal_gas_vane_deviation(mycontroller.M);//�����ȼ����ƫ
		//aero_model.Rudder_math = Vector3d(0, 0, 0);//�������ȸ�0
		//myguidar.SM_3.M = M_sum(gas_vane_math, aero_model.Rudder_math);
		//cout << "������ǵ��� " << myguidar.SM_3.theta_dot << endl;
	   //cout << "������� " << myguidar.SM_3.theta << endl;
	   //cout << "������� " << myguidar.SM_3.theta << endl;
	   //cout << "����ƫ�ǵ��� " << myguidar.SM_3.sigma_dot << endl;
	   //cout << "����ƫ�� " << myguidar.SM_3.sigma << endl;

	   //����3s���һ��,ȡ�����о������
	   //cout << "������: " << myguidar.SM_3.A << endl;
		myguidar.SM_3.rotate(mycontroller.M);
		//todo��Ҫ������������Ķ���ܵ�ƫ�ƣ�Ȼ���ټ����ʵ�ʵ������˴�����ԭ���������ŵ��壬�������ص������ǣ�
		//����ִ�л���ģ�飬����->������ܽǶ�->ʵ�ʵ�����(ʵ�ʵ��Ƶ�ϵͳ�Ϳ���ϵͳ�����)

		myguidar.SM_3.dynamic(Expected_Acc* myguidar.SM_3.m);
		//myguidar.SM_3.dynamic(F_sum());
		//cout <<"��̬��: " << myguidar.SM_3.Eular_angle << endl;
		//if (int(myguidar.SM_3.time * 100) % 100 == 0)//�����������е�ƫ��
		if (count % 100 == 0 && (count >= 100 || count == 0)) {

			cout << myguidar.SM_3.time << "ʱ�̣�������ǰλ��: " << myguidar.SM_3.pos.transpose() << endl;
			cout << "-----------------------------------------------" << endl;
			//cout << myguidar.SM_3.time << "ʱ�̣������Ƕ�: " << program_angle.transpose() * RAD2DEG << endl;
			//cout << "ʵ����̬��: " << myguidar.SM_3.Eular_angle.transpose() << endl;
			//cout << "�������: " << myguidar.SM_3.theta<< endl;
			//cout << myguidar.SM_3.time << "ʱ�̣�����: " << myguidar.SM_3.M.transpose() << endl;
			if (Ifdisplay == 1)
			{
			
	//			cout << "��������: " << myguidar.SM_3.m << endl;
		//		cout << "������ǵ���: " << myguidar.SM_3.theta_dot * RAD2DEG << endl;
				cout << "�������: " << myguidar.SM_3.theta * RAD2DEG << endl;
		//		cout << "����ƫ�ǵ���: " << myguidar.SM_3.sigma_dot * RAD2DEG << endl;
				cout << "����ƫ��: " << myguidar.SM_3.sigma * RAD2DEG << endl;

			/*	cout << "����: " << myguidar.SM_3.G.transpose() << endl;
				cout << "��������: " << myguidar.SM_3.A.transpose() << endl;
				cout << "����: " << myguidar.SM_3.P.transpose() << endl;
				cout << "������1: " << K1.transpose() << endl;
				cout << "������2: " << K2.transpose() << endl;*/
				cout << "�������ܺ���" << (Expected_Acc * myguidar.SM_3.m).transpose() << endl;//	cout << "�������ܺ���" << F_sum().transpose() << endl;
				//	cout << "��̬��: " << myguidar.SM_3.Eular_angle.transpose() * RAD2DEG << endl;
			//	cout << "������ٶ�: " << myguidar.SM_3.attitude_omega.transpose() * RAD2DEG << endl;
				cout << "���ٶ�: " << myguidar.SM_3.acc.transpose() << endl;
				cout << "�ٶ�: " << myguidar.SM_3.vel.transpose() << endl;
				cout << "-----------------------------------------------" << endl;
			}

		}

		count += 1;
		myguidar.SM_3.time += myguidar.SM_3.dt;
	}
	cout << myguidar.SM_3.time << "ʱ�̵ڶ������ƶν�����������ǰλ��" << myguidar.SM_3.pos.transpose() << endl;

	mymissile= myguidar.SM_3;
}


//powered_stage