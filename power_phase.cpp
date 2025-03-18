//主动段_程序转弯：标准弹道生成，推力矢量+pid控制
#include"powered_phase.h"
using namespace std;
powered_phase_stand_ballastic::powered_phase_stand_ballastic(missile m1,string aero_data, double P, Vector3d K1,\
	Vector3d K2, Vector3d PM, Vector3d AM, Vector3d point, double second_engine_worktime1):mymissile(m1),myguidar(point,0, second_engine_worktime1), mycontroller()
{
	//mymissile=m1;//有问题, 类的对象成员 必须用成员初始化列表来赋值，在构造函数里赋值调用的仍然是默认!!!
	this->G = mymissile.m*9.8;
	this->P = P;
	this->K1 = K1;
	this->K2= K2;
	this->x_g << 3.18, 0, 0;
	this->x_c << 0.1775, 0, 0;
	this->aero_model = Aero_F(aero_data);
	this->aero_model.status = m1.Extract_status();//不能用mymissile定义，容易出错,此处传递信息对了
	this->PM=PM;
	this->AM = AM;
	//mycontroller;//默认构造函数,类的构造函数最好在初始化列表
	
}

double powered_phase_stand_ballastic::Gravity(double h,double m) {

	double gn = earth_g0;
	double r = earth_r0;
	//cout << "高度" << h << endl;
	//cout << "R" <<r << endl;
	double g = gn * pow(r / (r + h), 2);
	G = g * m;
	//cout << "重力加速度" << g << endl;
	return G;
}


Vector3d powered_phase_stand_ballastic::thrust(Vector3d gas_vane_deviation) {

	//考虑实际物理 摆动发动机时再分配 ,逆时针为正,以下为X型
	//gas_vane_math[0] = ((gas_vane_physical[0] + gas_vane_physical[1]) - (gas_vane_physical[2] + gas_vane_physical[3]))/4;//俯仰舵
	//gas_vane_math[1] = ((gas_vane_physical[0] + gas_vane_physical[3]) - (gas_vane_physical[1] + gas_vane_physical[2])) / 4;//偏航舵
	//gas_vane_math[2] = ((gas_vane_physical[0] + gas_vane_physical[1] + gas_vane_physical[2] + gas_vane_physical[3])) / 4;//滚转舵
	//P_t = Vector3d(P, P * sqrt(2) / 2 * gas_vane_deviation[2], -1 * P * sqrt(2) / 2 * gas_vane_deviation[1]);
	
	
	//SM3的摆动喷管是+型（不是燃气舵）
	gas_vane_math= gas_vane_deviation;
	P = mass_flow * specific_impulse+S*(Pa- aero_model.atmos_m.g_CalAtomspherePressure(mymissile.pos[1]));
	P_t = Vector3d(P, P / 2 * gas_vane_deviation[2], -1 * P / 2 * gas_vane_deviation[1]);//力矩方向作为舵方向，为了不产生偏航运动，摆动舵[1]=0
	return P_t;
}


//--------------------------------------------变质量系/自转产生的非惯性力科氏力----------------------------------------------------------
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
//---------------------------------------------------------气动力与力矩计算------------------------------------------------------------
Vector3d powered_phase_stand_ballastic::Aero(Vector3d rudder_deviation) {
	//虽然调用自己的status结构体，但是参数要根据导弹进行更新
	aero_model.status= mymissile.Extract_status();//可删，用update代替
	aero_model.CalAero(aero_model.status, rudder_deviation);//气动模块嵌入了导弹的状态信息，后续可以拆开
	mymissile.A = aero_model.F_Aero;
	return mymissile.A;
}

Vector3d powered_phase_stand_ballastic::thrust_Moment(Vector3d gas_vane_deviation) {

	gas_vane_math = gas_vane_deviation;
	// X型
	//PM[0] = -1 * P* sqrt(2) * x_c[0] *  gas_vane_math[0];
	//PM[1] = -1 * P* sqrt(2) / 2 * x_g[0] *  gas_vane_math[1];
	//PM[2] = -1 * P* sqrt(2) / 2 * x_g[0] *  gas_vane_math[2];
	//十字型
	PM[0] = -1 * P  * x_c[0] * gas_vane_math[0];//btt控制，不产生滚转，且不产生偏航，控制 [0] [1]均为0
	PM[1] = -1 * P / 2 * x_g[0] * gas_vane_math[1];
	PM[2] = -1 * P / 2 * x_g[0] * gas_vane_math[2];
	//垂直发射段，只控制俯仰通道即可，[2]
	return PM;
}
Vector3d powered_phase_stand_ballastic::Aero_Moment(Vector3d rudder_deviation) {

	aero_model.status = mymissile.Extract_status();//可删
	aero_model.CalAero(aero_model.status, rudder_deviation);
	AM = aero_model.M_Aero;
	return AM;
}

//合力->弹道/速度坐标系
Vector3d powered_phase_stand_ballastic::F_sum() {
	//重力 地面->弹道坐标系
	mymissile.G =Vector3d(0,-1*Gravity(mymissile.pos[1], mymissile.m),0);
	Vector3d G_trajectory = G2T(mymissile.G, mymissile.theta, mymissile.sigma);
	//cout << "G_trajectory " << G_trajectory << endl;
	//推力 弹体->弹道坐标系
	mymissile.P = thrust(aero_model.Rudder_math);
	Vector3d P_trajectory = B2T(mymissile.P, mymissile.alpha, mymissile.beta);
	//cout << "P_trajectory" << P_trajectory << endl;
	mymissile.A = Aero(aero_model.Rudder_math);
	//cout << "mymissile.A" << mymissile.A << endl;
	K1 = K_1();// 科氏力
	K2 = K_2();//附加科氏力
	F=G_trajectory + P_trajectory + mymissile.A + K1 + K2;
	F[2]=0;//强制侧向力为0
	return F;
}
//合力矩->弹体坐标系
Vector3d powered_phase_stand_ballastic::M_sum(Vector3d gas_vane_deviation, Vector3d rudder_deviation) {
	//重力不产生力矩，空气、推力、姿控发动机才有力矩
	PM = thrust_Moment(gas_vane_deviation);
	AM = Aero_Moment(rudder_deviation);
	return PM + AM;
}
//更新状态和质量变换
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
//		turn_phase = 3;//分成三段转弯
//	switch (turn_phase)
//	{
//	case 1:
//		program_angle = Vector3d(PI / 2, 0, 0);
//		program_omega = Vector3d(0, 0, 0);
//		break;
//	case 2:
//
//		program_angle = Vector3d((PI / 2) * exp(-1 * (time_current - 3) * log(3.0 / 2.0) / 3), 0, 0);//log函数 以e为底 log10 以10为底
//		//cout << exp(-1.0 * (time_current - 1.0) * log(3.0 / 2.0) / 2.0) << endl;//注意要用小数，否则是int类型 一直是1
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
		turn_phase = 4;//分成四段转弯
	switch (turn_phase)
	{
	case 1:
		program_angle = Vector3d(PI / 2, 0, 0);
		program_omega = Vector3d(0, 0, 0);
		break;
	case 2:
		
		program_angle = Vector3d((PI/2)*exp(-1*(time_current-1)*log(3.0/2.0)/2.0), 0, 0);//log函数 以e为底 log10 以10为底
		//cout << exp(-1.0 * (time_current - 1.0) * log(3.0 / 2.0) / 2.0) << endl;//注意要用小数，否则是int类型 一直是1
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
	//根据期望力矩反算气
	gas_vane_math[0] = Expected_torque[0] / (-1 * P * x_c[0]);
	gas_vane_math[1] = Expected_torque[1] / (-1 * P / 2 * x_g[0]);
	gas_vane_math[2] = Expected_torque[2] / (-1 * P / 2 * x_g[0]);
	//限幅
	gas_vane_math = sature(gas_vane_math, 7.5);//todo 
	return gas_vane_math;
}
//MK72推力矢量控制
void powered_phase_stand_ballastic::boost_turning_phase(bool Ifdisplay) {
	int count = 0;
	while (mymissile.time < first_engine_worktime-0.01)
	{
		updata_parameter();
		Program_angle_turning(mymissile.time);
		Expected_M = mycontroller.control_instruct_eular(aero_model.status, program_angle, program_omega);
		//控制器产生期望力矩
		gas_vane_math= Cal_gas_vane_deviation(mycontroller.M);//反算出燃气舵偏
		aero_model.Rudder_math = Vector3d(0, 0, 0);//气动舵先给0
		mymissile.M = M_sum(gas_vane_math, aero_model.Rudder_math);
		 //cout << "弹道倾角导数 " << mymissile.theta_dot << endl;
		//cout << "弹道倾角 " << mymissile.theta << endl;
		//cout << "弹道倾角 " << mymissile.theta << endl;
		//cout << "弹道偏角导数 " << mymissile.sigma_dot << endl;
		//cout << "弹道偏角 " << mymissile.sigma << endl;

		//程序3s输出一次,取整，有精度误差
		//cout << "气动力: " << mymissile.A << endl;
		mymissile.rotate(mycontroller.M);
		//todo，要根据力矩算出拍动喷管的偏移，然后再计算出实际的力（此处出错原因，推力沿着弹体，但推力矩的力不是）
		//加入执行机构模块，力矩->解算喷管角度->实际的力。(实际的制导系统和控制系统有耦合)
	
		
		mymissile.dynamic(F_sum());

		//cout <<"姿态角: " << mymissile.Eular_angle << endl;
	
		//if (int(mymissile.time * 100) % 100 == 0)//浮点数精度有点偏差
		if (count%100 == 0&&(count>=100||count==0)) {
			
			//cout << mymissile.time << "时刻，导弹当前位置: " << mymissile.pos.transpose() << endl;
			cout << mymissile.time << "时刻，期望角度: " << program_angle.transpose() * RAD2DEG << endl;
			cout << "实际姿态角: " << mymissile.Eular_angle.transpose() * RAD2DEG << endl;

			cout << "弹道倾角 " << mymissile.theta * RAD2DEG << endl;
			cout << "弹道偏角 " << mymissile.sigma * RAD2DEG << endl;

			//cout << "弹道倾角: " << mymissile.theta<< endl;
			//cout << mymissile.time << "时刻，力矩: " << mymissile.M.transpose() << endl;
			if (Ifdisplay == 1)
			{
				cout << "-----------------------------------------------" << endl;
				cout << "导弹质量: " << mymissile.m << endl;
				cout << "弹道倾角导数: " << mymissile.theta_dot * RAD2DEG << endl;
				cout << "弹道倾角: " << mymissile.theta * RAD2DEG << endl;
				cout << "弹道偏角导数: " << mymissile.sigma_dot * RAD2DEG << endl;
				cout << "弹道偏角: " << mymissile.sigma * RAD2DEG << endl;

				cout << "重力: " << mymissile.G.transpose() << endl;
				cout << "空气动力: " << mymissile.A.transpose() << endl;
				cout << "推力: " << mymissile.P.transpose() << endl;
				cout << "科氏力1: " << K1.transpose() << endl;
				cout << "科氏力2: " << K2.transpose() << endl;
				cout << "导弹所受合力" << F_sum().transpose() << endl;
				cout << "姿态角: " << mymissile.Eular_angle.transpose() * RAD2DEG << endl;
				cout << "弹体角速度: " << mymissile.attitude_omega.transpose() * RAD2DEG << endl;
				cout << "加速度: " << mymissile.acc.transpose() << endl;
				cout << "速度: " << mymissile.vel.transpose() << endl;
				cout << "-----------------------------------------------" << endl;
			}

		}

		count+=1;
		mymissile.time += mymissile.dt;
	}
	cout <<"-----------"<< mymissile.time << "时刻第一级助推段结束，导弹当前位置" << "-----------" << mymissile.pos.transpose() << endl;
	cout <<"-------------------------------------------------------------------------------------" << endl;

}

//MK104双推力固体固体续航发动机，具有弹翼和控制舵面，气动力控制

Vector3d powered_phase_stand_ballastic::Calculate_Expected_control_angle(Vector3d F_sum) {

	//弹道坐标系
	mymissile.G = Vector3d(0, -1 * Gravity(mymissile.pos[1], mymissile.m), 0);
	Vector3d G_trajectory = G2T(mymissile.G, mymissile.theta, mymissile.sigma);
	
	mymissile.P = thrust(aero_model.Rudder_math);//推力 弹体->弹道坐标系
	Vector3d P_trajectory = B2T(mymissile.P, mymissile.alpha, mymissile.beta);
	K1 = K_1();// 科氏力
	K2 = K_2();//附加科氏力

	expected_Aero_F = F_sum- G_trajectory- P_trajectory-K1-K2;
	aero_model.update_Aero(mymissile.Extract_status(), aero_model.Rudder_math);//空气舵偏暂时认为是0
	Vector3d expected_Cx = expected_Aero_F / (aero_model.q_dyn * mymissile.S_ref);
	//根据 气动系数 反 插值到 攻角

}


void powered_phase_stand_ballastic::powered_stage(bool Ifdisplay) {
	//抛去第一级发动机MK72(750,507)，第二级发动机MK104(550 422)工作
	//更新弹体质量
	mymissile.m = 1501 - engine_m;//750
	//跟新发动机参数
	engine_m = 550,
	gas_m = 422,
	mass_flow = 9.5925, 
	specific_impulse = 2294;
	second_engine_worktime = 44;//44
	myguidar.SM_3 = mymissile;//把导弹参数传入到制导模块
	//气动参数更新和导弹质心的改变,之后再考虑todo


	//制导和控制方法也不同
	int count = 0;
	while (myguidar.SM_3.time - first_engine_worktime < second_engine_worktime-0.01)
	{
	
		updata_parameter();
		myguidar.expect_point = Vector3d{11000,25000,0};
		Expected_Acc = myguidar.powered_stage(9999, myguidar.SM_3.time - first_engine_worktime);
		myguidar.SM_3.M = Vector3d(0, 0, 0);
		//	Expected_M = mycontroller.control_instruct_eular(aero_model.status, program_angle, program_omega);
		//控制器产生期望力矩
		//gas_vane_math = Cal_gas_vane_deviation(mycontroller.M);//反算出燃气舵偏
		//aero_model.Rudder_math = Vector3d(0, 0, 0);//气动舵先给0
		//myguidar.SM_3.M = M_sum(gas_vane_math, aero_model.Rudder_math);
		//cout << "弹道倾角导数 " << myguidar.SM_3.theta_dot << endl;
	   //cout << "弹道倾角 " << myguidar.SM_3.theta << endl;
	   //cout << "弹道倾角 " << myguidar.SM_3.theta << endl;
	   //cout << "弹道偏角导数 " << myguidar.SM_3.sigma_dot << endl;
	   //cout << "弹道偏角 " << myguidar.SM_3.sigma << endl;

	   //程序3s输出一次,取整，有精度误差
	   //cout << "气动力: " << myguidar.SM_3.A << endl;
		myguidar.SM_3.rotate(mycontroller.M);
		//todo，要根据力矩算出拍动喷管的偏移，然后再计算出实际的力（此处出错原因，推力沿着弹体，但推力矩的力不是）
		//加入执行机构模块，力矩->解算喷管角度->实际的力。(实际的制导系统和控制系统有耦合)

		myguidar.SM_3.dynamic(Expected_Acc* myguidar.SM_3.m);
		//myguidar.SM_3.dynamic(F_sum());
		//cout <<"姿态角: " << myguidar.SM_3.Eular_angle << endl;
		//if (int(myguidar.SM_3.time * 100) % 100 == 0)//浮点数精度有点偏差
		if (count % 100 == 0 && (count >= 100 || count == 0)) {

			cout << myguidar.SM_3.time << "时刻，导弹当前位置: " << myguidar.SM_3.pos.transpose() << endl;
			cout << "-----------------------------------------------" << endl;
			//cout << myguidar.SM_3.time << "时刻，期望角度: " << program_angle.transpose() * RAD2DEG << endl;
			//cout << "实际姿态角: " << myguidar.SM_3.Eular_angle.transpose() << endl;
			//cout << "弹道倾角: " << myguidar.SM_3.theta<< endl;
			//cout << myguidar.SM_3.time << "时刻，力矩: " << myguidar.SM_3.M.transpose() << endl;
			if (Ifdisplay == 1)
			{
			
	//			cout << "导弹质量: " << myguidar.SM_3.m << endl;
		//		cout << "弹道倾角导数: " << myguidar.SM_3.theta_dot * RAD2DEG << endl;
				cout << "弹道倾角: " << myguidar.SM_3.theta * RAD2DEG << endl;
		//		cout << "弹道偏角导数: " << myguidar.SM_3.sigma_dot * RAD2DEG << endl;
				cout << "弹道偏角: " << myguidar.SM_3.sigma * RAD2DEG << endl;

			/*	cout << "重力: " << myguidar.SM_3.G.transpose() << endl;
				cout << "空气动力: " << myguidar.SM_3.A.transpose() << endl;
				cout << "推力: " << myguidar.SM_3.P.transpose() << endl;
				cout << "科氏力1: " << K1.transpose() << endl;
				cout << "科氏力2: " << K2.transpose() << endl;*/
				cout << "导弹所受合力" << (Expected_Acc * myguidar.SM_3.m).transpose() << endl;//	cout << "导弹所受合力" << F_sum().transpose() << endl;
				//	cout << "姿态角: " << myguidar.SM_3.Eular_angle.transpose() * RAD2DEG << endl;
			//	cout << "弹体角速度: " << myguidar.SM_3.attitude_omega.transpose() * RAD2DEG << endl;
				cout << "加速度: " << myguidar.SM_3.acc.transpose() << endl;
				cout << "速度: " << myguidar.SM_3.vel.transpose() << endl;
				cout << "-----------------------------------------------" << endl;
			}

		}

		count += 1;
		myguidar.SM_3.time += myguidar.SM_3.dt;
	}
	cout << myguidar.SM_3.time << "时刻第二级助推段结束，导弹当前位置" << myguidar.SM_3.pos.transpose() << endl;

	mymissile= myguidar.SM_3;
}


//powered_stage