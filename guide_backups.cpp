/*

#include<iostream>
#include<cmath>
#include <Eigen/Dense>
#include <Eigen/Core>



using namespace std;


terminal_phase_guidance::terminal_phase_guidance(missile m1, missile target1)//ֱ���к��в����ĵ�����Ŀ�����й������
{
	//ֻ�ǳ�ʼ�������ܸ��£�vel���ڵ�������ϵ���׸��£������ٶ�����ϵ
	this->SM_3 = m1;
	this->target = target1;
	this->r_dot = 0;
	this->theta_dot = 0;
	this->fai_dot = 0;
	this->r = sqrt(pow(target.pos[0] - SM_3 .pos[0], 2) + pow(target.pos[1] - SM_3 .pos[1], 2) + pow(target.pos[2] - SM_3 .pos[2], 2));
	//1.ע�⣺��������ϵ���ڵ����ƫ�ǡ���� ���㹫ʽ,���ݳ�ʼ�ٶ�ʸ�������ʼ�ǣ���Ŀ����ϵ�Ƕ������λ�Ʒ��򣬵����ĽǶ����ٶȷ���Ķ���
	this->theta = atan((this->target.y - this->SM_3 .y) / sqrt(pow(this->target.x - this->SM_3 .x, 2) + pow(this->target.z - this->SM_3 .z, 2)));
	this->fai = atan(-(this->target.z) / this->target.x);

}

terminal_phase_guidance::terminal_phase_guidance(missile m1, Eigen::Vector3d target)
{
	//ֻ�ǳ�ʼ�������ܸ��£�vel���ڵ�������ϵ���׸��£������ٶ�����ϵ
	this->SM_3 = m1;
	this->expect_point = target;
	this->r_dot = 0;
	this->theta_dot = 0;
	this->fai_dot = 0;
	this->r = sqrt(pow(expect_point[0] - expect_point[0], 2) + pow(expect_point[1] - SM_3.pos[1], 2) + pow(expect_point[2] - SM_3.pos[2], 2));
	//1.ע�⣺��������ϵ���ڵ����ƫ�ǡ���� ���㹫ʽ,���ݳ�ʼ�ٶ�ʸ�������ʼ�ǣ���Ŀ����ϵ�Ƕ������λ�Ʒ��򣬵����ĽǶ����ٶȷ���Ķ���
	this->theta = atan((this->expect_point[1] - this->SM_3.y) / sqrt (pow(this->expect_point[0] - this->SM_3.x, 2) + pow(this->expect_point[2] - this->SM_3.z, 2)));
	this->fai = atan(-(this->expect_point[2]) / this->expect_point[0]);

}

terminal_Guider::terminal_Guider(missile m1, missile target):terminal_guidance(m1, target) {

	//this->SM_3 = m1;
	//this->expect_point = target;
}
Data_Buffer::Data_Buffer(string root1)
{
	outfile.open(root1);

}
void Data_Buffer::record(string data)
{
	outfile << data << endl;
}
void Data_Buffer::save( )
{
	outfile.close();
}


missile Guider::proportional_guide(int Ky = 4, int Kz = 4, Data_Buffer dandao_data)
{
		//*******************�������****************
	//ofstream outfile;
	//outfile.open("./3D_proportional_guide.txt");
	string r1 = "./3D_proportional_guide.txt";
	//outfile << "����λ��X Y Z �ٶ�v  �������theta ����ƫ��fai  " << "Ŀ��λ��X Y Z �ٶ�v �������theta ����ƫ��fai  " << endl; //�Ƕ�theta eta
	string r2 = "����λ��X Y Z �ٶ�v  �������theta ����ƫ��fai Ŀ��λ��X Y Z �ٶ�v �������theta ����ƫ��fai ";
	//*****************�Ƶ�����ѭ��******************
	//		//�������ʴ�С����
	//this->SM_3.vel = Vector3d(this->SM_3.velocity, 0.0, 0.0);
	//this->expect_point.vel = Vector3d(this->expect_point.velocity, 0.0, 0.0);
	double dt = 0.01;
	double count = 0;
	while (this->r > 10)
	{
		count += dt;
		outfile << this->SM_3.x << " " << this->SM_3.y << " " << this->SM_3.z << " " << this->SM_3.vel[0] << " " << this->SM_3.theta * 180 / 3.1415926 << " " << this->SM_3.fai * 180 / 3.1415926
			<< " " << this->expect_point.x << " " << this->expect_point.y << " " << this->expect_point.z << " " << this->expect_point.vel[0] << " " << this->expect_point.theta * 180 / 3.1415926 << " " << this->expect_point.fai * 180 / 3.1415926 << endl;

		//����Ƕ�ʱ������ת�Ƕȣ�this->SM_3.theta * 180 / 3.1415926 << this->SM_3.fai * 180 / 3.1415926
		//*******������Ŀ��λ��״̬��ƫ�ǡ���ǡ�����*************

		this->r_dot = (G2T(T2G(this->expect_point.vel, this->expect_point.theta, this->expect_point.fai), this->theta, this->fai) - G2T(T2G(this->SM_3.vel, this->SM_3.theta, this->SM_3.fai), this->theta, this->fai))[0];
		this->theta_dot = (G2T(T2G(this->expect_point.vel, this->expect_point.theta, this->expect_point.fai), this->theta, this->fai) - G2T(T2G(this->SM_3.vel, this->SM_3.theta, this->SM_3.fai), this->theta, this->fai))[1] / this->r;
		this->fai_dot = (-1 * (G2T(T2G(this->expect_point.vel, this->expect_point.theta, this->expect_point.fai), this->theta, this->fai)[2]) + (G2T(T2G(this->SM_3.vel, this->SM_3.theta, this->SM_3.fai), this->theta, this->fai))[2]) / (this->r * cos(this->theta));

		//������������
		this->r = sqrt(pow(this->expect_point.x - this->SM_3.x, 2) + pow(this->expect_point.y - this->SM_3.y, 2) + pow(this->expect_point.z - this->SM_3.z, 2));//
		//std::cout << this->r << endl;
		//this->r += this->r_dot*dt;
		this->fai += this->fai_dot * dt;
		this->theta += this->theta_dot * dt;
		//ֻ������������"�����Ƕȱ仯��"��ֻ�ʺ��㷨��������û��"��������������"��ʵ�ʵĸı��ٶȷ���
		//(�����Ϣ-(�Ƶ���)->�����Ƕ�->�������->���õ�������-(������)->��ƫ-(����ѧģ�͡��˶�ѧ��ʽ)->ʵ������ʵ�ʹ���->λ�ñ任->�����Ϣ
		this->SM_3.ay = Ky * -1 * (this->r_dot) * this->theta_dot;
		this->SM_3.az = -Kz * -1 * (this->r_dot) * this->fai_dot * cos(this->theta);

		//this->SM_3.ay = Ky * (this->r_dot * this->theta_dot+this->r+theta_dot2);
		//this->SM_3.az = -Kz * abs(this->r_dot) * this->fai_dot * cos(this->theta);
		//***������***
		this->SM_3.Fx = 0;
		this->SM_3.Fy = this->SM_3.m * this->SM_3.ay;
		this->SM_3.Fz = this->SM_3.m * this->SM_3.az;
		//cout << this->SM_3.ay << "    " << this->SM_3.az << endl;//ay���ٶȴ���� 40   az -40
		this->SM_3.acc = Vector3d(this->SM_3.Fx / this->SM_3.m, 0, 0);
		this->SM_3.theta_dot = this->SM_3.Fy / (this->SM_3.m * this->SM_3.velocity);
		this->SM_3.fai_dot = -1 * this->SM_3.Fz / (this->SM_3.m * this->SM_3.velocity * cos(this->SM_3.theta));

		//ʵ���ǿ��Ʒ�����أ��ı��ٶȷ����ƫ��	//���ĳ�������任�ǳ��󣬿����ǻ���ʱ������û��dt
		this->SM_3.vel += this->SM_3.acc * dt;
		this->SM_3.theta += this->SM_3.theta_dot * dt;
		this->SM_3.fai += this->SM_3.fai_dot * dt;
		this->SM_3.velocity = this->SM_3.vel[0];
		//cout << this->SM_3.vel << endl;
		//vel���ٶ�����ϵ��ת������������ϵ�����ڵ�������ϵ��xyz���㣬������Ŀ��Ե�λ�ø���
		this->SM_3.x += this->SM_3.velocity * cos(this->SM_3.theta) * cos(this->SM_3.fai) * dt;
		this->SM_3.y += this->SM_3.velocity * sin(this->SM_3.theta) * dt;
		this->SM_3.z += -1 * this->SM_3.velocity * cos(this->SM_3.theta) * sin(this->SM_3.fai) * dt;
		this->expect_point.x += this->expect_point.velocity * cos(this->expect_point.theta) * cos(this->expect_point.fai) * dt;
		this->expect_point.y += this->expect_point.velocity * sin(this->expect_point.theta) * dt;
		this->expect_point.z += -1 * this->expect_point.velocity * cos(this->expect_point.theta) * sin(this->expect_point.fai) * dt;
		this->SM_3.pos = Vector3d(this->SM_3.x, this->SM_3.y, this->SM_3.z);
		this->expect_point.pos = Vector3d(this->expect_point.x, this->expect_point.y, this->expect_point.z);
		//xyzһֱ�ڸ��£�����û�и�pos��ֵ������pos���ǳ�ֵ
		//c++��ifֻ��һ����䣬û{}��ֱ���˳�
		if (this->r_dot > 0 || count >= 60)//|| count >= 60
		{
			std::cout << "δ�ܸ���" << endl;
			break;
		}
	}
	outfile.close();

	std::cout << this->expect_point.pos[0] - this->SM_3.pos[0] << std::endl;
	std::cout << this->expect_point.pos[1] - this->SM_3.pos[1] << std::endl;
	std::cout << this->expect_point.pos[2] - this->SM_3.pos[2] << std::endl;
	std::cout << "�Ѱ���:" << sqrt(pow(this->expect_point.pos[0] - this->SM_3.pos[0], 2) + pow(this->expect_point.pos[1] - this->SM_3.pos[1], 2) + pow(this->expect_point.pos[2] - this->SM_3.pos[2], 2));
	return SM_3;
}





missile terminal_phase_Guider::slidemode(int k1 = 6, int k2 = 0.1, int k3 = 6, int k4 = 0.1, double theta_e = PI / 6, double eta_e = PI / 6, string root = "./3D_slidemode_guide.txt")//��ģ���ܲ�����Ӱ�죬���·�ɢ,һ��������һ�������� �Ϳ��ܵ��������ͷ�ɢ
{//���ɷ��ʱ�������Ҫ��public

	//string root = "./3D_slidemode_guide.txt";
	//------------------------�������--------------
	//ofstream outfile;  outfile.open("./3D_slidemode_guide.txt");

	//outfile << "����λ��X Y Z �ٶ�v  �������theta ����ƫ��fai  " << "Ŀ��λ��X Y Z �ٶ�v �������theta ����ƫ��fai  " << endl; //�Ƕ�theta eta
	//ofstream accfile;
	//accfile.open("./3D_acc_guide.txt");

	//*****************�Ƶ�����ѭ��******************
	//		//�������ʴ�С����
	//this->SM_3.vel = Vector3d(this->SM_3.velocity, 0.0, 0.0);
	//this->expect_point.vel = Vector3d(this->expect_point.velocity, 0.0, 0.0);
	double dt = 0.01;
	double count = 0;
	double s1 = 0;
	double s2 = 0;
	double s1_dot = 0;
	double s2_dot = 0;
	double ay_dm = 0;
	double az_dm = 0;
	double ep1 = 0.001;
	double ep2 = 0.001;
	bool flag = 1;
	double temp = this->r;
	while (this->r > 4)
	{
		count += dt;	//����Ƕ�ʱ������ת�Ƕȣ�this->SM_3.theta * 180 / 3.1415926 << this->SM_3.fai * 180 / 3.1415926
		outfile << this->SM_3.x << " " << this->SM_3.y << " " << this->SM_3.z << " " << this->SM_3.vel[0] << " " << this->SM_3.theta * 180 / 3.1415926 << " " << this->SM_3.fai * 180 / 3.1415926
			<< " " << this->expect_point.x << " " << this->expect_point.y << " " << this->expect_point.z << " " << this->expect_point.vel[0] << " " << this->expect_point.theta * 180 / 3.1415926 << " " << this->expect_point.fai * 180 / 3.1415926 << endl;
		//outfile << this->SM_3.x << " " << this->SM_3.y << " " << this->SM_3.z << "   " << this->expect_point.x << " " << this->expect_point.y << " " << this->expect_point.z << " " << endl;
		//����˶�ѧģ�ͣ�r��theta��eta������ M��T��V��theta��eta����
		this->r_dot = (this->expect_point.velocity * cos(this->expect_point.theta) * cos(this->theta) * cos(this->expect_point.fai - this->fai) + this->expect_point.velocity * sin(this->expect_point.theta) * sin(this->theta)) - (this->SM_3.velocity * cos(this->SM_3.theta) * cos(this->theta) * cos(this->SM_3.fai - this->fai) + this->SM_3.velocity * sin(this->SM_3.theta) * sin(this->theta));
		this->theta_dot = ((this->expect_point.velocity * cos(this->theta) * sin(this->expect_point.theta) - this->expect_point.velocity * sin(this->theta) * cos(this->expect_point.theta) * cos((this->expect_point.fai - this->fai))) - (this->SM_3.velocity * cos(this->theta) * sin(this->SM_3.theta) - this->SM_3.velocity * sin(this->theta) * cos(this->SM_3.theta) * cos((this->SM_3.fai - this->fai)))) / this->r;
		this->fai_dot = (this->SM_3.velocity * cos(this->SM_3.theta) * sin(this->fai - this->SM_3.fai) - this->expect_point.velocity * cos(this->expect_point.theta) * sin(this->fai - this->expect_point.fai)) / (this->r * cos(this->theta));
		//������������
		this->r = sqrt(pow(this->expect_point.x - this->SM_3.x, 2) + pow(this->expect_point.y - this->SM_3.y, 2) + pow(this->expect_point.z - this->SM_3.z, 2));
		this->fai += this->fai_dot * dt;
		this->theta += this->theta_dot * dt;
		//���������ɣ�����ay��az  
		//epsilion = y0-(this->expect_point.y - this->SM_3.y);
		//t_go = this->r / this->r_dot + 0.0000001;
		s1 = this->theta_dot + k1 * (this->theta - theta_e);
		s2 = this->fai_dot + k3 * (this->fai - eta_e);
		ay_dm = sat(-2 * this->r_dot * this->theta_dot - this->r * cos(this->theta) * sin(this->theta) * this->fai_dot * this->fai_dot + this->r * k1 * this->theta_dot + this->r * ep1 * sign(s1) + this->r * k2 * s1, 250.0);//ģ�庯�����Զ�������ע������������ģ���������
		az_dm = sat(2 * this->r_dot * this->fai_dot * cos(this->theta) - 2 * this->r * this->fai_dot * this->theta_dot * sin(this->theta) - this->r * cos(this->theta) * k3 * this->fai_dot - this->r * cos(this->theta) * (ep2 * sign(s2) + k4 * s2), 250.0);
		//ay_dm = -2 * this->r_dot * this->theta_dot - this->r * cos(this->theta) * sin(this->theta) * this->fai_dot * this->fai_dot + this->r * k1 * this->theta_dot + this->r * ep1 * sign(s1) + this->r * k2 * s1;
		//az_dm = 2 * this->r_dot * this->fai_dot * cos(this->theta) - 2 * this->r * this->fai_dot * this->theta_dot * sin(this->theta) - this->r * cos(this->theta) * k3 * this->fai_dot - this->r * cos(this->theta) * (ep2 * sign(s2) + k4 * s2);

		//ay_dm = -this->r * this->fai_dot * this->fai_dot * sin(this->theta) * cos(this->theta) + 6 * (this->r_dot * this->r_dot * this->theta) / this->r - 6 * this->r_dot * this->theta_dot;
		//az_dm = 2 * this->r * this->theta_dot * this->fai_dot * sin(this->theta) + 6 * (this->r_dot * this->r_dot * this->fai * cos(this->theta)) / this->r - 6 * this->fai_dot * this->r_dot * cos(this->theta_dot);
		this->SM_3.acc = G2T(T2G(Vector3d(0, ay_dm, az_dm), this->theta, this->fai), this->SM_3.theta, this->SM_3.fai);
		accfile << this->SM_3.acc[0] << " " << this->SM_3.acc[1] << " " << this->SM_3.acc[2] << endl;
		//ay az����Թ�������ϵ�ļ��ٶ�   ���ǵ�Ŀ��������ϵ �ǵ���ϵ
		//this->SM_3.acc = G2T(Vector3d(0, ay_dm, az_dm),this->SM_3.theta,this->SM_3.fai);//ay az����ڵ��� �ǶԵ�
		//cout << this->SM_3.acc[0]<< "  "<<this->SM_3.acc[1] << "  " << this->SM_3.acc[2] << endl;
		//this->SM_3.acc =Vector3d(0, ay_dm, az_dm);
		this->SM_3.F = this->SM_3.m * this->SM_3.acc;
		this->SM_3.Fx = 0;// this->SM_3.F[0];//ʵ�ʵ���û�취���ٶ�
		this->SM_3.Fy = this->SM_3.F[1];
		this->SM_3.Fz = this->SM_3.F[2];
		//cout << this->SM_3.acc << "    " << endl;
		//ay��az->Fy Fz
		//this->SM_3.Fy = this->SM_3.m * this->SM_3.ay;// ;ay_dm
		//this->SM_3.Fz = this->SM_3.m * this->SM_3.az;// ;az_dm
		//�������Ķ���ѧģ��
		this->SM_3.theta_dot = this->SM_3.Fy / (this->SM_3.m * this->SM_3.velocity);
		this->SM_3.fai_dot = -1 * this->SM_3.Fz / (this->SM_3.m * this->SM_3.velocity * cos(this->SM_3.theta));
		//this->SM_3.vel+=Vector3d(this->SM_3.Fx / this->SM_3.m, 0, 0);
		//�ٰѼ��ٶ�ת�����ٶ�����ϵ��
		this->SM_3.acc = Vector3d(this->SM_3.Fx / this->SM_3.m, 0, 0);

		//�����˶�ѧģ��
		//if (this->SM_3.velocity >=800)
		//	this->SM_3.acc = Vector3d(0, 0, 0);
		this->SM_3.vel += this->SM_3.acc * dt;
		//cout << this->SM_3.vel << endl;
		this->SM_3.theta += this->SM_3.theta_dot * dt;
		this->SM_3.fai += this->SM_3.fai_dot * dt;
		this->SM_3.velocity = this->SM_3.vel[0];
		//�ٶ�����ϵת������������ϵ
		this->SM_3.x += this->SM_3.velocity * cos(this->SM_3.theta) * cos(this->SM_3.fai) * dt;
		this->SM_3.y += this->SM_3.velocity * sin(this->SM_3.theta) * dt;
		this->SM_3.z += -1 * this->SM_3.velocity * cos(this->SM_3.theta) * sin(this->SM_3.fai) * dt;
		this->expect_point.x += this->expect_point.velocity * cos(this->expect_point.theta) * cos(this->expect_point.fai) * dt;
		this->expect_point.y += this->expect_point.velocity * sin(this->expect_point.theta) * dt;
		this->expect_point.z += -1 * this->expect_point.velocity * cos(this->expect_point.theta) * sin(this->expect_point.fai) * dt;
		this->SM_3.pos = Vector3d(this->SM_3.x, this->SM_3.y, this->SM_3.z);
		this->expect_point.pos = Vector3d(this->expect_point.x, this->expect_point.y, this->expect_point.z);
		//���������˳����
		if (this->r > temp)//r��С�ˣ�����r_dot����0��   todo
		{
			flag = 1;
			if (temp > 15)
				flag = 0;
			break;
		}
		temp = this->r;
	}
	outfile.close();
	accfile.close();
	if (!flag)
		std::cout << "�Ѱ���>15��δ�ܸ���" << endl;
	std::cout << this->expect_point.pos[0] - this->SM_3.pos[0] << std::endl;
	std::cout << this->expect_point.pos[1] - this->SM_3.pos[1] << std::endl;
	std::cout << this->expect_point.pos[2] - this->SM_3.pos[2] << std::endl;
	std::cout << "�Ѱ���:" << (temp + this->r) / 2;

	return SM_3;
}


*/