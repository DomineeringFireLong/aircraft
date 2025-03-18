#include<iostream>
#include<cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include<graphics.h>
#include<fstream>
#include<random>
#include"controller.h"
#include"guide.h"
#include"powered_phase.h"
#include"aerodynamic.h"

using namespace std;
using namespace Eigen;


//ָ�뺯����������һ������������ֵ��ָ��

int* func_sum(int n)
{
    if (n < 0)
    {
        cout<<"error:n must be > 0\n";
        exit(-1);
    }
    static int sum = 0;//һ��ľֲ������Ǵ����ջ���ģ������������󣬾ͻᱻ�ͷŵ���stastic�������ӱ����������������������������ڼ�
    int* p = &sum;//���� ���ؾֲ�����ָ��  ��1.stastic(����ڴ����) 2.ȫ�ֱ���(��������ݶ�)
    for (int i = 0; i < n; i++)
    {
        sum += i;
    }
    return p;
}
//����ָ�룬��������һ��ָ�룬ָ������ָ��//�����Ķ����ڴ���Σ�ÿ�������ڴ���ζ�����ڵ�ַ�����亯��ָ��

int max1(int a, int b)
{
    return a > b ? a : b;
}
/*

int main(void)
{
    //int num = 0;
    //cout << "please input one number:";
    //cin>>num;
    //int* p = func_sum(num);
    //cout <<* p;

    int (*p)(int, int); //����ָ��Ķ���
    p = max1;
    int ret = p(10, 15);
    //int ret = (*max1)(10,15);
    //int ret = (*p)(10,15);//(*p)������ľ��Ǻ���ָ����ָ���ֵ��Ҳ���Ǻ�������
    
    cout << ret << endl;
    return 0;
}
*/

void main() {


	////Integration_test();
	//missile SM3(1501,3.4,Vector3d(23.64544063, 5427.480754, 5427.480754),Vector3d::Zero(), Vector3d(0.01,0,0), Vector3d::Zero(), Vector3d(90.0*DEG2RAD, 0, 0), \
	//	90.0 * DEG2RAD, 0, Vector3d::Zero(),Vector3d::Zero(),Vector3d::Zero(), 0, 0,0,0,0,0,0,0);

	//bool Ifdisplay = 0;	//���SM3����Ϣ��ֵ�� �����ζ���ĳ�Ա�������󣬱��������û��ߵ�ַ�������޷����ݳ���Ϣ�����õ���Ĭ�Ϲ��캯��
	//
	//powered_phase_stand_ballastic SM3_phase_one(SM3, "F:/SM_3/",174e3,Vector3d(0, 0, 0), \
	//	Vector3d(0, 0, 0),Vector3d(0, 0, 0), Vector3d(0, 0, 0));

	//SM3_phase_one.mymissile.time = 0.0;
	//SM3_phase_one.boost_turning_phase(Ifdisplay);
	//Ifdisplay = 1;
	//SM3_phase_one.powered_stage(Ifdisplay);

    Data_Processor t1;
    t1.test();
	

}


//int main()
//{
//    plt::plot({ 1,3,2,4 });
//    plt::show();
//    std::cout << "Hello World!\n";
//}