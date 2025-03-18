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


//指针函数：本质是一个函数，返回值是指针

int* func_sum(int n)
{
    if (n < 0)
    {
        cout<<"error:n must be > 0\n";
        exit(-1);
    }
    static int sum = 0;//一般的局部变量是存放再栈区的，当函数结束后，就会被释放掉，stastic可以增加变量的生命周期至整个程序运行期间
    int* p = &sum;//避免 返回局部变量指针  ：1.stastic(存放在代码段) 2.全局变量(存放在数据段)
    for (int i = 0; i < n; i++)
    {
        sum += i;
    }
    return p;
}
//函数指针，本质上是一个指针，指向函数的指针//函数的定义在代码段，每个函数在代码段都有入口地址，即其函数指针

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

    int (*p)(int, int); //函数指针的定义
    p = max1;
    int ret = p(10, 15);
    //int ret = (*max1)(10,15);
    //int ret = (*p)(10,15);//(*p)所代表的就是函数指针所指向的值，也就是函数本身
    
    cout << ret << endl;
    return 0;
}
*/

void main() {


	////Integration_test();
	//missile SM3(1501,3.4,Vector3d(23.64544063, 5427.480754, 5427.480754),Vector3d::Zero(), Vector3d(0.01,0,0), Vector3d::Zero(), Vector3d(90.0*DEG2RAD, 0, 0), \
	//	90.0 * DEG2RAD, 0, Vector3d::Zero(),Vector3d::Zero(),Vector3d::Zero(), 0, 0,0,0,0,0,0,0);

	//bool Ifdisplay = 0;	//想把SM3的信息赋值给 主动段对象的成员导弹对象，必须用引用或者地址，否则无法传递出信息，调用的是默认构造函数
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