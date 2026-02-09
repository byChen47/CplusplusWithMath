#include <iostream>
#include <cmath>
#include <iomanip>
// #include <vector>

using namespace std;

// 定义函数f(x) = x^3 - 2x -5
double f(double x)
{
    return x * x * x - 2 * x -5;    
}

// 二分法求解函数，最终得到的是根
// 也就是区间中点 c 的值
double bisectionMethod(double a ,double b ,double tolerance , int max_iter,bool showIteration)
{
    if(f(a ) * f( b) > 0)
    {
        cout <<"错误：函数在区间端点必须异号！ " << endl;
        return NAN; // 返回非数值表示错误
    }

    double c = 0; // 中点
    int iterations = 0;
    if(showIteration)
    {
        cout<<"=====================================================" << endl;
        cout<<setw(8) <<"迭代次数" <<setw(10) <<"a"<<setw(16) <<"b" 
            <<setw(13)<<"c"<<setw(18) <<"f(c)" << setw(28) <<"区间长度" <<endl;
    }

    // 二分法迭代求解
    for(iterations = 1; iterations < max_iter ; iterations++)
    {
        c = (a + b )/ 2;
        if(showIteration)
        {
            cout <<setw(8) <<fixed << setprecision(6) << iterations
                 <<setw(15)<<fixed << setprecision(6) << a
                 <<setw(15)<<fixed << setprecision(6) << b
                 <<setw(15)<<fixed << setprecision(6) << c
                 <<setw(15)<<fixed << setprecision(6) <<f(c)
                 <<setw(20)<<fixed << setprecision(6) <<(b -a)
                 << endl;
        }
        // 检查是否找到根
        if(f(c) == 0 || (b - a) / 2 <tolerance)
        {
            break;
        }

        // 确定新的区间
        if( f(a) * f(c) < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
    }
    if(showIteration)
    {
        cout<<"=====================================" << endl;
    }
    cout<<"最终迭代次数是: " << iterations <<endl;
    return c;
}

int main()
{
    // 设置输出格式
    cout <<fixed <<setprecision(6);

    // 定义参数
    double a = 2;
    double b = 3;
    double tolerance = 1e-7;
    int max_iter = 100;
    cout<<"二分法求解方程根 - C++实现" << endl;
    // cout<<"=========================================" <<endl;
    double root = bisectionMethod(a,b,tolerance ,max_iter,true);

    cout <<"方程的根约为: " <<root <<endl;
    cout <<"函数值f(" <<root<<") = " << f(root) << endl;

}