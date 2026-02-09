#include <iostream>
#include <cmath>
#include <iomanip>

// 定义迭代的函数

double function(double x)
{
    return x*x*x -2*x -5;
}

// 定义导数
double derivative(double x)
{
    return 3*x*x -2;
}

// 牛顿法实现，牛顿法是通过迭代减小真实值与理论值之间的误差来实现求解
// 所以需要进行迭代，需要一个初始的值，一个迭代次数，一个误差精度

double newtonMethod(double initiaGuess ,double tolerance ,int maxIterations, double & result)
{
    double x = initiaGuess;   // 随便给出一个数值
    double x_prev;            // 上一步的值，主要是为了测量两个值之间的误差
    int iteration = 0;        // 设置初始迭代次数为0
    bool success = false; 

    std::cout << "迭代过程" << std::endl;
    std::cout << "n\t x_n\t\t f(x_n)\t\t 误差" << std::endl;
    std::cout <<"-------------------------------------" <<std::endl;

    do
    {
        /* code */
        x_prev = x;
        
        // 检查导数是否为0 
        double df = derivative(x);
        if(fabs(df) < 1e-15)
        {
            std::cout << " 错误，导数为零，算法无法继续"  <<std::endl;
            return false;
        }

        // 牛顿迭代法求解

        x = x - function(x)/df;
        iteration++;
        
        // 输出相关信息
        std::cout << iteration << "\t"<< x
        << "\t" << function(x)
        << "\t" << fabs(x-x_prev) << std::endl;

        // 检查是否收敛
        double error  = fabs(x-x_prev);
        double func_val = fabs(function(x));

        // std::cout << "error = " << error << ", func_val = " << func_val << std::endl;

        if(error <= tolerance && func_val <= tolerance )
        {
            success  = true;
            break;
        }

        //如果达到最大迭代次数，退出
        if(iteration >maxIterations)
        {
            std::cout << "达到最大迭代次数，可能为收敛" << std::endl;
            success  = true;
            break;
        }

    } while (true);
    
    result  =  x;
    return success;
}

// 主函数实现

int main()
{
    // 参数设置
    double initiaGuess = 15;
    double tolerance = 1e-15;
    double maxIterations = 20;

    std::cout << " 求解方程： x^3 - 2x -5 =0 " << std::endl;
    std::cout << "初始猜测值：  " << initiaGuess << std::endl;
    std::cout << "容差： " <<tolerance << std::endl;

    // 使用牛顿迭代法
    double root;
    double success = newtonMethod(initiaGuess,tolerance, maxIterations,root);

    if(success)
    {
        std::cout << "\n最终结果为: " << std:: endl;
        std::cout << "\n方程的根为: " << root << std::endl;
        std::cout << "\nf(根) =  " << function(root) << std::endl;

        // 验证结果
        std::cout << "\n验证 \t\n"  << std::endl;
        std::cout << root << "^3 - 2*"<<root << " -5" << " = "<< function(root) <<std::endl; 
        std::cout << "\n满足方程等于零" << std::endl;

    }
    return 0;
}