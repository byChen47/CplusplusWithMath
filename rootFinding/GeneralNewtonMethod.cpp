// File: rootFinding/GeneralNewtonMethod.cpp
// Category: Nonlinear equation solver
// Method: General Newton-type method
// Purpose: Explore a more general Newton-style iteration workflow in C++.

#include <iostream>
#include <cmath>
#include <functional>

using namespace std;
// 通用牛顿迭代法
double newtowMethod(function<double(double)> f,
                    function<double(double)> df,
                    double initialGuess = 0.0,
                    double tolerance = 1e-10,
                    double max_iter = 100)
// f是一个方程，数值是double 类型，里边的变量是  double类型
// df同样也是
{
    double x = initialGuess;
    for(int i = 0 ; i < max_iter; i++)
    {
        double fx = f(x);
        double dfx = df(x);

        // 防止除零错误
        if(fabs(dfx) < 1e-15)
        {
            cout<<"警告:导数为0,停止迭代"<< endl;
            break;
        }

        double x_new = x - fx/ dfx;
        
        if(fabs(x_new - x) < tolerance)
        {
            return x_new;
        }
        x = x_new;
    }
    return x;
    
}

int main()
{
    double a = 7.0;

    auto f_sqrt = [a](double x){return x * x -a ;};
    auto df_sqrt = [a](double x ){return 2 * x ;};

    double sqrt_result = newtowMethod(f_sqrt,df_sqrt,a/2);
    std::cout << "求 " << a << " 的平方根：" << std::endl;
    std::cout << "牛顿法结果: " << sqrt_result << std::endl;
    std::cout << "标准库结果: " << std::sqrt(a) << std::endl;
    std::cout << "误差: " << std::fabs(sqrt_result - std::sqrt(a)) << std::endl;
    std::cout << "\n" << std::string(40, '-') << "\n" << std::endl;
}
