// File: ode/EulerMethod.cpp
// Category: Ordinary differential equation solver
// Method: Explicit Euler
// Purpose: Demonstrate forward Euler time integration for simple ODE problems.

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// using namespace std;
// 显示欧拉方法

class ExplicitEulerSolver 
{
    private:

    double (*f)(double,double);  // 函数，
    double t0;                   // 初始化值
    double y0;                   // 初始化值
    double h;                    // 时间步长
    int steps;                   // 时间步数



    public:
    // 构造函数
    
    ExplicitEulerSolver(double(*f)(double ,double),double initail_t,double initial_y,double steps_size)
    :f(f),t0(initail_t),y0(initial_y),h(steps_size)
    {

    }

    void setStepsSize(double steps_size)
    {
        h = steps_size;
    }
    
    // 单步计算dy/dx

    double step(double t ,double y) const
    {
        return y + h * f( t , y );
    }

    
    // 计算指定步数的结果
    std::vector<std::pair<double,double>> solverSteps(int num_steps) const
    {
        std::vector<std::pair<double ,double>> results;
        double t = t0;
        double y = y0;
        results.push_back({t,y});
        for(int i ; i < num_steps ; i++)
        {
            y = y + h * f(t,y);
            t = t + h;
            results.push_back({t,y});
        }
        
        return results;

    }

    // 计算直到指定时间的结果

    std::vector<std::pair<double,double>> solveUntil(double end_time) const
    {
        std::vector<std::pair<double ,double>> results;
        double t = t0;
        double y = y0;
        results.push_back({t,y});
        while (t < end_time)
        {
            double dt = std::min(h,end_time - t);
            y = y * dt + f( t, y );
            t += dt;
            results.push_back({t , y});
        }
        
        
        return results;

    }

        // 计算指定时间点的近似值
    double evaluateAt(double target_t) const 
    {
        double t = t0;
        double y = y0;
        
        while (t < target_t) {
            double dt = std::min(h, target_t - t);
            y = y + dt * f(t, y); // f(t,y) 就是泰勒公式的值，泰勒展开
            t += dt;
        }
        
        return y;
    }

};

int main() 
{
    std::cout << "=== 面向对象的显式欧拉法实现 ===" << std::endl;
    
    // 示例1：基本方程
    std::cout << "\n1. 求解 dy/dt = y - t^2 + 1" << std::endl;
    auto func1 = [](double t, double y) { return y - t*t + 1; };
    
    ExplicitEulerSolver solver1(func1, 1, 0.7, 0.1);
    auto results1 = solver1.solveUntil(3.0);
    
    std::cout << "时间区间 [0, 3], 步长 h = 0.1" << std::endl;
    std::cout << std::setw(10) << "t" << std::setw(15) << "y" << std::endl;
    for (const auto& point : results1) 
    {
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << point.first
                  << std::setw(15) << std::fixed << std::setprecision(6) << point.second << std::endl;
    }

    // 计算待定时间点的值
    std::cout << "\n计算特定时间点的近似值 " <<std::endl;
    auto func3 = []( double t, double y )
    {
        (void) t;
        return -0.5 * y;
    
    };  // 定义函数以及返回值

    ExplicitEulerSolver solver3(func3 ,0.0,100,0.01) ;// 里边的数值就是定义函数的初始值

    double target_t = 2.0 ; // 目标时间
    double approx_value = solver3.evaluateAt(target_t);

    double exact_value = 100.0 * exp(-0.5 * target_t);

    std::cout << "\n方程: dy/dt = -0.5y,y(0) = 100 " << std::endl;
    std::cout << "\n在t     =      " << target_t << " 时："<<std::endl;
    std::cout << "\n显示欧拉近似值: " << std::fixed << std::setprecision(6) << approx_value << std::endl;
    std::cout << "\n精确值:        " << std::fixed << std::setprecision(6) << exact_value << std::endl;
    std::cout << "\n绝对误差:      " << std::scientific << std::abs(approx_value - exact_value) << std::endl;

    return 0;



}
