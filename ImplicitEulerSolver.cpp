//**************隐式欧拉法 ****************//
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>

using namespace std;
class ImplicitEulerSolver
{
private:
    // 按照初始化顺序声明成员变量
    double t0, y0; //定义初始条件
    double h; // 定义时间步长
    double tolerance; // 定义牛顿迭代法的收敛误差
    int max_iterations; // 定义牛顿迭代法的最大迭代次数
    
    // 定义微分方程 f(t, y)
    std::function<double(double, double)> f;
    
    // 定义f对y的偏导数函数
    std::function<double(double, double)> df_dy;

public:

    // 构造函数2：df_dy采用解析导数求解
    ImplicitEulerSolver(
        std::function<double(double, double)> func,
        std::function<double(double, double)> derivative,
        double initial_t, double initial_y,
        double step_size, double tol = 1e-8,
        int max_iter = 100
    ) : 
        t0(initial_t), y0(initial_y), h(step_size),  // 按声明顺序初始化
        tolerance(tol), max_iterations(max_iter),
        f(func), df_dy(derivative)  // 使用传入的导数
    {}

    // 使用牛顿法求解隐式方程
    double solveImplicitEquation(double t_n, double y_n, double t_next)
    {
        // 初始猜测，使用显式欧拉法预测
        double y_guess = y_n + h * f(t_n, y_n);  // 修正：添加步长h
        
        for (int iter = 0; iter < max_iterations; ++iter)
        {
            // 计算 F(y) = y - y_n - h * f(t_next, y)
            double F = y_guess - y_n - h * f(t_next, y_guess);
            
            // 计算导数 dF/dy = 1 - h * df/dy
            double dF = 1.0 - h * df_dy(t_next, y_guess);
            
            // 避免除零错误
            if (fabs(dF) < 1e-15)
            {
                std::cerr << "警告: 导数太小, 使用简单迭代法" << std::endl;
                // 使用简单迭代法
                y_guess = y_n + h * f(t_next, y_guess);
                continue;
            }
            
            // 牛顿法更新
            double delta = F / dF;
            y_guess = y_guess - delta;
            
            // 检查收敛
            if (fabs(delta) < tolerance)
            {
                return y_guess;
            }
        }
        
        std::cerr << "警告: 达到最大迭代步数" << std::endl;
        return y_guess;  // 返回最后一次迭代的结果
    }



    // 指定求解时间
    std::vector<std::pair<double, double>> solverUntil(double end_time)
    {
        std::vector<std::pair<double, double>> results;
        double t = t0;
        double y = y0;
        double current_h = h;

        results.push_back({t, y});

        while (t < end_time)
        {
            // 调整最后一步的步长
            if (t + current_h > end_time)
            {
                current_h = end_time - t;
            }

            double t_next = t + current_h;
            y = solveImplicitEquation(t, y, t_next);
            t = t_next;
            results.push_back({t, y});

            // 恢复原始步长
            current_h = h;
        }
        
        return results;  // 添加返回语句
    }
};

// 示例1：线性方程 dy/dt = -λy
double linearEquation(double t, double y)
{
    (void)t;  // 标记t未使用
    double lambda = 1000.0;  // 大参数
    return -lambda * y;
}

double linearEquationDerivative(double t, double y)
{
    (void)t;  // 标记t未使用
    (void)y;  // 标记y未使用
    double lambda = 1000.0;
    return -lambda;  // ∂f/∂y = -λ
}

// 示例2：非线性方程 dy/dt = y - y^3
double nonLinearEquation(double t, double y)
{
    (void)t;  // 标记t未使用
    return y - y * y * y;
}

double nonLinearEquationDerivative(double t, double y)
{
    (void)t;  // 标记t未使用
    return 1.0 - 3.0 * y * y;  // ∂(y - y^3)/∂y = 1 - 3y^2
}

// 测试简单的指数衰减方程，便于调试
void testSimpleExponential()
{
    std::cout << "\n=== 简单指数衰减方程测试 ===" << std::endl;
    std::cout << "方程: dy/dt = -y, y(0) = 5" << std::endl;
    
    auto simpleEq = [](double t, double y) -> double {
        (void)t;
        return -y;
    };
    
    auto simpleDeriv = [](double t, double y) -> double {
        (void)t;
        (void)y;
        return -1.0;
    };
    
    ImplicitEulerSolver solver(simpleEq, simpleDeriv, 0.0, 1.0, 0.05);
    auto results = solver.solverUntil(2.0);
    
    std::cout << std::setw(10) << "t" << std::setw(15) << "数值解"
              << std::setw(15) << "精确解" << std::setw(15) << "误差" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    
    for (const auto &point : results)
    {
        double t = point.first;
        double y_num = point.second;
        double y_exact = exp(-t);
        double error = fabs(y_num - y_exact);
        
        std::cout << std::setw(10) << std::fixed << std::setprecision(4) << t;
        std::cout << std::setw(15) << std::fixed << std::setprecision(6) << y_num;
        std::cout << std::setw(15) << std::scientific << std::setprecision(6) << y_exact;
        std::cout << std::setw(15) << std::scientific << std::setprecision(6) << error << std::endl;
    }
}

// 求解线性方程
void solverLinearEquation()
{
    std::cout << "=== 隐式欧拉法稳定性测试 ===" << std::endl;
    std::cout << "方程: dy/dt = -1000y, y(0) = 1" << std::endl;

    double h = 0.05;
    
    // 使用带导数的构造函数
    ImplicitEulerSolver solver(linearEquation, linearEquationDerivative, 0.0, 1.0, h);
    
    auto results = solver.solverUntil(2.0);  // 先计算较短时间，便于调试
    
    std::cout << std::setw(10) << "t" << std::setw(15) << "数值解"
              << std::setw(15) << "精确解" << std::setw(15) << "误差" << std::endl;
    std::cout << std::string(55, '-') << std::endl;

    for (const auto &point : results)
    {
        double t = point.first;
        double y_num = point.second;
        double y_exact = exp(-1000.0 * t);
        double error = fabs(y_num - y_exact);
        
        std::cout << std::setw(10) << std::fixed << std::setprecision(4) << t;
        std::cout << std::setw(15) << std::scientific << std::setprecision(6) << y_num;
        std::cout << std::setw(15) << std::scientific << std::setprecision(6) << y_exact;
        std::cout << std::setw(15) << std::scientific << std::setprecision(6) << error << std::endl;
    }
}

// 求解非线性方程
void solverNonLinearEquation()
{
    std::cout << "\n=== 非线性方程求解 ===" << std::endl;
    std::cout << "方程: dy/dt = y - y^3, y(0) = 1.0" << std::endl;

    // 使用带导数的构造函数
    ImplicitEulerSolver solver(nonLinearEquation, nonLinearEquationDerivative, 0.0, 0.5, 0.05);
    
    auto results = solver.solverUntil(3.0);

    std::cout << std::setw(10) << "t" << std::setw(15) << "y(t)" << std::endl;
    std::cout << std::string(25, '-') << std::endl;

    for (const auto &point : results)
    {
        // 每0.2显示一次，避免输出太多
        if (fmod(point.first, 0.1) < 1e-10) 
        {
            std::cout << std::setw(10) << std::fixed << std::setprecision(1) << point.first;
            std::cout << std::setw(15) << std::fixed << std::setprecision(6) << point.second << std::endl;
        }
    }
}

int main()
{

    std::cout << "隐式欧拉法(Implicit Euler Method)实现" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // 先测试简单的方程
    testSimpleExponential();
    
    // 测试刚性方程（使用较小的结束时间）
    solverLinearEquation();
    
    // 求解非线性方程
    solverNonLinearEquation();

    return 0;
}