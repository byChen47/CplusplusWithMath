#include <iostream>
#include <cmath>
#include <functional>

// 通用的牛顿迭代法
double newtonMethod(std::function<double(double)> f, 
                   std::function<double(double)> df,
                   double initialGuess,
                   double tolerance = 1e-10,
                   int maxIterations = 100) {
    
    double x = initialGuess;
    
    for (int i = 0; i < maxIterations; i++) {
        double fx = f(x);
        double dfx = df(x);
        
        // 防止除零错误
        if (std::fabs(dfx) < 1e-15) {
            std::cout << "警告：导数为零，停止迭代。" << std::endl;
            break;
        }
        
        double x_new = x - fx / dfx;
        
        // 检查收敛
        if (std::fabs(x_new - x) < tolerance) {
            return x_new;
        }
        
        x = x_new;
    }
    
    return x;
}

int main() {
    // 示例1：求平方根 f(x) = x² - a
    double a = 7.0;
    
    // 使用lambda表达式定义函数和导数
    auto f_sqrt = [a](double x) { return x * x - a; };
    auto df_sqrt = [](double x) { return 2 * x; };
    
    double sqrt_result = newtonMethod(f_sqrt, df_sqrt, a/2);
    
    std::cout << "求 " << a << " 的平方根：" << std::endl;
    std::cout << "牛顿法结果: " << sqrt_result << std::endl;
    std::cout << "标准库结果: " << std::sqrt(a) << std::endl;
    std::cout << "误差: " << std::fabs(sqrt_result - std::sqrt(a)) << std::endl;
    
    std::cout << "\n" << std::string(40, '-') << "\n" << std::endl;
    
    // 示例2：求解方程 f(x) = x³ - 2x - 5 = 0
    auto f_cubic = [](double x) { return x * x * x - 2 * x - 5; };
    auto df_cubic = [](double x) { return 3 * x * x - 2; };
    
    double cubic_root = newtonMethod(f_cubic, df_cubic, 2.0);
    
    std::cout << "求解方程 x³ - 2x - 5 = 0：" << std::endl;
    std::cout << "近似根: " << cubic_root << std::endl;
    std::cout << "验证 f(" << cubic_root << ") = " << f_cubic(cubic_root) << std::endl;
    
    return 0;
}