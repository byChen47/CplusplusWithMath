// PBiCGsolver.h
#ifndef PBICGSOLVER_H
#define PBICGSOLVER_H

#include <vector>
#include <functional>
#include <memory>

class PBiCGsolver
{
public:
    
    using Vector = std::vector<double>;
    using MatrixVectorFunc = std::function<void(const Vector&,Vector&)>;

    // 求解器配置
    struct SolverConfig
    {
        double tolerance;
        int max_iterations;
        bool verbose;
        int print_interval;

        // 构造函数提供默认值
        SolverConfig():
            tolerance(1e-8),
            max_iterations(1000),
            verbose(true),
            print_interval(100)
            {

            }
    };

    // 求解结果
    struct SolverResult
    {
        bool converged;
        int iterations;
        double residual_norm;//残差向量范数
        double initial_residual;
        double final_residual;
        long long compute_time_ms;

        // 构造函数提供默认值
        SolverResult():
            converged(false),
            iterations(0),
            residual_norm(0.0),
            initial_residual(0.0),
            final_residual(0.0),
            compute_time_ms(0)
            {

            }
    };
    
    PBiCGsolver() = default;

    // 使用矩阵向量乘积函数求解
    SolverResult solve(
        const MatrixVectorFunc& A,
        const MatrixVectorFunc& A_T,
        const Vector& b,
        Vector & x,
        const MatrixVectorFunc& M,
        const MatrixVectorFunc& M_T,
        const SolverConfig& config
    );
    
    // 便捷接口（对称预处理）
    SolverResult solve(
        const MatrixVectorFunc& A,
        const Vector& b,
        Vector & x,
        const MatrixVectorFunc& M,
        const SolverConfig& config

    );

    // 无预处理版本
    SolverResult solve(
        const MatrixVectorFunc& A,
        const Vector& b,
        Vector& x,
        const SolverConfig&  config
    );

private:
    double dotProduct(const Vector& a,const Vector&b);
    double norm(const Vector& v);
    void axpy(double alpha,const Vector & x,Vector & y);
    void copyVector(const Vector& src,Vector& dst);

    // 辅助函数，创建默认配置
    static SolverConfig defaultConfig()
    {
        return SolverConfig();
    }
};

#endif

