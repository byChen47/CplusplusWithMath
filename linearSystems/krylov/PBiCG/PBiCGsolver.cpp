// File: linearSystems/krylov/PBiCG/PBiCGsolver.cpp
// Category: Krylov / preconditioned iterative solver
// Method: PBiCG
// Purpose: Alternative PBiCG solver implementation retained for experimentation.

#include "PBiCGsolver.h"
#include <iostream>
#include <chrono>
#include <cmath>

double PBiCGsolver::dotProduct(const Vector & a, const Vector & b)
{
    double result = 0.0;
    for(size_t i = 0; i < b.size() ;i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

double PBiCGsolver::norm(const Vector & v)
{
    return std::sqrt(dotProduct(v,v));
}

void PBiCGsolver::axpy(double alpha,const Vector& x, Vector & y)
{
    for(size_t i = 0 ; i < y.size() ; i++)
    {
        y[i] += alpha * x[i];
    }
}

void PBiCGsolver::copyVector(const Vector & src, Vector & dst)
{
    dst = src;
}

PBiCGsolver::SolverResult PBiCGsolver::solve(
    const MatrixVectorFunc& A,
    const MatrixVectorFunc& A_T,
    const Vector& b,
    Vector& x,
    const MatrixVectorFunc& M,
    const MatrixVectorFunc& M_T,
    const SolverConfig& config
)
{
    SolverResult result;
    auto star_time = std::chrono::high_resolution_clock::now();

    int n = static_cast<int>(b.size());
    // 初始值
    Vector r(n),r_tilde(n),p(n),p_tilde(n);
    Vector Ap(n),A_T_p_tilde(n);
    Vector z(n),z_tilde(n);

    A(x,Ap);// 计算初始残差 r = b - A * x
    for(int i = 0 ; i < n; ++i)
    {
        r[i] = b[i] - Ap[i];
    }

    copyVector(r,r_tilde);
    result.initial_residual = norm(r);
    if(result.initial_residual < config.tolerance)
    {
        result.converged = true;
        result.iterations = 0;
        result.final_residual = result.initial_residual;
        auto end_time = std::chrono::high_resolution_clock::now();
        result.compute_time_ms = std::chrono::duration_cast<std::chrono::microseconds>
        (end_time - star_time).count();
    }
    if(config.verbose)
    {
        std::cout<<"PBiCG求解开始,初始残差: "<<result.initial_residual<<std::endl;
    }

    double rho_old = 1.0,rho_new,alpha,beta; // 初始化值
    result.converged = false;

    for(int k = 0 ; k < config.max_iterations; ++k)
    {
        M(r,z); //预处理残差
        M_T(r_tilde,z_tilde);

        // 计算内积检查稳定性
        rho_new = dotProduct(z,r_tilde);
        if(std::fabs(rho_new) < 1e-30)
        {
            if(config.verbose)
            {
                std::cout<<"PBiCG中断:rho接近于零" << std::endl;
            }
            break;
        }
        // 更新搜索方向
        if(k == 0) //第一次迭代
        {
            copyVector(z,p);
            copyVector(z_tilde,p_tilde);
        }
        else //后续迭代 
        {
            beta = rho_new/ rho_old;
            for(int i = 0 ; i < n ; ++i)
            {
                p[i] = z[i] + beta * p[i];
                p_tilde[i] = z_tilde[i] + beta * p_tilde[i];
            }
        }
        A(p,Ap);
        A_T(p_tilde,A_T_p_tilde);

        alpha = rho_new/dotProduct(p_tilde,Ap);

        // 更新残差和解
        axpy(alpha,p,x);
        for(int i = 0 ; i < n; ++i)
        {
            r[i] -=alpha * Ap[i];
            r_tilde[i] -=alpha * A_T_p_tilde[i];
        }
        double residual = norm(r);
        double rel_residual = residual/result.initial_residual;
        if(rel_residual < config.tolerance)
        {
            result.converged = true;
            result.iterations = k + 1;
            result.final_residual = residual;

            auto end_time = std::chrono::high_resolution_clock::now();
            result.compute_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>
            (end_time - star_time).count();

            if(config.verbose)
            {
                std::cout<<"\nPBiCG收敛于迭代 " <<result.iterations
                         <<",最终残差: " <<residual
                         <<",相对残差: " << rel_residual<<std::endl;

            }
            return result;
        }

        if(config.verbose&&(k + 1)% config.print_interval == 0)
        {
            std::cout<<"迭代"<<k+1
                     <<",相对残差 "<<rel_residual<<std::endl;
        }

        rho_old = rho_new;
    }

    result.converged = false;
    result.iterations = config.max_iterations;
    result.final_residual = norm(r);
    auto end_time = std::chrono::high_resolution_clock::now();
    result.compute_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>
    (end_time - star_time).count();

    if(config.verbose)
    {
        std::cout<<"PBiCG未在最大迭代次数内收敛" << std::endl;

    }
    return result;
}


PBiCGsolver::SolverResult PBiCGsolver::solve(
    const MatrixVectorFunc& A,
    const Vector& b,
    Vector& x,
    const MatrixVectorFunc& M,
    const SolverConfig& config
)
{
    // 对称预处理:M_T = M
    return solve(A,A,b,x,M,M,config);

}

PBiCGsolver::SolverResult PBiCGsolver::solve(
    const MatrixVectorFunc& A,
    const Vector& b,
    Vector& x,
    const SolverConfig& config
)
{
    // 无预处理:M = I (单位矩阵)
    auto identity = [](const Vector& v ,Vector& result)
    {
        result = v;
    };

    return solve(A,A,b,x,identity,identity,config);   
}
