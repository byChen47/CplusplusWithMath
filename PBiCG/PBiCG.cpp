#include "PBiCG.h"
#include <iostream>
#include <chrono>
#include <cmath>

double PBiCG::dotProduct(const Vector & a, const Vector & b)
{
    double result = 0.0;
    for(size_t i = 0; i < b.size() ; i++)
    {
        result += a[i] * b[i];
    }
    return result;

}
double PBiCG::norm(const Vector& v)
{
    return std::sqrt(dotProduct(v,v));
}
void PBiCG::axpy(double alpha,const Vector& x,Vector& y )
{
    for(size_t i = 0 ;i < y.size() ;i++ )
    {
        y[i] += alpha * x[i]; 
    }

}
void PBiCG::copyVector(const Vector& src, Vector& dst)
{
    dst = src;
}

PBiCG::SolverResult PBiCG::solve
(
    const MatrixVectorFunc& A,
    const MatrixVectorFunc& A_T,
    const Vector& b,
    Vector& x,
    const MatrixVectorFunc& M,
    const MatrixVectorFunc& M_T,
    const SolverConfig &config
)
{
    SolverResult result;
    auto star_time = std::chrono::high_resolution_clock::now();

    int n = static_cast<int>(b.size());

    // inital some velue we need
    Vector r(n,0.0),r_tilde(n,0.0),p(n,0.0),p_tilde(n,0.0);
    Vector Ap(n,0.0),A_T_P_tilde(n,0.0);
    Vector z(n,0.0),z_tilde(n,0.0);

    // calculate the residual
    A(x,Ap);
    for(int i = 0 ; i < n ; i++)
    {
        r[i] = b[i] - Ap[i];
    }
    
    // initail to solver the matrix 
    copyVector(r,r_tilde);
    result.initail_residual = norm(r);
    if(result.initail_residual < config.tolerance)
    {
        result.converged = true;
        result.iterations = 0;
        result.final_residual = result.initail_residual;
        auto end_time = std::chrono::high_resolution_clock::now();
        result.computer_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>
        (end_time - star_time).count();
    }
    if(config.verbose)
    {
        std::cout<<"PBiCG求解开始,初始残差: " << result.initail_residual<<std::endl;
    }

    double rho_old = 1.0,rho_new,alpha,beta;

    result.converged = false;
    for(int k = 0 ; k < config.max_iterations; k++)
    {
        // solver PBICG in precondition
        M(r,z); // mean M * r = z
        M_T(r_tilde,z_tilde); // mean M_T * r_tilde = z_tilde

        rho_new = dotProduct(z,r_tilde);
        if(std::fabs(rho_new) < 1e-30)
        {
            if(config.verbose)
            {
                std::cout<<"PBiCG中断, rho接近于0" << std::endl;
            }
            break;
        }

        // reflesh the find direction
        if( k == 0)
        {
            copyVector(z,p);
            copyVector(z_tilde,p_tilde);
        }
        else
        {
            beta = rho_new/rho_old;
            for(int i = 0 ; i < n ; i++)
            {
                p[i] = z[i] + beta * p[i];
                p_tilde[i] = z_tilde[i] + beta * p_tilde[i];
            }
        }
        // Matrix product
        A(p,Ap);
        A_T(p_tilde,A_T_P_tilde);

        // reflesh step
        alpha = rho_new/ dotProduct(p_tilde,Ap);

        // reflesh result and residual
        axpy(alpha,p,x);
        for(int i = 0 ; i < n ; i++)
        {
            r[i] -= alpha * Ap[i];
            r_tilde[i] -= alpha * A_T_P_tilde[i];
        }

        double residual = norm(r);
        double rel_residual = residual / result.initail_residual;
        if(rel_residual< config.tolerance)
        {
            result.converged = true;
            result.iterations = k + 1;
            result.final_residual = residual;

            auto end_time = std::chrono::high_resolution_clock::now();
            result.computer_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>
            (end_time - star_time).count();

            if(config.verbose)
            {
                std::cout<<"\nPBiCG收敛于迭代: " << result.iterations
                         <<",最终残差: " <<residual
                         <<",相对残差: " <<rel_residual
                         <<std::endl;
            }
            return result;
        }
        if(config.verbose &&(k + 1)%config.print_interval == 0)
        {
            std::cout<<"迭代:" << k +1 
                     <<",相对残差:" << rel_residual<< std::endl;
        }
        rho_old = rho_new;

    }
    result.converged = true;
    result.iterations = config.max_iterations;
    result.final_residual = norm(r);
    auto end_time = std::chrono::high_resolution_clock::now();
    result.computer_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>
    (end_time - star_time).count();

    if(config.verbose)
    {
        std::cout<<"PBiCG未在最大迭代次数收敛"<< std::endl;
    }
    return result;
    
}

// symmetry precondition deal with
PBiCG::SolverResult PBiCG::solve
(
    const MatrixVectorFunc& A,
    const Vector& b,
    Vector & x,
    const MatrixVectorFunc& M,
    const SolverConfig& config
)
{
    return solve(A,A,b,x,M,M,config);

}

// no precondition deal with
PBiCG::SolverResult PBiCG::solve
(
    const MatrixVectorFunc& A,
    const Vector& b,
    Vector& x,
    const SolverConfig& config
)
{
    auto identity = [](const Vector & v,Vector& result)
    {
        result = v;
    };
    return solve(A,A,b,x,identity,identity,config);
}