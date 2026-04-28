// File: linearSystems/krylov/PBiCG/pbigc_solver.cpp
// Category: Krylov / preconditioned iterative solver
// Method: PBiCG
// Purpose: Additional PBiCG-style solver implementation for comparison and learning.

// pbigc_solver.cpp
#include "pbigc_solver.h"
#include <iostream>
#include <chrono>
#include <cmath>

double PBiCGSolver::dotProduct(const Vector& a, const Vector& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

double PBiCGSolver::norm(const Vector& v) {
    return std::sqrt(dotProduct(v, v));
}

void PBiCGSolver::axpy(double alpha, const Vector& x, Vector& y) {
    for (size_t i = 0; i < y.size(); ++i) {
        y[i] += alpha * x[i];
    }
}

void PBiCGSolver::copyVector(const Vector& src, Vector& dst) {
    dst = src;
}

PBiCGSolver::SolverResult PBiCGSolver::solve(
    const MatrixVectorFunc& A,
    const MatrixVectorFunc& A_T,
    const Vector& b,
    Vector& x,
    const MatrixVectorFunc& M,
    const MatrixVectorFunc& M_T,
    const SolverConfig& config) 
    {
    
    SolverResult result;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int n = static_cast<int>(b.size());
    Vector r(n), r_tilde(n), p(n), p_tilde(n);
    Vector Ap(n), A_T_p_tilde(n);
    Vector z(n), z_tilde(n);
    
    // 计算初始残差 r = b - A*x
    A(x, Ap);
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - Ap[i];
    }
    copyVector(r, r_tilde);
    
    result.initial_residual = norm(r);
    if (result.initial_residual < config.tolerance) {
        result.converged = true;
        result.iterations = 0;
        result.final_residual = result.initial_residual;
        auto end_time = std::chrono::high_resolution_clock::now();
        result.compute_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time).count();
        return result;
    }
    
    if (config.verbose) {
        std::cout << "PBiCG求解开始，初始残差: " << result.initial_residual << std::endl;
    }
    
    double rho_old = 1.0, rho_new, alpha, beta;
    result.converged = false;
    
    for (int k = 0; k < config.max_iterations; ++k) {
        // 预处理步骤
        M(r, z);
        M_T(r_tilde, z_tilde);
        
        rho_new = dotProduct(z, r_tilde);
        
        if (std::fabs(rho_new) < 1e-30) {
            if (config.verbose) {
                std::cout << "PBiCG中断: rho接近零" << std::endl;
            }
            break;
        }
        
        if (k == 0) {
            copyVector(z, p);
            copyVector(z_tilde, p_tilde);
        } else {
            beta = rho_new / rho_old;
            for (int i = 0; i < n; ++i) {
                p[i] = z[i] + beta * p[i];
                p_tilde[i] = z_tilde[i] + beta * p_tilde[i];
            }
        }
        
        A(p, Ap);
        A_T(p_tilde, A_T_p_tilde);
        
        alpha = rho_new / dotProduct(p_tilde, Ap);
        
        // 更新解和残差
        axpy(alpha, p, x);
        for (int i = 0; i < n; ++i) {
            r[i] -= alpha * Ap[i];
            r_tilde[i] -= alpha * A_T_p_tilde[i];
        }
        
        double residual = norm(r);
        double rel_residual = residual / result.initial_residual;
        
        if (rel_residual < config.tolerance) {
            result.converged = true;
            result.iterations = k + 1;
            result.final_residual = residual;
            
            auto end_time = std::chrono::high_resolution_clock::now();
            result.compute_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time).count();
            
            if (config.verbose) {
                std::cout << "\nPBiCG收敛于迭代 " << result.iterations 
                         << ", 最终残差: " << residual 
                         << ", 相对残差: " << rel_residual << std::endl;
            }
            return result;
        }
        
        if (config.verbose && (k + 1) % config.print_interval == 0) {
            std::cout << "迭代 " << k + 1 
                     << ", 相对残差: " << rel_residual << std::endl;
        }
        
        rho_old = rho_new;
    }
    
    result.converged = false;
    result.iterations = config.max_iterations;
    result.final_residual = norm(r);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.compute_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time).count();
    
    if (config.verbose) {
        std::cout << "PBiCG未在最大迭代次数内收敛" << std::endl;
    }
    
    return result;
}

PBiCGSolver::SolverResult PBiCGSolver::solve(
    const MatrixVectorFunc& A,
    const Vector& b,
    Vector& x,
    const MatrixVectorFunc& M,
    const SolverConfig& config) {
    
    // 对称预处理：M_T = M
    return solve(A, A, b, x, M, M, config);
}

PBiCGSolver::SolverResult PBiCGSolver::solve(
    const MatrixVectorFunc& A,
    const Vector& b,
    Vector& x,
    const SolverConfig& config) {
    
    // 无预处理：M = I (单位矩阵)
    auto identity = [](const Vector& v, Vector& result) {
        result = v;
    };
    
    return solve(A, A, b, x, identity, identity, config);
}
