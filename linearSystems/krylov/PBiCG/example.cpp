// example.cpp
// File: linearSystems/krylov/PBiCG/example.cpp
// Category: Krylov / preconditioned iterative solver
// Method: PBiCG example
// Purpose: Demonstrate usage of an alternate PBiCG solver interface on a test matrix.

#include "PBiCGsolver.h"
#include <iostream>
#include <vector>
#include <cmath>

// 创建非对称稀疏矩阵函数
class NonSymmetricMatrix {
private:
    int n;
    std::vector<std::vector<double>> data;
    
public:
    NonSymmetricMatrix(int size) : n(size), data(size, std::vector<double>(size, 0.0)) {
        // 创建非对称三对角矩阵
        for (int i = 0; i < n; ++i) {
            // 主对角线
            data[i][i] = 4.0;
            
            // 上对角线和下对角线（非对称）
            if (i > 0) {
                data[i][i-1] = -1.0 - 0.1 * i;  // 下对角线
            }
            if (i < n - 1) {
                data[i][i+1] = -1.0 + 0.1 * i;  // 上对角线
            }
            
            // 添加一些非对称元素
            if (i > 1) {
                data[i][i-2] = 0.2;
            }
            if (i < n - 2) {
                data[i][i+2] = -0.2;
            }
        }
    }
    
    void operator()(const std::vector<double>& x, std::vector<double>& result) const {
        result.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (std::fabs(data[i][j]) > 1e-12) {
                    result[i] += data[i][j] * x[j];
                }
            }
        }
    }
    
    // 转置矩阵向量乘积
    void transposeMul(const std::vector<double>& x, std::vector<double>& result) const {
        result.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (std::fabs(data[i][j]) > 1e-12) {
                    result[j] += data[i][j] * x[i];  // 注意索引交换
                }
            }
        }
    }
    
    // 获取矩阵大小
    int size() const { return n; }
    
    // 获取对角线元素（用于Jacobi预处理器）
    std::vector<double> getDiagonal() const {
        std::vector<double> diag(n);
        for (int i = 0; i < n; ++i) {
            diag[i] = data[i][i];
        }
        return diag;
    }
};

// Jacobi预处理器
class JacobiPreconditioner {
private:
    std::vector<double> inv_diag;
    
public:
    JacobiPreconditioner(const NonSymmetricMatrix& A) {
        inv_diag = A.getDiagonal();
        for (size_t i = 0; i < inv_diag.size(); ++i) {
            if (std::fabs(inv_diag[i]) > 1e-12) {
                inv_diag[i] = 1.0 / inv_diag[i];
            } else {
                inv_diag[i] = 1.0;  // 避免除零
            }
        }
    }
    
    void operator()(const std::vector<double>& r, std::vector<double>& z) const {
        z.resize(r.size());
        for (size_t i = 0; i < r.size(); ++i) {
            z[i] = inv_diag[i] * r[i];
        }
    }
};

int main() {
    const int n = 1000;
    
    // 创建矩阵
    NonSymmetricMatrix A(n);
    
    // 创建右侧向量
    std::vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        b[i] = std::sin(i * 0.1);
    }
    
    // 初始解
    std::vector<double> x(n, 0.0);
    
    // 创建预处理器
    JacobiPreconditioner jacobi_precond(A);
    
    // 配置求解器
    // PBiCGSolver solver;
    // PBiCGSolver::SolverConfig config;

    PBiCGsolver solver;
    PBiCGsolver::SolverConfig config;  
    config.tolerance = 1e-8;
    config.max_iterations = 10000;
    config.verbose = true;
    config.print_interval = 10;
    
    // 使用lambda函数包装矩阵向量乘积
    auto matVec = [&A](const std::vector<double>& v, std::vector<double>& result) {
        A(v, result);
    };
    
    auto matVecTranspose = [&A](const std::vector<double>& v, std::vector<double>& result) {
        A.transposeMul(v, result);
    };
    
    // 求解
    auto result = solver.solve(matVec, matVecTranspose, b, x, 
                               jacobi_precond, jacobi_precond, config);
    
    // 输出结果
    std::cout << "\n=== 求解结果 ===" << std::endl;
    std::cout << "收敛: " << (result.converged ? "是" : "否") << std::endl;
    std::cout << "迭代次数: " << result.iterations << std::endl;
    std::cout << "初始残差: " << result.initial_residual << std::endl;
    std::cout << "最终残差: " << result.final_residual << std::endl;
    std::cout << "计算时间: " << result.compute_time_ms << " ms" << std::endl;
    
    // 验证解
    std::vector<double> Ax(n);
    A(x, Ax);
    
    double max_error = 0.0;
    for (int i = 0; i < n; ++i) {
        double error = std::fabs(Ax[i] - b[i]);
        if (error > max_error) {
            max_error = error;
        }
    }
    
    std::cout << "\n验证最大误差: " << max_error << std::endl;
    
    // 显示部分解
    std::cout << "\n前10个解的值:" << std::endl;
    for (int i = 0; i < std::min(20, n); ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    
    return 0;
}
