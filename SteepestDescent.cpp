#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>

using namespace std;
/*==========最速下降法求解对称正定矩阵==========*/
class SteepestDescent
{
private:
    vector<vector<double>> A;   // 系数矩阵
    vector<double> b;           // 右端向量
    vector<double> x;           // 解向量
    vector<double> r;           // 残差向量
    int max_iter = 1000;               // 最大迭代次数
    double tolerance = 1e-12;           // 收敛容差
    int n;                      // 矩阵的维度
    
public:
    // 构造函数
    SteepestDescent(const vector<vector<double>> & A_matrix,
                    const vector<double> & b_vector,
                    const vector<double> & x_initial = {},
                    int max_iterations = 1000,
                    double tol = 1e-12)
    : A(A_matrix),b(b_vector),max_iter(max_iterations),tolerance(tol)
    {
        n = b.size();

        // 初始化解向量
        if(x_initial.empty())
        {
            x = vector<double>(n,0.0);
        }
        else
        {
            x = x_initial;
        }

        // 检查矩阵的维度
        if(A.size() != n || A[0].size() != n)
        {
            throw invalid_argument("Matrix dimension don't match");
        }

    }

    /*==== 计算向量的点积 ===*/
    double dot_product(const vector<double> & v1,const vector<double> & v2)
    {
        double result = 0.0;
        for(int i = 0 ; i < n ; i++)
        {
            result += v1[i] * v2[i];
        }
        return result;
    }

    /*=== 计算矩阵向量的乘法 ===*/
    vector<double> matrix_product(const vector<vector<double>> & mat,
                                  const vector<double> & vec)
    {
        vector<double> result(n,0.0);
        for(int i = 0 ; i < n ; i++)
        {
            for(int j = 0; j < n ; j++)
            {
                result[i] += mat[i][j] * vec[j];
            }
        }
        return result;
    }

    /*=== 计算向量范数：即模 ===*/
    double vector_norm(const vector<double> vec)
    {
        return sqrt(dot_product(vec,vec));
    }

    /*=== 计算残差 ===*/
    void compute_residual()
    {
        vector<double> Ax = matrix_product(A , x);
        r.resize(n,0.0);
        for(int i = 0; i < n; i++)
        {
            r[i] = b[i] - Ax[i];
        }
    }

    /*=== 主求解函数 ===*/
    void solve()
    {
        compute_residual();
        double residual_norm = vector_norm(r);//计算残差的向量范数

        cout<<"初始残差范数: " << residual_norm << endl;
        for(int iter = 0; iter < max_iter ; iter++)
        {
            // 计算Ar = A * r
            vector<double>  Ar = matrix_product(A,r);

            // 计算步长 alpha = (r^T * r)/(r^T * A *r)
            double r_dot_r = dot_product(r,r);
            double r_dot_Ar = dot_product(r,Ar);
            if (fabs(r_dot_Ar) < 1e-15) {
                std::cout << "警告: 分母接近零，迭代终止" << std::endl;
                break;
            }

            double alpha = r_dot_r /r_dot_Ar;

            // 更新解 x = x + alpha * r
            for(int i = 0 ; i < n ; i++)
            {
                x[i] += alpha * r[i];
            }

            // 更新残差 r = r - alpha * A * r
            // 使用的是递推公式，而不是 r = r - alpha * A * x
            for(int i = 0 ; i < n; i++)
            {
                r[i] -= alpha * Ar[i];
            }
            residual_norm = vector_norm(r);

            // 每100次迭代输出一次进度
            if(iter %100 == 0)
            {
                cout<<"迭代" <<iter<<": 残差范数 = " <<residual_norm<<endl;
            }

            // 检查收敛
            if(residual_norm < tolerance)
            {
                cout<<"收敛于迭代" <<iter <<endl;
                break;
            }
        }

        cout<<"最终残差范数: " << residual_norm<<endl;
    }

    // 获取解
    vector<double>get_solution() const
    {
        return x;
    }

    // 计算误差
    double compute_error() const
    {
        vector<double> Ax(n,0.0);
        for(int i = 0 ; i < n; i++)
        {
            for(int j = 0 ; j < n ;j++)
            {
                Ax[i] += A[i][j] * x[j];
            }
        }
        double error = 0;
        for(int i = 0 ; i < n;i ++)
        {
            error += (Ax[i]  - b[i]) * (Ax[i]  - b[i]);
        }
        return sqrt(error);
    }

};

int main()
{
        // 示例：求解小型对称正定线性方程组
    int n = 3;
    
    // 创建对称正定矩阵
    std::vector<std::vector<double>> A = {
        {4.0, 1.0, 1.0},
        {1.0, 3.0, -1.0},
        {1.0, -1.0, 2.0}
    };
    
    // 创建右端向量
    std::vector<double> b = {6.0, 4.0, 2.0};
    
    // 真实解（用于验证）
    std::vector<double> true_solution = {1.0, 1.0, 1.0};
    
    // 创建最速下降法求解器
    SteepestDescent solver(A, b, std::vector<double>(n, 0.0), 1000, 1e-12);
    
    // 记录开始时间
    auto start = std::chrono::high_resolution_clock::now();
    
    // 求解
    solver.solve();
    
    // 记录结束时间
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // 获取解
    std::vector<double> solution = solver.get_solution();
    
    // 输出结果
    std::cout << "\n解向量: ";
    for (double val : solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    std::cout << "真实解: ";
    for (double val : true_solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    std::cout << "误差: " << solver.compute_error() << std::endl;
    std::cout << "计算时间: " << duration.count() << " 毫秒" << std::endl;
    
    return 0;

}