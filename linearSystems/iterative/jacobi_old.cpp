// File: linearSystems/iterative/jacobi_old.cpp
// Category: Linear system iterative method
// Method: Jacobi iteration
// Purpose: Preserve an older Jacobi implementation for comparison with newer revisions.

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

class JacobiSolver
{
private:
    /* data */
    vector<vector<double>> A;
    vector<double> b;
    vector<double> x;
    vector<double>x_old;
    double tolerance;
    int max_iterations;
    int n;              //矩阵的维度
public:
    // 构造函数
    JacobiSolver(const vector<vector<double>>& coeff_matrix, 
                 const vector<double>& right_hand_side,
                 int max_iter = 1000, 
                 double tol = 1e-6)
        : A(coeff_matrix), b(right_hand_side), max_iterations(max_iter), tolerance(tol) {
        
        n = A.size();
        x.resize(n, 0.0);      // 初始解设为0
        x_old.resize(n, 0.0);
        
        // 检查矩阵是否为方阵
        if (A.size() != b.size()) {
            throw invalid_argument("Matrix dimensions do not match!");
        }
        
        for (int i = 0; i < n; i++) {
            if (A[i].size() != n) {
                throw invalid_argument("Matrix is not square!");
            }
        }
    }
    // 检查对角线元素是否为0（雅可比迭代的必要条件）
    bool checkDiagonalDominance()
    {
        for (int i = 0; i < n; i++)
        {
            if(fabs(A[i][i]) < 1e-10) //检查对角线元素小于某一个值可认为是0
            {
                return false;
            }
        }
        return true;                    
    }
    // 执行一次雅可比迭代
    void iterate()
    {
        x_old = x;
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            // 计算（L+U * x
            for (int j = 0; j < n; j++)
            {
                if(j != i)
                {
                    sum += A[i][j] * x_old[j];
                }
                
            }
            x[i] = (b[i] - sum) / A[i][i];
            
        }
        
    }
    // 计算当前误差
    double calculateError()
    {
        double max_error = 0.0;
        for (int i = 0; i < n ; i++)
        {
            double error = fabs(x[i] - x_old[i]);
            if (error > max_error)
            {
                error = max_error;
            }
            
        }
        return max_error;
    }

    // 求解线性方程组
    pair<vector<double>,int> solver(const vector<double>& initial_guass = {})
    {
        // 如果提供了初始的猜测值，就shiyongta
        if(!initial_guass.empty())
        {
            if(initial_guass.size() != n)
            {
                throw invalid_argument("Initial guass size does not match matrix dimensions");
            }
            x = initial_guass;
        }
        if(!checkDiagonalDominance())
        {
            cout<<"Waring:Diagonal elements contain zeros or vert small values !" << endl;
        }

        int iteration = 0;
        double error = tolerance ; //保证迭代一次
        cout<<fixed<<setprecision(6) << endl;;
        cout << "Start Jacobi Iteration ..." << endl;
        cout <<"Iteration\tMax Error" << endl;
        cout <<"----------------------"<< endl;
        while (iteration < max_iterations && error > tolerance)
        {
            iterate();
            error = calculateError();
            iteration++;

            if(iteration % 100 == 0)
            {
                cout <<iteration <<"\t\t" << error << endl;
            }
        }
        cout <<"---------------------" << endl;
        cout <<"Converged after " << iteration <<" iteration ." << endl;
        cout << "Final error: " << error <<endl;

        return make_pair(x,iteration);
        
    }

    // 打印解向量
    void printSolution() const
    {
        cout<<"Solution vector : " << endl;
        for (int i = 0; i < n; i++)
        {
            cout<<"x[" <<i << "]"  <<fixed<<setprecision(10) <<x[i] << endl;
        }
        
    }
                    
};


// 示例：求解线性方程组
int main()
{
    // 示例1：3x3方程组
    // 10x + 2y + z = 9
    // 2x + 20y - 3z = -44
    // -2x + 3y + 10z = 22
    vector<vector<double>> A1 = {
        {10, 2, 1},
        {2, 20, -3},
        {-2, 3, 10}
    };
    
    vector<double> b1 = {9, -44, 22};
    
    cout << "========== Example 1 ==========" << endl;
    JacobiSolver solver1(A1, b1, 1000, 1e-10);
    
    try {
        auto result1 = solver1.solver();
        solver1.printSolution();
        
        // 验证解
        cout << "\nVerification:" << endl;
        vector<double> x = result1.first;
        for (int i = 0; i < 3; i++) {
            double sum = 0;
            for (int j = 0; j < 3; j++) {
                sum += A1[i][j] * x[j];
            }
            cout << "Equation " << i+1 << ": " << sum << " ≈ " << b1[i] 
                 << " (error: " << fabs(sum - b1[i]) << ")" << endl;
        }
    } catch (const exception& e) {
        cout << "Error: " << e.what() << endl;
    }

}
