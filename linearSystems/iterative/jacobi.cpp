// File: linearSystems/iterative/jacobi.cpp
// Category: Linear system iterative method
// Method: Jacobi iteration
// Purpose: Demonstrate classical fixed-point iteration for linear systems.

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

class jacobiSolver
{
private:
    
    vector<vector<double>> A;
    vector<double> b;
    vector<double> x;
    vector<double> x_old;
    int n;
    int max_iterations;
    double tolerance = 1e-12;

public:
    // 雅可比矩阵基本设置和判断
    jacobiSolver(const vector<vector<double>> & coeff_matrix,
                 const vector<double>& right_hand_side,
                 int max_iter = 1000,
                 double tol = 1e-10)
    : A(coeff_matrix),b(right_hand_side),max_iterations(max_iter),tolerance(tol)
    {
        n = A.size();
        x.resize(n,0.0);
        x_old.resize(n,0.0);

        //检查矩阵是否为方阵
        if(A.size() != b.size())
        {
            throw invalid_argument("Matrir dimensions donot match");
        }
        for(int i = 0; i < n ; i++)
        {
            if(A[i].size() != n)
            throw invalid_argument("Matrix is no square !");
        }
    }

    // 检查对角线元素是否为零（雅可比迭代的必要条件）
    bool checkDDiagonalDominance()
    {
        for(int i = 0 ; i < n ; i++)
        {
            if(fabs(A[i][i]) < 1e-12)
            {
                return false;
            }
        }
        return true;
    }

    //执行一次雅可比迭代
    void iterate()
    {
        x_old = x;
        for(int i = 0;i <n ;i++)
        {
            double sum = 0.0;
            for(int j = 0 ; j < n ; j++)
            {
                if(j != i)
                {
                    sum += A[i][j] * x_old[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
    }

    // 计算当前误差（无穷范数)
    // 当前x 和 x_old之间的误差
    double calculateError()
    {
        double max_error = 0.0;
        for(int i = 0 ; i < n ; i++)
        {
            double error = fabs(x[i] - x_old[i]);
            if(error > max_error)
            {
                max_error = error;
            }
        }
        return max_error;
    }

    // 求解线性方程组
    pair<vector<double>,int> solver(const vector<double>& initial_guess = {})
    {
        // 看看有没有提供最初的猜测值，如果提供就使用它
        if(!initial_guess.empty())
        {
            if(initial_guess.size() != n)
            {
                throw invalid_argument("Initial size does not match matrix dimensions");
            }
            x = initial_guess;
        }
        if(!checkDDiagonalDominance())
        {
            cout<< "Warning : Diagonal elements contents zeros" << endl;
        }
        int iteration = 0;
        double error = tolerance + 1.0 ; // 确保至少迭代一次

        cout<< fixed << setprecision(6)<< endl;
        cout <<"Starting Jacobi iteration ..." << endl;
        cout << "Iteration\tMax Error" << endl;
        cout <<"-----------------------" << endl;

        while (iteration < max_iterations && error > tolerance)
        {
            iterate();
            error  = calculateError();
            iteration ++;

            if(iteration % 100 == 0)
            {
                cout<< iteration << "\t\t" << error<< endl;
            }
        }

        cout <<"----------------------" << endl;
        cout <<"Converged after :" << iteration <<"Iteration." << endl;
        cout <<"Final Error is : " << error << endl;

        return make_pair(x,iteration);
        
    }

    void printSolution()
    {
        cout << "Solition Vector: " << endl;
        for(int i = 0 ; i < n; i ++)
        {
            cout<<"x[" <<i <<"]" << fixed <<setprecision(10) <<x[i] << endl;
        }
    }


};

int main()
{
     // 示例1：3x3方程组
    // 10x + 2y + z = 9
    // 2x + 20y - 3z = -44
    // -2x + 3y + 10z = 22
    vector<vector<double>> A1 = {
        {15, 2, 1},
        {2, 6, -3},
        {-2, 3, 10}
    };

    vector<double> b1 = {9,-33,12};

    cout<<"==========Example 1=============" << endl;

    jacobiSolver solver1(A1,b1,1000,1e-10);

    try
    {
        auto result1 = solver1.solver();
        solver1.printSolution();

        // 验证解
        cout<<"\nVertification: " << endl;
        vector<double> x = result1.first;
        for(int i = 0; i < 3; i++)
        {
            double sum = 0;
            for(int j = 0 ; j < 3 ; j++)
            {
                sum += A1[i][j] * x[j];
            }
            cout <<"Equation " <<i +1 << ": " << sum << " = " << b1[i]
                 << "(Error: " << fabs(sum-b1[i]) << ")" << endl;
        }
    }
    catch(const std::exception& e)
    {
        std::cerr <<"Error "<< e.what() << '\n';
    }

    // 示例2 ：4×4方程组
    vector<vector<double>> A2 = 
    {
        {10 ,-1 , 2 ,0},
        {-1 ,11 , -1, 3},
        {2 ,-1 , 10 ,-2},
        {0 ,3 , -1 ,8}
    };

    vector<double> b2 = {6,25,-11,15};

    cout<<"=============Example 2===========" << endl;
    jacobiSolver solver2(A2,b2,1000,1e-12);
    try
    {
        // 使用初始猜测值
        // vector<double> initial_guess = {0.5,0.8,2.0,-1};
        auto result2 = solver2.solver();
        solver2.printSolution();

        cout<<"\nVertification: " << endl;
        vector<double> x = result2.first;
        for(int i = 0; i < 4; i++)
        {
            double sum = 0;
            for(int j = 0 ; j < 4 ; j++)
            {
                sum += A2[i][j] * x[j];
            }
            cout <<"Equation " <<i +1 << ": " << sum << " = " << b2[i]
                 << "(Error: " << fabs(sum-b2[i]) << ")" << endl;
        }
    }
    catch(const std::exception& e)
    {
        std::cerr<< "Error" << e.what() << '\n';
    }
    
    return 0;
    
}
