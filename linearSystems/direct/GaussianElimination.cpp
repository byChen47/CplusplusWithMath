// 高斯消元法

// File: linearSystems/direct/GaussianElimination.cpp
// Category: Linear system direct method
// Method: Gaussian elimination
// Purpose: Solve dense linear systems and verify the resulting solution vector.

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdexcept>
using namespace std;

int findPivotRow(const vector<vector<double>> & matrix,int col ,int n);
void printMatrix(const vector<vector<double>> & matrix);
bool verifySolution(const vector<vector<double>>& A,
                    const vector<double>& b,
                    const vector<double>& x);

/**
 * 列主元法：寻找第col列中绝对值最大的行（从第row行开始）
 */

int findPivotRow(const vector<vector<double>> & matrix ,int col ,int n)
{
    int maxRow = col;
    double maxValue = abs(matrix[col][col]);
    for (int i = col+1; i < n; i++)
    {
        if(abs(matrix[i][col]) > maxValue)
        {
            maxValue = abs(matrix[i][col]);
            maxRow = i;
        }
    }
    return maxRow;

}



bool verifySolution(const vector<vector<double>> & A,const vector<double> & b,const vector<double> & x)
{
    static constexpr double EPS = 1e-10;
    int n = A.size();
    for (int  i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j  < n ;j++)
        {
            sum += A[i][j] * x[j];
        }

        if(abs(sum - b[i]) > EPS)
        {
            cout << "第 " << i + 1 << " 个方程验证失败"
                 << sum  << " != " << b[i]
                 << " (差值： " << sum - b[i] << ")" << endl;
            return false;
        }

        
    }
    return true;
}



class GaussianElimination
{
private:
    /* data */

    static constexpr double EPS = 1e-10;  //容差，用于判断是否为0

public:

    /**
     * 高斯消元法求解线性方程组 Ax = b
     * @param A 系数矩阵 (n x n)
     * @param b 常数向量 (n)
     * @return 解向量 x
     */
    static void printMatrix(const vector<vector<double>> & matrix)
    {
        int n = matrix.size();
        int m = matrix[0].size();
        
        cout << "增广矩阵" <<endl;

        for (int i = 0; i < n; i++)
        {
            cout << "| " << endl;
            for (int  j = 0; j < m ; j++)
            {
                if(j == m-2)
                {
                    cout << setw(8) << fixed << setprecision(4) << matrix[i][j] << " |";
                }
                else
                {
                    cout << setw(8) << fixed << setprecision(4) << matrix[i][j] << " ";
                }
                cout << "|" << endl;
            }
            cout << endl;
            
        }
        
    }

     static vector<double> solver(vector<vector<double>> A,vector<double> b) // 定义矩阵A和b的值 
     {
        int n = A.size();

        // 1.0 构造增广矩阵
        vector<vector<double>> augmented(n,vector<double>(n+1.0,0.0)); // n行 n+1列
        for (int  i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                augmented[i][j] = A[i][j]; // 拿到A矩阵
            }
            augmented[i][n] = b[i];  // 添加b的值
            
        }

        // 2.0 前向消元（化为上三角矩阵）
        for (int  k = 0; k < n; k++)
        {
            // 2.1 列主元法，寻找到k列中绝对值最大的行
            int maxRow = findPivotRow(augmented,k,n);

            // 2.2 交换当前行与主元行
            if(maxRow != k)
            {
                swap(augmented[k] ,augmented[maxRow]);
            }
            // 2.3 检查主元是否为0（奇异矩阵，det = 0，不可逆; det != 0,可逆，非奇异矩阵）

            if(abs(augmented[k][k]) < EPS)
            {
                throw runtime_error(" 矩阵奇异,无唯一解");
            }

            // 2.4 归一化，将最大值，也就是主元化为1
            double pivot =  augmented[k][k];
            for (int j = k ; j < n ; j++)
            {
                augmented[k][j] /= pivot;
            }

            // 2.5 消元； 将下方行的第 k 列消为0
            for (int i = k + 1 ;  i < n ; i ++)
            {
                double factor = augmented[i][k];
                for (int  j = k +1 ;j <= n; j++)
                {
                    augmented[i][j] -= factor * augmented[k][j]; 
                }
                 // 注意：经过归一化后，augmented[i][i] = 1，所以不需要再除以它
            }
                       
        }

        // 3.0 回代求解
        vector<double> x(n,0.0); // 定义一个n哼1列的未知向量
        for (int i = n-1; i >= 0 ; i--)
        {
            x[i] =  augmented[i][n];
            for (int j =  i + 1; j < n ; j ++)
            {
                x[i] -= augmented[i][j] * x[j];
            }
            
        }

        return x;
        
     }

    static bool vertifySolution(const vector<vector<double>> & A, 
                                const vector<double> & b,
                                const vector<double> & x)
    {
        int n =  A.size();
        static constexpr double eps = 1e-10;
        for (int i = 0; i < n ; i ++)
        {
            double sum = 0;
            for (int j = 0 ; j < n ; j++)
            {
                sum += A[i][j] * x[j];
            }
            if (abs(sum - b[i]) > eps)
            {
                cout << "第 " << i + 1<< "个方程验证失败:"
                    << sum   << "!= " << b[i]
                    << "差值为： " << sum -b[i] << endl;
                return false;
            }
            
            
        }
        return true;
        
    };     
};



void testGaussianElimination()
{
    vector<vector<double>> A = 
    {
        {2 ,  1 ,-1},
        {-3 ,-1 , 2},
        {-2 , 1 , 2}
    };

    vector<double> b = {8 , -11 ,-3};
    cout << "========== 高斯消元法求解线性方程组 ==========" << endl;
    cout << "方程组：" << endl;
    cout << "  2x +  y -  z =  8" << endl;
    cout << " -3x -  y + 2z = -11" << endl;
    cout << " -2x +  y + 2z =  -3" << endl;
    cout << endl;

    vector<vector<double>> augmented(3,vector<double>(4, 0.0));
    for (int i = 0; i < 3; i ++)
    {
        for (int j = 0 ; j < 4; j++)
        {
            augmented[i][j] = A[i][j];
        }
        augmented[i][3] = b[i];
        
    }

    GaussianElimination::printMatrix(augmented);

    try
    {
        // 求解
        vector<double> solution = GaussianElimination::solver(A,b);

        // 打印结果

        cout << "解向量" <<endl;
        for (int i = 0; i < solution.size(); i++)
        {
            cout << "x" <<i +1 << "=" << solution[i] << endl;
            cout << endl;
        }
        
        //验证
        cout << "验证结果： " << endl;
        if( GaussianElimination::vertifySolution(A,b,solution))
        {
            cout << "所有方程验证通过" << endl;
            
        } 
        else
        {
            cout << "验证失败" << endl;
        }
        // 计算残差
        cout<<"残差向量" << endl;
        for (int i = 0; i< A.size(); i++)
        {
            double sum = 0.0;
            for (int j = 0; j < A[i].size(); j++)
            {
                sum += A[i][j] * solution[j];
            }
            cout << "方程 " << i + 1<< ": " << sum-b[i] << endl;
        }
        
    }
    catch(const exception & e)
    {
        cerr<< "错误" << e.what() << endl;
    }
    

}

int main()
{
    testGaussianElimination();
}
