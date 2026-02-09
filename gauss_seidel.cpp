#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

class Gauss_Seidel
{
private:

    vector<vector<double>> A;  // 系数矩阵
    vector<double> b;          // 常数向量
    vector<double> x;          // 解向量
    int n;                     // 矩阵得维度
    int maxIter = 1000;               // 最大迭代次数
    double tol = 1e-12;                   // 容差


public:
    // 析构函数
    Gauss_Seidel(const vector<vector<double>> & coeffMatrix,
                 const vector<double> & constVector,
                 double tolerance = 1e-12,int max_Iteration = 1000)
                : A(coeffMatrix),b(constVector),tol(tolerance),maxIter(max_Iteration)
    {
        n = A.size();
        x.resize(n,0.0);
    }

    // 检查矩阵是否对角线占优
    bool isDiagonallyDominant()
    {
        for(int i = 0; i < n ; i++)
        {
            double sum = 0.0;
            for(int j = 0 ; j < n ; j++)
            {
                if(j != i)
                {
                    sum += fabs(A[i][j]);
                }
            }
            if(fabs(A[i][i] <= sum))
            {
                return false;
            }
        }

        return true;
    }

    // 高斯赛德尔求解
    bool solve()
    {
        if(!isDiagonallyDominant())
        {
            cout <<"警告：矩阵可能不是对角占优得，可能不收敛 !" <<endl;
        }

        vector<double> xOld(n,0.0);
        int iter = 0;
        double error = 0.0;

        cout<<"迭代过程： " << endl;
        cout << "迭代次数\t解向量\t\t\t\t误差" << endl;
        cout <<"--------------------------------" << endl;

        do
        {
            xOld = x; // 保存上一次迭代
            
            // 更新每个分量
            for(int i = 0; i < n ; i ++)
            {
                double sum = 0.0;
                for(int j = 0; j < n ; j++)
                {
                    if(j != i)
                    {
                        sum += A[i][j] * x[j];
                    }
                }
                // 更新x[i],注意避免除以0
                if(A[i][i] != 0)
                {
                    x[i] = (b[i] - sum) / A[i][i];
                }
                else
                {
                    cout<<"错误,对角线上存在元素0 ! " << endl;
                }
            }

            // 计算最大误差
            error = 0.0;
            for(int i = 0; i < n ; i++)
            {
                error = max(error ,fabs(x[i] - xOld[i]));
            }

            // 输出当前迭代结果
            cout<<setw(4)<<iter <<"\t\t";
            for(int i = 0; i < n ;i++)
            {
                cout<< fixed<<setprecision(6) <<x[i] <<" ";
            }
            cout <<"\t" <<error <<endl;
            iter++;
            
        } while (error > tol && iter < maxIter);
        
        if(iter > maxIter)
        {
            cout<<"达到最大迭代次数,可能为收敛! " <<endl;
            return false;
        }
        cout<<"\n收敛于" <<iter <<"次迭代" << endl;
        return true;

    }

    // 获取解向量
    vector<double> getSolution() const
    {
        return x;
    }

    // 打印解
    void printSolution() const
    {
        cout<<"\n解向量: " << endl;
        for(int i = 0; i <n ; i++)
        {
            cout<<"x[" <<i<<"] = " <<fixed<<setprecision(10)<< x[i] <<endl; 
        }
    }

    // 验证解
    void verifySolution() const
    {
        cout<<"\n验证解(Ax  -  b) : " << endl;
        for(int i = 0; i <n ; i++)
        {
            double sum = 0;
            for(int j = 0; j <n ; j++)
            {
                sum += A[i][j] * x[j];
            }
            double residual = 0.0;
            residual = sum - b[i];
            cout <<"方程 " <<i <<": " << fixed <<setprecision(10) <<residual <<endl;
        }
    }

};

// 测试用例1：对角线占优得3×3系统
void testCase1()
{
    cout <<"======测试用例1: 3 x 3线性方程组====" << endl;
    vector<vector<double>> A = {
        {2 , 1 ,1},
        {1 , 6 ,1},
        {1 , -1 ,8}
    };
    vector<double> b = {7,9,5};
    Gauss_Seidel solver(A,b,1e-12,1000);
    if(solver.solve())
    {
        solver.printSolution();
        solver.verifySolution();
    }
}

// 测试用例2: 4 × 4线性系统
void testCase2()
{
    cout <<"======测试用例2: 4x4 线性系统" << endl;
    vector<vector<double>> A = {
        {10 ,-1 , 2 , 0},
        {-1 ,11 ,-1 , 3},
        { 2 ,-1 ,10 ,-1},
        { 0 , 3 , -1, 8}
    };
    vector<double> b = {6,25,-11,15};

    Gauss_Seidel solver(A,b,1e-12,1000);
    if (solver.solve())
    {
        solver.printSolution();
        solver.verifySolution();
    }
    
}

//用户输入自定义系统
void userInputSystem()
{
    cout <<"\n====自定义线性方程组====" << endl;

    int n;
    cout<<"请输入方程组得维度: " ;
    cin>> n;

    vector<vector<double>> A(n,vector<double>(n));
    vector<double> b(n);

    cout<<"请输入系数矩阵A: " << endl;
    for(int i = 0 ; i < n ; i++)
    {
        cout<<"第" <<i +1<<"行（"<< n <<"）个系数: ";
        for(int j = 0; j <n ;j++)
        {
            cin>>A[i][j];
        }
    }

    cout<<"请输入常数向量b (" << n <<")个值：";
    for(int i = 0; i < n ; i++)
    {
        cin>>b[i];
    }

    double tol ;
    int maxIter ;
    cout <<"请输入容差(默认1e-12): ";
    cin >>tol;
    cout << endl;

    cout<<"请输入最大迭代次数(默认1000): ";
    cin >> maxIter;
    cout << endl;
    Gauss_Seidel solver(A ,b ,tol ,maxIter);
    if(solver.solve())
    {
        solver.getSolution();
        solver.printSolution();
        solver.verifySolution();
    }
}

int main()
{
    int choice;
    cout<<"高斯-赛德尔迭代法求解线性方程组 " << endl;
    cout <<"=============================" << endl;
    cout<<" 1.运行测试用例1(3x3系统) " << endl;
    cout<<" 2.运行测试用例2(4x4系统) " << endl;
    cout<<" 3.输入自定义系统" << endl;
    cout<<" 请选择： ";
    cin >> choice;
    switch (choice)
    {
    case 1:
        testCase1();
        break;
    case 2:
        testCase2();
        break;
    case 3:
        userInputSystem();
        break;
    
    default:
        cout<<"无效选择! " <<endl;
        
    }

    return 0;
}