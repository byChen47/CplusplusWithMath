#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>


class LUDecomposition
{
private:
    std::vector<std::vector<double>> L;
    std::vector<std::vector<double>> U;


public:

    // 普通LU分解(要求矩阵顺序主子式不为0，即矩阵为非奇异阵，可逆)
    bool decompose(const std::vector<std::vector<double>> & A)
    {
        // 检查一下矩阵是否为空
        if(A.empty())
        {
            std::cerr<<"Error:Matrix is empty" << std::endl;
            return false;
        }

        size_t n = A.size();
        // 检查是否为方阵
        for (size_t i = 0; i < n; i++)
        {
            if (A[i].size() != n)
            {
                std::cerr<< "Error: Matrix must be square! Row" << i
                         << "has size " << A[i].size() << ",expected"
                         << n << std::endl;
                return false;
            }
            
        }

        // 检查完毕后应该要确定矩阵A没问题了
        // 可以开始进行LU分解了
        // 初始化L和U
        L = std::vector<std::vector<double>>(n, std::vector<double>(n,0.0));
        U = std::vector<std::vector<double>>(n, std::vector<double>(n,0.0));

        for (size_t i = 0; i < n; i++)
        {
            // 计算U的第i行
            for (size_t k = i; k < n; k++)
            {
                double sum = 0.0;
                for (size_t j = 0; j < i; j++)
                {
                    sum += L[i][j] * U[j][k];
                }
                U[i][k] = A[i][k] -sum;

                // 检查顺序主子式是否为0

                if(i == k && std::abs(U[i][i] < 1e-12))
                {
                    std::cerr<< "Error: Zero pivot encounterd at position <"
                             << i <<"," << i <<")" << std::endl;
                    std::cerr<<"This means a leading pricipal minor is zero ." 
                             << std::endl;
                    return false;
                }
                
            }
            
            // 计算L的第i列
            for (size_t k = i; k < n; k++)
            {
                if(i== k)
                {
                    L[i][i] = 1.0; // 对角线为1
                }
                else
                {
                    double sum = 0.0;
                    for (size_t j = 0; j < i; j++)
                    {
                        sum += L[k][j] * U[j][i];
                    }
                    L[k][i] = (A[k][i] - sum) /U[i][i];
                    
                }
            }
            
        }
        return true;
    }


    // 设置打印矩阵样式
    void printMatrix(const std::vector<std::vector<double>> & mat,
                     const std::string & name)
    {
        std::cout << name <<":" << std::endl;
        for (const auto& row: mat)
        {
            for (const auto & val: row)
            {
                std::cout << std::setw(10) << std::fixed
                          << std::setprecision(4) << val<< " ";
            }
            std::cout << std::endl;
            
        }

    }                     
    // 打印L和U矩阵
    void printResult()
    {
        printMatrix(L,"L matrix");
        printMatrix(U,"U matrix");
    }

    // 验证分解结果L和U相乘 - A应该小于1e-12，认为是准确的
    bool vertify(const std::vector<std::vector<double>> &A)
    {
        // 判断矩阵是否为空
        if (A.empty() || L.empty() || U.empty())
        {
            return false;
        }

        size_t n = A.size();
        double maxError = 0.0;

        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                double sum = 0.0;
                for (size_t k = 0; k < n; k++)
                {
                    sum += L[i][k] * U[k][j];
                }
                maxError = std::max(maxError,std::abs(A[i][j] - sum));
                
            }
            
        }
        std::cout << "Maximun reconstrution error: " << std::scientific
                  <<maxError << std::endl;
        return maxError < 1e-6;
    }

    // 解线性方程组 Ax = b
    std::vector<double> solver(const std::vector<double> & b)
    {
        //判断线性方程组的形式
        if (L.empty()|| U.empty() || L.size() != b.size())
        {
            std::cerr<< "Error : invalid input for solver" << std::endl;
        }
        size_t n = b.size();
        std::vector<double> y(n,0.0);
        std::vector<double> x(n,0.0);
        
        // 前向代入Ly = b;
        for (size_t i = 0; i < n ; i++)
        {
            double sum = 0.0;
            for (size_t j = 0 ; j < i ; j++)
            {
                sum += L[i][j] * y[j];
            }
            y[i] = b[i] - sum;
        }
        // 后向带入 Ux = y
        for(size_t i = n; i-- > 0;)
        {
            // 安全的向后循环
            double sum = 0.0;
            for (size_t j = i + 1; j < n ; j++)
            {
                sum +=U[i][j] * x[j]; 
            }
            x[i] = (y[i] - sum) / U[i][i];
            
        }

        return x;
    }

    // 计算行列式和获取顺序主子式函数主要是为了确定矩阵是正定的
    // A = LU |A| = det(L) * det(U)
    // det(L) = 1..0;det(U) = u[i][i] * u[i+1][i+1]...*u[n][n];
    // 计算行列式
    double determinat()
    {
        if(U.empty())
        {
            return 0.0;
        }
        double det = 1.0;
        for(size_t i = 0; i < U.size(); i++)
        {
            det *= U[i][i]; //U的对角线成绩几位行列式值
        }
        return det;
    }
    // 获取顺序主子式
    std::vector<double> getLeadingPrincipalMinors()
    {
        std::vector<double> minors;// 存储各阶顺序主子式的值
        if(U.empty())
        {
            return minors;
        }
        for(size_t k = 1; k <= U.size(); k++)
        {
            double minor = 1.0;
            for(size_t i = 0 ; i < k ;i++)
            {
                minor *= U[i][i];// 计算前k个对角线元素的乘积
            }
            minors.push_back(minor);
        }

        return minors;
    }
    
};
    // 辅助函数：创建测试矩阵
std::vector<std::vector<double>> createTestMatrix(int size)
{
    std::vector<std::vector<double>> A(size,std::vector<double>(size,0.0));
    // 创建一个对角占优矩阵，确保顺序主子式不为0
    for(int i = 0 ; i < size ; i++)
    {
        for(int j = 0; j < size; j++)
        {
            if(i == j)
            {
                A[i][j] = 2.0 + 0.1 * i; //
            }
            else if (std::abs(i - j) == 1)
            {
                A[i][j] = 0.5;// 次对角线元素
            }
            else
            {
                A[i][j] = 0.1/(1 + std::abs(i - j)); //衰减元素
            }
            
        }
    }

    return A;

}

// 示例用法
int main()
{
    std::cout << "=== LU Decomposition with Leading Principal  ===\n" << std::endl;
    // 测试矩阵
    std::vector<std::vector<double>> A = 
    {
        {2.0, 1.0, 5.0},
        {4.0, 4.0,-4.0},
        {1.0, 3.0, 1.0}
    };
    std::vector<double> b = {5.0, 0.0, 6.0};

    LUDecomposition lu; // lu就是这个类的表示
    std::cout << "Original matrix A :" << std::endl;
    for(const auto & row :A)
    {
        for(double val: row)
        {
            std::cout<< std::setw(8) << std::fixed<< std::setprecision(2) << val << " ";
            // setw 是设置字符宽度，这里设置8个，如果不满足8个，就会自动设置
        }
        std::cout << std::endl;
    }
    if(lu.decompose(A))
    {
        std::cout << "LU Decomposition successful" << std::endl;
        std::cout << "All leading principal minors are noe-zero" << std::endl;

        lu.printResult();

        // 验证分解
        if(lu.vertify(A))
        {
            std::cout <<"✔ LU decomposition is correct!" << std::endl;
        }
        else
        {
            std::cout << "✗ LU decompositon has large error!" << std::endl;
        }

        // 计算行列式
        std::cout<<"\nDeterminant of A: " << lu.determinat()<< std::endl;

        // 获取顺序主子式
        std::vector<double> minors = lu.getLeadingPrincipalMinors();
        for(size_t i = 0 ; i < minors.size(); i++)
        {
            std::cout<< "Δ_" << (i +1) << " = " << minors[i];
            if(std::abs(minors[i] < 1e-12))
            {
                std::cout <<"(Warning :close to zero !)";
            }
            std::cout << std::endl;
        }
        // 解方程
        std::vector<double> x = lu.solver(b);
        std::cout <<"\nSolution to Ax = b: " << std::endl;
        std::cout<<"b = [";
        for(size_t i = 0 ; i < b.size() ; i++)
        {
            std::cout<<b[i];
            if(i < b.size() - 1)
            {
                std::cout << ",";
            }
        } 
        std::cout << "]";

        std::cout<< "x = [" ;
        for(size_t i = 0;i < x.size(); i++)
        {
            std::cout <<x[i];
            if( i < x.size() - 1)
            {
                std::cout <<",";
            }
        }
        std::cout<<"]" << std::endl;

        // 验证解的准确性

        std::cout <<"\nVertification (Ax shold equal b):" << std::endl;
        for(size_t i = 0; i < A.size() ; i++)
        {
            double sum = 0.0;
            for(size_t j = 0; j < A[i].size() ; j++)
            {
                sum += A[i][j] * x[j];
            }
            std::cout<<" Row" <<i <<","<<sum<<"(expected " << b[i]
                     <<",error = " << std::abs(sum - b[i]) << ")" << std::endl;
        }

        // 测试更大的矩阵
        std::cout <<"\n====Testing with large matrix===" << std::endl;
        std::vector<std::vector<double>> A_large = createTestMatrix(7);
        LUDecomposition lu2;
        if (lu2.decompose(A_large))
        {
            std::cout<<"7x7 matrix LU decomposition successful" << std::endl;
            std::vector<double> minors_large = lu2.getLeadingPrincipalMinors();
            std::cout<<"Leading pricipal minors for 7x7 matrix:" << std::endl;
            for(size_t i = 0; i < minors_large.size();i++)
            {
                std::cout<<" Δ_" <<(i + 1) << " = " << minors_large[i] << std::endl;
            }
        }
        else
        {
            std::cout<<"LU decomposition failed !"<< std::endl;
            std::cout<<"This usually means one of the leading principal minors is zero "<< std::endl;

        }

       return 0; 


    }
    
}