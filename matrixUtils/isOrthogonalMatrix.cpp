// File: matrixUtils/isOrthogonalMatrix.cpp
// Category: Matrix / algebra utility
// Method: Orthogonality check
// Purpose: Test whether a matrix satisfies an orthogonality condition.

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
// 判断矩阵是不是正交矩阵

bool isOrthogonalMatrix(const std::vector<std::vector<double>> &matrix )
{
    size_t n = matrix.size();
    if(n == 0)
    {
        return false; // 空矩阵
    }

    for (size_t i = 0; i < n; i++)
    {
        if(matrix[i].size() != n )
        {
            return false;
        }
    }
    double epsilon = 1e-8;

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            double result = 0.0;
            for (size_t k = 0; k < n; i++)
            {
                result += matrix[k][i] * matrix[k][j];
            }

            if (i==j )
            {
                if (std::abs(result - 1.0) > epsilon)
                {
                    return false;
                }
                
            }
            else
            {
                if(std::abs(result) >epsilon)
                {
                    return false;
                }
            }
            
            
        }
        
    }
    
    return true;
};

int main()
{
    std::vector<std::vector<double>> matrix = 
    {
        {0.6 ,0.8 ,0.0},
        {-0.8,0.6 ,0.0},
        {0.0 ,0.0 ,1.0}
    };

    std::cout << "矩阵" << std::endl;
    for (size_t i = 0; i < matrix.size(); i++)
    {
        for (size_t j = 0; j < matrix.size(); j++)
        {
            std::cout<< std::setw(3)<< matrix[i][j] << "\t";
        }
        std::cout << std::endl;
        
        
    }

    //判断是不是正交矩阵
    if(isOrthogonalMatrix(matrix))
    {
        std::cout<< "this is a Orthogonal Matrix" << std::endl;
    }
    else
    {
        std::cout<< "this isn't a Orthogonal Matrix" << std::endl;
    }

    return 0;
    
}
