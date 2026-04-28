// File: linearSystems/krylov/PCGSolver.cpp
// Category: Krylov / preconditioned iterative solver
// Method: Preconditioned Conjugate Gradient
// Purpose: Explore PCG solvers, preconditioners, and residual behavior on test systems.

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <chrono>
#include <random>
#include <iomanip>
#include <memory>

using namespace std;

//===============数学工具函数================

// 向量点积
double dot(const vector<double> & a,const vector<double> & b)
{
    double result = 0.0;
    for(size_t i = 0; i < a.size();++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

// 向量范数
double norm(const vector<double> & v)
{
    return sqrt(dot(v,v));
}

// 向量减法
vector<double> subtract(const vector<double>& a, const vector<double> &b)
{
    vector<double> result(a.size());
    for(size_t i = 0; i < a.size();++i)
    {
        result[i]  = a[i] -b[i];
    }
    return result;
}

// 向量加法
vector<double> add(const vector<double> & a,const vector<double> & b)
{
    vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++)
    {
        /* code */
        result[i] = a[i] + b[i];
    }
    return result;
    
}

// 向量标量乘法
vector<double> scalar_multiply(const vector<double> & v, double scalar)
{
    vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
        /* code */
        result[i] = v[i] * scalar;
    }
    return result;
    
}

// 向量复制
vector<double> copy_vector(const vector<double> & v)
{
    return vector<double>(v);
}

// ===========预条件子类接口===========
class Preconditioner
{
private:
    /* data */
public:
    virtual ~Preconditioner() = default;
    virtual vector<double>apply(const vector<double>& r) const = 0;
    virtual string name() const = 0;
};

// =========具体预条件子实现===========
//1.对角预条件子（jacobi预条件子）
class DiagonalPreconditioner :public Preconditioner
{
private:
    vector<double> diag_inv;

public:
    DiagonalPreconditioner(const vector<vector<double>> & A)
    {
        int n = A.size();
        diag_inv.resize(n);
        for(int i = 0; i < n; i++)
        {
            if(abs(A[i][i])> 1e-12)
            {
                diag_inv[i]  = 1/A[i][i];
            }
            else
            {
                diag_inv[i] = 1.0;
            }
        }
    }
    // 这是一个对角预条件子（Jacobi预条件子） 的应用函数
    vector<double> apply(const vector<double>& r) const override
    {
        vector<double> z(r.size());
        for (size_t i = 0; i < r.size(); i++)
        {
            /* code */
            z[i] = diag_inv[i] * r[i];
        }
        return z;
        
    }

    string name() const override
    {
        return "Diagonal(Jacobi) Preconditioner";
    }

};

// 2.单位预条件子（无预条件）
class IdentityPreconditioner:public Preconditioner
{
private:
    /* data */
public:
    vector<double> apply(const vector<double> & r) const override
    {
        return copy_vector(r);
    }

    string name() const override
    {
        return "Identity(No) Preconditioner";
    }
};

// 3.对称超松弛预条件按子(SSOR)
class SSORPreconditioner : public Preconditioner
{
private:
    const vector<vector<double>> & A;
    double omega; // 高斯赛德尔超松弛因子，通常去(0，2)
    
    //前向代入解 (D + ωL)v = r
    vector<double> forward_substitution(const vector<double>& r) const
    {
        int n = r.size();
        vector<double>  z = r;

        for(int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for(int j = 0 ; j < i ;j++)
            {
                sum += A[i][j] * z[j];
            }
            z[i] = (r[i] - omega * sum) / A[i][i];
        }
        return z;
    }

    // 后向代入解 (D + ωU)z = r
    vector<double> backward_substitution(const vector<double> & r) const
    {
        int n = r.size();
        vector<double> z =r;
        for(int i = n  -1 ; i >=0 ; --i)
        {
            double sum = 0.0;
            for(int j = i+1; j < n ; j++)
            {
                sum += A[i][j] * z[j];
            }
            z[i] = (r[i] - omega  * sum ) /A[i][i];
        }
        return z;
    }


public:

    SSORPreconditioner(const vector<vector<double>> & matrix ,double w = 1.0)
    :A(matrix),omega(w){}

    vector<double> apply(const vector<double> & r) const override
    {
        // SSOR预条件子: M = (1/(ω(2-ω))) * (D + ωL)D^{-1}(D + ωU)
        // 简化实现: 先前向代入，然后后向代入
        vector<double> z1 = forward_substitution(r);
        vector<double> z2 = backward_substitution(z1);

        // 缩放因子
        double scale = omega * ( 2.0 -omega);
        if(abs(scale ) > 1e-12)
        {
            return scalar_multiply(z2 ,1.0/scale);
        }
        return z2;
    }

    string name() const override
    {
        return "SSOR Precondition (omega = " + to_string(omega) + ")";
    }

};

// ===================线性系统类=================
class LinearSystem
{
private:
    vector<vector<double>> A ;// 系数矩阵
    vector<double> b;         // 右侧向量

public:
    LinearSystem(const vector<vector<double>> & matrix, const vector<double> & rhs)
    : A(matrix),b(rhs)
    {
        if(matrix.size() != rhs.size())
        {
            throw invalid_argument("Matrix and vector dimensions do not match");
        }
    }

    // 矩阵向量乘法
    vector<double> mat_vec (const vector<double> & x) const
    {
        int n = A.size();
        vector<double> result(n,0.0);
        for(int i = 0 ; i < n ; i++)
        {
            for(int j = 0 ; j < n ; j++)
            {
                result[i] += A[i][j] * x[j];
            }
        }
        return result;
    }

    // 获取矩阵维度
    int size() const
    {
        return A.size();
    }

    // 获取右侧向量
    const vector<double> & get_rhs() const
    {
        return b;
    }

    // 检查矩阵是否对称
    bool is_symetric(double tolerance = 1e-10) const
    {
        int n = A.size();
        for(int i = 0 ; i < n; i++)
        {
            for(int j = 0 ; j< n ;j++)
            {
                if(abs(A[i][j] - A[j][i]) >tolerance)
                {
                    return false;
                }
            }
        }
        return true;
    }

    // 检查矩阵是否正定(通过检查所有顺序主子式)
    bool is_positiver_definite() const
    {
        int n = static_cast<int>(A.size());
        vector<vector<double>>A_copy = A;
        // 检查所有的顺序主子式
        for(int k = 1; k <= n ; k++)
        {
            // 提取k行k列子矩阵
            double det = 1.0;
            for(int i = 0; i < k ; i++)
            {
                // 简单的高斯消元求行列式
                if(abs(A_copy[i][i]) <1e-12)
                {
                    return false;
                }
                det *= A_copy[i][i];

                for(int j = i + 1; j < k ; j++)
                {
                    double factor = A_copy[j][i] /A_copy[i][i];
                    for(int m = i ; m < k ;m++)
                    {
                        A_copy[j][m]  -= factor * A_copy[i][m];
                    }
                }
            }
            if(det <= 0)
            {
                return false;
            }
        }
        return true;
    }

    // 获取矩阵
    const vector<vector<double>> & get_matrix() const
    {
        return A;
    }

};

// =======================PCG求解器类==============
class PCGSolver
{
private:
    int max_iterations;
    double tolerance;
    bool verbose;
    /* data */
public:
    PCGSolver(int max_iter = 1000, double tol = 1e-10 ,bool verb = false)
    : max_iterations(max_iter),tolerance(tol),verbose(verb){}


    // PCG求解主函数
    vector<double> solver(
        const LinearSystem & system,
        const Preconditioner & Preconditioner,
        const vector<double>& initial_guess = vector<double>()
    ) const
    {
        auto start_time = chrono::high_resolution_clock::now();

        int n = system.size();
        const vector<double> & b = system.get_rhs();

        // 初始解
        vector<double> x;
        if(initial_guess.empty() || initial_guess.size() != static_cast<size_t>(n))
        {
            x = vector<double>(n,0.0);
        }
        else
        {
            x = copy_vector(initial_guess);
        }

        // 计算初始残差 r = b -Ax
        vector<double> r = subtract(b,system.mat_vec(x));
        vector<double> z = Preconditioner.apply(r);
        vector<double> p = copy_vector(z);

        double rz_old = dot(r,z);
        double r_norm = norm(r);
        double b_norm = norm(b);
        double relative_residual = r_norm / (b_norm + 1e-16);

        if(verbose)
        {
            cout <<"Using "<<Preconditioner.name()<< endl;
            cout <<"Initial relative residual : " << relative_residual<<endl;
        }

        int iter;
        // PCG主迭代
        for(iter = 0; iter< max_iterations; iter++)
        {
            // 检查收敛
            if(relative_residual < tolerance)
            {
                break;
            }

            // 计算Ap
            vector<double> Ap = system.mat_vec(p);

            // 计算步长 alpha
            double pAp = dot(p,Ap);
            if(abs(pAp) < 1e-16)
            {
                if(verbose)
                {
                    cout<<"Warning:P^1 A p is too small at iterations" << iter <<endl;
                }
                break;
            }

            double alpha = rz_old / pAp;

            // 更新解和残差
            x = add(x,scalar_multiply(p,alpha));
            r = subtract(r,scalar_multiply(Ap,alpha));

            // 预条件残差
            z = Preconditioner.apply(r);

            // 计算新的内积
            double rz_new = dot(r,z);

            // 更新搜索方向
            double beta = rz_new /rz_old;
            p = add(z,scalar_multiply(p,beta));

            // 更新内积
            rz_old = rz_new;

            // 计算相对残差
            r_norm = norm(r);
            relative_residual = r_norm /(b_norm + 1e-16);

            // 输出迭代信息
            if(verbose&&iter % 10 == 0)
            {
                cout<<"Iterations " << setw(4) <<iter
                <<": relativer residual  = " << scientific << setprecision(3)
                <<relative_residual<< endl;
            }


        }
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

        if(verbose)
        {
            cout<<"\nPCG converged in " << iter<<" iterations" << endl;
            cout<<"Final relative residual : " << scientific<<relative_residual<<endl;
            cout<<"Time elapsed : " << duration.count() <<" ms" << endl;

            // 验证解
            vector<double> Ax = system.mat_vec(x);
            vector<double> final_residual = subtract(b,Ax);
            double final_norm = norm(final_residual);
            cout<<"Final absolute residual : " << final_norm << endl;
        }
        return x;

    }

    // 设置参数
    void set_max_iterations(int max_iter)
    {
        max_iterations = max_iter;
    }

    void set_tolerance(double tol)
    {
        tolerance = tol;
    }
    
    void set_verbose(bool verb)
    {
        verbose = verb;
    }
    
};

// ==================矩阵生成工具================

// 生成随机对称正定矩阵
vector<vector<double>> generate_spd_matrix(int n)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0,1.0);

    // 生成随机矩阵B
    vector<vector<double>> B(n,vector<double>(n));
    for(int i = 0 ; i < n ; i++)
    {
        for(int j = 0; j< n ; j++)
        {
            B[i][j] = dist(gen);
        }
    }

    // 计算 B^T B + nI，确保对称正定
    vector<vector<double>> A(n ,vector<double>(n,0.0));

    for(int i = 0 ; i < n ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            // 计算 B^T B
            for(int k = 0; k < n ; k++)
            {
                A[i][j] += B[k][i] * B[k][j];
            }
            // 架对角增强
            if(i == j)
            {
                A[i][j] += n * 2.0;
            }
        }
    }

    return A;
}

// 生成对角占优矩阵
vector<vector<double>> generate_diagonally_dominant_matrix(int n)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double>dist(0.,1.0);

    vector<vector<double>> A(n,vector<double>(n,0.0));
    for(int i = 0; i < n ; i++)
    {
        double row_sum = 0.0;
        for(int j = 0; j < n ; j++)
        {
            if(i != j)
            {
                A[i][j]  = dist(gen) * 0.1;//非对角元较小
                row_sum += abs(A[i][j]);
            }
        }
        A[i][i] = row_sum + 1.0 +dist(gen);
    }

    // 确保对称
    for(int i = 0; i < n ; i++)
    {
        for(int j = 0 ; j < n ; j++)
        {
            A[j][i]  = A[i][j];
        }
    }

    return A;
}

// 生成右侧向量
vector<double> generate_rhs(int n)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0,10.0);

    vector<double> b(n);
    for(int i = 0 ; i < n ; i++)
    {
        b[i] = dist(gen);
    }

    return b;
}

// =========测试和示例=========

void test_basic_pcg()
{
    cout<<"===========基本PCG测试=========="<<endl;

    // 创建一个小型测试问题
    vector<vector<double>> A = 
    {
        {4, 1 , 0},
        {1, 3 ,-1},
        {0, -1, 2}
    };
    vector<double> b = {6 , 4, 1};
    LinearSystem system(A ,b);
    
    // 测试无预条件
    IdentityPreconditioner identity_prec;
    PCGSolver solver1(100,1e-10,true);

    cout<<"\n1.无预条件子： " << endl;
    vector<double> x1 = solver1.solver(system,identity_prec);

    cout<<"SOlution: ";
    for(double xi :x1)
    {
        cout<<xi<<" ";
    }
    cout <<endl;

    // 测试对角预条件
    DiagonalPreconditioner diag_prec(A);
    PCGSolver solver2(100,1e-10,true);

    cout<<"\n2.对角预条件子： "<< endl;
    vector<double>x2 = solver2.solver(system,diag_prec);

    cout<<"Solution: ";
    for(auto xi : x2)
    {
        cout<<xi<<" ";
    }
    cout<<endl;

    // 验证解的正确性
    vector<double>Ax = system.mat_vec(x2);
    cout<<"Ax: ";
    for(double val :Ax)
    {
        cout<<val<<" ";
    }
    cout <<endl;

    cout<<"b: ";
    for(double val : b)
    {
        cout <<val<<" ";
    }
    cout<< endl;

}

void test_performance()
{
    cout<<"\n============PCG性能比较=========="<< endl;

    int n = 200;
    cout<<"Problem size : " << n << "x" << n<< endl;

    // 生成测试问题
    auto A =generate_diagonally_dominant_matrix(n);
    auto b = generate_rhs(n);
    LinearSystem system(A,b);

    // 不同预条件子测试
    vector<unique_ptr<Preconditioner>> preconditions;
    preconditions.push_back(unique_ptr<Preconditioner>(new IdentityPreconditioner()));
    preconditions.push_back(unique_ptr<Preconditioner>(new DiagonalPreconditioner(A)));
    preconditions.push_back(unique_ptr<Preconditioner>(new SSORPreconditioner(A,1.0)));

    PCGSolver solver(500,1e-10,true);
    for(const auto& prec :preconditions)
    {
        cout<<"\n---------------------------" << endl;

        auto start = chrono::high_resolution_clock::now();
        vector<double> x = solver.solver(system,*prec);
        auto end = chrono::high_resolution_clock::now();

        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

        // 计算最终残差
        vector<double> Ax = system.mat_vec(x);
        vector<double> residual = subtract(b,Ax);
        double res_norm = norm(residual);

        cout<< prec-> name()<<": " << endl;
        cout<<"Time: " <<duration.count() <<"ms"<<endl;
        cout<<"Final residual: "<< scientific<<res_norm<<endl;

    }

}

void test_large_problem()
{
    cout<<"\n===========大规模问题测试=========="<< endl;
    int sizes[] = {100,500,1000};
    for(int n : sizes)
    {
        cout<<"\n===Size: " <<n <<"x "<<n<<"==="<< endl;
        auto A =generate_spd_matrix(n);
        auto b = generate_rhs(n);
        LinearSystem system(A,b);

        // 验证矩阵性质
        cout<<"Matrix is symetric : " <<(system.is_symetric() ? "Yes": "No") << endl;
        cout<<"Matrix is positive definite: " <<(system.is_positiver_definite() ?"Yes" :"No")<< endl;

        // 测试对角预条件
        DiagonalPreconditioner diag_prec(A);
        PCGSolver solver(1000,1e-10,false);

        auto start = chrono::high_resolution_clock::now();
        vector<double> x = solver.solver(system,diag_prec);
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

        // 验证解的质量
        vector<double> Ax = system.mat_vec(x);
        vector<double> residual = subtract(b,Ax);
        double res_norm = norm(residual);
        double b_norm = norm(b);
        double relative_residual = res_norm / b_norm;

        cout<<"Solution quelatity: "<<endl;
        cout<<"Time: " <<duration.count()<<"ms"<<endl;
        cout<<"Absolution residual : " <<scientific<<res_norm<<endl;
        cout<<"Relative redidual: "<<scientific<<relative_residual<<endl;

    }
}

//============主函数=============
int main()
{
    cout<<fixed<<setprecision(6);
    try
    {
        //运行各种测试
        test_basic_pcg();
        test_performance();
        test_large_problem();
        cout<<"\n===========所有测试完成==========="<< endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
    
}
