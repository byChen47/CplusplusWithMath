// PBiCG solver .h
#ifndef PBICG_H
#define PBICG_H

#include <vector>
#include <functional>
#include<memory>

class PBiCG
{
public:
    using Vector = std::vector<double>;
    
    // const Vector is input number,Vector is output number
    // Vector =  MatrixVectorFunc * const Vector
    using MatrixVectorFunc = std::function<void(const Vector & ,Vector&)>;


    // solver config
    struct SolverConfig
    {
        double tolerance;
        int max_iterations;
        bool verbose;
        int print_interval;

        // initail config
        SolverConfig():
        tolerance(1e-8),
        max_iterations(1000),
        verbose(true),
        print_interval(100)
        {

        }
    };
    
    struct SolverResult
    {
        bool converged;
        int iterations;
        double residual_norm;
        double initail_residual;
        double final_residual;
        long long computer_time_ms;

        // solver result initial
        SolverResult():
        converged(false),
        iterations(0),
        residual_norm(0.0),
        initail_residual(0.0),
        final_residual(0.0),
        computer_time_ms(0)
        {

        }
    };
    
    PBiCG() = default;

    // use matrix vector to solver function
    // solve functiono return to solverResult
    SolverResult solve
    (
        const MatrixVectorFunc& A,
        const MatrixVectorFunc& A_T,
        const Vector& b,
        Vector& x,
        const MatrixVectorFunc& M,
        const MatrixVectorFunc& M_T,
        const SolverConfig &config
    );

    // symmetry precondition deal with
    SolverResult solve
    (
        const MatrixVectorFunc& A,
        const Vector& b,
        Vector & x,
        const MatrixVectorFunc& M,
        const SolverConfig& config
    );

    // no precondition deal with
    SolverResult solve
    (
        const MatrixVectorFunc& A,
        const Vector& b,
        Vector& x,
        const SolverConfig& config
    );

private:
    double dotProduct(const Vector & a, const Vector & b);
    double norm(const Vector& v);
    void axpy(double alpha,const Vector& x,Vector& y );
    void copyVector(const Vector& src, Vector& dst);


};



#endif