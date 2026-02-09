#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// define function
double f(double x)
{
    return x * x * x -2 * x * x - 5 * x + 3;
}

// bisection Method to solver function ,the result is root
// root is a number
double bisectionMeshod(double a,double b,int max_iter, double tol,bool showIteration)
{
    if(f(a) * f(b) > 0)
    {
        cout << "Error : function must in range [a , b],and lable is different"<< endl;
        return NAN;
    }

    double c = 0;
    int iteration = 0;

    if(showIteration)
    {
        cout<<"======================================" << endl;
        cout<<setw(8)<<"迭代次数" <<setw(10) <<"a" <<setw(16)<<"b"
        <<setw(13)<<"c" <<setw(18)<<"f(c)" <<setw(28) <<"区间长度"<<endl;

    }
    
    // bisection Method
    for(iteration = 1 ; iteration < max_iter; iteration++)
    {
        c = (a + b) / 2;
        if(showIteration)
        {
            cout <<setw(8) <<fixed<<setprecision(6) <<iteration
                 <<setw(15)<<fixed<<setprecision(6)<<a
                 <<setw(15)<<fixed<<setprecision(6)<<b
                 <<setw(15)<<fixed<<setprecision(6)<<c
                 <<setw(15)<<fixed<<setprecision(6)<<f(c)
                 <<setw(20)<<fixed<<setprecision(6)<<(b - a )/2
                 <<endl;
        }
        if(f(c) == 0 || (b - a) / 2 < tol)
        {
            break;
        }
        if(f(a) * f(c) < 0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
    }
    if(showIteration)
    {
        cout<<"=================================="<<endl;
    }
    cout<<"最终迭代次数是: " << iteration<<endl;
    return c;
}

int main()
{
    cout<<fixed<<setprecision(6);
    double a = 2;
    double b = 6;
    int max_iter = 100;
    double tolerance = 1e-4;
    cout<<"Bisection solver function -- C++ "<<endl;
    double root1 = bisectionMeshod(a,b,max_iter,tolerance,true);
    cout<<"方程的根约为: " <<root1<<endl;
    cout<<"函数值f("<<root1<<") = " <<f(root1) <<endl; 
}