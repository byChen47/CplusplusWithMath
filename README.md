# CplusplusWithMath

`CplusplusWithMath` is a personal C++ practice repository for advanced numerical computation methods.

`CplusplusWithMath` 是一个用于学习和实现高等数值计算方法的 C++ 代码仓库。

This repository is not a single packaged library.  
It is a collection of small numerical-method programs, examples, and solver experiments.

这个仓库不是一个完整封装好的统一库，  
而是一个以“数值方法练习”和“算法实现”为主的代码集合。

## Repository Goal

The main purpose of this repository is:

- use C++ to implement common numerical computation methods
- understand the algorithmic ideas behind those methods
- verify convergence, residuals, and numerical behavior through simple examples
- build a personal code base for mathematical computing practice

这个仓库主要想做的是：

- 用 C++ 实现常见的数值计算方法
- 理解这些方法背后的算法思想
- 用简单例子验证收敛性、残差和数值表现
- 作为个人高等数值计算练习代码库不断积累

## Topics Covered

From the current source files, the repository mainly involves the following categories.

### 1. Nonlinear Equation Solvers

- `bisectionMethod.cpp`
- `bisectionMethod2.cpp`
- `newtonMethod.cpp`
- `GeneralNewtonMethod.cpp`

These files focus on root-finding methods for nonlinear equations, such as:

- bisection method
- Newton method
- generalized Newton-type implementation ideas

### 2. Linear System Direct Methods

- `GaussianElimination.cpp`
- `LUDecomposition.cpp`
- `LU_LUDecomposition.cpp`

These files mainly correspond to direct methods for solving linear systems:

- Gaussian elimination
- pivoting / elimination ideas
- LU decomposition

### 3. Linear System Iterative Methods

- `jacobi.cpp`
- `jacobi_old.cpp`
- `gauss_seidel.cpp`
- `SOR_gauss_seidel.cpp`
- `SteepestDescent.cpp`

These files are related to iterative solvers for linear algebra problems:

- Jacobi iteration
- Gauss-Seidel iteration
- SOR iteration
- steepest descent method

### 4. Krylov-Type / Preconditioned Iterative Solvers

- `PCGSolver.cpp`
- `PBiCG/PBiCG.cpp`
- `PBiCG/PBiCGexample.cpp`
- `PBiCG/PBiCGsolver.cpp`
- `PBiCG/pbigc_solver.cpp`

These parts are closer to more advanced iterative linear algebra:

- PCG
- preconditioning ideas
- PBiCG / BiCG-style solvers
- simple test problems for symmetric and non-symmetric systems

### 5. Ordinary Differential Equation Solvers

- `EulerMethod.cpp`
- `ImplicitEulerSolver.cpp`

These files demonstrate ODE time-integration ideas, including:

- explicit Euler method
- implicit Euler method

### 6. Matrix / Algebra Utility Experiments

- `isOrthogonalMatrix.cpp`
- `a.cpp`

These files are smaller mathematical or matrix-related experiments.

## Repository Structure

```text
.
├── README.md
├── Makefile
├── *.cpp                    standalone numerical-method programs
├── PBiCG/                   BiCG / preconditioned iterative solver experiments
└── .vscode/                 editor configuration
```

Most files in the root directory are standalone programs.  
In practice, many of them can be compiled and tested independently.

仓库根目录下的大部分 `.cpp` 文件都更接近“独立小程序”，  
通常可以单独编译、单独运行、单独学习。

## Build and Run

### Option 1: compile a single file directly

For most standalone examples, you can compile them with `g++`:

```bash
g++ -std=c++11 -O2 -Wall -Wextra -o bisectionMethod bisectionMethod.cpp
./bisectionMethod
```

Another example:

```bash
g++ -std=c++11 -O2 -Wall -Wextra -o GaussianElimination GaussianElimination.cpp
./GaussianElimination
```

### Option 2: use the provided `Makefile`

The root `Makefile` is currently set up for one target at a time.

At the moment it contains entries like:

```makefile
TARGET = GeneralNewtonMethod
SOURCE = GeneralNewtonMethod.cpp
```

So if you want to build another file with `make`, you usually need to edit `TARGET` and `SOURCE` first.

Then run:

```bash
make
make run
make clean
```

### Option 3: build the `PBiCG` examples

Inside `PBiCG/` there is another `Makefile`.

Typical usage:

```bash
cd PBiCG
make
make run
```

## What You Can Learn from This Repository

This repository is useful if you are learning:

- numerical methods
- computational mathematics
- iterative linear solvers
- matrix computation
- root-finding methods
- basic ODE numerical integration
- how to express mathematical algorithms in C++

If you are studying advanced numerical computation, this repository can serve as:

- a homework-style reference
- a personal experiment collection
- a starting point for building cleaner reusable solvers later

## Current Characteristics of the Code

From the current codebase, some practical characteristics are clear:

- many programs include embedded test examples in `main()`
- the repository is oriented more toward learning than toward library packaging
- several files are iterative revisions of similar ideas
- naming style is not yet fully unified
- some files are more polished than others

这也说明这个仓库更像“学习中的算法实现记录”，  
而不是一个已经完全工程化的数值计算框架。

## Suggested Future Improvements

If you want to keep improving this repository, the most valuable next steps would be:

- add a short description comment at the top of each source file
- group files into folders such as `rootFinding/`, `linearSystem/`, `ODE/`, `Krylov/`
- add input/output examples for each method
- unify file naming style
- separate algorithm classes from example `main()` programs
- add a table comparing methods, convergence conditions, and typical use cases

## Recommended Audience

This repository is especially suitable for:

- students learning numerical analysis
- beginners practicing C++ with mathematical algorithms
- people who want simple executable examples instead of abstract formulas only

## Short Chinese Summary

这个仓库主要是在做：

- 非线性方程求根方法
- 线性方程组直接法
- 线性方程组迭代法
- PCG / PBiCG 等进阶迭代求解器
- 常微分方程的欧拉类方法
- 一些矩阵与代数计算小实验

整体定位更偏向：

- 高等数值计算课程练习
- 数学算法的 C++ 表达
- 个人学习记录与后续扩展基础
