# CplusplusWithMath

`CplusplusWithMath` is a personal C++ repository for advanced numerical computation methods.

`CplusplusWithMath` 是一个以高等数值计算方法为主题的 C++ 学习与实现仓库。

This repository is not a single packaged library.  
It is a collection of small numerical-method programs, algorithm experiments, and solver prototypes.

这个仓库不是一个统一封装好的数值计算库，  
而是一个围绕“数值方法实现”和“数学算法练习”逐步积累起来的代码集合。

## Repository Goal

The main goals of this repository are:

- implement common numerical methods in C++
- understand the algorithmic structure behind those methods
- verify convergence, residuals, and numerical behavior with executable examples
- build a reusable personal code base for computational mathematics practice

这个仓库的核心目标是：

- 用 C++ 实现常见数值计算方法
- 理解这些方法背后的算法结构
- 通过可运行的小程序验证收敛性、残差和数值表现
- 逐步形成一个可复用的个人数值计算练习代码库

## Method Overview

| Category | Representative methods | Source location |
| --- | --- | --- |
| Nonlinear equations | Bisection, Newton, General Newton | `rootFinding/` |
| Linear systems: direct methods | Gaussian elimination, LU decomposition | `linearSystems/direct/` |
| Linear systems: iterative methods | Jacobi, Gauss-Seidel, SOR, Steepest Descent | `linearSystems/iterative/` |
| Linear systems: Krylov methods | PCG, PBiCG variants | `linearSystems/krylov/` |
| ODE solvers | Explicit Euler, Implicit Euler | `ode/` |
| Matrix / algebra utilities | Orthogonality checks, small experiments | `matrixUtils/` |

## Repository Structure

```text
.
├── README.md
├── Makefile
├── rootFinding/
├── linearSystems/
│   ├── direct/
│   ├── iterative/
│   └── krylov/
│       ├── PCGSolver.cpp
│       └── PBiCG/
├── ode/
└── matrixUtils/
```

Most `.cpp` files are standalone learning-oriented programs.  
They usually contain both the method implementation and a small `main()` example.

仓库中的大多数 `.cpp` 文件都还是“独立可运行”的练习程序，  
通常既包含算法实现，也包含一个用于测试的 `main()` 示例。

## Topics Covered

### 1. Nonlinear Equation Solvers

- `rootFinding/bisectionMethod.cpp`
- `rootFinding/bisectionMethod2.cpp`
- `rootFinding/newtonMethod.cpp`
- `rootFinding/GeneralNewtonMethod.cpp`

These programs focus on scalar nonlinear root-finding and Newton-type methods.

### 2. Linear System Direct Methods

- `linearSystems/direct/GaussianElimination.cpp`
- `linearSystems/direct/LUDecomposition.cpp`
- `linearSystems/direct/LU_LUDecomposition.cpp`

These files cover direct solution strategies for linear systems.

### 3. Linear System Iterative Methods

- `linearSystems/iterative/jacobi.cpp`
- `linearSystems/iterative/jacobi_old.cpp`
- `linearSystems/iterative/gauss_seidel.cpp`
- `linearSystems/iterative/SOR_gauss_seidel.cpp`
- `linearSystems/iterative/SteepestDescent.cpp`

These files focus on classical iterative solvers and related experiments.

### 4. Krylov and Preconditioned Solvers

- `linearSystems/krylov/PCGSolver.cpp`
- `linearSystems/krylov/PBiCG/`

This part is closer to more advanced numerical linear algebra, including:

- conjugate-gradient style solvers
- preconditioning ideas
- BiCG / PBiCG experiments
- symmetric and non-symmetric test systems

### 5. Ordinary Differential Equation Solvers

- `ode/EulerMethod.cpp`
- `ode/ImplicitEulerSolver.cpp`

These programs demonstrate time integration for ODEs.

### 6. Matrix / Algebra Utility Experiments

- `matrixUtils/isOrthogonalMatrix.cpp`
- `matrixUtils/a.cpp`

These are small matrix-related or algebra-related experiments.

## Build and Run

### Option 1: compile a single file directly

For any standalone source file, you can compile it directly with `g++`:

```bash
mkdir -p build
g++ -std=c++11 -O2 -Wall -Wextra -o build/bisectionMethod rootFinding/bisectionMethod.cpp
./build/bisectionMethod
```

Another example:

```bash
mkdir -p build
g++ -std=c++11 -O2 -Wall -Wextra -o build/GaussianElimination linearSystems/direct/GaussianElimination.cpp
./build/GaussianElimination
```

### Option 2: use the root `Makefile`

The root `Makefile` is now designed for quick single-file builds through `SOURCE=...`.

Examples:

```bash
make SOURCE=rootFinding/bisectionMethod.cpp
make SOURCE=linearSystems/direct/GaussianElimination.cpp
make run SOURCE=ode/EulerMethod.cpp
make clean
```

Build products are written to `build/`.

### Option 3: use the `PBiCG` subdirectory `Makefile`

Inside `linearSystems/krylov/PBiCG/` there is a dedicated `Makefile`.

Typical usage:

```bash
cd linearSystems/krylov/PBiCG
make
make run
```

## What You Can Learn from This Repository

This repository is useful if you are learning:

- numerical methods
- computational mathematics
- linear system solvers
- direct and iterative matrix methods
- root-finding methods
- basic ODE numerical integration
- how to express mathematical algorithms in C++

If you are studying advanced numerical computation, this repository can serve as:

- a homework-style reference
- a personal experiment collection
- a prototype base for later refactoring into reusable solver classes

## Current Characteristics of the Code

From the current codebase, several practical characteristics are clear:

- many programs contain embedded test cases in `main()`
- the code is more learning-oriented than library-oriented
- some methods have multiple draft or experimental versions
- naming style and architecture are still evolving
- some files are polished demonstrations, while others are exploratory experiments

这也说明这个仓库更像“数值算法学习记录”，  
而不是一个已经完全工程化、接口统一的数学计算框架。

## Suggested Future Improvements

The most valuable next steps for this repository would be:

- unify naming conventions more strictly
- separate reusable solver classes from demo `main()` programs
- add example input/output snapshots for each method
- add convergence-condition notes for each algorithm
- add a comparison table for accuracy, convergence, and typical use cases
- gradually refactor repeated helper functions into shared utilities

## Recommended Audience

This repository is especially suitable for:

- students learning numerical analysis
- beginners using C++ to understand mathematical algorithms
- readers who prefer small executable examples over formula-only explanations

## Short Chinese Summary

这个仓库目前主要覆盖：

- 非线性方程求根
- 线性方程组直接法
- 线性方程组迭代法
- PCG / PBiCG 等 Krylov 类方法
- 常微分方程的欧拉类时间推进方法
- 一些矩阵与代数计算小实验

整体定位偏向：

- 高等数值计算课程练习
- 数学算法的 C++ 表达
- 个人学习记录与后续整理基础
