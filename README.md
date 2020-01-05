# Two-dimensional Unsteady Heat Conduction

![Compiler](https://img.shields.io/badge/GNU-pass%20(v8.1.0+)-brightgreen.svg)
![Compiler](https://img.shields.io/badge/Intel-not%20tested-yellow.svg)
![Compiler](https://img.shields.io/badge/IBM%20XL-not%20tested-yellow.svg)
![License](https://img.shields.io/badge/License-MIT-blue.svg)

This repository gives Fortran 90 codes to solve two-dimensional unsteady heat conduction problem:

- Numeric solutions which are programed in both **explicit** and **implicit** discrete method are included.
- Analytic solution to this problem (Laplace equation) is derived, through separate variable method.

## Contents

- [Problem Definition](#problem-definition)
- [Dimensionless Laplace Equation](#dimensionless-laplace-equation)
- [Numerical Solution](#numerical-solution)
    + [Explicit Method](#explicit-method)
    + [Implicit Method](#implicit-method)
        - [Jacobi Iteration (Point)](#jacobi-iteration-point)
        - [Jacobi Iteration (Line)](#jacobi-iteration-line)
        - [Gauss-Seidel Iteration (Point)](#gauss-seidel-iteration-point)
        - [Gauss-Seidel Iteration (Line)](#gauss-seidel-iteration-line)
- [Analytical Solution](#analytical-solution)
- [Reference](#reference)
- [License](#license)


## Problem Definition

<img width="300px" align="right" src="./doc/problem_def.svg"></img>

The governing equation of the two-dimensional unsteady heat conduction problem defined in a square region is

<p align="center"><img alt="$$&#10;\rho c \frac{\partial T}{\partial t} = \lambda\left( \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} \right),&#10;$$" src="./doc/formula/9828042524f60a32db2ae6e910a98c87.svg" align="middle" width="192.63089999999997pt" height="40.118265pt"/></p>

and the boundary conditions are

<p align="center"><img alt="$$&#10;\left\{&#10;\begin{array}{ll}&#10;\lambda\left.\dfrac{\partial T}{\partial x}\right\vert_{x=0} = -q_w, &amp; T(x,0,t) = T_1,\\[1.5em]&#10;\lambda\left.\dfrac{\partial T}{\partial x}\right\vert_{x=H} = q_e, &amp; T(x,H,t) = T_0.&#10;\end{array}&#10;\right.&#10;$$" src="./doc/formula/51ed6b93b4033818ed145b03f1624afb.svg" align="middle" width="278.8632pt" height="94.68557999999999pt"/></p>

Where, <img alt="$\rho$" src="./doc/formula/6dec54c48a0438a5fcde6053bdb9d712.svg" align="middle" width="8.498985000000003pt" height="14.155350000000013pt"/>, <img alt="$c$" src="./doc/formula/3e18a4a28fdee1744e5e3f79d13b9ff6.svg" align="middle" width="7.113876000000004pt" height="14.155350000000013pt"/> and <img alt="$\lambda$" src="./doc/formula/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg" align="middle" width="9.589140000000002pt" height="22.831379999999992pt"/> are respectively the density, the specific heat capacity and the thermal conductivity. 

## Dimensionless Laplace Equation 

Define <img alt="$(X, Y)=(x, y)/H$" src="./doc/formula/4ea8f2877e54d3167d1bcaf3dfffa2fe.svg" align="middle" width="130.55559000000002pt" height="24.65759999999998pt"/>, <img alt="$\theta = (T-T_0)/(T_1-T_0)$" src="./doc/formula/cda463e98fb3236a421f10f3df44b694.svg" align="middle" width="166.89535499999997pt" height="24.65759999999998pt"/>, <img alt="$\tau=\lambda t/(\rho cH^2)$" src="./doc/formula/0fe3cd381b962f488e01df3d80fcdd50.svg" align="middle" width="105.48152999999999pt" height="26.76201000000001pt"/>, thus the terms of the Laplace equation can be transformed into

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\dfrac{\partial T}{\partial t}&#10;= \dfrac{\partial \left[ (T_1-T_0)\theta + T_0 \right]}{\partial \left[\rho cH^2\tau/\lambda \right]}&#10;= \dfrac{\lambda(T_1-T_0)}{\rho c H^2} \dfrac{\partial\theta}{\partial\tau}, \\[1.5em]&#10;\dfrac{\partial T}{\partial x}&#10;= \dfrac{\partial\left[ (T_1-T_0)\theta + T_0 \right]}{\partial \left[HX \right]}&#10;= \dfrac{T_1-T_0}{H} \dfrac{\partial \theta}{\partial X}, \\[1.5em]&#10;\dfrac{\partial T}{\partial y}&#10;= \dfrac{\partial\left[ (T_1-T_0)\theta + T_0 \right]}{\partial \left[HY \right]}&#10;= \dfrac{T_1-T_0}{H} \dfrac{\partial \theta}{\partial Y}, \\[1.5em]&#10;\dfrac{\partial^2 T}{\partial x^2}&#10;= \dfrac{\partial}{\partial \left[HX \right]}\left(\dfrac{T_1-T_0}{H} \dfrac{\partial\theta}{\partial X}\right)&#10;= \dfrac{T_1-T_0}{H^2} \dfrac{\partial^2 \theta}{\partial X^2}, \\[1.5em]&#10;\dfrac{\partial^2 T}{\partial y^2}&#10;= \dfrac{\partial}{\partial \left[HY \right]}\left(\dfrac{T_1-T_0}{H} \dfrac{\partial\theta}{\partial Y}\right)&#10;= \dfrac{T_1-T_0}{H^2} \dfrac{\partial^2 \theta}{\partial Y^2}.&#10;\end{array}&#10;$$" src="./doc/formula/3f4e0389264538b270864e6047595f33.svg" align="middle" width="341.9427pt" height="257.26964999999996pt"/></p>

Then a dimensionless governing equation can be derived:

<p align="center"><img alt="$$&#10;\frac{\partial\theta}{\partial\tau} = \frac{\partial^2 \theta}{\partial X^2} + \frac{\partial^2 \theta}{\partial Y^2},&#10;$$" src="./doc/formula/86b144d436c0b1a2fd1cbeb62e985139.svg" align="middle" width="137.260035pt" height="35.777445pt"/></p>

and the dimensionless boundary conditions are

<p align="center"><img alt="$$&#10;\left\{&#10;\begin{array}{ll}&#10;\left. \dfrac{\partial \theta}{\partial X}\right\vert_{X=0} = \dfrac{-H q_w}{\lambda (T_1-T_0)}, &#10;&amp; \theta(X,0,\tau) = 1, \\[1.5em]&#10;\left. \dfrac{\partial \theta}{\partial X}\right\vert_{X=1} = \dfrac{H q_e}{\lambda (T_1-T_0)}, &#10;&amp; \theta(X,1,\tau) = 0.&#10;\end{array}&#10;\right.&#10;$$" src="./doc/formula/c7b92fb48b84980ad77be4a5a5ac46f4.svg" align="middle" width="312.7476pt" height="94.68557999999999pt"/></p>


## Numerical Solution

Physical parameters listed below will be used in  following simulations and analyses.

<p align="center"><img alt="$$&#10;\begin{matrix}&#10;H = 3\text{m}, \quad \rho = 7820\text{kg/m}^\text{3}, \quad c = 460\text{J/(kg}\cdot\text{K)}, \quad \lambda = 15\text{W/(m}\cdot\text{K)}, \\[0.5em]&#10;T_0 = 250^\circ\text{C}, \quad T_1 = 400^\circ\text{C}, \quad q_w = q_e = 750\text{W/m}^\text{2}.&#10;\end{matrix}&#10;$$" src="./doc/formula/6d705c71d9db000246d645818e83a386.svg" align="middle" width="493.01009999999997pt" height="49.607414999999996pt"/></p>

### Explicit Method

Soon...

<p align="center">
  <img width="260px" src="./doc/explicit_ctr.png"></img>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img width="269px" src="./doc/explicit_dat.svg"></img>
</p>

### Implicit Method

Soon...

#### Jacobi Iteration (Point)

Soon...

#### Jacobi Iteration (Line)

Soon...

#### Gauss-Seidel Iteration (Point)

Soon...

#### Gauss-Seidel Iteration (Line)

Soon...


## Analytical Solution

Soon...

<p align="center">
  <img width="260px" src="./doc/analytic_ctr.png"></img>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img width="269px" src="./doc/analytic_dat.svg"></img>
</p>

## Reference

- 吴清松. 计算热物理引论[M]. 合肥: 中国科学技术大学出版社, 2009.
- Peaceman D W, Rachford H H. The numerical solution of parabolic and elliptic differential equations[J]. Journal of the Society for industrial and Applied Mathematics, 1955, 3: 28-41.
- 吴崇试. 数学物理方法[M]. 北京: 北京大学出版社, 2003.
- 顾樵. 数学物理方法[M]. 北京: 科学出版社, 2015.


## License

[MIT](LICENSE) © Hawk Shaw, University of Science and Technology of China.