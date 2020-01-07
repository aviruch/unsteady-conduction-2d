# Two-dimensional Unsteady Heat Conduction

![Compiler](https://img.shields.io/badge/GNU-pass%20(v8.1.0+)-brightgreen.svg)
![Compiler](https://img.shields.io/badge/Intel-not%20tested-yellow.svg)
![Compiler](https://img.shields.io/badge/IBM%20XL-not%20tested-yellow.svg)
![License](https://img.shields.io/badge/License-MIT-blue.svg)

This repository gives Fortran 90 codes to solve two-dimensional unsteady heat conduction problem:

- Numeric solutions which are programed in both **explicit** and **implicit** discrete method are included.
- **Analytic** solution to this problem (Laplace equation) is derived, through separate variable method.

## Contents

- [Installation](#installation)
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

## Installation

- Clone this repository using [Git](https://git-scm.com/).

```
$ git clone https://github.com/Arsennnic/unsteady-conduction-2d.git
```

- Compile with [CMake](https://cmake.org/) and [MinGW](http://www.mingw.org/).
```
$ cd unsteady-conduction-2d/
$ cmake -G "MinGW Makefiles" -B build/
$ cd build/
$ mingw32-make
```

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

<p align="center"><img alt="$$&#10;\left\{&#10;\begin{array}{ll}&#10;\left. \dfrac{\partial \theta}{\partial X}\right\vert_{X=0} = \dfrac{-H q_w}{\lambda (T_1-T_0)} = -q_w^\ast, &#10;&amp; \theta(X,0,\tau) = 1, \\[1.5em]&#10;\left. \dfrac{\partial \theta}{\partial X}\right\vert_{X=1} = \dfrac{H q_e}{\lambda (T_1-T_0)} = q_e^\ast, &#10;&amp; \theta(X,1,\tau) = 0.&#10;\end{array}&#10;\right.&#10;$$" src="./doc/formula/134d9e93196d6e61c77c6b514371c603.svg" align="middle" width="365.43045pt" height="94.68557999999999pt"/></p>


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

Considering that the Laplace equation is a kind of linear equation, the original equation and boundary conditions can be isolated to three parts <img alt="$\text{E}_1$" src="./doc/formula/278b90d99558726b26fb43f8def8b042.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/>, <img alt="$\text{E}_2$" src="./doc/formula/975aeec5482df93394e935220136a624.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/> and <img alt="$\text{E}_3$" src="./doc/formula/f4983d931a0acb7db268de7692c60aff.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/>. Each part has only one non-homogeneous item in boundary conditions.

<p align="center"><img alt="$$&#10;\text{E}_1\left\{&#10;\begin{aligned}&#10;    &amp; \dfrac{\partial^2\theta_1}{\partial X^2} + \dfrac{\partial^2\theta_1}{\partial Y^2} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_1}{\partial X}\right\vert_{X=0} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_1}{\partial X}\right\vert_{X=1} = 0, \\&#10;    &amp; \left. \theta_1\right\vert_{Y=0} = 1, \\&#10;    &amp; \left. \theta_1\right\vert_{Y=1} = 0.&#10;\end{aligned}&#10;\right. \quad &#10;\text{E}_2\left\{&#10;\begin{aligned}&#10;    &amp; \dfrac{\partial^2\theta_2}{\partial X^2} + \dfrac{\partial^2\theta_2}{\partial Y^2} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_2}{\partial X}\right\vert_{X=0} = -q_w^\ast, \\&#10;    &amp; \left. \dfrac{\partial \theta_2}{\partial X}\right\vert_{X=1} = 0, \\&#10;    &amp; \left. \theta_2\right\vert_{Y=0} = 0, \\&#10;    &amp; \left. \theta_2\right\vert_{Y=1} = 0.&#10;\end{aligned}&#10;\right. \quad &#10;\text{E}_3\left\{&#10;\begin{aligned}&#10;    &amp; \dfrac{\partial^2\theta_3}{\partial X^2} + \dfrac{\partial^2\theta_3}{\partial Y^2} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_3}{\partial X}\right\vert_{X=0} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_3}{\partial X}\right\vert_{X=1} = q_e^\ast, \\&#10;    &amp; \left. \theta_3\right\vert_{Y=0} = 0, \\&#10;    &amp; \left. \theta_3\right\vert_{Y=1} = 0.&#10;\end{aligned}&#10;\right.&#10;$$" src="./doc/formula/efaf119f534a3ef7930139cf65f507e1.svg" align="middle" width="530.7324pt" height="179.44905pt"/></p>

Dimensionless temperature <img alt="$\theta$" src="./doc/formula/27e556cf3caa0673ac49a8f0de3c73ca.svg" align="middle" width="8.173588500000005pt" height="22.831379999999992pt"/> will be the addition of <img alt="$\theta_1$" src="./doc/formula/edcbf8dd6dd9743cceeee21183bbc3b6.svg" align="middle" width="14.269530000000003pt" height="22.831379999999992pt"/>, <img alt="$\theta_2$" src="./doc/formula/f1fe0aebb1c952f09cdbfd83af41f50e.svg" align="middle" width="14.269530000000003pt" height="22.831379999999992pt"/> and <img alt="$\theta_3$" src="./doc/formula/ef3e4ae43ab69ed7bc41775203af5d03.svg" align="middle" width="14.269530000000003pt" height="22.831379999999992pt"/>. With **separate variable method**, <img alt="$\theta_i$" src="./doc/formula/f166369f3ef0a7ff052f1e9bbf57d2e2.svg" align="middle" width="12.367905000000004pt" height="22.831379999999992pt"/> (<img alt="$i \in\{1,2,3\}$" src="./doc/formula/347b6e009f0f84af572622332ca74e39.svg" align="middle" width="81.462315pt" height="24.65759999999998pt"/>) can be written as
<p align="center"><img alt="$$&#10;\theta_i(X,Y) = F_i(X)G_i(Y),&#10;$$" src="./doc/formula/dc256a6a3214ac1cde8f0cffac3b5548.svg" align="middle" width="175.07325pt" height="16.438356pt"/></p>

where <img alt="$F_i(X)$" src="./doc/formula/c990059c5755182724a5d30891979abb.svg" align="middle" width="43.737705000000005pt" height="24.65759999999998pt"/> and <img alt="$G_i(X)$" src="./doc/formula/9865f3644044214150e99633b4063575.svg" align="middle" width="46.091595pt" height="24.65759999999998pt"/> are not always 0. Substituting the above expression into dimensionless Laplace's equation, 
<p align="center"><img alt="$$&#10;\frac{F_i^{\prime\prime}(X)}{F_i(X)} = -\frac{G_i^{\prime\prime}(Y)}{G_i(Y)}.&#10;$$" src="./doc/formula/5cfb5830d8d9c2107a583bfecbdd02d4.svg" align="middle" width="141.44509499999998pt" height="38.864264999999996pt"/></p>

The left side of this equation is a function of <img alt="$X$" src="./doc/formula/cbfb1b2a33b28eab8a3e59464768e810.svg" align="middle" width="14.908740000000003pt" height="22.46574pt"/> which has nothing to do with <img alt="$Y$" src="./doc/formula/91aac9730317276af725abd8cef04ca9.svg" align="middle" width="13.196370000000005pt" height="22.46574pt"/>, whereas the right side of this equation is a function of <img alt="$Y$" src="./doc/formula/91aac9730317276af725abd8cef04ca9.svg" align="middle" width="13.196370000000005pt" height="22.46574pt"/> which has nothing to do with <img alt="$X$" src="./doc/formula/cbfb1b2a33b28eab8a3e59464768e810.svg" align="middle" width="14.908740000000003pt" height="22.46574pt"/>. Therefore, the both sides must be equal to one constant <img alt="$-k_i$" src="./doc/formula/deff2863eef31b4c81b8466899cd9535.svg" align="middle" width="25.994265000000002pt" height="22.831379999999992pt"/>.
<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;    F_i^{\prime\prime}(X) + k_iF_i(X) = 0, \\[0.5em]&#10;    G_i^{\prime\prime}(Y) - k_iG_i(Y) = 0.&#10;\end{array}&#10;$$" src="./doc/formula/a6d61fb103680d14d1f6adedb2cf1afd.svg" align="middle" width="161.512395pt" height="44.58201pt"/></p>

The general solutions of <img alt="$F_i(X)$" src="./doc/formula/c990059c5755182724a5d30891979abb.svg" align="middle" width="43.737705000000005pt" height="24.65759999999998pt"/> and <img alt="$G_i(Y)$" src="./doc/formula/d4660a5b2a4edfb974de111adec12717.svg" align="middle" width="44.37939pt" height="24.65759999999998pt"/> are
<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;    F_i(X) = \left\{\begin{aligned}&#10;        &amp;A_i^0X + B_i^0, &amp; k_i=0, \\&#10;        &amp;A_i^+\sin\sqrt{k_i}X + B_i^+\cos\sqrt{k_i}X, &amp; k_i&gt; 0, \\&#10;        &amp;A_i^-\sinh\sqrt{-k_i}X + B_i^-\cosh\sqrt{-k_i}X, &amp; k_i&lt; 0; \\&#10;    \end{aligned}\right. \\[2.5em]&#10;    G_i(Y) = \left\{\begin{aligned}&#10;        &amp;C_i^0Y + D_i^0, &amp; k_i=0, \\&#10;        &amp;C_i^+\sinh\sqrt{k_i}Y + D_i^+\cosh\sqrt{k_i}Y, &amp; k_i&gt; 0, \\&#10;        &amp;C_i^-\sin\sqrt{-k_i}Y + D_i^-\cos\sqrt{-k_i}Y, &amp; k_i&lt; 0.&#10;    \end{aligned}\right. &#10;\end{array}&#10;$$" src="./doc/formula/6fbb78257d19a7967abfec230e9113c9.svg" align="middle" width="404.73839999999996pt" height="162.12982499999998pt"/></p>






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