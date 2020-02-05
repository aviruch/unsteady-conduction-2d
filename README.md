# Two-dimensional Unsteady Heat Conduction

![Compiler](https://img.shields.io/badge/GNU-pass%20(v6.3.0+)-brightgreen.svg)
![Compiler](https://img.shields.io/badge/Intel-not%20tested-yellow.svg)
![Compiler](https://img.shields.io/badge/IBM%20XL-not%20tested-yellow.svg)
![License](https://img.shields.io/badge/License-MIT-blue.svg)

This repository gives Fortran 90 codes to solve two-dimensional unsteady heat conduction problem:

- Numerical solutions which are programed in both **explicit** and **implicit** discrete method are included.
- **Analytical** solution to this problem (Laplace equation) is derived, through separate variable method.

## Contents

- [Installation](#installation)
- [Problem Definition](#problem-definition)
- [Dimensionless Laplace Equation](#dimensionless-laplace-equation)
- [Numerical Solution](#numerical-solution)
    + [Explicit Method](#explicit-method)
    + [Implicit Method](#implicit-method)
        - [Jacobi Iteration (Point)](#jacobi-iteration-point)
        - [Gauss-Seidel Iteration (Point)](#gauss-seidel-iteration-point)
        - [Jacobi Iteration (Line)](#jacobi-iteration-line)
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

The governing equation of the two-dimensional unsteady heat conduction problem defined in a rectangle region is

<p align="center"><img alt="$$&#10;\rho c \frac{\partial T}{\partial t} = \lambda\left( \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} \right),&#10;$$" src="./doc/formula/9828042524f60a32db2ae6e910a98c87.svg" align="middle" width="192.63089999999997pt" height="40.118265pt"/></p>

and the boundary conditions are

<p align="center"><img alt="$$&#10;\left\{&#10;\begin{array}{ll}&#10;\lambda\left.\dfrac{\partial T}{\partial x}\right\vert_{x=0} = -q_w, &amp; \left.T\right\vert_{y=0} = T_s,\\[1.5em]&#10;\lambda\left.\dfrac{\partial T}{\partial x}\right\vert_{x=L_x} = q_e, &amp; \left.T\right\vert_{y=L_y} = T_n.&#10;\end{array}&#10;\right.&#10;$$" src="./doc/formula/d0aea2a408cb1b6de1c4ce1e0836e853.svg" align="middle" width="264.29204999999996pt" height="96.329475pt"/></p>

Where, <img alt="$\rho$" src="./doc/formula/6dec54c48a0438a5fcde6053bdb9d712.svg" align="middle" width="8.498985000000003pt" height="14.155350000000013pt"/>, <img alt="$c$" src="./doc/formula/3e18a4a28fdee1744e5e3f79d13b9ff6.svg" align="middle" width="7.113876000000004pt" height="14.155350000000013pt"/> and <img alt="$\lambda$" src="./doc/formula/fd8be73b54f5436a5cd2e73ba9b6bfa9.svg" align="middle" width="9.589140000000002pt" height="22.831379999999992pt"/> are respectively the density, the specific heat capacity and the thermal conductivity. 

## Dimensionless Laplace Equation 

Define <img alt="$L=(L_x+L_y)/2$" src="./doc/formula/02798c679271dd38ed11e84fe9887f9e.svg" align="middle" width="120.97222499999998pt" height="24.65759999999998pt"/>, <img alt="$X=x/L$" src="./doc/formula/713305657e680d89a41300143615ed28.svg" align="middle" width="65.62776pt" height="24.65759999999998pt"/>, <img alt="$Y=y/L$" src="./doc/formula/4fb66c2d606e2d16f422ccdad20efb18.svg" align="middle" width="63.169754999999995pt" height="24.65759999999998pt"/>, <img alt="$\theta = (T-T_s)/(T_n-T_s)$" src="./doc/formula/aea328b6657ae3c3938017026bb90517.svg" align="middle" width="167.771505pt" height="24.65759999999998pt"/>, <img alt="$\tau=\lambda t/(\rho cL^2)$" src="./doc/formula/958e60ad13023f0e91ce2b1584a68174.svg" align="middle" width="101.66870999999999pt" height="26.76201000000001pt"/>, thus the terms of the Laplace equation can be transformed into

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\dfrac{\partial T}{\partial t}&#10;= \dfrac{\partial \left[ (T_n-T_s)\theta + T_s \right]}{\partial \left[\rho cL^2\tau/\lambda \right]}&#10;= \dfrac{\lambda(T_n-T_s)}{\rho c L^2} \dfrac{\partial\theta}{\partial\tau}, \\[1.5em]&#10;\dfrac{\partial T}{\partial x}&#10;= \dfrac{\partial\left[ (T_n-T_s)\theta + T_s \right]}{\partial \left[LX \right]}&#10;= \dfrac{T_n-T_s}{L} \dfrac{\partial \theta}{\partial X}, \\[1.5em]&#10;\dfrac{\partial T}{\partial y}&#10;= \dfrac{\partial\left[ (T_n-T_s)\theta + T_s \right]}{\partial \left[LY \right]}&#10;= \dfrac{T_n-T_s}{L} \dfrac{\partial \theta}{\partial Y}, \\[1.5em]&#10;\dfrac{\partial^2 T}{\partial x^2}&#10;= \dfrac{\partial}{\partial \left[LX \right]}\left(\dfrac{T_n-T_s}{L} \dfrac{\partial\theta}{\partial X}\right)&#10;= \dfrac{T_n-T_s}{L^2} \dfrac{\partial^2 \theta}{\partial X^2}, \\[1.5em]&#10;\dfrac{\partial^2 T}{\partial y^2}&#10;= \dfrac{\partial}{\partial \left[LY \right]}\left(\dfrac{T_n-T_s}{L} \dfrac{\partial\theta}{\partial Y}\right)&#10;= \dfrac{T_n-T_s}{L^2} \dfrac{\partial^2 \theta}{\partial Y^2}.&#10;\end{array}&#10;$$" src="./doc/formula/e404ec5ef5d8db15b2c0f3a0d341e4cd.svg" align="middle" width="340.58145pt" height="257.26964999999996pt"/></p>

Then a dimensionless governing equation can be derived:

<p align="center"><img alt="$$&#10;\frac{\partial\theta}{\partial\tau} = \frac{\partial^2 \theta}{\partial X^2} + \frac{\partial^2 \theta}{\partial Y^2},&#10;$$" src="./doc/formula/86b144d436c0b1a2fd1cbeb62e985139.svg" align="middle" width="137.260035pt" height="35.777445pt"/></p>

and the dimensionless boundary conditions are

<p align="center"><img alt="$$&#10;\left\{&#10;\begin{array}{ll}&#10;\left. \dfrac{\partial \theta}{\partial X}\right\vert_{X=0} = \dfrac{-L q_w}{\lambda (T_n-T_s)} = -\phi_w, &#10;&amp; \left.\theta\right\vert_{Y=0} = 0, \\[1.5em]&#10;\left. \dfrac{\partial \theta}{\partial X}\right\vert_{X=L_x/L} = \dfrac{L q_e}{\lambda (T_n-T_s)} = \phi_e, &#10;&amp; \left.\theta\right\vert_{Y=L_y/L} = 1.&#10;\end{array}&#10;\right.&#10;$$" src="./doc/formula/6dfc7e5ce2592f6568c0e270f4a1ca15.svg" align="middle" width="379.5363pt" height="98.63106pt"/></p>


## Numerical Solution

Physical parameters listed below will be used in  following simulations and analyses.

<p align="center"><img alt="$$&#10;\begin{matrix}&#10;L_x = 3\text{m}, \quad L_y = 4.5\text{m}, \quad \rho = 7820\text{kg/m}^\text{3}, \quad c = 460\text{J/(kg}\cdot\text{K)}, \quad \lambda = 15\text{W/(m}\cdot\text{K)}, \\[0.5em]&#10;T_s = 500\text{K}, \quad T_n = 300\text{K}, \quad q_w = 800\text{W/m}^\text{2}, \quad q_e = 800\text{W/m}^\text{2}.&#10;\end{matrix}&#10;$$" src="./doc/formula/580cf32af458febe628e01bedaa074ee.svg" align="middle" width="596.9271pt" height="49.607414999999996pt"/></p>

### Explicit Method

Integral numerical formulation can be written in explicit format:

<p align="center"><img alt="$$&#10;\left( \theta_{i, j} - \theta_{i, j}^0 \right)\Delta X \Delta Y&#10;= \left( \frac{\theta_{i - 1, j}^0 - \theta_{i, j}^0}{\delta X_e} - \frac{\theta_{i, j}^0 - \theta_{i + 1, j}^0}{\delta X_w} \right) \Delta Y \Delta\tau + \left( \frac{\theta_{i, j - 1}^0 - \theta_{i, j}^0}{\delta Y_n} - \frac{\theta_{i, j}^0 - \theta_{i, j + 1}^0}{\delta Y_s} \right) \Delta X \Delta\tau.&#10;$$" src="./doc/formula/2d8661453f15fa8d39c0218bbea8a002.svg" align="middle" width="725.53635pt" height="49.31553pt"/></p>

Assumed that <img alt="$\delta X_w = \delta X_e = \Delta X$" src="./doc/formula/6ed26e67ceebdffe6e87a74961a9b621.svg" align="middle" width="133.236345pt" height="22.831379999999992pt"/> and <img alt="$\delta Y_s = \delta Y_n = \Delta Y$" src="./doc/formula/061f8e328602728638dfaa50b0fc6920.svg" align="middle" width="121.64756999999999pt" height="22.831379999999992pt"/>, the formulation above can be simplified as

<p align="center"><img alt="$$&#10;\theta_{i, j} = (1 - 2\alpha_X - 2\alpha_Y) \theta_{i, j}^0 + \alpha_X \theta_{i - 1, j}^0 + \alpha_X \theta_{i + 1, j}^0 + \alpha_Y \theta_{i, j - 1}^0 + \alpha_Y \theta_{i, j + 1}^0,&#10;$$" src="./doc/formula/b1f945ab0f9fb6a0fde91c558eac9281.svg" align="middle" width="525.3286499999999pt" height="20.504055pt"/></p>

where <img alt="$\alpha_X = \Delta\tau / (\Delta X)^2$" src="./doc/formula/040ad3b1a03349e62ce33b2b0a7acc3e.svg" align="middle" width="123.83992499999998pt" height="26.76201000000001pt"/> and <img alt="$\alpha_Y = \Delta\tau / (\Delta Y)^2$" src="./doc/formula/744f942adb324c4f8861a86248fb4aa7.svg" align="middle" width="121.01133000000002pt" height="26.76201000000001pt"/>, subscript <img alt="$i$" src="./doc/formula/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663295000000005pt" height="21.683310000000006pt"/> ranges from <img alt="$2$" src="./doc/formula/76c5792347bb90ef71cfbace628572cf.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> to <img alt="$N_X - 1$" src="./doc/formula/532b7605c4b8992fb438154234441fc6.svg" align="middle" width="54.014894999999996pt" height="22.46574pt"/> and subscript <img alt="$j$" src="./doc/formula/36b5afebdba34564d884d347484ac0c7.svg" align="middle" width="7.710483000000004pt" height="21.683310000000006pt"/> ranges from <img alt="$2$" src="./doc/formula/76c5792347bb90ef71cfbace628572cf.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> to <img alt="$N_Y - 1$" src="./doc/formula/aa0f44a25a38cafd470df7c14f76c6a8.svg" align="middle" width="52.89850499999999pt" height="22.46574pt"/>. <img alt="$N_X$" src="./doc/formula/4855e1a0b6043df68024a06e6eb9d355.svg" align="middle" width="24.882495000000002pt" height="22.46574pt"/> and <img alt="$N_Y$" src="./doc/formula/451a9ba1c48b04fdcd537308328e8b66.svg" align="middle" width="23.766105000000007pt" height="22.46574pt"/> are respectively number of mesh nodes. 

Bottom and top sides are Dirichlet's boundaries, where non-dimensional temperature of boundaries <img alt="$\theta_B$" src="./doc/formula/85adc7353e07e5bc03447f24cf3e589b.svg" align="middle" width="18.209565pt" height="22.831379999999992pt"/> is introduced: 

- When <img alt="$i \neq 1$" src="./doc/formula/573091b03bafd592cd7acd78b9fbd635.svg" align="middle" width="35.800050000000006pt" height="22.831379999999992pt"/> or <img alt="$N_X$" src="./doc/formula/4855e1a0b6043df68024a06e6eb9d355.svg" align="middle" width="24.882495000000002pt" height="22.46574pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, <img alt="$\theta_B = \theta_s = 0$" src="./doc/formula/3d8edabcee55caa0a40be072b5f54525.svg" align="middle" width="85.82903999999999pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\frac{\theta_{i, 1} - \theta_B}{\Delta Y / 2} = \frac{\theta_{i, 2} - \theta_B}{3\Delta Y / 2} &#10;\quad \Rightarrow \quad &#10;\theta_{i, 1} = \frac{1}{3} \theta_{i, 2}; &#10;$$" src="./doc/formula/3bcc796d212ac5350e7537dfe75b7838.svg" align="middle" width="297.7623pt" height="37.92162pt"/></p>

- When <img alt="$i \neq 1$" src="./doc/formula/573091b03bafd592cd7acd78b9fbd635.svg" align="middle" width="35.800050000000006pt" height="22.831379999999992pt"/> or <img alt="$N_X$" src="./doc/formula/4855e1a0b6043df68024a06e6eb9d355.svg" align="middle" width="24.882495000000002pt" height="22.46574pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, <img alt="$\theta_B = \theta_n = 1$" src="./doc/formula/4aa03f5fa64a9e288ee61dd524b16e5a.svg" align="middle" width="87.750795pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\frac{\theta_{i, N_Y} - \theta_B}{\Delta Y / 2} = \frac{\theta_{i, N_Y - 1} - \theta_B}{3 \Delta Y / 2} &#10;\quad \Rightarrow \quad&#10;\theta_{i, N_Y} = \frac{1}{3} \theta_{i, N_Y - 1} + \frac{2}{3}.&#10;$$" src="./doc/formula/2d9daca071c442ee33c5508c2b924585.svg" align="middle" width="418.62975pt" height="37.92162pt"/></p>

Left and right sides are Neumann's boundaries, where non-dimensional heat flux <img alt="$\phi_B$" src="./doc/formula/8216311d8816a6198c9522ae9d7abbc0.svg" align="middle" width="20.287245000000002pt" height="22.831379999999992pt"/> is introduced: 

- When <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j \neq 1$" src="./doc/formula/ba27f6e293029d5bbf4188faa7b0598c.svg" align="middle" width="37.847370000000005pt" height="22.831379999999992pt"/> or <img alt="$N_Y$" src="./doc/formula/451a9ba1c48b04fdcd537308328e8b66.svg" align="middle" width="23.766105000000007pt" height="22.46574pt"/>, <img alt="$\phi_B = 2(\theta_B - \theta_{1, j}) / \Delta X = \phi_w$" src="./doc/formula/e41d1f1c2de897c472fd2e43375cce22.svg" align="middle" width="206.61250499999997pt" height="24.65759999999998pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{1, j} = (1 - \alpha_X - 2\alpha_Y) \theta_{1, j}^0 + \frac{\alpha_X}{2} \Delta X \phi_w + \alpha_X \theta_{2, j}^0 + \alpha_Y \theta_{1, j - 1}^0 + \alpha_Y \theta_{1, j + 1}^0; &#10;$$" src="./doc/formula/361d069d8354e412cd701e03939b6315.svg" align="middle" width="522.9378pt" height="29.474114999999998pt"/></p>

- When <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j \neq 1$" src="./doc/formula/ba27f6e293029d5bbf4188faa7b0598c.svg" align="middle" width="37.847370000000005pt" height="22.831379999999992pt"/> or <img alt="$N_Y$" src="./doc/formula/451a9ba1c48b04fdcd537308328e8b66.svg" align="middle" width="23.766105000000007pt" height="22.46574pt"/>, <img alt="$\phi_B = 2(\theta_B - \theta_{N_X, j}) / \Delta X = \phi_e$" src="./doc/formula/3a40dc2713a8f763af89d563dcd422f9.svg" align="middle" width="217.489305pt" height="24.65759999999998pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{N_X, j} = (1 - \alpha_X - 2\alpha_Y) \theta_{N_X, j}^0 + \alpha_X \theta_{N_X - 1, j}^0 + \frac{\alpha_X}{2} \Delta X \phi_e + \alpha_Y \theta_{N_X, j - 1}^0 + \alpha_Y \theta_{N_X, j + 1}^0. &#10;$$" src="./doc/formula/8e38d1e95ad608aae7ef3b508647c2a7.svg" align="middle" width="608.4771pt" height="29.474114999999998pt"/></p>

Four corners are the combination of Dirichlet's and Neumann's boundaries: 

- When <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, <img alt="$\theta_B = \theta_s = 0$" src="./doc/formula/3d8edabcee55caa0a40be072b5f54525.svg" align="middle" width="85.82903999999999pt" height="22.831379999999992pt"/>, <img alt="$\phi_B = \phi_w$" src="./doc/formula/6eb8e82173e2cfec99643fc2f2a8c895.svg" align="middle" width="62.640600000000006pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\theta_{1, 1}&#10;= (1 - 2\alpha_X - 2\alpha_Y) \theta_{1, 1}^0 + \alpha_X \Delta X \phi_w + \alpha_X \theta_{2, 1}^0 + \alpha_Y \theta_{1, 2}^0;&#10;$$" src="./doc/formula/f3daeb6c33b0806fcf29cf6d0e08b651.svg" align="middle" width="428.44725pt" height="20.504055pt"/></p>

- When <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, <img alt="$\theta_B = \theta_n = 1$" src="./doc/formula/4aa03f5fa64a9e288ee61dd524b16e5a.svg" align="middle" width="87.750795pt" height="22.831379999999992pt"/>, <img alt="$\phi_B = \phi_w$" src="./doc/formula/6eb8e82173e2cfec99643fc2f2a8c895.svg" align="middle" width="62.640600000000006pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\theta_{1, N_Y}&#10;= (1 - 2\alpha_X - 2\alpha_Y) \theta_{1, N_Y}^0 + \alpha_X \Delta X \phi_w + \alpha_X \theta_{2, N_Y}^0 + \alpha_Y \theta_{1, N_Y - 1}^0 + 2\alpha_Y;&#10;$$" src="./doc/formula/55addc3670d500ee319fe41f78f99376.svg" align="middle" width="550.4383499999999pt" height="20.504055pt"/></p>

- When <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, <img alt="$\theta_B = \theta_s = 0$" src="./doc/formula/3d8edabcee55caa0a40be072b5f54525.svg" align="middle" width="85.82903999999999pt" height="22.831379999999992pt"/>, <img alt="$\phi_B = \phi_e$" src="./doc/formula/8ca4a9e038afebda2391d1a006780e78.svg" align="middle" width="59.058285000000005pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\theta_{N_X, 1}&#10;= (1 - 2\alpha_X - 2\alpha_Y) \theta_{N_X, 1}^0 + \alpha_X \theta_{N_X - 1, 1}^0 + \alpha_X \Delta X \phi_e + \alpha_Y \theta_{N_X, 2}^0;&#10;$$" src="./doc/formula/7e70bea11bda7f7102e4a07ea591673a.svg" align="middle" width="499.5276pt" height="20.504055pt"/></p>

- When <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, <img alt="$\theta_B = \theta_n = 1$" src="./doc/formula/4aa03f5fa64a9e288ee61dd524b16e5a.svg" align="middle" width="87.750795pt" height="22.831379999999992pt"/>, <img alt="$\phi_B = \phi_e$" src="./doc/formula/8ca4a9e038afebda2391d1a006780e78.svg" align="middle" width="59.058285000000005pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\theta_{N_X, N_Y}&#10;= (1 - 2\alpha_X - 2\alpha_Y) \theta_{N_X, N_Y}^0 + \alpha_X \theta_{N_X - 1, N_Y}^0 + \alpha_X \Delta X \phi_e + \alpha_Y \theta_{N_X, N_Y - 1}^0 + 2\alpha_Y.&#10;$$" src="./doc/formula/a73b2a5a2723051cd0e9b6c792317130.svg" align="middle" width="621.5187pt" height="20.504055pt"/></p>

Note that the explicit method is available only if the diffusion number <img alt="$\alpha_X + \alpha_Y \leq 0.5$" src="./doc/formula/32ecff65e2c1e6e5ca59f966453e56a3.svg" align="middle" width="107.92171500000002pt" height="21.18732pt"/>. 

Relative codes can be referred to [explicit.f90](explicit.f90). Non-dimensional temperature on each nodes as well as the contour is shown here. Compared with the analytical solution, the average value and the standard deviation of the relative error are respectively <img alt="$-1.256 \times 10^{-2}$" src="./doc/formula/4e566f6d17876819b64190b9d7d7cb29.svg" align="middle" width="103.58469pt" height="26.76201000000001pt"/> and <img alt="$5.119 \times 10^{-4}$" src="./doc/formula/688127a718ec21bfae60e8b0a6d4b77d.svg" align="middle" width="90.799335pt" height="26.76201000000001pt"/>.

<p align="center">
  <img height="520px" src="./doc/explicit_ctr.png"></img>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img height="520px" src="./doc/explicit_dat.png"></img>
</p> 

### Implicit Method

Integral numerical formulation can be written in implicit format:

<p align="center"><img alt="$$&#10;\left( \theta_{i, j} - \theta_{i, j}^0 \right)\Delta X \Delta Y&#10;= \left( \frac{\theta_{i - 1, j} - \theta_{i, j}}{\delta X_e} - \frac{\theta_{i, j} - \theta_{i + 1, j}}{\delta X_w} \right) \Delta Y \Delta\tau + \left( \frac{\theta_{i, j - 1} - \theta_{i, j}}{\delta Y_n} - \frac{\theta_{i, j} - \theta_{i, j + 1}}{\delta Y_s} \right) \Delta X \Delta\tau.&#10;$$" src="./doc/formula/cb0b18464cd3edcc071fb8bf900137ba.svg" align="middle" width="721.8848999999999pt" height="39.45249pt"/></p>

Assumed that <img alt="$\delta X_w = \delta X_e = \Delta X$" src="./doc/formula/6ed26e67ceebdffe6e87a74961a9b621.svg" align="middle" width="133.236345pt" height="22.831379999999992pt"/> and <img alt="$\delta Y_s = \delta Y_n = \Delta Y$" src="./doc/formula/061f8e328602728638dfaa50b0fc6920.svg" align="middle" width="121.64756999999999pt" height="22.831379999999992pt"/>, the formulation above can be simplified as

<p align="center"><img alt="$$&#10;(1 + 2\alpha_X + 2\alpha_Y) \theta_{i, j} - \alpha_X \theta_{i - 1, j} - \alpha_X \theta_{i + 1, j} - \alpha_Y \theta_{i, j - 1} - \alpha_Y \theta_{i, j + 1} = \theta_{i, j}^0,&#10;$$" src="./doc/formula/ab49c522a7916fae93f4c9a6f18c5fa2.svg" align="middle" width="525.3286499999999pt" height="20.504055pt"/></p>

where subscript <img alt="$i$" src="./doc/formula/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663295000000005pt" height="21.683310000000006pt"/> ranges from <img alt="$2$" src="./doc/formula/76c5792347bb90ef71cfbace628572cf.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> to <img alt="$N_X - 1$" src="./doc/formula/532b7605c4b8992fb438154234441fc6.svg" align="middle" width="54.014894999999996pt" height="22.46574pt"/> and subscript <img alt="$j$" src="./doc/formula/36b5afebdba34564d884d347484ac0c7.svg" align="middle" width="7.710483000000004pt" height="21.683310000000006pt"/> ranges from <img alt="$2$" src="./doc/formula/76c5792347bb90ef71cfbace628572cf.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> to <img alt="$N_Y - 1$" src="./doc/formula/aa0f44a25a38cafd470df7c14f76c6a8.svg" align="middle" width="52.89850499999999pt" height="22.46574pt"/>. 

Bottom and top sides are Dirichlet's boundaries, where non-dimensional temperature of boundaries <img alt="$\theta_B$" src="./doc/formula/85adc7353e07e5bc03447f24cf3e589b.svg" align="middle" width="18.209565pt" height="22.831379999999992pt"/> is introduced: 

- When <img alt="$i \neq 1$" src="./doc/formula/573091b03bafd592cd7acd78b9fbd635.svg" align="middle" width="35.800050000000006pt" height="22.831379999999992pt"/> or <img alt="$N_X$" src="./doc/formula/4855e1a0b6043df68024a06e6eb9d355.svg" align="middle" width="24.882495000000002pt" height="22.46574pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, <img alt="$\theta_B = \theta_s = 0$" src="./doc/formula/3d8edabcee55caa0a40be072b5f54525.svg" align="middle" width="85.82903999999999pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\frac{\theta_{i, 1} - \theta_B}{\Delta Y / 2} = \frac{\theta_{i, 2} - \theta_B}{3\Delta Y / 2} &#10;\quad \Rightarrow \quad &#10;\theta_{i, 1} - \frac{1}{3} \theta_{i, 2} = 0; &#10;$$" src="./doc/formula/3291f561b2e56c5e289857c558e0c450.svg" align="middle" width="326.07135pt" height="37.92162pt"/></p>

- When <img alt="$i \neq 1$" src="./doc/formula/573091b03bafd592cd7acd78b9fbd635.svg" align="middle" width="35.800050000000006pt" height="22.831379999999992pt"/> or <img alt="$N_X$" src="./doc/formula/4855e1a0b6043df68024a06e6eb9d355.svg" align="middle" width="24.882495000000002pt" height="22.46574pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, <img alt="$\theta_B = \theta_n = 1$" src="./doc/formula/4aa03f5fa64a9e288ee61dd524b16e5a.svg" align="middle" width="87.750795pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\frac{\theta_{i, N_Y} - \theta_B}{\Delta Y / 2} = \frac{\theta_{i, N_Y - 1} - \theta_B}{3 \Delta Y / 2} &#10;\quad \Rightarrow \quad&#10;\theta_{i, N_Y} - \frac{1}{3} \theta_{i, N_Y - 1} = \frac{2}{3}.&#10;$$" src="./doc/formula/f5d1e9527aa66d619069b5ddfce2f3d3.svg" align="middle" width="418.62975pt" height="37.92162pt"/></p>

Left and right sides are Neumann's boundaries, where non-dimensional heat flux <img alt="$\phi_B$" src="./doc/formula/8216311d8816a6198c9522ae9d7abbc0.svg" align="middle" width="20.287245000000002pt" height="22.831379999999992pt"/> is introduced: 

- When <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j \neq 1$" src="./doc/formula/ba27f6e293029d5bbf4188faa7b0598c.svg" align="middle" width="37.847370000000005pt" height="22.831379999999992pt"/> or <img alt="$N_Y$" src="./doc/formula/451a9ba1c48b04fdcd537308328e8b66.svg" align="middle" width="23.766105000000007pt" height="22.46574pt"/>, <img alt="$\phi_B = 2(\theta_B - \theta_{1, j}) / \Delta X = \phi_w$" src="./doc/formula/e41d1f1c2de897c472fd2e43375cce22.svg" align="middle" width="206.61250499999997pt" height="24.65759999999998pt"/>, 

<p align="center"><img alt="$$&#10;(1 + \alpha_X + 2\alpha_Y) \theta_{1, j} - \alpha_X \theta_{2, j} - \alpha_Y \theta_{1, j - 1} - \alpha_Y \theta_{1, j + 1} = \theta_{1, j}^0 + \frac{\alpha_X}{2} \Delta X \phi_w; &#10;$$" src="./doc/formula/2ab918a51cf9a73df082106a2831ad5a.svg" align="middle" width="522.9378pt" height="29.474114999999998pt"/></p>

- When <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j \neq 1$" src="./doc/formula/ba27f6e293029d5bbf4188faa7b0598c.svg" align="middle" width="37.847370000000005pt" height="22.831379999999992pt"/> or <img alt="$N_Y$" src="./doc/formula/451a9ba1c48b04fdcd537308328e8b66.svg" align="middle" width="23.766105000000007pt" height="22.46574pt"/>, <img alt="$\phi_B = 2(\theta_B - \theta_{N_X, j}) / \Delta X = \phi_e$" src="./doc/formula/3a40dc2713a8f763af89d563dcd422f9.svg" align="middle" width="217.489305pt" height="24.65759999999998pt"/>, 

<p align="center"><img alt="$$&#10;(1 + \alpha_X + 2\alpha_Y) \theta_{N_X, j} - \alpha_X \theta_{N_X - 1, j} - \alpha_Y \theta_{N_X, j - 1} - \alpha_Y \theta_{N_X, j + 1} = \theta_{N_X, j}^0 + \frac{\alpha_X}{2} \Delta X \phi_e. &#10;$$" src="./doc/formula/0e27e03e6bffa335435fa43e18b36b0c.svg" align="middle" width="608.4771pt" height="29.474114999999998pt"/></p>

Four corners are the combination of Dirichlet's and Neumann's boundaries: 

- When <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, <img alt="$\theta_B = \theta_s = 0$" src="./doc/formula/3d8edabcee55caa0a40be072b5f54525.svg" align="middle" width="85.82903999999999pt" height="22.831379999999992pt"/>, <img alt="$\phi_B = \phi_w$" src="./doc/formula/6eb8e82173e2cfec99643fc2f2a8c895.svg" align="middle" width="62.640600000000006pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;(1 + 2\alpha_X + 2\alpha_Y) \theta_{1, 1} - \alpha_X \theta_{2, 1} - \alpha_Y \theta_{1, 2}&#10;= \theta_{1, 1}^0 + \alpha_X \Delta X \phi_w;&#10;$$" src="./doc/formula/b5da43e0b1cb60f05bc1634018da39d4.svg" align="middle" width="428.44725pt" height="20.504055pt"/></p>

- When <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, <img alt="$\theta_B = \theta_n = 1$" src="./doc/formula/4aa03f5fa64a9e288ee61dd524b16e5a.svg" align="middle" width="87.750795pt" height="22.831379999999992pt"/>, <img alt="$\phi_B = \phi_w$" src="./doc/formula/6eb8e82173e2cfec99643fc2f2a8c895.svg" align="middle" width="62.640600000000006pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;(1 + 2\alpha_X + 2\alpha_Y) \theta_{1, N_Y} - \alpha_X \theta_{2, N_Y} - \alpha_Y \theta_{1, N_Y - 1} &#10;= \theta_{1, N_Y}^0 + \alpha_X \Delta X \phi_w+ 2\alpha_Y;&#10;$$" src="./doc/formula/41eb053d451ef939f72bcd2b74ad1aa0.svg" align="middle" width="550.4383499999999pt" height="20.504055pt"/></p>

- When <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, <img alt="$\theta_B = \theta_s = 0$" src="./doc/formula/3d8edabcee55caa0a40be072b5f54525.svg" align="middle" width="85.82903999999999pt" height="22.831379999999992pt"/>, <img alt="$\phi_B = \phi_e$" src="./doc/formula/8ca4a9e038afebda2391d1a006780e78.svg" align="middle" width="59.058285000000005pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;(1 + 2\alpha_X + 2\alpha_Y) \theta_{N_X, 1} - \alpha_X \theta_{N_X - 1, 1} - \alpha_Y \theta_{N_X, 2}&#10;= \theta_{N_X, 1}^0 + \alpha_X \Delta X \phi_e;&#10;$$" src="./doc/formula/a8348b93512e385aea5c1184e6e63d58.svg" align="middle" width="499.5276pt" height="20.504055pt"/></p>

- When <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, <img alt="$\theta_B = \theta_n = 1$" src="./doc/formula/4aa03f5fa64a9e288ee61dd524b16e5a.svg" align="middle" width="87.750795pt" height="22.831379999999992pt"/>, <img alt="$\phi_B = \phi_e$" src="./doc/formula/8ca4a9e038afebda2391d1a006780e78.svg" align="middle" width="59.058285000000005pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;(1 + 2\alpha_X + 2\alpha_Y) \theta_{N_X, N_Y} - \alpha_X \theta_{N_X - 1, N_Y} - \alpha_Y \theta_{N_X, N_Y - 1}&#10;= \theta_{N_X, N_Y}^0 + \alpha_X \Delta X \phi_e + 2\alpha_Y.&#10;$$" src="./doc/formula/b1d821235c45730994672501a6bcc6b6.svg" align="middle" width="621.5187pt" height="20.504055pt"/></p> 

Therefore, <img alt="$n$" src="./doc/formula/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.867000000000003pt" height="14.155350000000013pt"/> equations form a system of linear equations, which can also be represented in matrix

<p align="center"><img alt="$$&#10;\mathbf{A} \theta = b,&#10;$$" src="./doc/formula/702d7396b8bf81c5d145c0f0659f547f.svg" align="middle" width="56.0043pt" height="14.611871999999998pt"/></p>

where <img alt="$\mathbf{A}$" src="./doc/formula/96458543dc5abd380904d95cae6aa2bc.svg" align="middle" width="14.292300000000003pt" height="22.557149999999986pt"/> is a <img alt="$n$" src="./doc/formula/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.867000000000003pt" height="14.155350000000013pt"/>-by-<img alt="$n$" src="./doc/formula/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.867000000000003pt" height="14.155350000000013pt"/> square matrix, <img alt="$\theta$" src="./doc/formula/27e556cf3caa0673ac49a8f0de3c73ca.svg" align="middle" width="8.173588500000005pt" height="22.831379999999992pt"/> and <img alt="$b$" src="./doc/formula/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg" align="middle" width="7.054855500000005pt" height="22.831379999999992pt"/> are <img alt="$n$" src="./doc/formula/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.867000000000003pt" height="14.155350000000013pt"/>-by-<img alt="$1$" src="./doc/formula/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> vectors. <img alt="$n = N_X \cdot N_Y$" src="./doc/formula/7b01443604322acb73658eeace1e43cc.svg" align="middle" width="93.12699pt" height="22.46574pt"/>.

The task to be done, is to solve this system of linear equations. Different iterative methods are supported.

#### Jacobi Iteration (Point)

Square matrix <img alt="$\mathbf{A}$" src="./doc/formula/96458543dc5abd380904d95cae6aa2bc.svg" align="middle" width="14.292300000000003pt" height="22.557149999999986pt"/> can be decomposed into a diagonal component <img alt="$\mathbf{D}$" src="./doc/formula/17104becada06c6cda0447c33ec6c846.svg" align="middle" width="14.497725000000003pt" height="22.557149999999986pt"/>, and the remainder <img alt="$\mathbf{R}$" src="./doc/formula/6423e0d54c2545769ad013e5f6a4cf94.svg" align="middle" width="14.178120000000003pt" height="22.557149999999986pt"/>: 

<p align="center"><img alt="$$&#10;\mathbf{A} = \mathbf{D} + \mathbf{R},&#10;$$" src="./doc/formula/c12185a4a7eed7f90aba49fc8665e00c.svg" align="middle" width="89.54285999999999pt" height="14.474889pt"/></p>

where,

<p align="center"><img alt="$$&#10;\mathbf{D} = \left[&#10;\begin{matrix}&#10;a_{11} &amp; 0 &amp; \cdots &amp; 0 \\&#10;0 &amp; a_{22} &amp; \cdots &amp; 0 \\&#10;\vdots &amp; \vdots &amp; \ddots &amp; \vdots \\&#10;0 &amp; 0 &amp; \cdots &amp; a_{nn}&#10;\end{matrix}\right], \quad&#10;\mathbf{R} = \left[&#10;\begin{matrix}&#10;0 &amp; a_{12} &amp; \cdots &amp; a_{1n} \\&#10;a_{21} &amp; 0 &amp; \cdots &amp; a_{2n} \\&#10;\vdots &amp; \vdots &amp; \ddots &amp; \vdots \\&#10;a_{n1} &amp; a_{n2} &amp; \cdots &amp; 0&#10;\end{matrix}\right].&#10;$$" src="./doc/formula/0ff96f71b5e06dfecb087eb3d157022f.svg" align="middle" width="436.16595pt" height="88.76801999999999pt"/></p>

The solution is then obtained iteratively via 

<p align="center"><img alt="$$&#10;\theta^{(k+1)} = \mathbf{D}^{-1} \left( b - \mathbf{R} \theta^{(k)} \right),&#10;$$" src="./doc/formula/62bd4c2be3557dd47eb6c882ec5f2795.svg" align="middle" width="194.7825pt" height="29.589285pt"/></p>

where <img alt="$\theta^{(k)}$" src="./doc/formula/919cf38c70a1f2cf163309261dcb6d3e.svg" align="middle" width="25.713600000000003pt" height="29.19113999999999pt"/> is the <img alt="$k$" src="./doc/formula/63bb9849783d01d91403bc9a5fea12a2.svg" align="middle" width="9.075495000000004pt" height="22.831379999999992pt"/>-th approximation or iteration of vector <img alt="$\theta$" src="./doc/formula/27e556cf3caa0673ac49a8f0de3c73ca.svg" align="middle" width="8.173588500000005pt" height="22.831379999999992pt"/>. Thus the element-based formula is 

<p align="center"><img alt="$$&#10;\theta_p^{(k+1)} = \frac{1}{a_{pp}} \left( b_p - \sum_{q \neq p}a_{pq} \theta_q^{(k)} \right), \quad p, q \in \left\{ 1, 2, \ldots , n \right\}.&#10;$$" src="./doc/formula/fc6137e37e6cbd64290689fa30b1bb93.svg" align="middle" width="398.06744999999995pt" height="59.178735pt"/></p>

#### Gauss-Seidel Iteration (Point)

Soon...

#### Jacobi Iteration (Line)

Soon...

#### Gauss-Seidel Iteration (Line)

Soon...


## Analytical Solution

Considering that the Laplace equation is a kind of linear equation, the original equation and boundary conditions can be isolated to three parts <img alt="$\text{E}_1$" src="./doc/formula/278b90d99558726b26fb43f8def8b042.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/>, <img alt="$\text{E}_2$" src="./doc/formula/975aeec5482df93394e935220136a624.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/> and <img alt="$\text{E}_3$" src="./doc/formula/f4983d931a0acb7db268de7692c60aff.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/>. Each part has only one non-homogeneous item in boundary conditions.

<p align="center"><img alt="$$&#10;\text{E}_1\left\{&#10;\begin{aligned}&#10;    &amp; \dfrac{\partial^2\theta_1}{\partial X^2} + \dfrac{\partial^2\theta_1}{\partial Y^2} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_1}{\partial X}\right\vert_{X=0} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_1}{\partial X}\right\vert_{X=L_x/L} = 0, \\&#10;    &amp; \left. \theta_1\right\vert_{Y=0} = 0, \\&#10;    &amp; \left. \theta_1\right\vert_{Y=L_y/L} = 1.&#10;\end{aligned}&#10;\right. \quad &#10;\text{E}_2\left\{&#10;\begin{aligned}&#10;    &amp; \dfrac{\partial^2\theta_2}{\partial X^2} + \dfrac{\partial^2\theta_2}{\partial Y^2} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_2}{\partial X}\right\vert_{X=0} = -\phi_w, \\&#10;    &amp; \left. \dfrac{\partial \theta_2}{\partial X}\right\vert_{X=L_x/L} = 0, \\&#10;    &amp; \left. \theta_2\right\vert_{Y=0} = 0, \\&#10;    &amp; \left. \theta_2\right\vert_{Y=L_y/L} = 0.&#10;\end{aligned}&#10;\right. \quad &#10;\text{E}_3\left\{&#10;\begin{aligned}&#10;    &amp; \dfrac{\partial^2\theta_3}{\partial X^2} + \dfrac{\partial^2\theta_3}{\partial Y^2} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_3}{\partial X}\right\vert_{X=0} = 0, \\&#10;    &amp; \left. \dfrac{\partial \theta_3}{\partial X}\right\vert_{X=L_x/L} = \phi_e, \\&#10;    &amp; \left. \theta_3\right\vert_{Y=0} = 0, \\&#10;    &amp; \left. \theta_3\right\vert_{Y=L_y/L} = 0.&#10;\end{aligned}&#10;\right.&#10;$$" src="./doc/formula/88f2fba908c9b1e141a21455b18dc50f.svg" align="middle" width="543.45555pt" height="185.56889999999999pt"/></p>

Dimensionless temperature <img alt="$\theta$" src="./doc/formula/27e556cf3caa0673ac49a8f0de3c73ca.svg" align="middle" width="8.173588500000005pt" height="22.831379999999992pt"/> will be the addition of <img alt="$\theta_n$" src="./doc/formula/6198455ff8721b0169e94091580d971b.svg" align="middle" width="15.842970000000003pt" height="22.831379999999992pt"/>, <img alt="$\theta_2$" src="./doc/formula/f1fe0aebb1c952f09cdbfd83af41f50e.svg" align="middle" width="14.269530000000003pt" height="22.831379999999992pt"/> and <img alt="$\theta_3$" src="./doc/formula/ef3e4ae43ab69ed7bc41775203af5d03.svg" align="middle" width="14.269530000000003pt" height="22.831379999999992pt"/>. With **separate variable method**, <img alt="$\theta_i$" src="./doc/formula/f166369f3ef0a7ff052f1e9bbf57d2e2.svg" align="middle" width="12.367905000000004pt" height="22.831379999999992pt"/> (<img alt="$i \in\{1,2,3\}$" src="./doc/formula/347b6e009f0f84af572622332ca74e39.svg" align="middle" width="81.462315pt" height="24.65759999999998pt"/>) can be written as
<p align="center"><img alt="$$&#10;\theta_i(X,Y) = F_i(X)G_i(Y),&#10;$$" src="./doc/formula/dc256a6a3214ac1cde8f0cffac3b5548.svg" align="middle" width="175.07325pt" height="16.438356pt"/></p>

where <img alt="$F_i(X)$" src="./doc/formula/c990059c5755182724a5d30891979abb.svg" align="middle" width="43.737705000000005pt" height="24.65759999999998pt"/> and <img alt="$G_i(X)$" src="./doc/formula/9865f3644044214150e99633b4063575.svg" align="middle" width="46.091595pt" height="24.65759999999998pt"/> are not always 0. Substituting the above expression into dimensionless Laplace's equation, 
<p align="center"><img alt="$$&#10;\frac{F_i^{\prime\prime}(X)}{F_i(X)} = -\frac{G_i^{\prime\prime}(Y)}{G_i(Y)}.&#10;$$" src="./doc/formula/5cfb5830d8d9c2107a583bfecbdd02d4.svg" align="middle" width="141.44509499999998pt" height="38.864264999999996pt"/></p>

The left side of this equation is a function of <img alt="$X$" src="./doc/formula/cbfb1b2a33b28eab8a3e59464768e810.svg" align="middle" width="14.908740000000003pt" height="22.46574pt"/> which has nothing to do with <img alt="$Y$" src="./doc/formula/91aac9730317276af725abd8cef04ca9.svg" align="middle" width="13.196370000000005pt" height="22.46574pt"/>, whereas the right side of this equation is a function of <img alt="$Y$" src="./doc/formula/91aac9730317276af725abd8cef04ca9.svg" align="middle" width="13.196370000000005pt" height="22.46574pt"/> which has nothing to do with <img alt="$X$" src="./doc/formula/cbfb1b2a33b28eab8a3e59464768e810.svg" align="middle" width="14.908740000000003pt" height="22.46574pt"/>. Therefore, the both sides must be equal to one constant <img alt="$-k_i$" src="./doc/formula/deff2863eef31b4c81b8466899cd9535.svg" align="middle" width="25.994265000000002pt" height="22.831379999999992pt"/>.
<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;    F_i^{\prime\prime}(X) + k_iF_i(X) = 0, \\[0.5em]&#10;    G_i^{\prime\prime}(Y) - k_iG_i(Y) = 0.&#10;\end{array}&#10;$$" src="./doc/formula/a6d61fb103680d14d1f6adedb2cf1afd.svg" align="middle" width="161.512395pt" height="44.58201pt"/></p>

The general solutions of <img alt="$F_i(X)$" src="./doc/formula/c990059c5755182724a5d30891979abb.svg" align="middle" width="43.737705000000005pt" height="24.65759999999998pt"/> and <img alt="$G_i(Y)$" src="./doc/formula/d4660a5b2a4edfb974de111adec12717.svg" align="middle" width="44.37939pt" height="24.65759999999998pt"/> are
<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;    F_i(X) = \left\{\begin{aligned}&#10;        &amp;A_i^0X + B_i^0, &amp; k_i=0, \\&#10;        &amp;A_i^+\sin\sqrt{k_i}X + B_i^+\cos\sqrt{k_i}X, &amp; k_i&gt; 0, \\&#10;        &amp;A_i^-\sinh\sqrt{-k_i}X + B_i^-\cosh\sqrt{-k_i}X, &amp; k_i&lt; 0; \\&#10;    \end{aligned}\right. \\[2.5em]&#10;    G_i(Y) = \left\{\begin{aligned}&#10;        &amp;C_i^0Y + D_i^0, &amp; k_i=0, \\&#10;        &amp;C_i^+\sinh\sqrt{k_i}Y + D_i^+\cosh\sqrt{k_i}Y, &amp; k_i&gt; 0, \\&#10;        &amp;C_i^-\sin\sqrt{-k_i}Y + D_i^-\cos\sqrt{-k_i}Y, &amp; k_i&lt; 0.&#10;    \end{aligned}\right. &#10;\end{array}&#10;$$" src="./doc/formula/6fbb78257d19a7967abfec230e9113c9.svg" align="middle" width="404.73839999999996pt" height="162.12982499999998pt"/></p>

When <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, substitute <img alt="$\theta_1 = F_1 G_1$" src="./doc/formula/3397215080e13de3dbc66868da1cfc8a.svg" align="middle" width="74.4315pt" height="22.831379999999992pt"/> into boundary conditions of <img alt="$\text{E}_1$" src="./doc/formula/278b90d99558726b26fb43f8def8b042.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/>: 

- If <img alt="$k_1 = 0$" src="./doc/formula/448308e9f699f5729392556e4bd0f0ba.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>, 

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\left.\dfrac{\partial \theta_1}{\partial X}\right\vert_{X=0} = \left.\dfrac{\partial \theta_1}{\partial X}\right\vert_{X=L_x/L} = A_1^0 G_1 = 0 &#10;\quad \Rightarrow \quad&#10;A_1^0 = 0; \\[1.5em]&#10;\left.\theta_1\right\vert_{Y=0} = D_1^0F_1 = 0&#10;\quad \Rightarrow \quad&#10;D_1^0 = 0; \\[1em]&#10;\left.\theta_1\right\vert_{Y=L_y/L} = \dfrac{C_1^0 L_y}{L} \cdot B_1^0 = 1&#10;\quad \Rightarrow \quad&#10;B_1^0 C_1^0 = \dfrac{L}{L_y}.&#10;\end{array}&#10;$$" src="./doc/formula/0900b4fcffd285f98bc792e2322ff801.svg" align="middle" width="378.97035pt" height="131.0562pt"/></p>

&ensp;&ensp;&ensp;&ensp;Therefore, 

<p align="center"><img alt="$$&#10;\theta_1(X,Y) = B_1^0 C_1^0 Y = \dfrac{L}{L_y} Y.&#10;$$" src="./doc/formula/3465678e047a736569646ee9dfce42c6.svg" align="middle" width="198.43064999999999pt" height="38.332634999999996pt"/></p>

- If <img alt="$k_1 &gt; 0$" src="./doc/formula/b59a2e7b013ab49bd37368a762f50327.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>,

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\left.\dfrac{\partial \theta_1}{\partial X}\right\vert_{X=0} = A_1^+ \sqrt{k_1} G_1 = 0 &#10;\quad \Rightarrow \quad&#10;A_1^+ = 0; \\[1.5em]&#10;\left.\dfrac{\partial \theta_1}{\partial X}\right\vert_{X=L_x/L} = -B_1^+ \sqrt{k_1} \sin\dfrac{\sqrt{k_1} L_x}{L} = 0 &#10;\quad \Rightarrow \quad&#10;B_1^+ = 0.&#10;\end{array}&#10;$$" src="./doc/formula/09dd11370feba0968a955736d21798fc.svg" align="middle" width="401.60835pt" height="98.93845499999999pt"/></p> 

&ensp;&ensp;&ensp;&ensp;Considered that <img alt="$F_1(X) \equiv 0$" src="./doc/formula/ebb529f149ef135059765a311a077842.svg" align="middle" width="75.77625pt" height="24.65759999999998pt"/> contradicts the assumption, the equation of <img alt="$\text{E}_1$" src="./doc/formula/278b90d99558726b26fb43f8def8b042.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/> has no solution when <img alt="$k_1 &gt; 0$" src="./doc/formula/b59a2e7b013ab49bd37368a762f50327.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>.

- If <img alt="$k_1 &lt; 0$" src="./doc/formula/ce08b0ba2b0ada65699b871e1122450d.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>, 

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\left.\dfrac{\partial \theta_1}{\partial X}\right\vert_{X=0} = A_1^- \sqrt{-k_1} G_1 = 0 &#10;\quad \Rightarrow \quad&#10;A_1^- = 0; \\[1.5em]&#10;\left.\dfrac{\partial \theta_1}{\partial X}\right\vert_{X=L_x/L} = B_1^- \sqrt{-k_1} \sinh\dfrac{\sqrt{-k_1} L_x}{L} = 0 &#10;\quad \Rightarrow \quad&#10;B_1^- = 0.&#10;\end{array}&#10;$$" src="./doc/formula/5a8e2baf91515af2c77db2f09fd9b295.svg" align="middle" width="423.8916pt" height="98.93845499999999pt"/></p> 

&ensp;&ensp;&ensp;&ensp;Considered that <img alt="$F_1(X) \equiv 0$" src="./doc/formula/ebb529f149ef135059765a311a077842.svg" align="middle" width="75.77625pt" height="24.65759999999998pt"/> contradicts the assumption, the equation of <img alt="$\text{E}_1$" src="./doc/formula/278b90d99558726b26fb43f8def8b042.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/> has no solution when <img alt="$k_1 &lt; 0$" src="./doc/formula/ce08b0ba2b0ada65699b871e1122450d.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>.

When <img alt="$i = 2$" src="./doc/formula/85fb1bddc68e7fda92939c7099078e6b.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, substitute <img alt="$\theta_2 = F_2 G_2$" src="./doc/formula/5bf11b180b342bdfd2863ae4d6e575c7.svg" align="middle" width="74.4315pt" height="22.831379999999992pt"/> into boundary conditions of <img alt="$\text{E}_2$" src="./doc/formula/975aeec5482df93394e935220136a624.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/>: 

- If <img alt="$k_2 = 0$" src="./doc/formula/876182a1798ad933978e4954ac486cf4.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>, 

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\left.\dfrac{\partial \theta_2}{\partial X}\right\vert_{X=0} = A_2^0 G_2 = -\phi_w &#10;\quad \Rightarrow \quad&#10;A_2^0 = -\phi_w\dfrac{1}{G_2}; \\[1.5em]&#10;\left.\dfrac{\partial \theta_2}{\partial X}\right\vert_{X=L_x/L} = A_2^0 G_2 = 0 &#10;\quad \Rightarrow \quad&#10;A_2^0 = 0.&#10;\end{array}&#10;$$" src="./doc/formula/bce571583a2e40161f4117c1b457ec57.svg" align="middle" width="340.48245pt" height="97.562355pt"/></p>

&ensp;&ensp;&ensp;&ensp;Considered that <img alt="$\phi_w \neq 0$" src="./doc/formula/24f563b0a443a6cd2605ab5832da0342.svg" align="middle" width="50.572665pt" height="22.831379999999992pt"/>, the equation of <img alt="$\text{E}_2$" src="./doc/formula/975aeec5482df93394e935220136a624.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/> has no solution when <img alt="$k_2 = 0$" src="./doc/formula/876182a1798ad933978e4954ac486cf4.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>.

- If <img alt="$k_2 &gt; 0$" src="./doc/formula/e306bf4981c57add7aa42f5ae4370824.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>, 

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\left.\theta_2\right\vert_{Y=0} = D_2^+ F_2 = 0  &#10;\quad \Rightarrow \quad&#10;D_2^+ = 0; \\[1.5em]&#10;\left.\theta_2\right\vert_{X=L_y/L} = C_2^+ \sinh\dfrac{\sqrt{k_2} L_y}{L} = 0 &#10;\quad \Rightarrow \quad&#10;C_2^+ = 0.&#10;\end{array}&#10;$$" src="./doc/formula/890950661a993ee88d10e567f7a0cc88.svg" align="middle" width="354.4794pt" height="81.218445pt"/></p> 

&ensp;&ensp;&ensp;&ensp;Considered that <img alt="$G_2(Y) \equiv 0$" src="./doc/formula/a02c00ce03ec81509800ec9d58c5296f.svg" align="middle" width="76.41776999999999pt" height="24.65759999999998pt"/> contradicts the assumption, the equation of <img alt="$\text{E}_2$" src="./doc/formula/975aeec5482df93394e935220136a624.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/> has no solution when <img alt="$k_2 &gt; 0$" src="./doc/formula/e306bf4981c57add7aa42f5ae4370824.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>.

- If <img alt="$k_2 &lt; 0$" src="./doc/formula/03c066ea87dd03f3b39ca8978ce373d5.svg" align="middle" width="46.06932pt" height="22.831379999999992pt"/>, 

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\left.\theta_2\right\vert_{Y=0} = D_2^- F_2 = 0  &#10;\quad \Rightarrow \quad&#10;D_2^- = 0; \\[1.5em]&#10;\left.\theta_2\right\vert_{X=L_y/L} = C_2^- \sin\dfrac{\sqrt{-k_2} L_y}{L} = 0 &#10;\quad \Rightarrow \quad&#10;C_2^- \neq 0, \quad k_2 = -\left( \dfrac{n \pi L}{L_y} \right)^2.&#10;\end{array}&#10;$$" src="./doc/formula/f566f817c5b26a55f418ec685fa6bfed.svg" align="middle" width="504.89175pt" height="87.89715pt"/></p>

&ensp;&ensp;&ensp;&ensp;Therefore, the eigen equation with undetermined coefficients is 

<p align="center"><img alt="$$&#10;\theta_2(X,Y) = \sum_{n = 1}^\infty \left( A_{2n}^- \sinh \dfrac{n \pi L}{L_y} X + B_{2n}^- \cosh \dfrac{n \pi L}{L_y} X \right) \cdot C_{2n}^- \sin \dfrac{n \pi L}{L_y} Y,&#10;$$" src="./doc/formula/2e0fe678029e52c3ab384e8e50c024bd.svg" align="middle" width="497.80994999999996pt" height="44.698829999999994pt"/></p>

&ensp;&ensp;&ensp;&ensp;where <img alt="$n \in \mathbb{N}^+$" src="./doc/formula/65b9daa7d31517c726c6180d181d979d.svg" align="middle" width="51.921704999999996pt" height="26.177579999999978pt"/>. According to two other boundary conditions, 

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;\displaystyle \sum_{n = 1}^\infty A_{2n}^- C_{2n}^- \dfrac{n \pi L}{L_y} \sin \dfrac{n \pi L}{L_y} Y = -\phi_w ; \\[1.5em]&#10;A_{2n}^- C_{2n}^- \cosh \dfrac{n \pi L_x}{L_y} + B_{2n}^- C_{2n}^- \sinh \dfrac{n \pi L_x}{L_y} = 0.&#10;\end{array}&#10;$$" src="./doc/formula/c7e868ea210ddcc868d583857b71f69a.svg" align="middle" width="322.9776pt" height="94.41036pt"/></p> 

&ensp;&ensp;&ensp;&ensp;Moreover, based on the orthonormality of eigen equations, 

<p align="center"><img alt="$$&#10;\int_0^1 \sin \dfrac{m \pi L}{L_y} Y \cdot \sin \dfrac{n \pi L}{L_y} Y \text{d} Y = \delta_{mn},&#10;$$" src="./doc/formula/6acb0a31ad1e15976d9ef9026a1a9534.svg" align="middle" width="261.87809999999996pt" height="41.705234999999995pt"/></p>

&ensp;&ensp;&ensp;&ensp;where Kronecker <img alt="$\delta_{mn} = 1$" src="./doc/formula/5fd9ebd06338c08649bd35668daec549.svg" align="middle" width="58.05558pt" height="22.831379999999992pt"/> only if <img alt="$m = n$" src="./doc/formula/3e648cddc9684ebe5e7f9309f3a495fb.svg" align="middle" width="46.217655pt" height="14.155350000000013pt"/>, otherwise <img alt="$\delta_{mn} = 0$" src="./doc/formula/ffe3c3a74aaeaf8e07bff0cb12ed60bc.svg" align="middle" width="58.05558pt" height="22.831379999999992pt"/>. Then coefficients can be determined:

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;A_{2n}^- C_{2n}^- = -\dfrac{L_y \phi_w}{n \pi L} \displaystyle \int_0^1 \sin \dfrac{n \pi L}{L_y} Y \text{d} Y = &#10;\left\{ \begin{array}{ll}&#10;    -\dfrac{2 L_y^2 \phi_w}{(n \pi L)^2}, &amp; n = 2l - 1, \\[1em]&#10;    0, &amp; n=2l;&#10;\end{array}\right. \\[2em]&#10;B_{2n}^- C_{2n}^- = -A_{2n}^- C_{2n}^- \coth \dfrac{n \pi L_x}{L_y} = &#10;\left\{ \begin{array}{ll}&#10;\dfrac{2 L_y^2 \phi_w}{(n \pi L)^2} \cdot \coth \dfrac{n \pi L_x}{L_y}, &amp; n = 2l - 1, \\[1em]&#10;0, &amp; n = 2l.&#10;\end{array}\right.&#10;\end{array}&#10;$$" src="./doc/formula/e06c7152f4f0808713106aab38901132.svg" align="middle" width="512.6434499999999pt" height="146.46686999999997pt"/></p>

&ensp;&ensp;&ensp;&ensp;Where, <img alt="$l \in \mathbb{N}^+$" src="./doc/formula/6d5bfc1800766bd3e19d9b6e0495d639.svg" align="middle" width="47.283060000000006pt" height="26.177579999999978pt"/>. 

When <img alt="$i = 3$" src="./doc/formula/b8fe7a8f30cfdea91cd4bb99033d6e21.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, convert the coordinate of <img alt="$\text{E}_3$" src="./doc/formula/f4983d931a0acb7db268de7692c60aff.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/> from (<img alt="$X$" src="./doc/formula/cbfb1b2a33b28eab8a3e59464768e810.svg" align="middle" width="14.908740000000003pt" height="22.46574pt"/>, <img alt="$Y$" src="./doc/formula/91aac9730317276af725abd8cef04ca9.svg" align="middle" width="13.196370000000005pt" height="22.46574pt"/>) to (<img alt="$L_x/L - X$" src="./doc/formula/2a8d70646337836ac9144eaa3f9c405b.svg" align="middle" width="73.86984pt" height="24.65759999999998pt"/>, <img alt="$Y$" src="./doc/formula/91aac9730317276af725abd8cef04ca9.svg" align="middle" width="13.196370000000005pt" height="22.46574pt"/>), then a solution which is similar to <img alt="$\text{E}_2$" src="./doc/formula/975aeec5482df93394e935220136a624.svg" align="middle" width="17.739810000000002pt" height="22.46574pt"/> can be written as 

<p align="center"><img alt="$$&#10;\theta_3(X,Y) = \sum_{n = 1}^\infty \left[ A_{3n}^- \sinh \dfrac{n \pi L}{L_y} \left( \dfrac{L_x}{L} - X \right) + B_{3n}^- \cosh \dfrac{n \pi L}{L_y} \left( \dfrac{L_x}{L} - X \right) \right] \cdot C_{3n}^- \sin \dfrac{n \pi L}{L_y} Y,&#10;$$" src="./doc/formula/41dabbba777fac61ae228832408a52d8.svg" align="middle" width="631.8411pt" height="44.698829999999994pt"/></p>

as well as the coefficients of the equation are

<p align="center"><img alt="$$&#10;\begin{array}{c}&#10;A_{3n}^- C_{3n}^- = -\dfrac{L_y \phi_e}{n \pi L} \displaystyle \int_0^1 \sin \dfrac{n \pi L}{L_y} Y \text{d} Y = &#10;\left\{ \begin{array}{ll}&#10;    -\dfrac{2 L_y^2 \phi_e}{(n \pi L)^2}, &amp; n = 2l - 1, \\[1em]&#10;    0, &amp; n=2l;&#10;\end{array}\right. \\[2em]&#10;B_{3n}^- C_{3n}^- = -A_{3n}^- C_{3n}^- \coth \dfrac{n \pi L_x}{L_y} = &#10;\left\{ \begin{array}{ll}&#10;\dfrac{2 L_y^2 \phi_e}{(n \pi L)^2} \cdot \coth \dfrac{n \pi L_x}{L_y}, &amp; n = 2l - 1, \\[1em]&#10;0, &amp; n = 2l.&#10;\end{array}\right.&#10;\end{array}&#10;$$" src="./doc/formula/06808d19d795eecffab150b67ff8c244.svg" align="middle" width="512.6434499999999pt" height="146.46686999999997pt"/></p>

Where, <img alt="$l \in \mathbb{N}^+$" src="./doc/formula/6d5bfc1800766bd3e19d9b6e0495d639.svg" align="middle" width="47.283060000000006pt" height="26.177579999999978pt"/>. In conclusion, the analytical solution can be derived with **linear superposition principle**: 

<p align="center"><img alt="$$&#10;\theta(X, Y) = \theta_1(X, Y) + \theta_2(X, Y) + \theta_3(X, Y).&#10;$$" src="./doc/formula/ec18ca9570ca10e8d387bc64d5079784.svg" align="middle" width="309.24629999999996pt" height="16.438356pt"/></p>

Relative codes can be referred to [analytic.f90](analytic.f90). Non-dimensional temperature on each nodes as well as the contour is shown here.

<p align="center">
  <img height="520px" src="./doc/analytic_ctr.png"></img>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img height="520px" src="./doc/analytic_dat.png"></img>
</p> 


## Reference

- 吴清松. 计算热物理引论[M]. 合肥: 中国科学技术大学出版社, 2009.
- Saad Y. Iterative methods for sparse linear systems[M]. SIAM, 2003.
- Peaceman D W, Rachford H H. The numerical solution of parabolic and elliptic differential equations[J]. Journal of the Society for industrial and Applied Mathematics, 1955, 3: 28-41.
- 吴崇试. 数学物理方法[M]. 北京: 北京大学出版社, 2003.
- 顾樵. 数学物理方法[M]. 北京: 科学出版社, 2015.


## License

[MIT](LICENSE) © Hawk Shaw, University of Science and Technology of China.