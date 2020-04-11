# Two-dimensional Unsteady Heat Conduction

![System](https://img.shields.io/badge/System-Windows%2010-brightgreen.svg)
![Compiler](https://img.shields.io/badge/Build-IVF%202019-brightgreen.svg)
![License](https://img.shields.io/badge/License-MIT-blue.svg)

This repository gives Fortran 90 codes to solve two-dimensional unsteady heat conduction problem:

- Numerical solutions which are programed in both **explicit** and **implicit** discrete method are included.
- **Analytical** solution to this problem (Laplace equation) is given to verify the numerical solutions.
- All equations are derived by **[Wolfram Mathematica 11](https://www.wolfram.com/mathematica/)**.

## Contents

- [Problem Definition](#problem-definition)
- [Dimensionless Laplace Equation](#dimensionless-laplace-equation)
- [Numerical Solution](#numerical-solution)
    + [Explicit Method](#explicit-method)
    + [Implicit Method](#implicit-method)
        - [Jacobi Iteration](#jacobi-iteration)
        - [Gauss-Seidel Iteration](#gauss-seidel-iteration)
        - [Jacobi Iteration (Block) + TDMA Algorithm](#jacobi-iteration-block--tdma-algorithm)
        - [Gauss-Seidel Iteration (Block) + TDMA Algorithm](#gauss-seidel-iteration-block--tdma-algorithm)
- [Analytical Solution](#analytical-solution)
- [Reference](#reference)
- [License](#license)

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

<p align="center"><img alt="$$&#10;\left\{&#10;\begin{array}{ll}&#10;\left. \dfrac{\partial \theta}{\partial X}\right\vert_{X=0} = \dfrac{-L q_w}{\lambda (T_n-T_s)} = -\phi_w, &#10;&amp; \left.\theta\right\vert_{Y=0} = 0, \\[1.5em]&#10;\left. \dfrac{\partial \theta}{\partial X}\right\vert_{X=\Gamma_x} = \dfrac{L q_e}{\lambda (T_n-T_s)} = \phi_e, &#10;&amp; \left.\theta\right\vert_{Y=\Gamma_y} = 1,&#10;\end{array}&#10;\right.&#10;$$" src="./doc/formula/42d1ce9a93cd521dbd13313ae3732468.svg" align="middle" width="353.71215pt" height="96.329475pt"/></p>

where <img alt="$\Gamma_x = L_x/L$" src="./doc/formula/93ec08939bc5c8650581fd987813dbb0.svg" align="middle" width="79.33794pt" height="24.65759999999998pt"/> and <img alt="$\Gamma_y = L_y/L$" src="./doc/formula/2d2b1f7f27a47191088893f12c084050.svg" align="middle" width="78.588345pt" height="24.65759999999998pt"/>.


## Numerical Solution

Physical parameters listed below will be used in following simulations and analyses.

<p align="center"><img alt="$$&#10;\begin{matrix}&#10;L_x = 3\text{m}, \quad L_y = 4.5\text{m}, \quad \rho = 7820\text{kg/m}^\text{3}, \quad c = 460\text{J/(kg}\cdot\text{K)}, \quad \lambda = 15\text{W/(m}\cdot\text{K)}, \\[0.5em]&#10;T_s = 500\text{K}, \quad T_n = 300\text{K}, \quad q_w = 800\text{W/m}^\text{2}, \quad q_e = 800\text{W/m}^\text{2}.&#10;\end{matrix}&#10;$$" src="./doc/formula/580cf32af458febe628e01bedaa074ee.svg" align="middle" width="596.9271pt" height="49.607414999999996pt"/></p>

Besides, the number of meshes in <img alt="$x$" src="./doc/formula/332cc365a4987aacce0ead01b8bdcc0b.svg" align="middle" width="9.395100000000005pt" height="14.155350000000013pt"/> and <img alt="$y$" src="./doc/formula/deceeaf6940a8c7a5a02373728002b0f.svg" align="middle" width="8.649300000000004pt" height="14.155350000000013pt"/> direction are respectively chosen as

<p align="center"><img alt="$$&#10;N_x = 10, \quad N_y = 15.&#10;$$" src="./doc/formula/261bde4714d17a33f46f6f4a08752aef.svg" align="middle" width="147.61593pt" height="15.93603pt"/></p>

You can change the values of all these variables in `params.f90`.


### Explicit Method

Integral numerical formulation can be written in explicit format:

<p align="center"><img alt="$$&#10;\left( \theta_{i, j} - \theta_{i, j}^0 \right)\Delta X \Delta Y&#10;= \left( \frac{\theta_{i + 1, j}^0 - \theta_{i, j}^0}{\delta X_e} - \frac{\theta_{i, j}^0 - \theta_{i - 1, j}^0}{\delta X_w} \right) \Delta Y \Delta\tau + \left( \frac{\theta_{i, j + 1}^0 - \theta_{i, j}^0}{\delta Y_n} - \frac{\theta_{i, j}^0 - \theta_{i, j - 1}^0}{\delta Y_s} \right) \Delta X \Delta\tau.&#10;$$" src="./doc/formula/735af3bb40332bae55bd21be0793bd2b.svg" align="middle" width="725.53635pt" height="49.31553pt"/></p>

Assumed that <img alt="$\delta X_w = \delta X_e = \Delta X$" src="./doc/formula/6ed26e67ceebdffe6e87a74961a9b621.svg" align="middle" width="133.236345pt" height="22.831379999999992pt"/> and <img alt="$\delta Y_s = \delta Y_n = \Delta Y$" src="./doc/formula/061f8e328602728638dfaa50b0fc6920.svg" align="middle" width="121.64756999999999pt" height="22.831379999999992pt"/>, the formulation above can be simplified as

<p align="center"><img alt="$$&#10;\theta_{i, j} = (1 - 2\alpha_X - 2\alpha_Y) \theta_{i, j}^0 + \alpha_X \theta_{i - 1, j}^0 + \alpha_X \theta_{i + 1, j}^0 + \alpha_Y \theta_{i, j - 1}^0 + \alpha_Y \theta_{i, j + 1}^0,&#10;$$" src="./doc/formula/b1f945ab0f9fb6a0fde91c558eac9281.svg" align="middle" width="525.3286499999999pt" height="20.504055pt"/></p>

where <img alt="$\alpha_X = \Delta\tau / (\Delta X)^2$" src="./doc/formula/040ad3b1a03349e62ce33b2b0a7acc3e.svg" align="middle" width="123.83992499999998pt" height="26.76201000000001pt"/>, <img alt="$\alpha_Y = \Delta\tau / (\Delta Y)^2$" src="./doc/formula/744f942adb324c4f8861a86248fb4aa7.svg" align="middle" width="121.01133000000002pt" height="26.76201000000001pt"/>, <img alt="$1 &lt; i &lt; N_X$" src="./doc/formula/a15d0c5928c1260bb2908eb8e089d036.svg" align="middle" width="82.60032pt" height="22.46574pt"/>, <img alt="$1 &lt; j &lt; N_Y$" src="./doc/formula/2c1321a8fb87e5cedd1cf8e3a2c1bc54.svg" align="middle" width="83.531085pt" height="22.46574pt"/>. <img alt="$N_X$" src="./doc/formula/4855e1a0b6043df68024a06e6eb9d355.svg" align="middle" width="24.882495000000002pt" height="22.46574pt"/> and <img alt="$N_Y$" src="./doc/formula/451a9ba1c48b04fdcd537308328e8b66.svg" align="middle" width="23.766105000000007pt" height="22.46574pt"/> are respectively number of mesh nodes. 

Boundary conditions:

- **South**: when <img alt="$1 &lt; i &lt; N_X$" src="./doc/formula/a15d0c5928c1260bb2908eb8e089d036.svg" align="middle" width="82.60032pt" height="22.46574pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{i, 1} = (1 - 2 \alpha_X - 3 \alpha_Y) \theta_{i, 1}^0 + \alpha_X \theta_{i - 1, 1}^0 + \alpha_X \theta_{i + 1, 1}^0 + \alpha_Y \theta_{i, 2}^0&#10;$$" src="./doc/formula/3d5ef3944057c6c88df7ebac5a1fcbfb.svg" align="middle" width="423.5253pt" height="20.504055pt"/></p>

- **North**: when <img alt="$1 &lt; i &lt; N_X$" src="./doc/formula/a15d0c5928c1260bb2908eb8e089d036.svg" align="middle" width="82.60032pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{i, N_Y} = (1 - 2 \alpha_X - 3 \alpha_Y) \theta_{i, N_Y}^0 + \alpha_X \theta_{i - 1, N_Y}^0 + \alpha_X \theta_{i + 1, N_Y}^0 + \alpha_Y \theta_{i, N_Y - 1}^0 + 2 \alpha_Y&#10;$$" src="./doc/formula/3e8a36950869d4ae5321772316c56812.svg" align="middle" width="559.2576pt" height="20.504055pt"/></p>

- **West**: when <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$1&lt; j &lt; N_Y$" src="./doc/formula/24175cc2e42c89c528adfd33bc75800f.svg" align="middle" width="83.531085pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{1, j} = (1 - \alpha_X - 2\alpha_Y) \theta_{1, j}^0 + \alpha_X \Delta X \phi_w + \alpha_X \theta_{2, j}^0 + \alpha_Y \theta_{1, j - 1}^0 + \alpha_Y \theta_{1, j + 1}^0&#10;$$" src="./doc/formula/3031a3994f8f0e80790b970e30c9a3de.svg" align="middle" width="513.60375pt" height="20.504055pt"/></p>

- **East**: when <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$1&lt; j &lt; N_Y$" src="./doc/formula/24175cc2e42c89c528adfd33bc75800f.svg" align="middle" width="83.531085pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{N_X, j} = (1 - \alpha_X - 2\alpha_Y) \theta_{N_X, j}^0 + \alpha_X \theta_{N_X - 1, j}^0 + \alpha_X \Delta X \phi_e + \alpha_Y \theta_{N_X, j - 1}^0 + \alpha_Y \theta_{N_X, j + 1}^0&#10;$$" src="./doc/formula/8f0fd0013f1c444b277e12f7b2f17345.svg" align="middle" width="599.14305pt" height="20.504055pt"/></p>

- **Southwest**: when <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{1, 1}&#10;= (1 - \alpha_X - 3 \alpha_Y) \theta_{1, 1}^0 + \alpha_X \Delta X \phi_w + \alpha_X \theta_{2, 1}^0 + \alpha_Y \theta_{1, 2}^0&#10;$$" src="./doc/formula/9db66dcbdebeac8eb809a18448ea3e61.svg" align="middle" width="414.8397pt" height="20.504055pt"/></p>

- **Northwest**: when <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{1, N_Y}&#10;= (1 - \alpha_X - 3 \alpha_Y) \theta_{1, N_Y}^0 + \alpha_X \Delta X \phi_w + \alpha_X \theta_{2, N_Y}^0 + \alpha_Y \theta_{1, N_Y - 1}^0 + 2 \alpha_Y&#10;$$" src="./doc/formula/793dd709c5341d4993e3a65fd6a8b439.svg" align="middle" width="536.8308pt" height="20.504055pt"/></p>

- **Southeast**: when <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{N_X, 1}&#10;= (1 - \alpha_X - 3 \alpha_Y) \theta_{N_X, 1}^0 + \alpha_X \theta_{N_X - 1, 1}^0 + \alpha_X \Delta X \phi_e + \alpha_Y \theta_{N_X, 2}^0&#10;$$" src="./doc/formula/5fbd80eb967a27e174e89f8baccef4fa.svg" align="middle" width="485.92005pt" height="20.504055pt"/></p>

- **Northeast**: when <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;\theta_{N_X, N_Y}&#10;= (1 - \alpha_X - 3 \alpha_Y) \theta_{N_X, N_Y}^0 + \alpha_X \theta_{N_X - 1, N_Y}^0 + \alpha_X \Delta X \phi_e + \alpha_Y \theta_{N_X, N_Y - 1}^0 + 2 \alpha_Y&#10;$$" src="./doc/formula/5113fdccadaaba215b5aaf6a173a7a36.svg" align="middle" width="607.9111499999999pt" height="20.504055pt"/></p>

>  Note: the explicit method is available only if the diffusion number <img alt="$\alpha_X + \alpha_Y \leq 0.5$" src="./doc/formula/32ecff65e2c1e6e5ca59f966453e56a3.svg" align="middle" width="107.92171500000002pt" height="21.18732pt"/>. 

<!-- 

Non-dimensional temperature on each nodes as well as the contour is shown here. Compared with the analytical solution, the average value and the standard deviation of the relative error are respectively <img alt="$-1.256 \times 10^{-2}$" src="./doc/formula/4e566f6d17876819b64190b9d7d7cb29.svg" align="middle" width="103.58469pt" height="26.76201000000001pt"/> and <img alt="$5.119 \times 10^{-4}$" src="./doc/formula/688127a718ec21bfae60e8b0a6d4b77d.svg" align="middle" width="90.799335pt" height="26.76201000000001pt"/>.

<p align="center">
  <img height="520px" src="./doc/explicit_ctr.png"></img>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img height="520px" src="./doc/explicit_dat.png"></img>
</p>  
-->

### Implicit Method

Integral numerical formulation can be written in implicit format:

<p align="center"><img alt="$$&#10;\left( \theta_{i, j} - \theta_{i, j}^0 \right)\Delta X \Delta Y&#10;= \left( \frac{\theta_{i - 1, j} - \theta_{i, j}}{\delta X_e} - \frac{\theta_{i, j} - \theta_{i + 1, j}}{\delta X_w} \right) \Delta Y \Delta\tau + \left( \frac{\theta_{i, j - 1} - \theta_{i, j}}{\delta Y_n} - \frac{\theta_{i, j} - \theta_{i, j + 1}}{\delta Y_s} \right) \Delta X \Delta\tau.&#10;$$" src="./doc/formula/cb0b18464cd3edcc071fb8bf900137ba.svg" align="middle" width="721.8848999999999pt" height="39.45249pt"/></p>

Assumed that <img alt="$\delta X_w = \delta X_e = \Delta X$" src="./doc/formula/6ed26e67ceebdffe6e87a74961a9b621.svg" align="middle" width="133.236345pt" height="22.831379999999992pt"/> and <img alt="$\delta Y_s = \delta Y_n = \Delta Y$" src="./doc/formula/061f8e328602728638dfaa50b0fc6920.svg" align="middle" width="121.64756999999999pt" height="22.831379999999992pt"/>, the formulation above can be simplified as

<p align="center"><img alt="$$&#10;(1 + 2\alpha_X + 2\alpha_Y) \theta_{i, j} - \alpha_X \theta_{i - 1, j} - \alpha_X \theta_{i + 1, j} - \alpha_Y \theta_{i, j - 1} - \alpha_Y \theta_{i, j + 1} = \theta_{i, j}^0,&#10;$$" src="./doc/formula/ab49c522a7916fae93f4c9a6f18c5fa2.svg" align="middle" width="525.3286499999999pt" height="20.504055pt"/></p>

where <img alt="$1 &lt; i &lt; N_X$" src="./doc/formula/a15d0c5928c1260bb2908eb8e089d036.svg" align="middle" width="82.60032pt" height="22.46574pt"/>, <img alt="$1 &lt; j &lt; N_Y$" src="./doc/formula/2c1321a8fb87e5cedd1cf8e3a2c1bc54.svg" align="middle" width="83.531085pt" height="22.46574pt"/>.

Boundary conditions:

- **South**: when <img alt="$1 &lt; i &lt; N_X$" src="./doc/formula/a15d0c5928c1260bb2908eb8e089d036.svg" align="middle" width="82.60032pt" height="22.46574pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, 

<p align="center"><img alt="$$&#10;(1 + 2 \alpha_X + 3 \alpha_Y) \theta_{i, 1} - \alpha_X \theta_{i - 1, 1} - \alpha_X \theta_{i + 1, 1} - \alpha_Y \theta_{i, 2} = \theta_{i, 1}^0&#10;$$" src="./doc/formula/4871d3eb69795dd95df590b82a2f587c.svg" align="middle" width="423.5253pt" height="20.504055pt"/></p>

- **North**: when <img alt="$1 &lt; i &lt; N_X$" src="./doc/formula/a15d0c5928c1260bb2908eb8e089d036.svg" align="middle" width="82.60032pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;(1 + 2 \alpha_X + 3 \alpha_Y) \theta_{i, N_Y} - \alpha_X \theta_{i - 1, N_Y} - \alpha_X \theta_{i + 1, N_Y} - \alpha_Y \theta_{i, N_Y - 1} = \theta_{i, N_Y}^0 + 2 \alpha_Y&#10;$$" src="./doc/formula/157ed91884274fc1eef469adc4da8aa7.svg" align="middle" width="559.2576pt" height="20.504055pt"/></p>

- **West**: when <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$1&lt; j &lt; N_Y$" src="./doc/formula/24175cc2e42c89c528adfd33bc75800f.svg" align="middle" width="83.531085pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;(1 + \alpha_X + 2\alpha_Y) \theta_{1, j} - \alpha_X \theta_{2, j} - \alpha_Y \theta_{1, j - 1} - \alpha_Y \theta_{1, j + 1} = \theta_{1, j}^0 + \alpha_X \Delta X \phi_w&#10;$$" src="./doc/formula/3953eac30b605c29fc0aaaf32810b894.svg" align="middle" width="513.60375pt" height="20.504055pt"/></p>

- **East**: when <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$1&lt; j &lt; N_Y$" src="./doc/formula/24175cc2e42c89c528adfd33bc75800f.svg" align="middle" width="83.531085pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;(1 + \alpha_X + 2\alpha_Y) \theta_{N_X, j} - \alpha_X \theta_{N_X - 1, j} - \alpha_Y \theta_{N_X, j - 1} - \alpha_Y \theta_{N_X, j + 1} = \theta_{N_X, j}^0 + \alpha_X \Delta X \phi_e&#10;$$" src="./doc/formula/c24e6f0911e05ba36cb9f9f7834a3326.svg" align="middle" width="599.14305pt" height="20.504055pt"/></p>

- **Southwest**: when <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, 

<p align="center"><img alt="$$&#10;(1 + \alpha_X + 3 \alpha_Y) \theta_{1, 1} - \alpha_X \theta_{2, 1} - \alpha_Y \theta_{1, 2}&#10;= \theta_{1, 1}^0 + \alpha_X \Delta X \phi_w&#10;$$" src="./doc/formula/5fff7e46b782e389f2f804b81231dc7b.svg" align="middle" width="414.8397pt" height="20.504055pt"/></p>

- **Northwest**: when <img alt="$i = 1$" src="./doc/formula/0ac75c805f5e7bf3181cb114d8ac5ae4.svg" align="middle" width="35.800050000000006pt" height="21.683310000000006pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;(1 + \alpha_X + 3 \alpha_Y) \theta_{1, N_Y} - \alpha_X \theta_{2, N_Y} - \alpha_Y \theta_{1, N_Y - 1} &#10;= \theta_{1, N_Y}^0 + \alpha_X \Delta X \phi_w + 2 \alpha_Y&#10;$$" src="./doc/formula/5d9536d8bbddbeb2eb76e3c9a1f33f89.svg" align="middle" width="536.8308pt" height="20.504055pt"/></p>

- **Southeast**: when <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j = 1$" src="./doc/formula/99c0c04664ba70db8077b2b6415b2d85.svg" align="middle" width="37.847370000000005pt" height="21.683310000000006pt"/>, 

<p align="center"><img alt="$$&#10;(1 + \alpha_X + 3 \alpha_Y) \theta_{N_X, 1} - \alpha_X \theta_{N_X - 1, 1} - \alpha_Y \theta_{N_X, 2}&#10;= \theta_{N_X, 1}^0 + \alpha_X \Delta X \phi_e&#10;$$" src="./doc/formula/f863743427b7201bafe6243e389d6200.svg" align="middle" width="485.92005pt" height="20.504055pt"/></p>

- **Northeast**: when <img alt="$i = N_X$" src="./doc/formula/f485c01952121672df0977bd2fd855f1.svg" align="middle" width="52.46339999999999pt" height="22.46574pt"/>, <img alt="$j = N_Y$" src="./doc/formula/7b78010d1a89d4144b8b19092fa058c0.svg" align="middle" width="53.394165pt" height="22.46574pt"/>, 

<p align="center"><img alt="$$&#10;(1 + \alpha_X + 3 \alpha_Y) \theta_{N_X, N_Y} - \alpha_X \theta_{N_X - 1, N_Y} - \alpha_Y \theta_{N_X, N_Y - 1}&#10;= \theta_{N_X, N_Y}^0 + \alpha_X \Delta X \phi_e + 2 \alpha_Y&#10;$$" src="./doc/formula/8a98d6e68fe86a9b5aeb7beac1f5042c.svg" align="middle" width="607.9111499999999pt" height="20.504055pt"/></p>

Therefore, <img alt="$n$" src="./doc/formula/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.867000000000003pt" height="14.155350000000013pt"/> equations form a system of linear equations, which can also be represented in matrix

<p align="center"><img alt="$$&#10;\mathbf{A} \theta = b,&#10;$$" src="./doc/formula/702d7396b8bf81c5d145c0f0659f547f.svg" align="middle" width="56.0043pt" height="14.611871999999998pt"/></p>

where <img alt="$\mathbf{A}$" src="./doc/formula/96458543dc5abd380904d95cae6aa2bc.svg" align="middle" width="14.292300000000003pt" height="22.557149999999986pt"/> is a <img alt="$n$" src="./doc/formula/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.867000000000003pt" height="14.155350000000013pt"/>-by-<img alt="$n$" src="./doc/formula/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.867000000000003pt" height="14.155350000000013pt"/> square matrix, <img alt="$\theta$" src="./doc/formula/27e556cf3caa0673ac49a8f0de3c73ca.svg" align="middle" width="8.173588500000005pt" height="22.831379999999992pt"/> and <img alt="$b$" src="./doc/formula/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg" align="middle" width="7.054855500000005pt" height="22.831379999999992pt"/> are <img alt="$n$" src="./doc/formula/55a049b8f161ae7cfeb0197d75aff967.svg" align="middle" width="9.867000000000003pt" height="14.155350000000013pt"/>-by-<img alt="$1$" src="./doc/formula/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> vectors. <img alt="$n = N_X \cdot N_Y$" src="./doc/formula/7b01443604322acb73658eeace1e43cc.svg" align="middle" width="93.12699pt" height="22.46574pt"/>.

The task to be done, is to solve this system of linear equations. Different iterative methods are supported.

#### Jacobi Iteration

Square matrix <img alt="$\mathbf{A}$" src="./doc/formula/96458543dc5abd380904d95cae6aa2bc.svg" align="middle" width="14.292300000000003pt" height="22.557149999999986pt"/> can be decomposed into a diagonal component <img alt="$\mathbf{D}$" src="./doc/formula/17104becada06c6cda0447c33ec6c846.svg" align="middle" width="14.497725000000003pt" height="22.557149999999986pt"/>, and the remainder <img alt="$\mathbf{R}$" src="./doc/formula/6423e0d54c2545769ad013e5f6a4cf94.svg" align="middle" width="14.178120000000003pt" height="22.557149999999986pt"/>: 

<p align="center"><img alt="$$&#10;\mathbf{A} = \mathbf{D} + \mathbf{R},&#10;$$" src="./doc/formula/c12185a4a7eed7f90aba49fc8665e00c.svg" align="middle" width="89.54285999999999pt" height="14.474889pt"/></p>

where,

<p align="center"><img alt="$$&#10;\mathbf{D} = \left[&#10;\begin{matrix}&#10;a_{11} &amp; 0 &amp; \cdots &amp; 0 \\&#10;0 &amp; a_{22} &amp; \cdots &amp; 0 \\&#10;\vdots &amp; \vdots &amp; \ddots &amp; \vdots \\&#10;0 &amp; 0 &amp; \cdots &amp; a_{nn}&#10;\end{matrix}\right], \quad&#10;\mathbf{R} = \left[&#10;\begin{matrix}&#10;0 &amp; a_{12} &amp; \cdots &amp; a_{1n} \\&#10;a_{21} &amp; 0 &amp; \cdots &amp; a_{2n} \\&#10;\vdots &amp; \vdots &amp; \ddots &amp; \vdots \\&#10;a_{n1} &amp; a_{n2} &amp; \cdots &amp; 0&#10;\end{matrix}\right].&#10;$$" src="./doc/formula/0ff96f71b5e06dfecb087eb3d157022f.svg" align="middle" width="436.16595pt" height="88.76801999999999pt"/></p>

The solution is then obtained iteratively via 

<p align="center"><img alt="$$&#10;\theta^{(k+1)} = \mathbf{D}^{-1} \left( b - \mathbf{R} \theta^{(k)} \right),&#10;$$" src="./doc/formula/62bd4c2be3557dd47eb6c882ec5f2795.svg" align="middle" width="194.7825pt" height="29.589285pt"/></p>

where <img alt="$\theta^{(k)}$" src="./doc/formula/919cf38c70a1f2cf163309261dcb6d3e.svg" align="middle" width="25.713600000000003pt" height="29.19113999999999pt"/> is the <img alt="$k$" src="./doc/formula/63bb9849783d01d91403bc9a5fea12a2.svg" align="middle" width="9.075495000000004pt" height="22.831379999999992pt"/>-th approximation or iteration of vector <img alt="$\theta$" src="./doc/formula/27e556cf3caa0673ac49a8f0de3c73ca.svg" align="middle" width="8.173588500000005pt" height="22.831379999999992pt"/>. Thus the element-based formula is 

<p align="center"><img alt="$$&#10;\theta_p^{(k+1)} = \frac{1}{a_{pp}} \left( b_p - \sum_{\substack{q = 1 \\ q \neq p}}^{n} a_{pq} \theta_q^{(k)} \right), \quad p \in \left\{ 1, 2, \ldots , n \right\}.&#10;$$" src="./doc/formula/146bcecd57131c95a147f8472a179c45.svg" align="middle" width="382.83299999999997pt" height="78.904815pt"/></p>

#### Gauss-Seidel Iteration

Square matrix <img alt="$\mathbf{A}$" src="./doc/formula/96458543dc5abd380904d95cae6aa2bc.svg" align="middle" width="14.292300000000003pt" height="22.557149999999986pt"/> can be decomposed into its lower triangular component <img alt="$\mathbf{L_\ast}$" src="./doc/formula/0bc7c34ea269f7f68091dcab926d65da.svg" align="middle" width="18.105120000000003pt" height="22.557149999999986pt"/>, and its strictly upper triangular component <img alt="$\mathbf{U}$" src="./doc/formula/35531be55273dc37ee90083451d089ff.svg" align="middle" width="14.543430000000004pt" height="22.557149999999986pt"/>:

<p align="center"><img alt="$$&#10;\mathbf{A} = \mathbf{L_\ast} + \mathbf{U},&#10;$$" src="./doc/formula/46269c08f5a760b6c772286d9a6f62e3.svg" align="middle" width="94.33743pt" height="14.474889pt"/></p>

where,

<p align="center"><img alt="$$&#10;\mathbf{L_\ast} = \left[&#10;\begin{matrix}&#10;a_{11} &amp; 0 &amp; \cdots &amp; 0 \\&#10;a_{21} &amp; a_{22} &amp; \cdots &amp; 0 \\&#10;\vdots &amp; \vdots &amp; \ddots &amp; \vdots \\&#10;a_{n1} &amp; a_{n2} &amp; \cdots &amp; a_{nn}&#10;\end{matrix}\right], \quad&#10;\mathbf{U} = \left[&#10;\begin{matrix}&#10;0 &amp; a_{12} &amp; \cdots &amp; a_{1n} \\&#10;0 &amp; 0 &amp; \cdots &amp; a_{2n} \\&#10;\vdots &amp; \vdots &amp; \ddots &amp; \vdots \\&#10;0 &amp; 0 &amp; \cdots &amp; 0&#10;\end{matrix}\right].&#10;$$" src="./doc/formula/239aba86f0d3c9d83a6ebf6d7fc2e7a6.svg" align="middle" width="426.56295pt" height="88.76801999999999pt"/></p>

The solution is then obtained iteratively via 

<p align="center"><img alt="$$&#10;\theta^{(k+1)} = \mathbf{L_\ast}^{-1} \left( b - \mathbf{U} \theta^{(k)} \right),&#10;$$" src="./doc/formula/4a87469744018e2d03d298c52a3442d5.svg" align="middle" width="199.57739999999998pt" height="29.589285pt"/></p>

Thus the element-based formula is 

<p align="center"><img alt="$$&#10;\theta_p^{(k+1)} = \frac{1}{a_{pp}} \left( b_p - \sum_{q = 1}^{p - 1} a_{pq} \theta_q^{(k + 1)} - \sum_{q = p + 1}^{n} a_{pq} \theta_q^{(k)} \right), \quad p \in \left\{ 1, 2, \ldots , n \right\}.&#10;$$" src="./doc/formula/9f43245120a3674f3c1ccf2eab3437c7.svg" align="middle" width="508.77750000000003pt" height="50.20125pt"/></p>

#### Jacobi Iteration (Block) + TDMA Algorithm

The methods of the previous sections are also referred to as *point* or *line* iterative methods, since they act on single entries of matrix <img alt="$\mathbf{A}$" src="./doc/formula/96458543dc5abd380904d95cae6aa2bc.svg" align="middle" width="14.292300000000003pt" height="22.557149999999986pt"/>. It is possible to devise block versions of the algorithms, provided that square matrix <img alt="$\mathbf{A}$" src="./doc/formula/96458543dc5abd380904d95cae6aa2bc.svg" align="middle" width="14.292300000000003pt" height="22.557149999999986pt"/> is divided to <img alt="$m$" src="./doc/formula/0e51a2dede42189d77627c4d742822c3.svg" align="middle" width="14.433210000000003pt" height="14.155350000000013pt"/>-by-<img alt="$m$" src="./doc/formula/0e51a2dede42189d77627c4d742822c3.svg" align="middle" width="14.433210000000003pt" height="14.155350000000013pt"/> diagonal blocks of matrix, while vector <img alt="$b$" src="./doc/formula/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg" align="middle" width="7.054855500000005pt" height="22.831379999999992pt"/> is divided to <img alt="$m$" src="./doc/formula/0e51a2dede42189d77627c4d742822c3.svg" align="middle" width="14.433210000000003pt" height="14.155350000000013pt"/>-by-<img alt="$1$" src="./doc/formula/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> blocks of vector:

<p align="center"><img alt="$$&#10;\mathbf{A} = \left[&#10;\begin{matrix}&#10;\mathbf{A}_{11} &amp; \mathbf{A}_{12} &amp; \cdots &amp; \mathbf{A}_{1m} \\&#10;\mathbf{A}_{21} &amp; \mathbf{A}_{22} &amp; \cdots &amp; \mathbf{A}_{2m} \\&#10;\vdots &amp; \vdots &amp; \ddots &amp; \vdots \\&#10;\mathbf{A}_{m1} &amp; \mathbf{A}_{m2} &amp; \cdots &amp; \mathbf{A}_{mm} \\&#10;\end{matrix}\right], \quad&#10;b = \left[&#10;\begin{matrix}&#10;b_{1} \\&#10;b_{2} \\&#10;\vdots \\&#10;b_{m} \\&#10;\end{matrix}\right].&#10;$$" src="./doc/formula/764003dff78ada7d75822fd1622b339a.svg" align="middle" width="338.6889pt" height="88.76801999999999pt"/></p>

The solution is then obtained iteratively via 

<p align="center"><img alt="$$&#10;\mathbf{A}_{pp} \theta_p^{(k + 1)} = b_p - \sum_{\substack{q = 1 \\ q \neq p}}^p \mathbf{A}_{pq} \theta_p^{(k)}, \quad p \in \left\{ 1, 2, \ldots , m \right\}.&#10;$$" src="./doc/formula/961536fc30194f28e63f2f64093c3baf.svg" align="middle" width="360.41279999999995pt" height="59.617305pt"/></p>

Where, subscript <img alt="$p$" src="./doc/formula/2ec6e630f199f589a2402fdf3e0289d5.svg" align="middle" width="8.270625000000004pt" height="14.155350000000013pt"/> is no longer the index of points, but the index of blocks.

In addition, TDMA algorithm introduces a method to directly solve *tridiagonal* systems of equations, which can be combined with block Jocobi iteration method to reduce the time consumption of iterative process. 

- When subscript <img alt="$i$" src="./doc/formula/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663295000000005pt" height="21.683310000000006pt"/> loops from <img alt="$1$" src="./doc/formula/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> to <img alt="$N_X$" src="./doc/formula/4855e1a0b6043df68024a06e6eb9d355.svg" align="middle" width="24.882495000000002pt" height="22.46574pt"/>, the <img alt="$N_Y$" src="./doc/formula/451a9ba1c48b04fdcd537308328e8b66.svg" align="middle" width="23.766105000000007pt" height="22.46574pt"/>-by-<img alt="$1$" src="./doc/formula/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> vector <img alt="$\theta_i$" src="./doc/formula/f166369f3ef0a7ff052f1e9bbf57d2e2.svg" align="middle" width="12.367905000000004pt" height="22.831379999999992pt"/> can be solved directly by TDMA algorithm:

<p align="center"><img alt="$$&#10;-a_{i, j - 1} \theta_{i, j - 1}^{(k + 1)} + a_{i, j} \theta_{i, j}^{(k + 1)} - a_{i, j + 1} \theta_{i, j + 1}^{(k + 1)} = a_{i - 1, j} \theta_{i - 1, j}^{(k)} + a_{i + 1, j} \theta_{i + 1, j}^{(k)} + b_{i, j}.&#10;$$" src="./doc/formula/7ffcd795adfa7aa2720d440bd65eb2f1.svg" align="middle" width="539.3684999999999pt" height="23.988359999999997pt"/></p>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Note that <img alt="$a_{i, 0} = a_{i, N_Y + 1} = 0$" src="./doc/formula/6886a07e1d0cd1a1c45b9f5f3ec2ea14.svg" align="middle" width="131.675445pt" height="21.18732pt"/> in this equation.

- When subscript <img alt="$j$" src="./doc/formula/36b5afebdba34564d884d347484ac0c7.svg" align="middle" width="7.710483000000004pt" height="21.683310000000006pt"/> loops from <img alt="$1$" src="./doc/formula/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> to <img alt="$N_Y$" src="./doc/formula/451a9ba1c48b04fdcd537308328e8b66.svg" align="middle" width="23.766105000000007pt" height="22.46574pt"/>, the <img alt="$N_X$" src="./doc/formula/4855e1a0b6043df68024a06e6eb9d355.svg" align="middle" width="24.882495000000002pt" height="22.46574pt"/>-by-<img alt="$1$" src="./doc/formula/034d0a6be0424bffe9a6e7ac9236c0f5.svg" align="middle" width="8.219277000000005pt" height="21.18732pt"/> vector <img alt="$\theta_j$" src="./doc/formula/455b7e5df6537b98819492ec6537494c.svg" align="middle" width="13.821390000000005pt" height="22.831379999999992pt"/> can be solved directly by TDMA algorithm:

<p align="center"><img alt="$$&#10;-a_{i - 1, j} \theta_{i - 1, j}^{(k + 1)} + a_{i, j} \theta_{i, j}^{(k + 1)} - a_{i + 1, j} \theta_{i + 1, j}^{(k + 1)} = a_{i, j - 1} \theta_{i, j - 1}^{(k)} + a_{i, j + 1} \theta_{i, j + 1}^{(k)} + b_{i, j}.&#10;$$" src="./doc/formula/c1eace5aa886d973b34d5ffa83b27af5.svg" align="middle" width="539.3684999999999pt" height="23.988359999999997pt"/></p>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Note that <img alt="$a_{0, j} = a_{N_X + 1, j} = 0$" src="./doc/formula/b2b75fc627963e9e50a9a7a6181aef2b.svg" align="middle" width="135.30198000000001pt" height="21.18732pt"/> in this equation.


#### Gauss-Seidel Iteration (Block) + TDMA Algorithm

Soon...


## Analytical Solution

Analytical solution is given by **[Wolfram Mathematica 11](https://www.wolfram.com/mathematica/)**:

<p align="center"><img alt="$$&#10;\theta(X,Y)\vert_{\tau\rightarrow\infty} = \frac{Y}{\Gamma_y} - \sum_{n=1}^\infty \frac{2\left( -1 + (-1)^n \right) \Gamma_y}{(n\pi)^2} \left[ \phi_e\cosh\frac{n\pi X}{\Gamma_y} + \phi_w\cosh\frac{n\pi (\Gamma_x - X)}{\Gamma_y} \right] \text{csch}\frac{n\pi\Gamma_x}{\Gamma_y} \sin\frac{n\pi Y}{\Gamma_y}&#10;$$" src="./doc/formula/23e304ea57b30f0e15dbf76a84fe2070.svg" align="middle" width="726.21285pt" height="44.698829999999994pt"/></p>

<!--
Non-dimensional temperature on each nodes as well as the contour is shown here. Relative codes can be referred to `analytic.f90`. 

<p align="center">
  <img height="520px" src="./doc/analytic_ctr.png"></img>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img height="520px" src="./doc/analytic_dat.png"></img>
</p>
-->


## Reference

- 吴清松. 计算热物理引论[M]. 合肥: 中国科学技术大学出版社, 2009.
- Saad Y. Iterative methods for sparse linear systems[M]. SIAM, 2003.
- Moukalled, F., Mangani, L., Darwish, M. The Finite Volume Method in Computational Fluid Dynamics: an Advanced Introduction with Openfoam® and Matlab®[M]. Springer, 2016.
- Peaceman D W, Rachford H H. The numerical solution of parabolic and elliptic differential equations[J]. Journal of the Society for industrial and Applied Mathematics, 1955, 3: 28-41.
- 吴崇试. 数学物理方法[M]. 北京: 北京大学出版社, 2003.
- 顾樵. 数学物理方法[M]. 北京: 科学出版社, 2015.


## License

[MIT](LICENSE) © Hawk Shaw, University of Science and Technology of China.