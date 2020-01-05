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

<img width="360px" align="right" src="./doc/problem_def.svg"></img>


Soon...


## Numerical Solution

Soon...

### Explicit Method

<table width="100%" border="0px" >
    <tr>
        <td><img width="300px" align="center" src="./doc/explicit_ctr.png"></img></td>
        <td><img width="310px" align="center" src="./doc/explicit_dat.svg"></img></td>
    </tr>
</table>

Soon...

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

<table width="100%" border="0px" >
    <tr>
        <td><img width="300px" align="center" src="./doc/analytic_ctr.png"></img></td>
        <td><img width="310px" align="center" src="./doc/analytic_dat.svg"></img></td>
    </tr>
</table>

## Reference

- 吴清松. 计算热物理引论[M]. 合肥: 中国科学技术大学出版社, 2009.
- Peaceman D W, Rachford H H. The numerical solution of parabolic and elliptic differential equations[J]. Journal of the Society for industrial and Applied Mathematics, 1955, 3: 28-41.
- 吴崇试. 数学物理方法[M]. 北京: 北京大学出版社, 2003.
- 顾樵. 数学物理方法[M]. 北京: 科学出版社, 2015.


## License

[MIT](LICENSE) © Hawk Shaw, University of Science and Technology of China.