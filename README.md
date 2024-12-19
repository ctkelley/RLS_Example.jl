# RLS_Example.jl
This repo has the Julia code for the example in the paper

__Min-Max Optimization for Robust Nonlinear Least Squares Problems__

by Xiaojun Chen and C. T. Kelley

You must know enough about Julia and GitHub to

   - Use the package manager
   - Know where to put projects so you can get to them with ```using XXX```
   - Clone a GitHub repo.

This repo is a Julia project __RLS_Example.jl__ . The codes in this project produce the tables and figures in the paper. To get started

  - Install the packages: LaTeXStrings, Optimization, OptimizationOptimJL, Polynomials, Printf, PyPlot, Zygote
      - [How to install packages](https://datatofish.com/install-package-julia/)

  - Clone this repo as RLS_Example.jl and either
  
    - store it as a subdirectory of any directory in your __DEPOT_PATH__ (optimal choice)
        
    - cd to the src subdirectory of __RLS_Example.jl__ and run it from in there
        
  - Start Julia and type ```using RLS_Example``` at the Julia prompt
  
To make the tables type ```tables()``` at the Julia prompt.

To make the figures type ```figures()``` at the Julia prompt.

Nothing in this repo pretends to be a general purpose project and its only mission is to make the figures/tables for the paper.

Here is an example of how one would use the code to print the tables.

```
    _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.11.2 (2024-12-01)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |


julia> using RLS_Example
[ Info: Precompiling RLS_Example [9d6a3767-96d8-5986-83e8-8cac8e739f00]

julia> tables()
Table 1
\begin{tabular}{llll} 
$\lambda$ &$\Psi_\lambda$ &$\| \nabla \Psi_\lambda \|$ &    slope \\ 
\hline 
0.00e+00 & 2.21e-02 & 1.14e-12 & 0.00e+00   \\ 
1.00e-01 & 8.30e-01 & 2.27e-12 & 7.27e-01   \\ 
1.00e+00 & 5.05e+00 & 2.49e-11 & 6.15e+00   \\ 
5.00e+00 & 7.72e+00 & 2.44e-09 & 8.16e+00   \\ 
1.00e+01 & 8.85e+00 & 6.24e-10 & 8.27e+00   \\ 
1.00e+02 & 1.25e+01 & 5.69e-08 & 8.44e+00   \\ 
\hline 
\end{tabular} 
Table 2
\begin{tabular}{llll} 
$\lambda$ &$\Psi_\lambda$ &$\| \nabla \Psi_\lambda \|$ &    slope \\ 
\hline 
0.00e+00 & 4.53e-02 & 1.03e-15 & 0.00e+00   \\ 
1.00e+00 & 1.17e+01 & 1.69e-11 & 8.29e-01   \\ 
5.00e-01 & 5.97e+00 & 2.02e-11 & 4.25e-01   \\ 
1.00e+01 & 8.02e+01 & 2.97e-12 & 8.09e+00   \\ 
1.00e+02 & 1.78e+02 & 2.80e-08 & 1.15e+01   \\ 
2.00e+02 & 2.22e+02 & 3.17e-08 & 1.17e+01   \\ 
\hline 
\end{tabular} 
Table 3
\begin{tabular}{llll} 
$\lambda$ &$\Psi_\lambda$ &$\| \nabla \Psi_\lambda \|$ &    slope \\ 
\hline 
0.00e+00 & 4.53e-02 & 1.03e-15 & 0.00e+00   \\ 
1.00e+00 & 1.17e+01 & 1.69e-11 & 8.76e+00   \\ 
5.00e-01 & 5.97e+00 & 2.02e-11 & 2.65e+00   \\ 
1.00e+01 & 8.02e+01 & 2.97e-12 & 9.47e+01   \\ 
1.00e+02 & 1.78e+02 & 2.80e-08 & 9.82e+01   \\ 
2.00e+02 & 2.22e+02 & 3.17e-08 & 9.83e+01   \\ 
\hline 
\end{tabular} 
```
