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
  
To make Table 5.1 type ```HistPrint()``` at the Julia prompt.

There are two figures. They differ only in the range of delta.  To make the figures type

  - ```PlotNew()``` at the prompt to get Figure 5.1 and
  - ```PlotNew(; longplot=true)``` to get Figure 5.2

Nothing in the repo pretends to be a general purpose project and its only mission is to make the figures/tables for the paper.

Here is an example of how one would use the code to print the table.

```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.10.4 (2024-06-04)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> using RLS_Example
[ Info: Precompiling RLS_Example [9d6a3767-96d8-5986-83e8-8cac8e739f00]

julia> HistPrint()
\begin{tabular}{llll} 
$\lambda$ &$\Psi_\lambda$ &$\| \nabla \Psi_\lambda \|$ &    slope \\ 
\hline 
0.00e+00 & 2.21e-02 & 1.14e-12 & 0.00e+00   \\ 
5.00e-04 & 2.65e-02 & 4.90e-13 & 4.68e-01   \\ 
1.00e-03 & 3.07e-02 & 4.84e-13 & 4.69e-01   \\ 
1.00e-02 & 1.07e-01 & 4.31e-11 & 4.78e-01   \\ 
1.00e-01 & 8.70e-01 & 8.63e-10 & 4.82e-01   \\ 
1.00e+00 & 8.49e+00 & 2.30e-08 & 4.82e-01   \\ 
\hline 
\end{tabular} 
```
