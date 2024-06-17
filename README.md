# RLS_Example.jl
This repo has the Julia code for the example in the paper

__Min-Max Optimization for Robust Nonlinear Least Squares Problems__

by Xiaojun Chen and C. T. Kelley

You must know enough about Julia to

   - Use the package manager
     
   - Know where to put projects so you can get to them with ```using XXX```

This repo is a Julia project __RLS_Example.jl__ . The codes in this project produce the tables and figures in the paper. To get started

  - Clone this repo as RLS_Example.jl and either
  
    - store it as a subdirectory of any directory in your __DEPOT_PATH__ (optimal choice)
        
    - cd to the src subdirectory of __RLS_Example.jl__ and run it from in there
        
  - Fire up Julia and type ```using RLS_Example``` at the Julia prompt
  
To make Table 5.1 type ```HistPrint()``` at the Julia prompt.

There are two figures. They differ only in the range of delta.  To make the figures type

  - ```PlotNew()``` at the prompt to get Figure 5.1 and
  - ```PlotNew(; longplot=true)``` to get Figure 5.2

Nothing in the repo pretends to be a general purpose project and its only mission is to make the figures/tables for the paper.
