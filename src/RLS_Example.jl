module RLS_Example

using LinearAlgebra
using Optimization
using OptimizationOptimJL
using Zygote
using Printf
using Polynomials
using PyPlot
using LaTeXStrings
using SparseArrays

include("Conv_Test.jl")
include("D1.jl")
include("Green_Example.jl")
include("Green_Utils.jl")
include("Generic_Utils.jl")
include("HistPrint.jl")
include("fprintTeX.jl")
include("Cmat.jl")
include("tables.jl")

export tables
export figures
export Conv_Test
export G2
export SetDEData
export HistPrint
export PlotNew
export Basic_Example
export Chen_Example
export PlotSimp
export T2
export ST2
export One_Shot2
#export Ex_Paper
export sL1
export do_opt
export znorm
export OptSimp

end
