module RLS_Example

using LinearAlgebra
using Optimization
using OptimizationOptimJL
using Zygote
using Printf
using Polynomials
using PyPlot
using LaTeXStrings

include("D1.jl")
include("G2.jl")
include("Green_Utils.jl")
include("Generic_Utils.jl")
include("HistPrint.jl")
include("fprintTeX.jl")

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
