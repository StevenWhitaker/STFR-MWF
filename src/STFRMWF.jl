module STFRMWF

# Function for referencing data/results files
modulepath(file) = joinpath(@__DIR__, file)

# Import dependencies
using SymPy
using ScanDesign, Distributions, JLD, Random
using ProgressMeter, Plots, LaTeXStrings
using ForwardDiff, Statistics
using StatsBase
using STFR, PERK
using LinearAlgebra
using BlochSim, Interpolations, NNLS
using MAT, ImageTransformations

# Make sure both PyPlot and GR are loaded
pyplot()
gr()

# Include source code
include("./scandesign/stfrscandesign.jl")
include("./scandesign/designA.jl")
include("./scandesign/designB.jl")
include("./scandesign/table2.jl")
include("./scandesign/figure2.jl")

include("./estimation/FitDist.jl")
include("./estimation/designA.jl")
include("./estimation/designAsorted.jl")
include("./estimation/designB.jl")
include("./estimation/stfrperk.jl")
include("./estimation/stfrblochsimperk.jl")
include("./estimation/mese.jl")
include("./estimation/brainweb.jl")
include("./estimation/testdata.jl")
include("./estimation/figure3.jl")
include("./estimation/figure4table3.jl")
include("./estimation/figure5.jl")
include("./estimation/figure6table4.jl")
include("./estimation/figureS1.jl")
include("./estimation/figureS2tableS1.jl")
include("./estimation/figureS3tableS2.jl")
include("./estimation/figureS4.jl")

include("./bias/stfrperk_bias.jl")
include("./bias/stfrblochsimperk_bias.jl")
include("./bias/crlb.jl")
include("./bias/designA.jl")
include("./bias/figureS5.jl")
include("./bias/figureS6.jl")

end
