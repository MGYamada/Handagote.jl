module Handagote

using LinearAlgebra
using OMEinsum
using Random
using ChainRulesCore
using Zygote, FillArrays

include("dual.jl")
include("hyperdual.jl")
include("rand.jl")
include("dualarray.jl")
include("hyperdualarray.jl")
include("einsum.jl")
include("autodiff.jl")
include("broadcast.jl")

export @dualein_str, dualein, @hyperdualein_str, hyperdualein
export AbstractDual, Dual, AbstractDualArray, DualArray, DualVector, DualMatrix
export HyperDual, HyperDualArray, HyperDualVector, HyperDualMatrix
export ɛ, ɛ₁, ɛ₂, ɛ₁ε₂, realpart, ɛpart, hyperrealpart, ɛ₁part, ɛ₂part, ɛ₁ε₂part
export isdual, isdualarray, ishyperdual, ishyperdualarray

end