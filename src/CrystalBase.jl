module CrystalBase

using LinearAlgebra
using StaticArrays
using OrderedCollections
using DocStringExtensions

include("type.jl")
include("lattice.jl")
include("atom.jl")
include("crystal.jl")
include("kpath.jl")

include("precompile.jl")

end
