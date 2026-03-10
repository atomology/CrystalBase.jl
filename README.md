# CrystalBase

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://atomology.github.io/CrystalBase.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://atomology.github.io/CrystalBase.jl/dev/)
[![Build Status](https://github.com/atomology/CrystalBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/atomology/CrystalBase.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/atomology/CrystalBase.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/atomology/CrystalBase.jl)

`CrystalBase.jl` is a tiny Julia package for working with crystal structures in
real and reciprocal space.

The package only defines a few container types and provides some utility functions.
It does not involve heavy computations but just to avoid rewriting boilerplate code
in other packages.

## Example

```julia
julia> using CrystalBase

# Define a matrix for lattice vectors
julia> a1, a2, a3 = [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]
julia> lattice = mat3(a1, a2, a3)
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0  0.0  0.0
 0.0  2.0  0.0
 0.0  0.0  3.0

# Compute the reciprocal lattice (as a matrix)
julia> recip_lattice = reciprocal_lattice(lattice)
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 6.28319  0.0      0.0
 0.0      3.14159  0.0
 0.0      0.0      2.0944

# Get the lattice vectors, which are the columns of the matrix
julia> lattice_vectors(recip_lattice)
3-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [6.283185307179586, 0.0, 0.0]
 [0.0, 3.141592653589793, 0.0]
 [0.0, 0.0, 2.0943951023931953]

# Compute the real-space lattice vectors
julia> real_lattice(recip_lattice)
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0  0.0  0.0
 0.0  2.0  0.0
 0.0  0.0  3.0

# Define a kpoint path
julia> points = [[i/10, 0.0, 0.0] for i in 0:5]
julia> indices = [1, 6]
julia> labels = ["Γ", "X"]
julia> KPath(recip_lattice, points, indices, labels)
KPath{Float64}([6.283185307179586 0.0 0.0; 0.0 3.141592653589793 0.0; 0.0 0.0 2.0943951023931953], StaticArraysCore.SVector{3, Float64}[[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0], [0.3, 0.0, 0.0], [0.4, 0.0, 0.0], [0.5, 0.0, 0.0]], [1, 6], ["Γ", "X"])
```
