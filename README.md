# CrystalUtils

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://atomology.github.io/CrystalUtils.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://atomology.github.io/CrystalUtils.jl/dev/)
[![Build Status](https://github.com/atomology/CrystalUtils.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/atomology/CrystalUtils.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/atomology/CrystalUtils.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/atomology/CrystalUtils.jl)

`CrystalUtils.jl` is a tiny Julia package for working with crystal structures in
real and reciprocal space.

The package only defines a few container types and provides some utility functions.
It does not involve heavy computations but just to avoid rewriting boilerplate code
in other packages.

## Example

```julia
julia> using CrystalUtils

# Define a matrix for lattice vectors
julia> a1, a2, a3 = [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]
julia> lattice = mat3(a1, a2, a3)
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0  0.0  0.0
 0.0  2.0  0.0
 0.0  0.0  3.0

# Compute the reciprocal lattice (as a matrix)
julia> recip_lattice = reciprocal_lattice(lattice)
3×3 SMatrix{3, 3, Float64, 9} with indices
    SOneTo(3)×SOneTo(3):
    6.283185307179586  0.0                0.0
    0.0                3.141592653589793  0.0
    0.0                0.0                2.0943951023931953

# Get the lattice vectors, which are the columns of the matrix
julia> lattice_vectors(recip_lattice)

# Compute the real-space lattice vectors
julia> real_lattice(recip_lattice)

# Define a kpoint path
julia> KPath()
```
