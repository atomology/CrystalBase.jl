# CrystalBase

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://atomology.github.io/CrystalBase.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://atomology.github.io/CrystalBase.jl/dev/)
[![Build Status](https://github.com/atomology/CrystalBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/atomology/CrystalBase.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/atomology/CrystalBase.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/atomology/CrystalBase.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

`CrystalBase.jl` is a tiny Julia package for working with crystal structures in
real and reciprocal space.

The package only defines a few container types and provides some utility functions.
It does not involve heavy computations but just to avoid rewriting boilerplate code
in other packages.

## Example

### Real and reciprocal lattice vectors

Lattices are 3 × 3 matrices whose **columns** are the lattice vectors.

```julia
julia> using CrystalBase

# Stack the lattice vectors as columns
julia> a1, a2, a3 = [1.0, 0.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 3.0]
julia> lattice = mat3(a1, a2, a3)
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0  1.0  0.0
 0.0  2.0  0.0
 0.0  0.0  3.0

julia> lattice[:, 2] == a2
true

# Compute the reciprocal lattice, again one vector per column
julia> recip_lattice = reciprocal_lattice(lattice)
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  6.28319  0.0      0.0
 -3.14159  3.14159  0.0
  0.0      0.0      2.0944

# Split the matrix into its columns, the reciprocal lattice vectors
julia> b1, b2, b3 = vec3(recip_lattice)
3-element StaticArraysCore.SVector{3, StaticArraysCore.SVector{3, Float64}} with indices SOneTo(3):
 [6.283185307179586, -3.141592653589793, 0.0]
 [0.0, 3.141592653589793, 0.0]
 [0.0, 0.0, 2.0943951023931953]

# aᵢ ⋅ bⱼ = 2π δᵢⱼ
julia> [a1 a2 a3]' * b1
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 6.283185307179586
 0.0
 0.0

# Compute the real-space lattice back
julia> real_lattice(recip_lattice) ≈ lattice
true
```

### Fractional to Cartesian coordinates interconversion

```julia
# Fractional coordinates weight the columns: 0.5 a1 + 0.5 a2 + 0.5 a3
julia> frac_coords = [0.5, 0.5, 0.5]
julia> cart_coords = frac_to_cart(lattice, frac_coords)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 1.0
 1.0
 1.5

# Also support multiple coordinates at once
julia> frac_to_cart(lattice, [frac_coords, [0.0, 1.0, 0.0]])
2-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [1.0, 1.0, 1.5]
 [1.0, 2.0, 0.0]

# Convert back to fractional coordinates
julia> cart_to_frac(lattice, cart_coords)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 0.5
 0.5
 0.5
```

### Crystal structures

```julia
julia> lattice = [0.0 2.715 2.715; 2.715 0.0 2.715; 2.715 2.715 0.0]
julia> si = Crystal(lattice, [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]], ["Si", "Si"])
Crystal{Float64}: Si2, 2 atoms
  lattice (Å):
    a₁ = 0.0    2.715  2.715
    a₂ = 2.715  0.0    2.715
    a₃ = 2.715  2.715  0.0
  atoms (fractional):
    Si  0.0   0.0   0.0
    Si  0.25  0.25  0.25

julia> n_atoms(si), formula(si), reciprocal_lattice(si) ≈ reciprocal_lattice(lattice)
(2, "Si2", true)

# Per-site labels refine the element: `Fe1`/`Fe2` are two Fe sites
julia> fe = Crystal(lattice, ["Fe1" => [0.0, 0.0, 0.0], "Fe2" => [0.5, 0.5, 0.5]])
julia> atom_symbols(fe), unique_species(fe)
(["Fe", "Fe"], ["Fe1", "Fe2"])

# With `using Spglib`: space group and symmetry-consistent k-point grids
julia> using Spglib
julia> spacegroup(si)
(symbol = "Fd-3m", number = 227)
julia> kgrid_from_density(si, 3.0)
(7, 7, 7)
```

### K-point paths in the Brillouin zone

A `KPath` stores polylines: connected `Subpath`s of labeled vertices with a
division count per segment. The dense k-points, the plot axis and the ticks
are derived on demand.

```julia
# From wannier90-style segments; divisions follow `bands_num_points`
julia> recip_lattice = reciprocal_lattice(lattice)
julia> kpoint_path = [["Γ" => [0.0, 0.0, 0.0], "X" => [0.5, 0.0, 0.5]],
                      ["X" => [0.5, 0.0, 0.5], "U" => [0.625, 0.25, 0.625]]]
julia> kp = KPath(recip_lattice, kpoint_path; n_points_first_segment = 5)
KPath{Float64}: 1 subpaths, 8 kpoints
  recip_lattice (Å⁻¹):
    b₁ = -1.15712437   1.15712437   1.15712437
    b₂ =  1.15712437  -1.15712437   1.15712437
    b₃ =  1.15712437   1.15712437  -1.15712437
  1: Γ—X—U  divisions 5 2

julia> kpoints(kp)[1:3]
3-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [0.0, 0.0, 0.0]
 [0.1, 0.0, 0.1]
 [0.2, 0.0, 0.2]

julia> axis(kp), tick_indices(kp), tick_labels(kp)
([0.0, 0.2314, …], [1, 6, 8], ["Γ", "X", "U"])

# Wrap an explicit k-point list verbatim (e.g. from wannier90 band.kpt + labelinfo)
julia> KPath(recip_lattice, kpoints(kp), [1, 6, 8], ["Γ", "X", "U"])

# Change the sampling without touching the vertices
julia> resample(kp; density = 50.0)

# With `using Spglib, Brillouin`: the standard path of a crystal
julia> KPath(si)
KPath{Float64}: 2 subpaths, 451 kpoints
  ...
  1: Γ—X—U  divisions 100 35
  2: K—Γ—L—W—X  divisions 106 87 71 50
```
