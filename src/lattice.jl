export reciprocal_lattice, real_lattice, lattice_vectors, frac2cart, cart2frac

"""
    reciprocal_lattice(lattice)
    reciprocal_lattice([a1, a2, a3])
    reciprocal_lattice(a1, a2, a3)

Compute reciprocal lattice vectors from lattice vectors.

# Arguments
lattice vectors, can be
  - a matrix (each column is a lattice vector)
  - a vector of lattice vectors
  - or anything [`mat3`](@ref) accepts

# Returns
Reciprocal lattice vectors as [`Mat3`](@ref) matrix.

# Examples
```jldoctest reciprocal_lattice; setup = :(using CrystalBase)
a1, a2, a3 = [0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0];
reciprocal_lattice(a1, a2, a3)
# output
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -2.69279   1.3464     0.448799
  2.24399  -1.122      0.673198
  2.0196    0.560999  -0.336599
```

```jldoctest reciprocal_lattice
lattice = [a1, a2, a3];
reciprocal_lattice(lattice)
# output
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -2.69279   1.3464     0.448799
  2.24399  -1.122      0.673198
  2.0196    0.560999  -0.336599
```

```jldoctest reciprocal_lattice
lattice = mat3(a1, a2, a3);
reciprocal_lattice(lattice)
# output
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -2.69279   1.3464     0.448799
  2.24399  -1.122      0.673198
  2.0196    0.560999  -0.336599
```
"""
function reciprocal_lattice end

function reciprocal_lattice(lattice::Mat3)
    M = mat3(lattice)
    # Always return a Mat3 as well
    return 2π * inv(lattice)'
end

function reciprocal_lattice(lattice...)
    return reciprocal_lattice(mat3(lattice...))
end

"""
    lattice_vectors(lattice)

Return (real or reciprocal) lattice vectors from [`Mat3`](@ref) matrix columns.

# Examples
```jldoctest lattice_vectors; setup = :(using CrystalBase)
lattice = mat3([0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]);
lattice_vectors(lattice)
# output
3-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [0.0, 1.0, 2.0]
 [3.0, 0.0, 4.0]
 [5.0, 6.0, 0.0]
```
"""
function lattice_vectors(lattice::Mat3)
    return collect(vec3(lattice))
end

"""
    real_lattice(recip_lattice)

Compute real-space lattice vectors from reciprocal lattice vectors.

# Arguments
Reciprocal lattice vectors, can be
  - a matrix (each column is a reciprocal lattice vector)
  - a vector of reciprocal lattice vectors
  - or anything [`mat3`](@ref) accepts

# Returns
Real-space lattice vectors as [`Mat3`](@ref) matrix.

# Examples
```jldoctest real_lattice; setup = :(using CrystalBase)
b1, b2, b3 = [0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]
real_lattice(b1, b2, b3)
# output
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -2.69279   1.3464     0.448799
  2.24399  -1.122      0.673198
  2.0196    0.560999  -0.336599
```

```jldoctest real_lattice
recip_lattice = [b1, b2, b3];
real_lattice(recip_lattice)
# output
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -2.69279   1.3464     0.448799
  2.24399  -1.122      0.673198
  2.0196    0.560999  -0.336599
```

```jldoctest real_lattice
recip_lattice = mat3(b1, b2, b3);
real_lattice(recip_lattice)
# output
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -2.69279   1.3464     0.448799
  2.24399  -1.122      0.673198
  2.0196    0.560999  -0.336599
```

```jldoctest real_lattice
reciprocal_lattice(real_lattice(recip_lattice)) ≈ recip_lattice
# output
true
```
"""
function real_lattice(recip_lattice...)
    return reciprocal_lattice(recip_lattice...)
end

"""
    frac2cart(lattice, vec)
    frac2cart(lattice, vecs)

Convert fractional to Cartesian coordinates based on lattice vectors.

# Arguments
- `lattice`: lattice vectors.
    - For lattice vectors, the unit is usually in angstrom.
    - For reciprocal lattice vectors, the unit is usually in 1/angstrom.
- `vec`: a vector or a list of vectors in fractional coordinates.

# Examples
```jldoctest frac2cart; setup = :(using CrystalBase: frac2cart)
lattice = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]];
positions = [[0.1, 0.2, 0.3], [1.0, 2.0, 3.0]];
frac2cart(lattice, positions[1])
# output
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 2.1
 1.9
 1.0
```

```jldoctest frac2cart
frac2cart(lattice, positions)
# output
2-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [2.1, 1.9, 1.0]
 [21.0, 19.0, 10.0]
```
"""
function frac2cart end

function frac2cart(lattice, vecs::AbstractVector{<:AbstractVector})
    mat = mat3(lattice)
    carts = Ref(mat) .* vecs
    return carts
end

function frac2cart(lattice, vec::AbstractVector{<:Real})
    return mat3(lattice) * vec
end

"""
    cart2frac(lattice, vec)
    cart2frac(lattice, vecs)

Convert Cartesian to fractional coordinates based on lattice vectors.

# Arguments
- `lattice`: lattice vectors.
    - For lattice vectors, the unit is usually in angstrom.
    - For reciprocal lattice vectors, the unit is usually in 1/angstrom.
- `vec`: a vector or a list of vectors in Cartesian coordinates.

# Examples
```jldoctest cart2frac; setup = :(using CrystalBase: cart2frac, frac2cart)
lattice = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]];
positions = [[2.1, 1.9, 1.0], [21.0, 19.0, 10.0]];
frac2cart(lattice, cart2frac(lattice, positions[1])) ≈ positions[1]
# output
true
```

```jldoctest cart2frac
frac2cart(lattice, cart2frac(lattice, positions)) ≈ positions
# output
true
```
"""
function cart2frac end

function cart2frac(lattice, vecs::AbstractVector{<:AbstractVector})
    mat = mat3(lattice)
    inv_mat = inv(mat)
    frac = Ref(inv_mat) .* vecs
    return frac
end

function cart2frac(lattice, vec::AbstractVector{<:Real})
    return inv(mat3(lattice)) * vec
end
