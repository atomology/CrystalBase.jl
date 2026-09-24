export reciprocal_lattice, real_lattice
export frac_to_cart, cart_to_frac

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

function reciprocal_lattice(lattice::AbstractMatrix)
    M = mat3(lattice)
    # Always return a Mat3 as well
    return 2π * inv(lattice)'
end

function reciprocal_lattice(lattice...)
    return reciprocal_lattice(mat3(lattice...))
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
    frac_to_cart(lattice, vec)
    frac_to_cart(lattice, vecs)

Convert fractional to Cartesian coordinates based on lattice vectors.

# Arguments
- `lattice`: lattice vectors.
    - For lattice vectors, the unit is usually in angstrom.
    - For reciprocal lattice vectors, the unit is usually in 1/angstrom.
- `vec`: a vector or a list of vectors in fractional coordinates.

# Examples
```jldoctest frac_to_cart; setup = :(using CrystalBase: frac_to_cart)
lattice = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]];
positions = [[0.1, 0.2, 0.3], [1.0, 2.0, 3.0]];
frac_to_cart(lattice, positions[1])
# output
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 2.1
 1.9
 1.0
```

```jldoctest frac_to_cart
frac_to_cart(lattice, positions)
# output
2-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [2.1, 1.9, 1.0]
 [21.0, 19.0, 10.0]
```
"""
function frac_to_cart end

function frac_to_cart(lattice, vecs::AbstractVector{<:AbstractVector})
    mat = mat3(lattice)
    carts = Ref(mat) .* vecs
    return carts
end

function frac_to_cart(lattice, vec::AbstractVector{<:Real})
    return mat3(lattice) * vec
end

"""
    cart_to_frac(lattice, vec)
    cart_to_frac(lattice, vecs)

Convert Cartesian to fractional coordinates based on lattice vectors.

# Arguments
- `lattice`: lattice vectors.
    - For lattice vectors, the unit is usually in angstrom.
    - For reciprocal lattice vectors, the unit is usually in 1/angstrom.
- `vec`: a vector or a list of vectors in Cartesian coordinates.

# Examples
```jldoctest cart_to_frac; setup = :(using CrystalBase: cart_to_frac, frac_to_cart)
lattice = [[0.0, 1.0, 2.0], [3.0, 0.0, 4.0], [5.0, 6.0, 0.0]];
positions = [[2.1, 1.9, 1.0], [21.0, 19.0, 10.0]];
frac_to_cart(lattice, cart_to_frac(lattice, positions[1])) ≈ positions[1]
# output
true
```

```jldoctest cart_to_frac
frac_to_cart(lattice, cart_to_frac(lattice, positions)) ≈ positions
# output
true
```
"""
function cart_to_frac end

function cart_to_frac(lattice, vecs::AbstractVector{<:AbstractVector})
    mat = mat3(lattice)
    inv_mat = inv(mat)
    frac = Ref(inv_mat) .* vecs
    return frac
end

function cart_to_frac(lattice, vec::AbstractVector{<:Real})
    return inv(mat3(lattice)) * vec
end

"""
    $(SIGNATURES)

Format a real number for display with up to 8 significant decimals.

8 decimals keeps a repeating fraction recognizable (`1/3` prints as
`0.33333333`, not an ambiguous `0.333333`) while staying short of float noise,
which a `Float64` sum only reaches around the 17th digit.

The format is always fixed point, never scientific, so a tiny residual prints
as `0.00000001` rather than an `1.0e-8` that would not line up with the column
around it. Trailing zeros are then dropped, so round values such as `0.25` stay
short. A negative zero prints as `0.0`, since its sign is an artifact of how the
value was computed, not a property of the position it describes.
"""
function _fmt(x::Real)
    y = float(x)
    # catches a tiny negative residual too, which would otherwise print as -0.0
    iszero(round(y; digits = 8)) && return "0.0"
    s = rstrip(@sprintf("%.8f", y), '0')
    return endswith(s, '.') ? s * "0" : s
end

"""
    $(SIGNATURES)

Format `values` as strings of equal width whose decimal points line up.

Unlike a fixed-precision format, trailing zeros are not padded onto the
fractional part, so `0.0` stays short while its dot still aligns with the dot
of `2.715265`.
"""
function _fmt_aligned(values)
    strs = map(_fmt, values)
    splits = map(strs) do s
        i = findfirst('.', s)
        isnothing(i) ? (s, "") : (s[1:(i - 1)], s[i:end])
    end
    width_int = maximum(length(first(p)) for p in splits; init = 0)
    width_frac = maximum(length(last(p)) for p in splits; init = 0)
    return map(p -> lpad(first(p), width_int) * rpad(last(p), width_frac), splits)
end

"""
    $(SIGNATURES)

Convert an integer `0 <= i <= 9` to its subscript character, e.g. `1 -> ₁`.
"""
function _subscript(i::Integer)
    @assert 0 <= i <= 9
    return Char(0x2080 + i)
end

"""
    $(SIGNATURES)

Print the columns of `mat` one per line, labelled `symbol` with a subscript.

All entries share one field width, so the decimal points align both down a
line and across lines.
"""
function _show_vectors(io::IO, mat::AbstractMatrix, symbol::Char; indent = "")
    fields = reshape(_fmt_aligned(vec(mat)), size(mat))
    for (i, col) in enumerate(eachcol(fields))
        println(io, rstrip(string(indent, symbol, _subscript(i), " = ", join(col, "  "))))
    end
    return
end

"""
    $(SIGNATURES)

Print `lattice` (columns are lattice vectors, in Å) as `a₁`, `a₂`, `a₃`.
"""
function show_lattice(io::IO, lattice::AbstractMatrix; indent = "")
    println(io, indent, "lattice (Å):")
    return _show_vectors(io, lattice, 'a'; indent = indent * "  ")
end

"""
    $(SIGNATURES)

Print `recip_lattice` (columns are reciprocal lattice vectors, in Å⁻¹) as
`b₁`, `b₂`, `b₃`.
"""
function show_recip_lattice(io::IO, recip_lattice::AbstractMatrix; indent = "")
    println(io, indent, "recip_lattice (Å⁻¹):")
    return _show_vectors(io, recip_lattice, 'b'; indent = indent * "  ")
end
