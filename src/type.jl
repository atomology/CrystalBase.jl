export vec3, mat3, stringvec3

"""
Length-3 vector type.

For atom positions, kpoints, etc.
"""
const Vec3{T} = SVector{3,T} where {T}

"""
    vec3(v)
    vec3(x, y, z)
    vec3(A)

Convert input to Vec3-compatible representation.

- For a vector input `v`, convert to `Vec3`.
- For separate `x`, `y`, `z` inputs, convert to `Vec3`.
- For a matrix input `A`, return a `Vec3` of column vectors.

!!! note

    This is not defined as a constructor of `Vec3` to avoid type piracy.

# Examples
```jldoctest vec3; setup = :(using CrystalBase, StaticArrays)
v = [1.0, 2.0, 3.0];
vec3(v)
# output
3-element SVector{3, Float64} with indices SOneTo(3):
 1.0
 2.0
 3.0
```

```jldoctest vec3
v = SVector(1.0, 2.0, 3.0);
vec3(v)
# output
3-element SVector{3, Float64} with indices SOneTo(3):
 1.0
 2.0
 3.0
```

```jldoctest vec3
vec3(1.0, 2.0, 3.0)
# output
3-element SVector{3, Float64} with indices SOneTo(3):
 1.0
 2.0
 3.0
```

```jldoctest vec3
A = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0];
vec3(A)
# output
3-element SVector{3, SVector{3, Float64}} with indices SOneTo(3):
 [1.0, 4.0, 7.0]
 [2.0, 5.0, 8.0]
 [3.0, 6.0, 9.0]
```
"""
function vec3 end

vec3(v::Vec3) = v
vec3(v::AbstractVector) = Vec3(v)
vec3(x, y, z) = Vec3(x, y, z)


"""
3 x 3 matrix type.

For lattice and reciprocal lattice.
"""
const Mat3{T} = SMatrix{3,3,T,9} where {T}

"""
    mat3(A)
    mat3(v1, v2, v3)
    mat3(cols)

Convert input to `Mat3`.

- For a matrix input `A`, convert directly to `Mat3`.
- For separate `v1`, `v2`, `v3` vector inputs, treat each as a column of the matrix.
- For a vector-of-vectors input `cols`, each inner vector is treated as one column.

!!! note

    This is not defined as a constructor of `Mat3` to avoid type piracy.

# Examples
```jldoctest mat3; setup = :(using CrystalBase, StaticArrays)
A = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0];
mat3(A)
# output
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0  2.0  3.0
 4.0  5.0  6.0
 7.0  8.0  9.0
```

```jldoctest mat3
cols = SMatrix{3,3}(A);
mat3(cols)
# output
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0  2.0  3.0
 4.0  5.0  6.0
 7.0  8.0  9.0
```

```jldoctest mat3
v1, v2, v3 = [1.0, 4.0, 7.0], [2.0, 5.0, 8.0], [3.0, 6.0, 9.0];
mat3(v1, v2, v3)
# output
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0  2.0  3.0
 4.0  5.0  6.0
 7.0  8.0  9.0
```

```jldoctest mat3
cols = [v1, v2, v3];
mat3(cols)
# output
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0  2.0  3.0
 4.0  5.0  6.0
 7.0  8.0  9.0
```
"""
function mat3 end

mat3(A::Mat3) = A
mat3(A::AbstractMatrix) = Mat3(A)
mat3(v1, v2, v3) = mat3([v1, v2, v3])

# Each vector is a column of the matrix.
mat3(cols::AbstractVector{<:AbstractVector}) = Mat3(reduce(hcat, cols))

# Each column of the matrix is a vector.
vec3(A::AbstractMatrix) = Vec3(Vec3.(eachcol(A)))


"""
Pair type associating a `String` with a `Vec3`.

E.g., for pair of atom name -> position, kpoint label -> coordinates, etc.
"""
const StringVec3{T} = Pair{String,Vec3{T}} where {T}

"""
    stringvec3(s, v)
    stringvec3(p)
    stringvec3(d)

Build a `StringVec3` pair from label/name and coordinates.

- `String` labels are stored as `String`.
- `Symbol` labels are accepted and converted to `String`.
- `Pair` and one-entry `Dict` inputs are also supported.

# Examples
```jldoctest stringvec3; setup = :(using CrystalBase, StaticArrays)
stringvec3("Si", [0.0, 0.5, 0.5])
# output
"Si" => [0.0, 0.5, 0.5]
```

```jldoctest stringvec3
stringvec3(:Γ, [0.0, 0.0, 0.0])
# output
"Γ" => [0.0, 0.0, 0.0]
```

```jldoctest stringvec3
stringvec3("K" => [0.25, 0.25, 0.25])
# output
"K" => [0.25, 0.25, 0.25]
```

```jldoctest stringvec3
stringvec3(Dict("K" => [0.25, 0.25, 0.25]))
# output
"K" => [0.25, 0.25, 0.25]
```
"""
function stringvec3 end

stringvec3(s::AbstractString, v) = StringVec3{eltype(v)}(String(s), vec3(v))
# opt-in Symbol input
stringvec3(s::Symbol, v) = stringvec3(string(s), v)
stringvec3(p::Pair) = stringvec3(p.first, p.second)
stringvec3(d::AbstractDict) = stringvec3(only(d))
