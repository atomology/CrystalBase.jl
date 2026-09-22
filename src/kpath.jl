export Subpath, KPath
export kpoints, kpoints_cart, axis, tick_indices, tick_labels, tick_positions, n_kpoints
export resample, unicode_kpoint_labels, unicode_kpoint_labels!

# A k-point path is stored as polylines: subpaths of vertices with
# per-segment division counts. The dense k-point list, the plot axis and the
# high-symmetry ticks are derived on access, never stored. Connectivity is
# structural: vertices inside one `Subpath` are connected, consecutive
# subpaths are separated by a discontinuity (the `|` of a band plot).

"""
    $(TYPEDEF)

One connected polyline of a [`KPath`](@ref).

Consecutive vertices define segments; `divisions[j]` is the number of
intervals the `j`-th segment is split into when sampled. A verbatim list of
k-points is a subpath whose divisions are all `1`.

# Fields
$(FIELDS)
"""
struct Subpath{T <: Real}
    """fractional k-point coordinates of the vertices, length `m ≥ 1`"""
    vertices::Vector{Vec3{T}}

    """one label per vertex; `""` marks an unlabeled vertex (no plot tick)"""
    labels::Vector{String}

    """intervals per segment, length `m - 1`, each `≥ 1`"""
    divisions::Vector{Int}

    function Subpath{T}(
            vertices::Vector{Vec3{T}}, labels::Vector{String}, divisions::Vector{Int}
        ) where {T <: Real}
        m = length(vertices)
        m >= 1 || throw(ArgumentError("a Subpath needs at least one vertex"))
        length(labels) == m ||
            throw(DimensionMismatch("labels has $(length(labels)) entries, expected $m"))
        length(divisions) == m - 1 ||
            throw(DimensionMismatch("divisions has $(length(divisions)) entries, expected $(m - 1)"))
        all(>=(1), divisions) || throw(ArgumentError("every division count must be ≥ 1"))
        return new{T}(vertices, labels, divisions)
    end
end

"""
    Subpath(vertices, labels = fill("", length(vertices)); divisions = 1)

Construct a [`Subpath`](@ref). `divisions` is an integer broadcast to every
segment, or one integer per segment.

# Examples
```jldoctest subpath; setup = :(using CrystalBase)
sp = Subpath([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0]], ["Γ", "X", "M"]; divisions = [4, 2]);
sp.divisions
# output
2-element Vector{Int64}:
 4
 2
```
"""
function Subpath(
        vertices::AbstractVector, labels::AbstractVector = fill("", length(vertices));
        divisions::Union{Integer, AbstractVector{<:Integer}} = 1,
    )
    T = _vertex_eltype(vertices)
    return Subpath{T}(
        Vector{Vec3{T}}(vec3.(vertices)),
        String.(labels),
        _divisions(divisions, length(vertices)),
    )
end

_divisions(d::Integer, m::Integer) = fill(Int(d), max(m - 1, 0))
_divisions(d::AbstractVector{<:Integer}, ::Integer) = Vector{Int}(d)

function _vertex_eltype(vertices::AbstractVector)
    isempty(vertices) && return Float64
    return float(promote_type(map(eltype, vertices)...))
end

n_segments(sp::Subpath) = length(sp.divisions)
n_kpoints(sp::Subpath) = 1 + sum(sp.divisions; init = 0)

"""
    $(TYPEDEF)

A k-point path in the Brillouin zone, stored as polylines, one per [`Subpath`](@ref).

Vertices within a subpath are connected; consecutive subpaths are separated
by a discontinuity. The dense k-point list ([`kpoints`](@ref)), the
cumulative plot axis ([`axis`](@ref)) and the high-symmetry ticks
([`tick_indices`](@ref), [`tick_labels`](@ref)) are derived on access.

# Fields
$(FIELDS)
"""
struct KPath{T <: Real}
    """reciprocal lattice, 3 × 3, each column is a reciprocal lattice vector (Å⁻¹)"""
    recip_lattice::Mat3{T}

    """connected pieces of the path"""
    subpaths::Vector{Subpath{T}}
end

"""
    KPath(recip_lattice, subpaths)

Construct a [`KPath`](@ref) from a vector of [`Subpath`](@ref)s.
"""
function KPath(recip_lattice::AbstractMatrix, subpaths::AbstractVector{<:Subpath})
    rlatt = mat3(recip_lattice)
    T = float(promote_type(eltype(rlatt), map(sp -> eltype(eltype(sp.vertices)), subpaths)...))
    return KPath{T}(Mat3{T}(rlatt), [_convert(Subpath{T}, sp) for sp in subpaths])
end

_convert(::Type{Subpath{T}}, sp::Subpath{T}) where {T} = sp
_convert(::Type{Subpath{T}}, sp::Subpath) where {T} =
    Subpath{T}(Vector{Vec3{T}}(sp.vertices), sp.labels, sp.divisions)

reciprocal_lattice(kpath::KPath) = kpath.recip_lattice
real_lattice(kpath::KPath) = real_lattice(kpath.recip_lattice)

# ---------------------------------------------------------------------------
# Constructors from explicit k-point lists (verbatim)
# ---------------------------------------------------------------------------

"""
    KPath(recip_lattice, kpoints; labels = fill("", length(kpoints)), break_tol = 3.0)

Wrap an explicit list of fractional `kpoints` verbatim: every point becomes a
vertex, every segment has one division, so [`kpoints`](@ref) round-trips
the input exactly.

`labels` gives one label per k-point (`""` for none). Discontinuities are
inferred from the sampling: a step longer than `break_tol` times the median
step starts a new [`Subpath`](@ref). Pass `break_tol = 0` to disable the
heuristic and keep one connected subpath.
"""
function KPath(
        recip_lattice::AbstractMatrix, kpoints::AbstractVector{<:AbstractVector{<:Real}};
        labels::AbstractVector = fill("", length(kpoints)), break_tol::Real = 3.0,
    )
    n = length(kpoints)
    n >= 1 || throw(ArgumentError("kpoints must not be empty"))
    length(labels) == n ||
        throw(DimensionMismatch("labels has $(length(labels)) entries, expected $n"))
    starts = _geometric_breaks(recip_lattice, kpoints, break_tol)
    return _verbatim(recip_lattice, kpoints, labels, starts)
end

"""
    KPath(recip_lattice, kpoints, indices, labels; break_tol = 0.0)

Wrap an explicit list of fractional `kpoints` verbatim, labeling
`kpoints[indices[i]]` with `labels[i]`. This is the layout of wannier90
`band.kpt` and `band.labelinfo.dat` files.

Two labeled points at consecutive indices mark a discontinuity, following
the wannier90 convention that a shared corner is written once. `break_tol`
additionally enables the geometric heuristic of the `labels`-keyword method.
"""
function KPath(
        recip_lattice::AbstractMatrix,
        kpoints::AbstractVector{<:AbstractVector{<:Real}},
        indices::AbstractVector{<:Integer},
        labels::AbstractVector;
        break_tol::Real = 0.0,
    )
    n = length(kpoints)
    length(indices) == length(labels) ||
        throw(DimensionMismatch("indices and labels must have the same length"))
    all(i -> 1 <= i <= n, indices) || throw(ArgumentError("indices out of range 1:$n"))
    full_labels = fill("", n)
    for (i, l) in zip(indices, labels)
        full_labels[i] = string(l)
    end
    starts = Set(_geometric_breaks(recip_lattice, kpoints, break_tol))
    for (a, b) in zip(indices[1:(end - 1)], indices[2:end])
        b == a + 1 && push!(starts, b)
    end
    return _verbatim(recip_lattice, kpoints, full_labels, sort!(collect(starts)))
end

# indices (into `kpoints`) that start a new subpath, always including 1
function _geometric_breaks(recip_lattice, kpoints, break_tol)
    starts = [1]
    (break_tol > 0 && length(kpoints) > 2) || return starts
    cart = frac_to_cart(recip_lattice, kpoints)
    steps = norm.(diff(cart))
    moved = filter(>(0), steps)
    isempty(moved) && return starts
    typical = _median(moved)
    for (i, s) in enumerate(steps)
        s > break_tol * typical && push!(starts, i + 1)
    end
    return starts
end

function _median(x::AbstractVector)
    s = sort(x)
    n = length(s)
    return isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

function _verbatim(recip_lattice, kpoints, labels, starts::AbstractVector{<:Integer})
    n = length(kpoints)
    bounds = vcat(starts, n + 1)
    subpaths = map(zip(bounds[1:(end - 1)], bounds[2:end])) do (a, b)
        Subpath(kpoints[a:(b - 1)], labels[a:(b - 1)]; divisions = 1)
    end
    return KPath(recip_lattice, subpaths)
end

# ---------------------------------------------------------------------------
# Constructor from wannier90 `kpoint_path` block
# ---------------------------------------------------------------------------

"""
    KPath(recip_lattice, kpoint_path; n_points_first_segment = 100)

Construct the polylines of a wannier90-style `kpoint_path`: a vector of
two-point segments, each a pair of `label => fractional coordinate`, as
returned for the `kpoint_path` block by `WannierIO.read_win`:

```julia
kpoint_path = [
    [:Γ => [0.0, 0.0, 0.0], :M => [0.5, 0.5, 0.0]],
    [:M => [0.5, 0.5, 0.0], :R => [0.5, 0.5, 0.5]],
]
```

Consecutive segments are joined into one [`Subpath`](@ref) when the end of one
and the start of the next have the same label and coordinate; otherwise a
discontinuity is inserted. Divisions follow wannier90: the first segment gets
`n_points_first_segment` intervals and every other segment the same spacing
(see [`resample`](@ref)), so the result reproduces `bands_num_points`.
"""
function KPath(
        recip_lattice::AbstractMatrix, kpoint_path::AbstractVector{<:AbstractVector{<:Pair}};
        n_points_first_segment::Integer = 100,
    )
    isempty(kpoint_path) && throw(ArgumentError("kpoint_path must not be empty"))
    T = float(eltype(mat3(recip_lattice)))
    vertices = Vector{Vector{Vec3{T}}}()
    labels = Vector{Vector{String}}()
    for segment in kpoint_path
        length(segment) == 2 || throw(ArgumentError("each kpoint_path entry needs exactly 2 kpoints"))
        (l1, k1), (l2, k2) = segment
        l1, l2 = string(l1), string(l2)
        v1, v2 = Vec3{T}(k1), Vec3{T}(k2)
        if !isempty(vertices) && labels[end][end] == l1 && isapprox(vertices[end][end], v1; atol = 1.0e-6)
            push!(vertices[end], v2)
            push!(labels[end], l2)
        else
            push!(vertices, [v1, v2])
            push!(labels, [l1, l2])
        end
    end
    subpaths = [Subpath(v, l; divisions = 1) for (v, l) in zip(vertices, labels)]
    return resample(KPath(recip_lattice, subpaths); n_points_first_segment)
end

"""
    resample(kpath; n_points_first_segment)
    resample(kpath; density)

Return a [`KPath`](@ref) with the same vertices and labels but new divisions.

- `n_points_first_segment`: wannier90 rule. The first segment of the first
  subpath gets this many intervals; every other segment gets the same
  spacing, rounded half up, at least 1.
- `density`: intervals per unit reciprocal length (Å); each segment gets
  `round(length · density)`, at least 1.
"""
function resample(
        kpath::KPath;
        n_points_first_segment::Union{Nothing, Integer} = nothing,
        density::Union{Nothing, Real} = nothing,
    )
    if !isnothing(n_points_first_segment) && isnothing(density)
        first_sp = first(kpath.subpaths)
        n_segments(first_sp) >= 1 ||
            throw(ArgumentError("the first subpath needs at least one segment"))
        first_length = _segment_lengths(kpath.recip_lattice, first_sp)[1]
        dk = first_length / n_points_first_segment
        return _resample(kpath, len -> len / dk)
    elseif isnothing(n_points_first_segment) && !isnothing(density)
        return _resample(kpath, len -> len * density)
    end
    throw(ArgumentError("pass exactly one of n_points_first_segment or density"))
end

function _resample(kpath::KPath{T}, intervals) where {T}
    subpaths = map(kpath.subpaths) do sp
        lengths = _segment_lengths(kpath.recip_lattice, sp)
        # Round half up to reproduce wannier90; Julia rounds half to even by default.
        divisions = [max(round(Int, intervals(len), RoundNearestTiesUp), 1) for len in lengths]
        Subpath{T}(sp.vertices, sp.labels, divisions)
    end
    return KPath{T}(kpath.recip_lattice, subpaths)
end

function _segment_lengths(recip_lattice, sp::Subpath)
    cart = frac_to_cart(recip_lattice, sp.vertices)
    return norm.(diff(cart))
end

# ---------------------------------------------------------------------------
# Derived views
# ---------------------------------------------------------------------------

# Dense points of one subpath, plus (local index, label) ticks.
function _sample(sp::Subpath{T}) where {T}
    points = Vector{Vec3{T}}(undef, n_kpoints(sp))
    ticks = Tuple{Int, String}[]
    points[1] = sp.vertices[1]
    isempty(sp.labels[1]) || push!(ticks, (1, sp.labels[1]))
    n = 1
    for j in 1:n_segments(sp)
        a, b = sp.vertices[j], sp.vertices[j + 1]
        d = sp.divisions[j]
        # Interpolate in fractional coordinates so vertices are hit exactly.
        for t in range(0, 1, d + 1)[2:end]
            n += 1
            points[n] = a + (b - a) * t
        end
        isempty(sp.labels[j + 1]) || push!(ticks, (n, sp.labels[j + 1]))
    end
    return points, ticks
end

# Dense points of the whole path, global (index, label) ticks and the set of
# global indices that start a new subpath (excluding the first).
function _sample(kpath::KPath{T}) where {T}
    points = Vec3{T}[]
    ticks = Tuple{Int, String}[]
    breaks = Set{Int}()
    for (i, sp) in enumerate(kpath.subpaths)
        offset = length(points)
        i > 1 && push!(breaks, offset + 1)
        sp_points, sp_ticks = _sample(sp)
        append!(points, sp_points)
        append!(ticks, [(idx + offset, l) for (idx, l) in sp_ticks])
    end
    return points, ticks, breaks
end

"""
    $(SIGNATURES)

Number of k-points along the path.
"""
n_kpoints(kpath::KPath) = sum(n_kpoints, kpath.subpaths; init = 0)

"""
    $(SIGNATURES)

Dense fractional k-points along the path, one segment sampled into
`divisions` intervals.
"""
kpoints(kpath::KPath) = _sample(kpath)[1]

"""
    $(SIGNATURES)

Dense k-points along the path in Cartesian coordinates (Å⁻¹).
"""
kpoints_cart(kpath::KPath) = frac_to_cart(kpath.recip_lattice, kpoints(kpath))

"""
    $(SIGNATURES)

Cumulative Cartesian distance (Å⁻¹) along the path, one value per k-point,
held flat across discontinuities between subpaths. This is the x axis of a
band-structure plot.
"""
function axis(kpath::KPath{T}) where {T}
    points, _, breaks = _sample(kpath)
    cart = frac_to_cart(kpath.recip_lattice, points)
    x = zeros(T, length(cart))
    for i in 2:length(cart)
        x[i] = x[i - 1] + (i in breaks ? zero(T) : norm(cart[i] - cart[i - 1]))
    end
    return x
end

# Ticks with labels straddling a break merged into "A|B".
function _merged_ticks(kpath::KPath)
    _, ticks, breaks = _sample(kpath)
    indices = Int[]
    labels = String[]
    for (idx, label) in ticks
        if !isempty(indices) && idx == indices[end] + 1 && idx in breaks
            labels[end] = labels[end] * "|" * label
        else
            push!(indices, idx)
            push!(labels, label)
        end
    end
    return indices, labels
end

"""
    tick_indices(kpath; merge = true)

Indices into [`kpoints`](@ref) of the labeled vertices.

With `merge = true` two labeled points straddling a discontinuity share one
tick, labeled `A|B`, as a band plot draws them. With `merge = false` every
labeled vertex is its own tick, the layout of wannier90 `labelinfo.dat`.
"""
tick_indices(kpath::KPath; merge::Bool = true) = _ticks(kpath, merge)[1]

"""
    tick_labels(kpath; merge = true)

Labels of the ticks, aligned with [`tick_indices`](@ref).
"""
tick_labels(kpath::KPath; merge::Bool = true) = _ticks(kpath, merge)[2]

function _ticks(kpath::KPath, merge::Bool)
    merge && return _merged_ticks(kpath)
    ticks = _sample(kpath)[2]
    return first.(ticks), last.(ticks)
end

"""
    tick_positions(kpath; merge = true)

Position of each tick on [`axis`](@ref).
"""
tick_positions(kpath::KPath; merge::Bool = true) = axis(kpath)[tick_indices(kpath; merge)]

# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------

const _UNICODE_LABELS = Dict(
    "GAMMA" => "Γ",
    "DELTA" => "Δ",
    "LAMBDA" => "Λ",
    "SIGMA" => "Σ",
    "0" => "₀",
    "1" => "₁",
    "2" => "₂",
    "3" => "₃",
    "4" => "₄",
    "5" => "₅",
    "6" => "₆",
    "7" => "₇",
    "8" => "₈",
    "9" => "₉",
)

"""
    unicode_kpoint_labels(label)
    unicode_kpoint_labels(labels)
    unicode_kpoint_labels(kpath)

Convert high-symmetry k-point labels to Unicode: `GAMMA` → `Γ`, and a
`_n` suffix to a subscript. For a [`KPath`](@ref) a new path is returned;
see [`unicode_kpoint_labels!`](@ref) to convert in place.

# Examples
```jldoctest unicode_kpoint_labels; setup = :(using CrystalBase)
unicode_kpoint_labels(["GAMMA", "DELTA_0", "LAMBDA_1", "SIGMA_2", "X"])
# output
5-element Vector{String}:
 "Γ"
 "Δ₀"
 "Λ₁"
 "Σ₂"
 "X"
```
"""
function unicode_kpoint_labels(label::AbstractString)
    if occursin("_", label)
        base, sub = split(label, "_"; limit = 2)
        return get(_UNICODE_LABELS, base, base) * get(_UNICODE_LABELS, sub, sub)
    end
    return get(_UNICODE_LABELS, label, String(label))
end

unicode_kpoint_labels(labels::AbstractVector{<:AbstractString}) = map(unicode_kpoint_labels, labels)

function unicode_kpoint_labels(kpath::KPath{T}) where {T}
    subpaths = [
        Subpath{T}(sp.vertices, unicode_kpoint_labels(sp.labels), sp.divisions)
            for sp in kpath.subpaths
    ]
    return KPath{T}(kpath.recip_lattice, subpaths)
end

"""
    $(SIGNATURES)

Convert the labels of `kpath` to Unicode in place, see [`unicode_kpoint_labels`](@ref).
"""
function unicode_kpoint_labels!(kpath::KPath)
    for sp in kpath.subpaths
        map!(unicode_kpoint_labels, sp.labels, sp.labels)
    end
    return kpath
end

# ---------------------------------------------------------------------------
# Comparison and printing
# ---------------------------------------------------------------------------

function Base.:(==)(a::Subpath, b::Subpath)
    return a.vertices == b.vertices && a.labels == b.labels && a.divisions == b.divisions
end

function Base.isapprox(a::Subpath, b::Subpath; kwargs...)
    a.labels == b.labels || return false
    a.divisions == b.divisions || return false
    length(a.vertices) == length(b.vertices) || return false
    return all(isapprox.(a.vertices, b.vertices; kwargs...))
end

Base.:(==)(a::KPath, b::KPath) = a.recip_lattice == b.recip_lattice && a.subpaths == b.subpaths

function Base.isapprox(a::KPath, b::KPath; kwargs...)
    isapprox(a.recip_lattice, b.recip_lattice; kwargs...) || return false
    length(a.subpaths) == length(b.subpaths) || return false
    return all(isapprox(x, y; kwargs...) for (x, y) in zip(a.subpaths, b.subpaths))
end

# Number of vertices up to which a fully labeled subpath is printed as corners.
const _SHOW_CORNERS_LIMIT = 12

# Describe a subpath in one line: corner form when every vertex is labeled and
# the list is short, tick form otherwise. `offset` shifts tick indices to the
# global k-point numbering.
function _describe(io::IO, sp::Subpath, offset::Integer = 0)
    labels = get(io, :unicode, true) ? unicode_kpoint_labels(sp.labels) : sp.labels
    m = length(sp.vertices)
    if all(!isempty, labels) && m <= _SHOW_CORNERS_LIMIT
        print(io, join(labels, "—"))
        isempty(sp.divisions) || print(io, "  divisions ", join(sp.divisions, " "))
    else
        print(io, n_kpoints(sp), " kpoints")
        _, ticks = _sample(sp)
        if isempty(ticks)
            print(io, ", no ticks")
        else
            print(io, "  ticks ", join((string(unicode_kpoint_labels(l), "@", i + offset) for (i, l) in ticks), " "))
        end
    end
    return
end

function Base.show(io::IO, sp::Subpath)
    print(io, "Subpath(")
    _describe(io, sp)
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", sp::Subpath{T}) where {T}
    print(io, "Subpath{$T}: ", length(sp.vertices), " vertices, ", n_kpoints(sp), " kpoints\n  ")
    return _describe(io, sp)
end

Base.summary(io::IO, kpath::KPath) =
    print(io, "KPath(", length(kpath.subpaths), " subpaths, ", n_kpoints(kpath), " kpoints)")

Base.show(io::IO, kpath::KPath) = summary(io, kpath)

function Base.show(io::IO, ::MIME"text/plain", kpath::KPath{T}) where {T}
    println(io, "KPath{$T}: ", length(kpath.subpaths), " subpaths, ", n_kpoints(kpath), " kpoints")
    show_recip_lattice(io, kpath.recip_lattice; indent = "  ")
    limit = get(io, :limit, false) ? 20 : typemax(Int)
    offset = 0
    for (i, sp) in enumerate(kpath.subpaths)
        if i > limit
            print(io, "  ⋮ (", length(kpath.subpaths) - limit, " more subpaths)")
            break
        end
        i > 1 && println(io)
        print(io, "  ", i, ": ")
        _describe(io, sp, offset)
        offset += n_kpoints(sp)
    end
    return
end
