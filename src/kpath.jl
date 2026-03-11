export KPath, KSegment, linear_path

"""
    $(TYPEDEF)

Kpoint path in the Brillouin zone.

Storing explicitly the kpoint coordinates along the path, as well as the indices
and labels of high-symmetry kpoints.

See also [`KSegment`](@ref) for an alternative representation of kpoint path,
which only stores the segments of high-symmetry kpoint path but not the explicit
kpoint coordinates along the path.

# Fields
$(FIELDS)
"""
struct KPath{T<:Real}
    "Reciprocal lattice vectors (in units of 1/L, where L is unit of length)"
    recip_lattice::Mat3{T}

    "Fractional kpoint coordinates along the kpath"
    points::Vector{Vec3{T}}

    "Indices of high-symmetry kpoints along the kpath"
    indices::Vector{Int}

    "Labels of high-symmetry kpoints"
    labels::Vector{String}
end

function KPath(recip_lattice, points, indices, labels)
    rlatt = mat3(recip_lattice)
    T = eltype(rlatt)
    return KPath{T}(rlatt, Vector{Vec3{T}}(points), Vector{Int}(indices), string.(labels))
end

function Base.length(kpath::KPath)
    return length(kpath.points)
end

"""
    $(SIGNATURES)

Convert labels of high-symmetry kpoints in `kpath` to unicode string.
"""
function unicode_kpoint_labels!(kpath::KPath)
    kpath.labels = unicode_kpoint_labels(kpath.labels)
end

function Base.show(io::IO, kpath::KPath)
    n_kpts = length(kpath.points)
    return print(io, "KPath($(n_kpts) kpoints, [$(join(labels, " → "))])")
end

function Base.show(io::IO, ::MIME"text/plain", kpath::KPath{T}) where {T}
    n_kpts = length(kpath.points)
    println(io, "KPath{$T} with $(n_kpts) k-points and $(length(labels)) high-symmetry labels:")
    println(io, "  Path: $(join(labels, " → "))")
    println(io, "  High-symmetry k-points:")
    for (i, (idx, lab)) in enumerate(zip(kpath.indices, labels))
        k = kpath.points[idx]
        print(io, "    $(lpad(idx, ndigits(n_kpts))): $(rpad(lab, maximum(length, labels)))  ")
        print(io, "($(join(round.(k, sigdigits=4), ", ")))")
        i < length(kpath.indices) && println(io)
    end
end

function reciprocal_lattice(kpath::KPath)
    return kpath.recip_lattice
end

"""
    $(TYPEDEF)

Segments of high-symmetry kpoint path, each segment is a continuous path
between high-symmetry kpoints. Different segments are disconnected in the
Brillouin zone.

See also [`KPath`](@ref) for an alternative representation of kpoint path.

# Fields
$(FIELDS)
"""
struct KSegment{T<:Real}
    "Reciprocal lattice vectors (in units of 1/L, where L is unit of length)"
    recip_lattice::Mat3{T}

    """
    Segments of high-symmetry kpoint path, each segment is a continuous path
    between high-symmetry kpoints. Different segments are disconnected in the
    Brillouin zone.
    """
    segments::Vector{Vector{String}}

    "Coordinates of high-symmetry kpoints"
    coords::OrderedDict{String,Vec3{T}}
end

function KSegment(recip_lattice, segments, coords)
    rlatt = mat3(recip_lattice)
    T = eltype(rlatt)
    return KSegment{T}(
        rlatt, Vector{Vector{String}}(segments), OrderedDict{String,Vec3{T}}(coords)
    )
end

function reciprocal_lattice(kseg::KSegment)
    return kseg.recip_lattice
end

function unicode_kpoint_labels!(kseg::KSegment)
    kseg.segments = map(kseg.segments) do seg
        unicode_kpoint_labels(seg)
    end
    kseg.coords = Dict(
        unicode_kpoint_labels(k) => v for (k, v) in kseg.coords
    )
end

function Base.show(io::IO, kseg::KSegment)
    segs = [join(s, " → ") for s in kseg.segments]
    return print(io, "KSegment([$(join(segs, " | "))])")
end

function Base.show(io::IO, ::MIME"text/plain", kseg::KSegment{T}) where {T}
    n_seg = length(kseg.segments)
    n_pts = length(kseg.coords)
    println(io, "KSegment{$T} with $(n_seg) segment$(n_seg == 1 ? "" : "s") and $(n_pts) high-symmetry k-point$(n_pts == 1 ? "" : "s"):")
    for (i, seg) in enumerate(kseg.segments)
        println(io, "  Segment $i: $(join(seg, " → "))")
    end
    print(io, "  Coordinates:")
    for (label, coord) in sort(collect(kseg.coords); by=first)
        k = round.(coord; sigdigits=4)
        print(io, "\n    $(rpad(label, maximum(length, keys(kseg.coords))))  ($(join(k, ", ")))")
    end
end

"""
    $(SIGNATURES)

Convert labels of high-symmetry kpoints to unicode string.

# Examples
```jldoctest unicode_kpoint_labels; setup = :(using CrystalBase)
CrystalBase.unicode_kpoint_labels(["GAMMA", "DELTA_0", "LAMBDA_1", "SIGMA_2", "X"])
# output
5-element Vector{String}:
 "Γ"
 "Δ₀"
 "Λ₁"
 "Σ₂"
 "X"
```
"""
function unicode_kpoint_labels(labels::AbstractVector{<:AbstractString})
    label_maps = Dict(
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
    return map(labels) do l
        if occursin("_", l)
            base, sub = split(l, "_"; limit=2)
            return get(label_maps, base, base) * get(label_maps, sub, sub)
        end
        return get(label_maps, l, l)
    end
end

function unicode_kpoint_labels(label::AbstractString)
    return unicode_kpoint_labels([label])[1]
end

"""
    $(SIGNATURES)

Group consecutive high-symmetry points.

If two high-symmetry kpoints are neighbors, group them together.

# Arguments
- `indices`: indices of high-symmetry kpoints, start from 1

# Return
- `groups`: a vector of vectors, if the inner vector contains more than 1 index,
    it means those indices are neighbors and are grouped together.

# Example
```jldoctest group_nearby_indices; setup = :(using CrystalBase)
CrystalBase.group_nearby_indices([1, 2, 4, 5, 6])
# output
2-element Vector{Vector{Int64}}:
 [1, 2]
 [4, 5, 6]
```
"""
function group_nearby_indices(indices::AbstractVector{T}) where {T<:Integer}
    groups = Vector{Vector{T}}()
    isempty(indices) && return groups
    push!(groups, [indices[1]])

    counter = 2
    for i in eachindex(indices)[2:end]
        if indices[i] == indices[i - 1] + 1
            push!(groups[counter - 1], indices[i])
        else
            push!(groups, [indices[i]])
            counter += 1
        end
    end

    return groups
end

"""
    $(SIGNATURES)

Merge consecutive high-symmetry points.

If two high-symmetry kpoints are neighbors, merge them into one,
with label `X|Y`, where `X` and `Y` are the original labels of
the two kpoints, respectively.

# Arguments
- `indices`: indices of high-symmetry kpoints, start from 1
- `labels`: labels of high-symmetry kpoints
"""
function merge_nearby_labels(
    indices::AbstractVector{<:Integer}, labels::AbstractVector{<:AbstractString}
)
    grps = group_nearby_indices(indices)
    labs = map(grps) do idxs
        js = map(idxs) do i
            findfirst(==(i), indices)
        end
        join(labels[js], "|")
    end
    jdxs = isempty(grps) ? grps : first.(grps)
    return jdxs, labs
end

function linear_path(
    kpoints_cart::AbstractVector{<:AbstractVector{<:Real}},
    indices::AbstractVector{<:Integer}=[],
)
    grps = group_nearby_indices(indices)
    x = [0.0; accumulate(+, norm.(diff(kpoints_cart)))]

    for idxs in grps
        (length(idxs) > 1) || continue
        for j in idxs[2:end]
            δ = x[j] - x[j - 1]
            x[j:end] .-= δ
        end
    end

    return x
end

"""
    $(SIGNATURES)

Get a 1D vector of cumulative distance along the kpath.
"""
function linear_path(kpath::KPath)
    kpts_cart = frac2cart(kpath.recip_lattice, kpath.points)
    x = linear_path(kpts_cart, kpath.indices)
    return x
end
