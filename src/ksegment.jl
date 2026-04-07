export KSegment

"""
    $(TYPEDEF)

Segments of high-symmetry kpoint path, each segment is a continuous path
between high-symmetry kpoints. Different segments are disconnected in the
Brillouin zone.

See also [`KPath`](@ref) for an alternative representation of kpoint path.

# Fields
$(FIELDS)
"""
struct KSegment{T <: Real}
    "Reciprocal lattice vectors (in units of 1/L, where L is unit of length)"
    recip_lattice::Mat3{T}

    """
    Segments of high-symmetry kpoint path, each segment is a continuous path
    between high-symmetry kpoints. Different segments are disconnected in the
    Brillouin zone.
    """
    segments::Vector{Vector{String}}

    "Coordinates of high-symmetry kpoints"
    coords::OrderedDict{String, Vec3{T}}
end

function KSegment(recip_lattice::AbstractMatrix, segments::AbstractVector, coords::AbstractDict)
    rlatt = mat3(recip_lattice)
    T = eltype(rlatt)
    return KSegment{T}(
        rlatt, Vector{Vector{String}}(segments), OrderedDict{String, Vec3{T}}(coords)
    )
end

function reciprocal_lattice(kseg::KSegment)
    return kseg.recip_lattice
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
    for (label, coord) in sort(collect(kseg.coords); by = first)
        k = round.(coord; sigdigits = 4)
        print(io, "\n    $(rpad(label, maximum(length, keys(kseg.coords))))  ($(join(k, ", ")))")
    end
    return
end

function _new_label(label::AbstractString, existing_labels::Union{AbstractVector, AbstractSet})
    max_tries = 10
    i = 1
    while i < max_tries
        new_label = "$(label)_$(i)"
        if new_label ∉ existing_labels
            return new_label
        end
        i += 1
    end
    return error("Cannot find a new label?")
end

"""
    $(SIGNATURES)

# Arguments
- `recip_lattice`: each column is a reciprocal lattice vector
- `kpoint_path`: a vector of paths, each path is a vector of two kpoints, 
    each kpoint is a `Pair` of label and fractional coordinate, e.g.,
    ```julia
    kpoint_path = [
        [:Γ => [0.0, 0.0, 0.0], :M => [0.5, 0.5, 0.0]],
        [:M => [0.5, 0.5, 0.0], :R => [0.5, 0.5, 0.5]],
    ]
    ```
    This is the same as the returned `kpoint_path` of `WannierIO.read_win()`.
"""
function KSegment(
        recip_lattice::AbstractMatrix, kpoint_path::AbstractVector{<:AbstractVector{<:Pair}}
    )
    coords = OrderedDict{String, Vec3{Float64}}()
    segments = Vector{Vector{String}}()

    warn_str =
        "Two kpoints in `kpoint_path` have the same label but different " *
        "coordinates, appending a number to the label of the 2nd kpoint"

    for path in kpoint_path
        length(path) == 2 || error("each path should have 2 kpoints")
        k1, k2 = path
        # start kpoint
        label1 = String(k1.first)
        v1 = Vec3{Float64}(k1.second)
        if (label1 ∈ keys(coords)) && (coords[label1] ≉ v1)
            # convert to Tuple for better printing without type info
            @warn warn_str label = label1 k1 = Tuple(coords[label1]) k2 = Tuple(v1)
            label1 = _new_label(label1, keys(coords))
        end
        coords[label1] = v1
        # end kpoint
        label2 = String(k2.first)
        v2 = Vec3{Float64}(k2.second)
        if (label2 ∈ keys(coords)) && (coords[label2] ≉ v2)
            @warn warn_str label = label2 k1 = Tuple(coords[label2]) k2 = Tuple(v2)
            label2 = _new_label(label2, keys(coords))
        end
        coords[label2] = v2

        if (length(segments) > 0) && (label1 == segments[end][end])
            push!(segments[end], label2)
        else
            push!(segments, [label1, label2])
        end
    end

    return KSegment(recip_lattice, segments, coords)
end

function KSegment end

"""
    KSegment(lattice, atom_positions, atom_labels)

Get a `KSegment` for arbitrary cell (can be non-standard).

Requires `using Spglib, Brillouin`.

# Arguments
- `lattice`: `3 * 3`, each column is a lattice vector
- `atom_positions`: `n_atoms`-length vector of fractional coordinates
- `atom_labels` or `atom_numbers`: `n_atoms` of atomic numbers (integer) or atomic labels (string)
"""
function KSegment(lattice::AbstractMatrix, atom_positions::AbstractVector, atom_symbols::AbstractVector{<:AbstractString})
    atom_numbers = atomic_number(atom_symbols)
    return KSegment(lattice, atom_positions, atom_numbers)
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
            base, sub = split(l, "_"; limit = 2)
            return get(label_maps, base, base) * get(label_maps, sub, sub)
        end
        return get(label_maps, l, l)
    end
end

function unicode_kpoint_labels(label::AbstractString)
    return unicode_kpoint_labels([label])[1]
end

function unicode_kpoint_labels!(kseg::KSegment)
    kseg.segments = map(kseg.segments) do seg
        unicode_kpoint_labels(seg)
    end
    return kseg.coords = Dict(
        unicode_kpoint_labels(k) => v for (k, v) in kseg.coords
    )
end

function Base.isapprox(a::KSegment, b::KSegment; kwargs...)
    # segments and labels must match exactly
    if a.segments != b.segments
        return false
    end
    # compare reciprocal lattices approximately
    if !isapprox(a.recip_lattice, b.recip_lattice; kwargs...)
        return false
    end
    # compare coords: same keys and approximate numeric values
    if keys(a.coords) != keys(b.coords)
        return false
    end
    for k in keys(a.coords)
        if !isapprox(a.coords[k], b.coords[k]; kwargs...)
            return false
        end
    end
    return true
end
