export KPath, linear_path

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
struct KPath{T <: Real}
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

function Base.collect(kpath::KPath)
    return copy(kpath.points)
end

function Base.show(io::IO, kpath::KPath)
    n_kpts = length(kpath.points)
    return print(io, "KPath($(n_kpts) kpoints, [$(join(labels, " → "))])")
end

function Base.show(io::IO, ::MIME"text/plain", kpath::KPath{T}) where {T}
    n_kpts = length(kpath.points)
    println(io, "KPath{$T} with $(n_kpts) k-points and $(length(kpath.labels)) high-symmetry labels:")
    println(io, "  Path: $(join(kpath.labels, " → "))")
    println(io, "  High-symmetry k-points:")
    for (i, (idx, lab)) in enumerate(zip(kpath.indices, kpath.labels))
        k = kpath.points[idx]
        print(io, "    $(lpad(idx, ndigits(n_kpts))): $(rpad(lab, maximum(length, kpath.labels)))  ")
        print(io, "($(join(round.(k, sigdigits = 4), ", ")))")
        (i < length(kpath.indices)) && println(io)
    end
    return
end

function reciprocal_lattice(kpath::KPath)
    return kpath.recip_lattice
end

"""
    $(SIGNATURES)

Convert labels of high-symmetry kpoints in `kpath` to unicode string.
"""
function unicode_kpoint_labels!(kpath::KPath)
    kpath.labels = unicode_kpoint_labels(kpath.labels)
    return nothing
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
function group_nearby_indices(indices::AbstractVector{T}) where {T <: Integer}
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

# Return
- `tick_indices`: indices of high-symmetry kpoints after merging
- `tick_labels`: labels of high-symmetry kpoints after merging
"""
function merge_nearby_labels(
        indices::AbstractVector{<:Integer}, labels::AbstractVector{<:AbstractString}
    )
    grps = group_nearby_indices(indices)
    tick_labs = map(grps) do idxs
        js = map(idxs) do i
            findfirst(==(i), indices)
        end
        join(labels[js], "|")
    end
    tick_idxs = isempty(grps) ? grps : first.(grps)
    return tick_idxs, tick_labs
end

"""
    $(SIGNATURES)

Get a 1D vector of cumulative distance along the kpath.
"""
function linear_path(
        kpoints_cart::AbstractVector{<:AbstractVector{<:Real}},
        indices::AbstractVector{<:Integer} = [],
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

function linear_path(kpoints_cart, indices, labels)
    x = linear_path(kpoints_cart, indices)
    tick_idxs, tick_labs = merge_nearby_labels(indices, labels)
    return x, tick_idxs, tick_labs
end

"""
    $(SIGNATURES)

Get a 1D vector of cumulative distance along the kpath, and the corresponding
tick indices and labels for high-symmetry kpoints.

# Return
- `x`: 1D vector of cumulative distance along the kpath
- `tick_indices`: indices of high-symmetry kpoints after merging
- `tick_labels`: labels of high-symmetry kpoints after merging
"""
function linear_path(kpath::KPath)
    kpts_cart = frac_to_cart(kpath.recip_lattice, kpath.points)
    return linear_path(kpts_cart, kpath.indices, kpath.labels)
end

"""
    $(SIGNATURES)

Generate a `KPath` containing kpoint coordinates that are exactly
the same as wannier90.

The kpoints are generated by the following criteria:
- the kpath spacing of remaining segments are kept the same as the first segment
- merge same high-symmetry labels at the corner between two segments; keep both
    labels if the two labels (ending of the 1st segment and starting point of the
    2nd segment) are different

# Arguments
- `kseg`: a `KSegment`
- `n_points_first_segment`: number of kpoints in the first segment, remaining
    segments will have the same spacing as the 1st segment. The default value
    is 100, which is the same as wannier90 default value of `bands_num_points`.

!!! note

    This reproduce exactly the wannier90 behavior, if
    - the `kseg` is generated from the `kpoint_path` obtained by
            `WannierIO.read_win` which parses the `kpoint_path` block of `win` file,
    - the `n_points` is the same as wannier90 `win` file input parameter `bands_num_points`,
        which again can be obtained by `WannierIO.read_win`.
"""
function KPath(kseg::KSegment, n_points_first_segment::Integer = 100)
    # Cartesian
    coords_cart = OrderedDict(k => frac_to_cart(kseg.recip_lattice, v) for (k, v) in kseg.coords)

    # kpath spacing from first two kpoints
    isempty(kseg.segments) && error("kseg should have at least one segment")
    (length(kseg.segments[1]) < 2) && error("the first segment should have at least two kpoints")
    k1, k2 = kseg.segments[1][1:2]
    v = coords_cart[k2] - coords_cart[k1]
    v_norm = norm(v)
    dk = v_norm / n_points_first_segment

    # kpoints along each segment
    kpaths = Vector{Vector{Vec3{Float64}}}()
    # symmetry points along each segment
    indices = Vector{Vector{Int}}()
    labels = Vector{Vector{String}}()

    for seg in kseg.segments
        kpoints_seg = Vector{Vec3{Float64}}()
        indices_seg = Vector{Int}()
        labels_seg = Vector{String}()

        n_seg = length(seg) - 1
        n_x_seg = 0
        for j in 1:n_seg
            k1 = seg[j]
            k2 = seg[j + 1]

            # One vector in each segment, from k1 to k2
            v = coords_cart[k2] - coords_cart[k1]
            v_norm = norm(v)

            # By default julia rounds to even when the value is exactly x.5,
            # but I want to round up to ensure the same behavior as wannier90.
            # See round(Int, 1.5) and round(Int, 2.5) in julia.
            n_v = round(Int, v_norm / dk, RoundNearestTiesUp)
            # ensure at least the two ending points are there, if the v_norm is too small
            n_v = max(n_v, 1)
            # Now operates with fractional coordinates, to avoid numerical
            # issues by doing frac -> cart -> frac transformation
            x_v = collect(range(0, 1, n_v + 1))
            # column vector * row vector = matrix
            kpt_v = (kseg.coords[k2] - kseg.coords[k1]) * x_v'
            kpt_v .+= kseg.coords[k1]

            if j == 1
                push!(indices_seg, 1)
                push!(labels_seg, k1)
            else
                # remove repeated points
                popfirst!(x_v)
                kpt_v = kpt_v[:, 2:end]
            end
            n_x_seg += length(x_v)
            push!(indices_seg, n_x_seg)
            push!(labels_seg, k2)

            append!(kpoints_seg, collect(eachcol(kpt_v)))
        end

        push!(kpaths, kpoints_seg)
        push!(indices, indices_seg)
        push!(labels, labels_seg)
    end

    points = reduce(vcat, kpaths)
    indices = reduce(vcat, indices)
    labels = reduce(vcat, labels)
    return KPath(kseg.recip_lattice, points, indices, labels)
end

function Base.isapprox(a::KPath, b::KPath; kwargs...)
    # indices and labels must match exactly
    if a.indices != b.indices || a.labels != b.labels
        return false
    end
    # reciprocal lattice and points compared approximately
    if !isapprox(a.recip_lattice, b.recip_lattice; kwargs...)
        return false
    end
    if length(a.points) != length(b.points)
        return false
    end
    return all(isapprox.(a.points, b.points; kwargs...))
end
