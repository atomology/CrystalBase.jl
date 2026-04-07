module CrystalBaseBrillouinExt

using LinearAlgebra
using OrderedCollections
using CrystalBase

# Use import instead of using to avoid name conflicts with CrystalBase.
import Brillouin

"""
Return the symmetry indices and labels.
"""
function symm_point_indices_labels(kpi::Brillouin.KPathInterpolant)
    kpi_frac = Brillouin.latticize(kpi)

    symm_point_indices = Vector{Int}()
    symm_point_labels = Vector{String}()
    ik0 = 0
    for lab in kpi_frac.labels
        for (ik, l) in lab
            push!(symm_point_indices, ik + ik0)
            push!(symm_point_labels, String(l))
        end
        ik0 += maximum(keys(lab))
    end

    perm = sortperm(symm_point_indices)
    symm_point_indices = symm_point_indices[perm]
    symm_point_labels = symm_point_labels[perm]

    return symm_point_indices, symm_point_labels
end

function CrystalBase.KPath(kpi::Brillouin.KPathInterpolant)
    kpi_frac = Brillouin.latticize(kpi)
    indices, labels = symm_point_indices_labels(kpi_frac)
    recip_lattice = reduce(hcat, kpi_frac.basis)
    return CrystalBase.KPath(recip_lattice, collect(kpi_frac), indices, labels)
end

function CrystalBase.KSegment(kpi::Brillouin.KPathInterpolant)
    kpi_frac = Brillouin.latticize(kpi)

    segments = Vector{Vector{String}}()
    coords = OrderedDict{String, Vec3{Float64}}()

    for (label, path) in zip(kpi_frac.labels, kpi_frac.kpaths)
        ik_labels = sort(pairs(label); by = first)
        push!(segments, collect(values(ik_labels)))
        for (ik, lab) in pairs(ik_labels)
            coords[string(lab)] = path[ik]
        end
    end

    return CrystalBase.KSegment(kpi.basis, segments, coords)
end

"""
The `Brillouin.KPath` is actually a container for the kpath segments,
while the `Brillouin.KPathInterpolant` is the one that contains the explicit list
of kpoints coordinates along the path.
"""
function CrystalBase.KSegment(kp::Brillouin.KPath)
    kp_frac = Brillouin.latticize(kp)
    recip_lattice = reduce(hcat, kp_frac.basis)
    segments = [string.(seg) for seg in kp_frac.paths]
    coords = OrderedDict(string(l) => vec3(k) for (l, k) in pairs(kp_frac.points))
    return CrystalBase.KSegment(recip_lattice, segments, coords)
end

end # module
