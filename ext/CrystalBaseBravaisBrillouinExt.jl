module CrystalBaseBravaisBrillouinExt

using CrystalBase

# Use import instead of using to avoid name conflicts with CrystalBase.
import Bravais
import Brillouin

"""
    Brillouin.KPathInterpolant(kpath::CrystalBase.KPath)

Sample `kpath` into a `KPathInterpolant`, one line per `Subpath`, in the
lattice (fractional) setting.
"""
function Brillouin.KPathInterpolant(kpath::CrystalBase.KPath)
    kpaths = Vector{Vector{Vec3{Float64}}}()
    labels = Vector{Dict{Int, Symbol}}()
    for sp in kpath.subpaths
        points, ticks = CrystalBase._sample(sp)
        push!(kpaths, Vector{Vec3{Float64}}(points))
        push!(labels, Dict(i => Symbol(l) for (i, l) in ticks))
    end
    basis = Bravais.ReciprocalBasis(lattice_vectors(kpath.recip_lattice))
    return Brillouin.KPathInterpolant(kpaths, labels, basis, Ref(Brillouin.LATTICE))
end

end # module
