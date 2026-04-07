module CrystalBaseBravaisBrillouinExt

using CrystalBase

# Use import instead of using to avoid name conflicts with CrystalBase.
import Bravais
import Brillouin

function Brillouin.KPathInterpolant(kpath::CrystalBase.KPath)
    # kpoints along path
    bkpaths = Vector{Vector{Vec3{Float64}}}()
    # symmetry points
    blabels = Vector{Dict{Int, Symbol}}()

    i0 = kpath.indices[1]  # 1st point
    lab = Dict{Int, Symbol}()  # label of each line
    push!(lab, i0 => Symbol(kpath.labels[1]))
    for (i, l) in zip(kpath.indices[2:end], kpath.labels[2:end])
        if i == i0 + 1
            push!(blabels, lab)
            lab = Dict{Int, Symbol}()
        end
        push!(lab, i => Symbol(l))
        i0 = i
    end
    push!(blabels, lab)

    for lab in blabels
        kp = Vector{Vec3{Float64}}()  # kpath of each line
        ik1 = minimum(keys(lab))
        ik2 = maximum(keys(lab))
        append!(kp, [v for v in kpath.points[ik1:ik2]])
        push!(bkpaths, kp)
    end

    for (i, lab) in enumerate(blabels)
        ik1 = minimum(keys(lab)) - 1
        blabels[i] = Dict(((k - ik1) => v) for (k, v) in lab)
    end

    basis = Bravais.ReciprocalBasis(collect(eachcol(kpath.recip_lattice)))
    setting = Ref(Brillouin.LATTICE)
    kpi = Brillouin.KPathInterpolant(bkpaths, blabels, basis, setting)
    return kpi
end

end # module
