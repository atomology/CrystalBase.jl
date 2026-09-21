module CrystalBaseBrillouinExt

using CrystalBase

# Use import instead of using to avoid name conflicts with CrystalBase.
import Brillouin

"""
    KPath(kpi::Brillouin.KPathInterpolant)

Wrap the explicit k-points of a `KPathInterpolant` verbatim, one `Subpath`
per connected line, keeping its labels.
"""
function CrystalBase.KPath(kpi::Brillouin.KPathInterpolant)
    kpi_frac = Brillouin.latticize(kpi)
    subpaths = map(zip(kpi_frac.kpaths, kpi_frac.labels)) do (points, labels)
        full_labels = fill("", length(points))
        for (i, l) in labels
            full_labels[i] = string(l)
        end
        Subpath(points, full_labels; divisions = 1)
    end
    return KPath(reduce(hcat, kpi_frac.basis), subpaths)
end

"""
    KPath(kp::Brillouin.KPath)

Route of a `Brillouin.KPath`: its high-symmetry points as vertices, with one
division per segment. Call [`resample`](@ref) to choose the sampling.
"""
function CrystalBase.KPath(kp::Brillouin.KPath)
    kp_frac = Brillouin.latticize(kp)
    subpaths = map(kp_frac.paths) do path
        Subpath([kp_frac.points[l] for l in path], string.(path); divisions = 1)
    end
    return KPath(reduce(hcat, kp_frac.basis), subpaths)
end

end # module
