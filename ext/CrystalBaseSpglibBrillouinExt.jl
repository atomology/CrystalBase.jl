module CrystalBaseSpglibBrillouinExt

using CrystalBase

# Use import instead of using to avoid name conflicts with CrystalBase.
import Spglib
import Brillouin

"""
    KPath(crystal::Crystal; n_points_first_segment = 100)

Standard high-symmetry k-point path of `crystal` (any setting, standard or
not), from `Brillouin.irrfbz_path`. Divisions follow the wannier90 rule, see
[`resample`](@ref).

Requires `using Spglib` and `import Brillouin`; not `using Brillouin`, which
exports its own `KPath` and so makes a bare `KPath` ambiguous.
"""
function CrystalBase.KPath(crystal::Crystal; n_points_first_segment::Integer = 100)
    bkpath = Brillouin.irrfbz_path(Spglib.Cell(crystal))
    return resample(KPath(bkpath); n_points_first_segment)
end

end # module
