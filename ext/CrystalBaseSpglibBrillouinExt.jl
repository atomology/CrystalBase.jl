module CrystalBaseSpglibBrillouinExt

using CrystalBase

# Use import instead of using to avoid name conflicts with CrystalBase.
import Spglib
import Brillouin

function CrystalBase.KSegment(
        lattice::AbstractMatrix, atom_positions::AbstractVector, atom_numbers::AbstractVector{<:Integer}
    )
    vecs = collect(eachcol(lattice))
    cell = Spglib.Cell(vecs, Vector.(atom_positions), atom_numbers)
    bkpath = Brillouin.irrfbz_path(cell)
    kseg = CrystalBase.KSegment(bkpath)
    return kseg
end

end # module
