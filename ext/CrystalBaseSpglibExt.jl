module CrystalBaseSpglibExt

using CrystalBase

# Use import instead of using to avoid name conflicts with CrystalBase.
import Spglib

"""
    Spglib.Cell(crystal::Crystal)

Convert a `Crystal` to a `Spglib.Cell`. Atoms are identified by atomic
number, so per-site labels (`Fe1`/`Fe2`) do not split the symmetry.
"""
function Spglib.Cell(crystal::Crystal)
    return Spglib.Cell(
        vec3(crystal.lattice), Vector.(crystal.atom_positions), crystal.atom_numbers
    )
end

"""
    Crystal(cell::Spglib.Cell)

Convert a `Spglib.Cell` to a `Crystal`; labels are the element symbols.
"""
function CrystalBase.Crystal(cell::Spglib.Cell)
    return Crystal(Matrix(cell.lattice), cell.positions, Int.(cell.atoms))
end

function CrystalBase.spacegroup(crystal::Crystal; symprec::Real = 1.0e-5)
    dataset = Spglib.get_dataset(Spglib.Cell(crystal), symprec)
    return (; symbol = dataset.international_symbol, number = Int(dataset.spacegroup_number))
end

function CrystalBase.coupled_axes(crystal::Crystal; symprec::Real = 1.0e-5)
    dataset = Spglib.get_dataset(Spglib.Cell(crystal), symprec)
    # union-find over the three axes
    parent = [1, 2, 3]
    find(a) = (
        while parent[a] != a
            a = parent[a]
        end; a
    )
    for rotation in dataset.rotations
        for i in 1:3, j in 1:3
            (i != j && rotation[i, j] != 0) || continue
            ri, rj = find(i), find(j)
            ri == rj || (parent[rj] = ri)
        end
    end
    return [find(i) for i in 1:3]
end

end # module
