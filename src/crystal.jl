export Crystal
export n_atoms, atom_symbols, unique_species, reduced_formula, atom_positions_cart
export spacegroup, kgrid_from_density

"""
    $(TYPEDEF)

A crystal structure: lattice plus atoms.

Element identity (`atom_numbers`) and per-site identity (`atom_labels`) are
stored separately. Labels default to the element symbol; a refined label such
as `Fe1`/`Fe2` distinguishes symmetry-inequivalent sites of the same element.
Every label must begin with the symbol of its element, see [`atomic_symbol`](@ref).

# Fields
$(FIELDS)
"""
struct Crystal{T <: Real}
    """unit cell, 3 × 3, each column is a lattice vector, in Å"""
    lattice::Mat3{T}

    """atomic positions, length-`n_atoms` vector of fractional coordinates"""
    atom_positions::Vector{Vec3{T}}

    """atomic numbers, length-`n_atoms`"""
    atom_numbers::Vector{Int}

    """per-site labels, length-`n_atoms`; each begins with the element symbol"""
    atom_labels::Vector{String}

    function Crystal{T}(
            lattice::Mat3{T},
            atom_positions::Vector{Vec3{T}},
            atom_numbers::Vector{Int},
            atom_labels::Vector{String},
        ) where {T <: Real}
        n = length(atom_positions)
        length(atom_numbers) == n ||
            throw(DimensionMismatch("atom_numbers has $(length(atom_numbers)) entries, expected $n"))
        length(atom_labels) == n ||
            throw(DimensionMismatch("atom_labels has $(length(atom_labels)) entries, expected $n"))
        for (z, label) in zip(atom_numbers, atom_labels)
            symbol = atomic_symbol(label)
            atomic_number(symbol) == z || throw(
                ArgumentError("label \"$label\" starts with $symbol, but the site is $(atomic_symbol(z))")
            )
        end
        return new{T}(lattice, atom_positions, atom_numbers, atom_labels)
    end
end

"""
    Crystal(lattice, atom_positions, atom_numbers, atom_labels)
    Crystal(lattice, atom_positions, atom_labels)
    Crystal(lattice, atom_positions, atom_numbers)
    Crystal(lattice, atoms)

Construct a [`Crystal`](@ref).

# Arguments
- `lattice`: anything [`mat3`](@ref) accepts; each column is a lattice vector in Å
- `atom_positions`: length-`n_atoms` vector of fractional coordinates
- `atom_numbers`: atomic numbers; derived from `atom_labels` when omitted
- `atom_labels`: per-site labels (strings or symbols); default to the element
  symbols when omitted
- `atoms`: length-`n_atoms` vector of `label => fractional position` pairs, as
  returned for the `atoms_frac` block of a wannier90 `win` file

# Examples
```jldoctest crystal; setup = :(using CrystalBase)
lattice = [0.0 2.715 2.715; 2.715 0.0 2.715; 2.715 2.715 0.0];
si = Crystal(lattice, [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]], ["Si", "Si"]);
n_atoms(si), formula(si), si.atom_numbers
# output
(2, "Si2", [14, 14])
```

```jldoctest crystal
fe = Crystal(lattice, ["Fe1" => [0.0, 0.0, 0.0], "Fe2" => [0.5, 0.5, 0.5]]);
atom_symbols(fe), unique_species(fe)
# output
(["Fe", "Fe"], ["Fe1", "Fe2"])
```
"""
function Crystal(
        lattice::AbstractMatrix,
        atom_positions::AbstractVector,
        atom_numbers::AbstractVector{<:Integer},
        atom_labels::AbstractVector,
    )
    latt = mat3(lattice)
    T = float(promote_type(eltype(latt), _position_eltype(atom_positions)))
    return Crystal{T}(
        Mat3{T}(latt),
        Vector{Vec3{T}}(vec3.(atom_positions)),
        Vector{Int}(atom_numbers),
        String.(atom_labels),
    )
end

function Crystal(
        lattice::AbstractMatrix,
        atom_positions::AbstractVector,
        atom_labels::AbstractVector{<:Union{AbstractString, Symbol}},
    )
    atom_numbers = Int[atomic_number(atomic_symbol(l)) for l in atom_labels]
    return Crystal(lattice, atom_positions, atom_numbers, atom_labels)
end

function Crystal(
        lattice::AbstractMatrix,
        atom_positions::AbstractVector,
        atom_numbers::AbstractVector{<:Integer},
    )
    return Crystal(lattice, atom_positions, atom_numbers, atomic_symbol(atom_numbers))
end

function Crystal(lattice::AbstractMatrix, atoms::AbstractVector{<:Pair})
    atom_labels = [string(a.first) for a in atoms]
    atom_positions = [a.second for a in atoms]
    return Crystal(lattice, atom_positions, atom_labels)
end

"""
    Crystal{T}(crystal::Crystal)

Convert the element type of `crystal`.
"""
function Crystal{T}(crystal::Crystal) where {T <: Real}
    return Crystal{T}(
        Mat3{T}(crystal.lattice),
        Vector{Vec3{T}}(crystal.atom_positions),
        crystal.atom_numbers,
        crystal.atom_labels,
    )
end
Crystal{T}(crystal::Crystal{T}) where {T <: Real} = crystal
Base.convert(::Type{Crystal{T}}, crystal::Crystal) where {T <: Real} = Crystal{T}(crystal)

_position_eltype(positions::AbstractVector) =
    isempty(positions) ? Float64 : promote_type(map(p -> eltype(p), positions)...)

"""
    $(SIGNATURES)

Number of atoms.
"""
n_atoms(crystal::Crystal) = length(crystal.atom_positions)

"""
    $(SIGNATURES)

Element symbols per site, derived from `atom_numbers`. Always the bare element
(`Fe`, never `Fe1`); see [`unique_species`](@ref) for labels.
"""
atom_symbols(crystal::Crystal) = atomic_symbol(crystal.atom_numbers)

"""
    $(SIGNATURES)

Distinct labels in first-appearance order, e.g. `["Fe1", "Fe2", "O"]`.
This is the `ATOMIC_SPECIES` set of a Quantum ESPRESSO input.
"""
function unique_species(crystal::Crystal)
    species = String[]
    for l in crystal.atom_labels
        l in species || push!(species, l)
    end
    return species
end

"""
    formula(crystal; order = :hill, reduced = false)

Chemical formula of the cell, see [`formula(::AbstractVector)`](@ref formula).
"""
formula(crystal::Crystal; kwargs...) = formula(atom_symbols(crystal); kwargs...)

"""
    $(SIGNATURES)

[`formula`](@ref) divided by the greatest common divisor of its counts, so a
primitive cell and its supercells share one `reduced_formula`.
"""
reduced_formula(crystal::Crystal; kwargs...) = formula(crystal; reduced = true, kwargs...)

real_lattice(crystal::Crystal) = crystal.lattice
reciprocal_lattice(crystal::Crystal) = reciprocal_lattice(crystal.lattice)

"""
    $(SIGNATURES)

Atomic positions in Cartesian coordinates (Å).
"""
atom_positions_cart(crystal::Crystal) = frac_to_cart(crystal.lattice, crystal.atom_positions)

function Base.:(==)(a::Crystal, b::Crystal)
    return a.lattice == b.lattice && a.atom_positions == b.atom_positions &&
        a.atom_numbers == b.atom_numbers && a.atom_labels == b.atom_labels
end

function Base.isapprox(a::Crystal, b::Crystal; kwargs...)
    a.atom_numbers == b.atom_numbers || return false
    a.atom_labels == b.atom_labels || return false
    isapprox(a.lattice, b.lattice; kwargs...) || return false
    return all(isapprox.(a.atom_positions, b.atom_positions; kwargs...))
end

function Base.show(io::IO, crystal::Crystal)
    return print(io, "Crystal(", formula(crystal), ", ", n_atoms(crystal), " atoms)")
end

function Base.show(io::IO, ::MIME"text/plain", crystal::Crystal{T}) where {T}
    println(io, "Crystal{$T}: ", formula(crystal), ", ", n_atoms(crystal), " atoms")
    show_lattice(io, crystal.lattice; indent = "  ")
    print(io, "  atoms (fractional):")
    width = maximum(length, crystal.atom_labels; init = 0)
    n = n_atoms(crystal)
    n_shown = min(n, get(io, :limit, false) ? 20 : n)
    # align the dots over all shown positions at once, hence the flat format pass
    fields = reshape(_fmt_aligned(reduce(vcat, @view(crystal.atom_positions[1:n_shown]); init = T[])), 3, :)
    for i in 1:n_shown
        label = rpad(crystal.atom_labels[i], width)
        print(io, "\n", rstrip(string("    ", label, "  ", join(view(fields, :, i), "  "))))
    end
    n_shown < n && print(io, "\n    ⋮ (", n - n_shown, " more)")
    return
end

"""
    spacegroup(crystal; symprec = 1e-5)

Space group of `crystal` as `(; symbol, number)`, e.g. `(symbol = "Fm-3m", number = 225)`.

Requires `using Spglib`. `symprec` is the spglib distance tolerance in Å.
"""
function spacegroup end

"""
    coupled_axes(crystal; symprec = 1e-5)

Group the three lattice axes into classes coupled by the crystal point group:
returns a length-3 vector of class ids. Two axes share a class when some
rotation (integer matrix in the lattice basis) has a nonzero off-diagonal
entry mixing them, so a Monkhorst--Pack grid must give them equal
subdivisions to stay symmetric. Requires `using Spglib`.
"""
function coupled_axes end

"""
    kgrid_from_density(recip_lattice, density)
    kgrid_from_density(crystal, density; symmetrize = true, symprec = 1e-5)

Monkhorst--Pack grid `(n1, n2, n3)` sized from the reciprocal lattice.

`density` is k-points per unit reciprocal length (Å); each axis gets
`ceil(|bᵢ| · density)`. Larger `density` gives a denser grid.

With a [`Crystal`](@ref) and `symmetrize = true` (requires `using Spglib`),
axes coupled by the point group (see [`coupled_axes`](@ref)) are raised to
their common maximum, so centered lattices get a mesh the point group leaves
invariant. Simple tetragonal, orthorhombic and hexagonal cells are unchanged.

# Examples
```jldoctest kgrid_from_density; setup = :(using CrystalBase)
lattice = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 8.0];
kgrid_from_density(reciprocal_lattice(lattice), 5.0)
# output
(8, 8, 4)
```
"""
function kgrid_from_density(recip_lattice::AbstractMatrix, density::Real)
    b = mat3(recip_lattice)
    return ntuple(i -> ceil(Int, norm(b[:, i]) * density), 3)
end

function kgrid_from_density(
        crystal::Crystal, density::Real; symmetrize::Bool = true, symprec::Real = 1.0e-5
    )
    grid = kgrid_from_density(reciprocal_lattice(crystal), density)
    symmetrize || return grid
    hasmethod(coupled_axes, Tuple{Crystal}) ||
        error("kgrid_from_density with symmetrize = true requires `using Spglib`")
    classes = coupled_axes(crystal; symprec)
    return ntuple(i -> maximum(grid[j] for j in 1:3 if classes[j] == classes[i]), 3)
end
