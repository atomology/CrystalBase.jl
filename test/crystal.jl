@testmodule CrystalEnv begin
    using CrystalBase
    lattice = [
        0.0       2.715265       2.715265
        2.715265       0.0       2.715265
        2.715265       2.715265       0.0
    ]
    positions = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]
    si = Crystal(lattice, positions, ["Si", "Si"])
end

@testitem "Crystal constructors" setup = [CrystalEnv] begin
    si = CrystalEnv.si
    @test si isa Crystal{Float64}
    @test si.lattice == CrystalEnv.lattice
    @test si.atom_positions == CrystalEnv.positions
    @test si.atom_numbers == [14, 14]
    @test si.atom_labels == ["Si", "Si"]

    # from atomic numbers: labels default to symbols
    @test Crystal(CrystalEnv.lattice, CrystalEnv.positions, [14, 14]) == si
    # from label => position pairs (wannier90 atoms_frac)
    @test Crystal(CrystalEnv.lattice, ["Si" => CrystalEnv.positions[1], :Si => CrystalEnv.positions[2]]) == si
    # integer positions and lattice are promoted to Float64
    @test Crystal([1 0 0; 0 1 0; 0 0 1], [[0, 0, 0]], ["H"]) isa Crystal{Float64}

    # refined labels
    fe = Crystal(CrystalEnv.lattice, ["Fe1" => [0.0, 0.0, 0.0], "Fe2" => [0.5, 0.5, 0.5], "O" => [0.25, 0.25, 0.25]])
    @test fe.atom_numbers == [26, 26, 8]
    @test atom_symbols(fe) == ["Fe", "Fe", "O"]
    @test unique_species(fe) == ["Fe1", "Fe2", "O"]

    # dummy sites of a synthetic model
    toy = Crystal(CrystalEnv.lattice, CrystalEnv.positions, ["X1", "X2"])
    @test toy.atom_numbers == [0, 0]
    @test atom_symbols(toy) == ["X", "X"]
    @test formula(toy) == "X2"

    @test_throws DimensionMismatch Crystal(CrystalEnv.lattice, CrystalEnv.positions, ["Si"])
    @test_throws DimensionMismatch Crystal(CrystalEnv.lattice, CrystalEnv.positions, [14], ["Si", "Si"])
    @test_throws ArgumentError Crystal(CrystalEnv.lattice, CrystalEnv.positions, [14, 14], ["Si", "O"])
    @test_throws ArgumentError Crystal(CrystalEnv.lattice, CrystalEnv.positions, ["Si", "Zz"])
end

@testitem "Crystal accessors" setup = [CrystalEnv] begin
    si = CrystalEnv.si
    @test n_atoms(si) == 2
    @test atom_symbols(si) == ["Si", "Si"]
    @test unique_species(si) == ["Si"]
    @test formula(si) == "Si2"
    @test reduced_formula(si) == "Si"
    @test real_lattice(si) == si.lattice
    @test reciprocal_lattice(si) ≈ reciprocal_lattice(CrystalEnv.lattice)
    @test atom_positions_cart(si) ≈ frac_to_cart(CrystalEnv.lattice, CrystalEnv.positions)

    shifted = Crystal(CrystalEnv.lattice, [p .+ 1.0e-7 for p in CrystalEnv.positions], ["Si", "Si"])
    @test isapprox(si, shifted; atol = 1.0e-6)
    @test !isapprox(si, shifted)
    @test si != shifted
    @test si == Crystal(CrystalEnv.lattice, CrystalEnv.positions, ["Si", "Si"])
end

@testitem "show Crystal" setup = [CrystalEnv] begin
    s = sprint(show, MIME("text/plain"), CrystalEnv.si)
    @test occursin("Crystal{Float64}: Si2, 2 atoms", s)
    @test occursin("Si          0.25        0.25        0.25", s)
    @test sprint(show, CrystalEnv.si) == "Crystal(Si2, 2 atoms)"
end

@testitem "kgrid_from_density" setup = [CrystalEnv] begin
    lattice = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 8.0]
    @test kgrid_from_density(reciprocal_lattice(lattice), 5.0) == (8, 8, 4)
    @test kgrid_from_density(CrystalEnv.si, 3.0; symmetrize = false) == (7, 7, 7)
end

@testitem "Spglib extension" setup = [CrystalEnv] begin
    using Spglib
    si = CrystalEnv.si
    cell = Spglib.Cell(si)
    @test cell.atoms == [14, 14]
    @test Crystal(cell) ≈ si

    @test spacegroup(si) == (symbol = "Fd-3m", number = 227)
    @test CrystalBase.coupled_axes(si) == [1, 1, 1]
    @test kgrid_from_density(si, 3.0) == (7, 7, 7)

    # body-centered tetragonal: the 4-fold rotation couples all three primitive
    # axes, so the anisotropic raw grid is raised to an isotropic one
    bct = Crystal([-2.0 2.0 2.0; 2.0 -2.0 2.0; 3.0 3.0 -3.0], [[0.0, 0.0, 0.0]], ["Fe"])
    @test spacegroup(bct).number == 139
    @test kgrid_from_density(bct, 3.0; symmetrize = false) == (6, 6, 7)
    @test kgrid_from_density(bct, 3.0) == (7, 7, 7)

    # simple tetragonal: only a and b are coupled
    st = Crystal([4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 8.0], [[0.0, 0.0, 0.0]], ["Fe"])
    @test CrystalBase.coupled_axes(st) == [1, 1, 3]
    @test kgrid_from_density(st, 5.0) == (8, 8, 4)

    # per-site labels do not split the symmetry
    fe = Crystal(CrystalEnv.lattice, ["Fe1" => [0.0, 0.0, 0.0], "Fe2" => [0.5, 0.5, 0.5]])
    @test spacegroup(fe) == spacegroup(Crystal(CrystalEnv.lattice, [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ["Fe", "Fe"]))
end
