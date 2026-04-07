@testitem "group_nearby_indices" begin
    using CrystalBase: group_nearby_indices
    indices = [1, 2, 4, 5, 6]
    grps = group_nearby_indices(indices)
    @test grps == [[1, 2], [4, 5, 6]]
end

@testitem "KSegment" begin
    using OrderedCollections: OrderedDict

    lattice = [
        -2.6988 0.0 -2.6988
        0.0 2.6988 2.6988
        2.6988 2.6988 0.0
    ]
    recip_lattice = reciprocal_lattice(lattice)
    kpoint_path = [
        ["L" => [0.5, 0.5, 0.5], "G" => [0.0, 0.0, 0.0]],
        ["G" => [0.0, 0.0, 0.0], "X" => [0.5, 0.0, 0.5]],
        ["X" => [0.5, -0.5, 0.0], "K" => [0.375, -0.375, 0.0]],
        ["K" => [0.375, -0.375, 0.0], "G" => [0.0, 0.0, 0.0]],
    ]
    kseg = KSegment(recip_lattice, kpoint_path)

    ref_kseg = (;
        recip_lattice = reduce(
            hcat, [
                [-1.164069, -1.164069, 1.164069],
                [1.164069, 1.164069, 1.164069],
                [-1.164069, 1.164069, -1.164069],
            ]
        ),
        segments = [["L", "G", "X"], ["X_1", "K", "G"]],
        coords = OrderedDict(
            "X_1" => [0.5, -0.5, 0.0],
            "G" => [0.0, 0.0, 0.0],
            "K" => [0.375, -0.375, 0.0],
            "L" => [0.5, 0.5, 0.5],
            "X" => [0.5, 0.0, 0.5],
        ),
    )
    @test isapprox(ref_kseg.recip_lattice, kseg.recip_lattice; atol = 1.0e-5)
    @test ref_kseg.segments == kseg.segments
    @test ref_kseg.coords == kseg.coords
end

@testitem "KSegment from structure" begin
    using Spglib, Brillouin
    lattice = [
        0.0       2.715265       2.715265
        2.715265       0.0       2.715265
        2.715265       2.715265       0.0
    ]
    atom_positions = [
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
    ]
    atom_symbols = ["Si", "Si"]
    kseg = KSegment(lattice, atom_positions, atom_symbols)

    @test kseg.recip_lattice ≈ reciprocal_lattice(lattice)
    ref_segments = [
        ["Γ", "X", "U"],
        ["K", "Γ", "L", "W", "X"],
    ]
    @test ref_segments == kseg.segments
    # Compare without order
    @test Dict(kseg.coords) == Dict(
        "K" => [0.375, 0.375, 0.75],
        "L" => [0.5, 0.5, 0.5],
        "U" => [0.625, 0.25, 0.625],
        "W" => [0.5, 0.25, 0.75],
        "X" => [0.5, 0.0, 0.5],
        "Γ" => [0.0, 0.0, 0.0],
    )
end

@testmodule KPathEnv begin
    using CrystalBase
    # From WannierDatasets/datasets/Si2/Si2.win, with `bands_num_points = 5`
    lattice = [
        0.0       2.715265       2.715265
        2.715265       0.0       2.715265
        2.715265       2.715265       0.0
    ]
    recip_lattice = reciprocal_lattice(lattice)
    kpoint_path = [
        ["G" => [0.0, 0.0, 0.0], "X" => [0.5, 0.0, 0.5]],
        ["X" => [0.5, 0.0, 0.5], "U" => [0.625, 0.25, 0.625]],
        ["U" => [0.625, 0.25, 0.625], "K" => [0.375, 0.375, 0.75]],
        ["K" => [0.375, 0.375, 0.75], "G" => [0.0, 0.0, 0.0]],
        ["G" => [0.0, 0.0, 0.0], "L" => [0.5, 0.5, 0.5]],
        ["L" => [0.5, 0.5, 0.5], "W" => [0.5, 0.25, 0.75]],
        ["W" => [0.5, 0.25, 0.75], "X" => [0.5, 0.0, 0.5]],
    ]

    # From `Si2_band.kpt`
    kpoints = [
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.1],
        [0.2, 0.0, 0.2],
        [0.3, 0.0, 0.3],
        [0.4, 0.0, 0.4],
        [0.5, 0.0, 0.5],
        [0.5625, 0.125, 0.5625],
        [0.625, 0.25, 0.625],
        [0.541667, 0.291667, 0.666667],
        [0.458333, 0.333333, 0.708333],
        [0.375, 0.375, 0.75],
        [0.3, 0.3, 0.6],
        [0.225, 0.225, 0.45],
        [0.15, 0.15, 0.3],
        [0.075, 0.075, 0.15],
        [0.0, 0.0, 0.0],
        [0.125, 0.125, 0.125],
        [0.25, 0.25, 0.25],
        [0.375, 0.375, 0.375],
        [0.5, 0.5, 0.5],
        [0.5, 0.4375, 0.5625],
        [0.5, 0.375, 0.625],
        [0.5, 0.3125, 0.6875],
        [0.5, 0.25, 0.75],
        [0.5, 0.166667, 0.666667],
        [0.5, 0.083333, 0.583333],
        [0.5, 0.0, 0.5],
    ]
    # From `Si2_band.labelinfo.dat`
    labels = [
        "G",
        "X",
        "U",
        "K",
        "G",
        "L",
        "W",
        "X",
    ]
    indices = [
        1,
        6,
        8,
        11,
        16,
        20,
        24,
        27,
    ]
    # From the 1st column of `Si2_band.dat`
    x = [
        0.0e+0,
        0.23140229e+0,
        0.46280457e+0,
        0.69420686e+0,
        0.92560915e+0,
        0.11570114e+1,
        0.13615441e+1,
        0.15660768e+1,
        0.18022507e+1,
        0.20384247e+1,
        0.22745987e+1,
        0.25200379e+1,
        0.2765477e+1,
        0.30109162e+1,
        0.32563554e+1,
        0.35017946e+1,
        0.37522949e+1,
        0.40027953e+1,
        0.42532956e+1,
        0.45037959e+1,
        0.47083286e+1,
        0.49128612e+1,
        0.51173939e+1,
        0.53219265e+1,
        0.55147618e+1,
        0.5707597e+1,
        0.59004323e+1,
    ]
end

@testitem "KPath" setup = [KPathEnv] begin
    kseg = KSegment(KPathEnv.recip_lattice, KPathEnv.kpoint_path)
    kp = KPath(kseg, 5)

    @test all(isapprox.(kp.points, KPathEnv.kpoints; atol = 1.0e-5))
    @test kp.labels == KPathEnv.labels
    @test kp.indices == KPathEnv.indices
end

@testitem "linear_path" setup = [KPathEnv] begin
    kseg = KSegment(KPathEnv.recip_lattice, KPathEnv.kpoint_path)
    kp = KPath(kseg, 5)
    x = linear_path(kp)[1]
    @test all(isapprox.(KPathEnv.x, x; atol = 1.0e-5))
end

@testitem "KPathInterpolant from KPath" setup = [KPathEnv] begin
    import Bravais, Brillouin

    kseg = KSegment(KPathEnv.recip_lattice, KPathEnv.kpoint_path)
    kp = KPath(kseg, 5)
    kpi = Brillouin.KPathInterpolant(kp)

    @test all(isapprox.(kpi, kp.points; atol = 1.0e-5))
    @test mat3(kpi.basis) == kp.recip_lattice
    @test kpi.labels[1] == Dict(i => Symbol(l) for (i, l) in zip(kp.indices, kp.labels))
end

@testitem "isapprox for KSegment and KPath" begin
    # simple path
    lattice = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    kpoint_path = [
        ["G" => [0.0, 0.0, 0.0], "X" => [0.5, 0.0, 0.0]],
        ["X" => [0.5, 0.0, 0.0], "L" => [1.0, 0.0, 0.0]],
    ]
    ks1 = KSegment(lattice, kpoint_path)
    # small perturbation in coordinates
    kpoint_path2 = [
        ["G" => [1.0e-7, 0.0, 0.0], "X" => [0.5 + 1.0e-7, 0.0, 0.0]],
        ["X" => [0.5 + 1.0e-7, 0.0, 0.0], "L" => [1.0 + 1.0e-7, 0.0, 0.0]],
    ]
    ks2 = KSegment(lattice, kpoint_path2)
    @test isapprox(ks1, ks2; atol = 1e-6)
    @test !isapprox(ks1, ks2)
    # different labels -> not approx
    kpoint_path3 = [
        ["G" => [0.0, 0.0, 0.0], "Y" => [0.5, 0.0, 0.0]],
    ]
    ks3 = KSegment(lattice, kpoint_path3)
    @test !isapprox(ks1, ks3; atol = 1e-6)

    # KPath via KSegment
    kp1 = KPath(ks1, 10)
    kp2 = KPath(ks2, 10)
    @test isapprox(kp1, kp2; atol = 1e-6)
    @test !isapprox(kp1, kp2)
end
