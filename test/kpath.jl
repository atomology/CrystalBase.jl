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
    labels = ["G", "X", "U", "K", "G", "L", "W", "X"]
    indices = [1, 6, 8, 11, 16, 20, 24, 27]
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

@testitem "Subpath" begin
    sp = Subpath([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0]], ["G", "X", "M"]; divisions = [4, 2])
    @test sp.divisions == [4, 2]
    @test CrystalBase.n_kpoints(sp) == 7
    @test length(Subpath([[0.0, 0.0, 0.0]]).divisions) == 0
    # integer divisions broadcast
    @test Subpath([[0, 0, 0], [1, 0, 0], [1, 1, 0]]; divisions = 3).divisions == [3, 3]
    @test eltype(Subpath([[0, 0, 0], [1, 0, 0]]).vertices) == Vec3{Float64}
    @test_throws ArgumentError Subpath(Vector{Vector{Float64}}())
    @test_throws DimensionMismatch Subpath([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], ["G"])
    @test_throws DimensionMismatch Subpath([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]; divisions = [1, 1])
    @test_throws ArgumentError Subpath([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]; divisions = 0)
end

@testitem "KPath from wannier90 kpoint_path" setup = [KPathEnv] begin
    kp = KPath(KPathEnv.recip_lattice, KPathEnv.kpoint_path; n_points_first_segment = 5)

    @test reciprocal_lattice(kp) ≈ KPathEnv.recip_lattice
    @test real_lattice(kp) ≈ KPathEnv.lattice
    # all segments share corners, so one connected subpath
    @test length(kp.subpaths) == 1
    @test kp.subpaths[1].labels == ["G", "X", "U", "K", "G", "L", "W", "X"]
    @test kp.subpaths[1].divisions == [5, 2, 3, 5, 4, 4, 3]
    @test n_kpoints(kp) == 27
    @test all(isapprox.(kpoints(kp), KPathEnv.kpoints; atol = 1.0e-5))
    @test kpoints_cart(kp) ≈ frac_to_cart(KPathEnv.recip_lattice, KPathEnv.kpoints) atol = 1.0e-5
    @test ticks(kp) == (; indices = KPathEnv.indices, labels = KPathEnv.labels)
    @test all(isapprox.(cumulative_distances(kp), KPathEnv.x; atol = 1.0e-5))

    # The default 100 points/segment should return 511 kpoints as in
    # `WannierDatasets/datasets/Si2/outputs/MDRS/Si2_band.kpt`
    @test n_kpoints(KPath(KPathEnv.recip_lattice, KPathEnv.kpoint_path)) == 511

    # Symbol labels are accepted
    kp_sym = KPath(KPathEnv.recip_lattice, [[:G => [0.0, 0.0, 0.0], :X => [0.5, 0.0, 0.5]]])
    @test kp_sym.subpaths[1].labels == ["G", "X"]
end

@testitem "KPath disconnected segments" begin
    using LinearAlgebra
    recip_lattice = Matrix{Float64}(I, 3, 3)
    # Same label, different coordinate: a break with both labels kept
    kpoint_path = [
        ["L" => [0.5, 0.5, 0.5], "G" => [0.0, 0.0, 0.0]],
        ["G" => [0.0, 0.0, 0.0], "X" => [0.5, 0.0, 0.5]],
        ["X" => [0.5, -0.5, 0.0], "K" => [0.375, -0.375, 0.0]],
        ["K" => [0.375, -0.375, 0.0], "G" => [0.0, 0.0, 0.0]],
    ]
    kp = KPath(recip_lattice, kpoint_path; n_points_first_segment = 4)
    @test length(kp.subpaths) == 2
    @test kp.subpaths[1].labels == ["L", "G", "X"]
    @test kp.subpaths[2].labels == ["X", "K", "G"]
    @test ticks(kp).labels == ["L", "G", "X|X", "K", "G"]

    # Two disconnected segments: A->B and C->D, distance is flat across the break
    kpoint_path = [
        ["A" => [0.0, 0.0, 0.0], "B" => [1.0, 0.0, 0.0]],
        ["C" => [0.0, 0.0, 0.0], "D" => [0.0, 2.0, 0.0]],
    ]
    kp = KPath(recip_lattice, kpoint_path; n_points_first_segment = 5)
    @test [sp.divisions for sp in kp.subpaths] == [[5], [10]]
    @test n_kpoints(kp) == 17
    @test ticks(kp) == (; indices = [1, 6, 17], labels = ["A", "B|C", "D"])
    unmerged = ticks(kp; merge = false)
    @test unmerged == (; indices = [1, 6, 7, 17], labels = ["A", "B", "C", "D"])
    x = cumulative_distances(kp)
    @test x[unmerged.indices] ≈ [0.0, 1.0, 1.0, 3.0]
    @test x[6] ≈ 1.0
    @test x[7] ≈ 1.0
    @test x[end] ≈ 3.0
end

@testitem "KPath verbatim from kpoints" setup = [KPathEnv] begin
    # wannier90 band.kpt + labelinfo layout
    kp = KPath(KPathEnv.recip_lattice, KPathEnv.kpoints, KPathEnv.indices, KPathEnv.labels)
    @test length(kp.subpaths) == 1
    @test all(sp -> all(==(1), sp.divisions), kp.subpaths)
    @test kpoints(kp) == KPathEnv.kpoints
    @test ticks(kp) == (; indices = KPathEnv.indices, labels = KPathEnv.labels)
    @test all(isapprox.(cumulative_distances(kp), KPathEnv.x; atol = 1.0e-5))

    # per-point labels keyword gives the same path
    full_labels = fill("", length(KPathEnv.kpoints))
    full_labels[KPathEnv.indices] .= KPathEnv.labels
    @test KPath(KPathEnv.recip_lattice, KPathEnv.kpoints; labels = full_labels) == kp

    # Consecutive labeled indices mark a break (wannier90 convention)
    route = KPath(
        KPathEnv.recip_lattice,
        [["A" => [0.0, 0.0, 0.0], "B" => [0.5, 0.0, 0.0]], ["C" => [0.0, 0.0, 0.0], "D" => [0.0, 0.5, 0.0]]];
        n_points_first_segment = 4,
    )
    dense = KPath(KPathEnv.recip_lattice, kpoints(route), [1, 5, 6, 10], ["A", "B", "C", "D"])
    @test length(dense.subpaths) == 2
    @test ticks(dense).labels == ["A", "B|C", "D"]
    @test cumulative_distances(dense) ≈ cumulative_distances(route)

    # Geometric break detection without labels
    geometric = KPath(KPathEnv.recip_lattice, kpoints(route))
    @test length(geometric.subpaths) == 2
    @test all(isempty, ticks(geometric).labels)
    @test length(KPath(KPathEnv.recip_lattice, kpoints(route); break_tol = 0).subpaths) == 1

    @test_throws DimensionMismatch KPath(KPathEnv.recip_lattice, KPathEnv.kpoints; labels = ["G"])
    @test_throws ArgumentError KPath(KPathEnv.recip_lattice, KPathEnv.kpoints, [100], ["G"])
end

@testitem "resample" setup = [KPathEnv] begin
    kp = KPath(KPathEnv.recip_lattice, KPathEnv.kpoint_path; n_points_first_segment = 5)
    kp100 = resample(kp; n_points_first_segment = 100)
    @test kp100 == KPath(KPathEnv.recip_lattice, KPathEnv.kpoint_path)
    @test kp100.subpaths[1].vertices == kp.subpaths[1].vertices
    @test kp100.subpaths[1].labels == kp.subpaths[1].labels

    dense = resample(kp; density = 10.0)
    lengths = CrystalBase._segment_lengths(KPathEnv.recip_lattice, kp.subpaths[1])
    @test dense.subpaths[1].divisions == max.(round.(Int, lengths .* 10.0, RoundNearestTiesUp), 1)

    @test_throws ArgumentError resample(kp)
    @test_throws ArgumentError resample(kp; density = 1.0, n_points_first_segment = 1)
end

@testitem "unicode_kpoint_labels" begin
    using LinearAlgebra
    @test CrystalBase.unicode_kpoint_labels(["GAMMA", "DELTA_0", "LAMBDA_1", "SIGMA_2", "X", ""]) ==
        ["Γ", "Δ₀", "Λ₁", "Σ₂", "X", ""]
    kp = KPath(Matrix{Float64}(I, 3, 3), [["GAMMA" => [0.0, 0.0, 0.0], "X_1" => [0.5, 0.0, 0.0]]])
    kp_unicode = unicode_kpoint_labels(kp)
    @test ticks(kp_unicode).labels == ["Γ", "X₁"]
    @test ticks(kp).labels == ["GAMMA", "X_1"]
    # multi-digit subscripts; other suffixes are kept verbatim
    @test unicode_kpoint_labels(["X_10", "GAMMA_12", "X_a", "X_"]) == ["X₁₀", "Γ₁₂", "X_a", "X_"]
    @test CrystalBase._subscript(3) == "₃"
    @test_throws ArgumentError CrystalBase._subscript("1a")
end

@testitem "isapprox for KPath" begin
    using LinearAlgebra
    lattice = Matrix{Float64}(I, 3, 3)
    kpoint_path = [
        ["G" => [0.0, 0.0, 0.0], "X" => [0.5, 0.0, 0.0]],
        ["X" => [0.5, 0.0, 0.0], "L" => [1.0, 0.0, 0.0]],
    ]
    kp1 = KPath(lattice, kpoint_path; n_points_first_segment = 10)
    kpoint_path2 = [
        ["G" => [1.0e-7, 0.0, 0.0], "X" => [0.5 + 1.0e-7, 0.0, 0.0]],
        ["X" => [0.5 + 1.0e-7, 0.0, 0.0], "L" => [1.0 + 1.0e-7, 0.0, 0.0]],
    ]
    kp2 = KPath(lattice, kpoint_path2; n_points_first_segment = 10)
    @test isapprox(kp1, kp2; atol = 1.0e-6)
    @test !isapprox(kp1, kp2)
    @test kp1 == KPath(lattice, kpoint_path; n_points_first_segment = 10)
    # different labels -> not approx
    kp3 = KPath(lattice, [["G" => [0.0, 0.0, 0.0], "Y" => [0.5, 0.0, 0.0]]]; n_points_first_segment = 10)
    @test !isapprox(kp1, kp3; atol = 1.0e-6)
end

@testitem "show KPath" setup = [KPathEnv] begin
    kp = KPath(KPathEnv.recip_lattice, KPathEnv.kpoint_path; n_points_first_segment = 5)
    s = sprint(show, MIME("text/plain"), kp)
    @test occursin("1 subpaths, 27 kpoints", s)
    @test occursin("recip_lattice (\u00c5\u207b\u00b9):", s)
    @test occursin("G—X—U—K—G—L—W—X  divisions 5 2 3 5 4 4 3", s)
    @test sprint(show, kp) == "KPath(1 subpaths, 27 kpoints)"

    dense = KPath(KPathEnv.recip_lattice, KPathEnv.kpoints, KPathEnv.indices, KPathEnv.labels)
    s = sprint(show, MIME("text/plain"), dense)
    @test occursin("27 kpoints  ticks G@1 X@6 U@8", s)

    unlabeled = KPath(KPathEnv.recip_lattice, KPathEnv.kpoints; break_tol = 0)
    @test occursin("27 kpoints, no ticks", sprint(show, MIME("text/plain"), unlabeled))

    # labels are displayed as stored; convert explicitly for unicode
    gamma = KPath(KPathEnv.recip_lattice, [["GAMMA" => [0.0, 0.0, 0.0], "X" => [0.5, 0.0, 0.5]]])
    @test occursin("GAMMA—X", sprint(show, MIME("text/plain"), gamma))
    @test occursin("Γ—X", sprint(show, MIME("text/plain"), unicode_kpoint_labels(gamma)))
end

@testitem "KPath from Crystal" begin
    using Spglib
    import Brillouin
    lattice = [
        0.0       2.715265       2.715265
        2.715265       0.0       2.715265
        2.715265       2.715265       0.0
    ]
    crystal = Crystal(lattice, [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]], ["Si", "Si"])
    kp = KPath(crystal)

    @test reciprocal_lattice(kp) ≈ reciprocal_lattice(lattice)
    @test real_lattice(kp) ≈ lattice
    @test [sp.labels for sp in kp.subpaths] == [["Γ", "X", "U"], ["K", "Γ", "L", "W", "X"]]
    @test kp.subpaths[1].divisions[1] == 100
    coords = Dict(l => v for sp in kp.subpaths for (l, v) in zip(sp.labels, sp.vertices))
    @test coords == Dict(
        "K" => [0.375, 0.375, 0.75],
        "L" => [0.5, 0.5, 0.5],
        "U" => [0.625, 0.25, 0.625],
        "W" => [0.5, 0.25, 0.75],
        "X" => [0.5, 0.0, 0.5],
        "Γ" => [0.0, 0.0, 0.0],
    )
    @test ticks(kp).labels == ["Γ", "X", "U|K", "Γ", "L", "W", "X"]
end

@testitem "Brillouin interop" setup = [KPathEnv] begin
    import Bravais, Brillouin

    kp = KPath(KPathEnv.recip_lattice, KPathEnv.kpoint_path; n_points_first_segment = 5)
    kpi = Brillouin.KPathInterpolant(kp)

    @test all(isapprox.(kpi, kpoints(kp); atol = 1.0e-5))
    @test mat3(kpi.basis) == kp.recip_lattice
    @test kpi.labels[1] == Dict(i => Symbol(l) for (i, l) in zip(ticks(kp)...))

    # round trip through the interpolant gives the verbatim path
    back = KPath(kpi)
    @test back ≈ KPath(KPathEnv.recip_lattice, kpoints(kp), ticks(kp)...)

    # two disconnected lines map to two subpaths
    kp2 = KPath(
        KPathEnv.recip_lattice,
        [["A" => [0.0, 0.0, 0.0], "B" => [0.5, 0.0, 0.0]], ["C" => [0.0, 0.0, 0.0], "D" => [0.0, 0.5, 0.0]]];
        n_points_first_segment = 4,
    )
    kpi2 = Brillouin.KPathInterpolant(kp2)
    @test length(kpi2.kpaths) == 2
    @test length(KPath(kpi2).subpaths) == 2
end
