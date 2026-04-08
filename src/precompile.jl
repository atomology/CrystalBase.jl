using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    lattice = mat3(
        [5.43, 0.0, 0.0],
        [0.0, 5.43, 0.0],
        [0.0, 0.0, 5.43],
    )
    recip_lattice = reciprocal_lattice(lattice)

    coords = OrderedDict(
        "GAMMA" => vec3(0.0, 0.0, 0.0),
        "X" => vec3(0.5, 0.0, 0.0),
        "M" => vec3(0.5, 0.5, 0.0),
        "R" => vec3(0.5, 0.5, 0.5),
    )
    segments = [["GAMMA", "X", "M"], ["M", "R"]]

    kseg = KSegment(recip_lattice, segments, coords)
    kpath = KPath(kseg, 80)
    frac_positions = [vec3(0.0, 0.0, 0.0), vec3(0.25, 0.25, 0.25)]

    @compile_workload begin
        vec3([1.0, 2.0, 3.0])
        mat3([1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0])
        mvec3([1.0, 2.0, 3.0])
        mmat3([1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0])
        stringvec3("Si", [0.0, 0.5, 0.5])

        reciprocal_lattice(lattice)
        real_lattice(recip_lattice)
        lattice_vectors(lattice)
        reciprocal_lattice_vectors(lattice)
        real_lattice_vectors(recip_lattice)
        frac_to_cart(lattice, frac_positions)
        cart_to_frac(lattice, frac_to_cart(lattice, frac_positions))

        atomic_number("Si")
        atomic_number(["Si", "O"])
        atomic_symbol(14)
        atomic_symbol([14, 8])

        unicode_kpoint_labels(["GAMMA", "DELTA_0", "SIGMA_2", "X"])
        linear_path(kpath)
        KPath(recip_lattice, frac_positions, [1, 2], ["GAMMA", "R"])
    end
end
