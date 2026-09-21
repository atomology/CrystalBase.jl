@testitem "atomic_number" begin
    @test atomic_number("O") == 8
    @test atomic_number(["H", "O"]) == [1, 8]
    @test atomic_number("X") == 0
    @test isnothing(atomic_number("Zz"))
end

@testitem "atomic_symbol" begin
    @test atomic_symbol(8) == "O"
    @test atomic_symbol([1, 8]) == ["H", "O"]
    @test atomic_symbol(0) == "X"
end

@testitem "label_symbol" begin
    @test label_symbol("Fe1") == "Fe"
    @test label_symbol("O") == "O"
    @test label_symbol("Co") == "Co"
    @test label_symbol("C1") == "C"
    @test label_symbol("Cl_dn") == "Cl"
    @test label_symbol("X1") == "X"
    @test_throws ArgumentError label_symbol("Zz")
    @test_throws ArgumentError label_symbol("")
end

@testitem "formula" begin
    @test formula(["Si", "O", "O"]) == "O2Si"
    @test formula(["Si", "O", "O"]; order = :alpha) == "O2Si"
    @test formula(["C", "O", "H", "H", "H", "H"]) == "CH4O"
    @test formula(["C", "C", "O", "O"]; order = :alpha) == "C2O2"
    @test formula(["Fe", "Fe", "O", "O", "O", "O"]; reduced = true) == "FeO2"
    @test formula(fill("Fe", 4) ∪ String[]; reduced = true) == "Fe"
    @test formula(vcat(fill("Fe", 4), fill("O", 6)); reduced = true) == "Fe2O3"
    @test formula(["Si", "Si"]; reduced = true) == "Si"
    @test formula(String[]) == ""
    @test_throws ArgumentError formula(["Si"]; order = :bogus)
end
