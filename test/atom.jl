@testitem "atomic_number" begin
    @test atomic_number("O") == 8
    @test atomic_number(["H", "O"]) == [1, 8]
end

@testitem "atomic_symbol" begin
    @test atomic_symbol(8) == "O"
    @test atomic_symbol([1, 8]) == ["H", "O"]
end
