@testitem "group_nearby_indices" begin
    using CrystalBase: group_nearby_indices
    indices = [1, 2, 4, 5, 6]
    grps = group_nearby_indices(indices)
    @test grps == [[1, 2], [4, 5, 6]]
end

@testitem "linear_path" begin
    # TODO
end
