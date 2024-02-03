using C1LinearElasticity: Mesh2D

@testset "Mesh2D tests" begin
    @test Mesh2D.renumber_nodes(3, 3) ==
          [1, 2, 5, 11, 3, 4, 6, 10, 8, 7, 9, 12, 15, 14, 13, 16]

    let mesh = Mesh2D.Mesh{Int,Float64}(1, 1, 1.0, 1.0, 1.0, 1.0)
        I, J = Mesh2D.sparsity(mesh)
        refJ = [j for i = 1:32 for j = i:32]
        @test J == refJ
        refI = [1, 33, 64, 94, 123, 151, 178, 204, 229, 253, 276, 298,
            319, 339, 358, 376, 393, 409, 424, 438, 451, 463, 474,
            484, 493, 501, 508, 514, 519, 523, 526, 528, 529]
        @test I == refI
    end
end