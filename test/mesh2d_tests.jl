using C1LinearElasticity: Mesh2D

@testset "Mesh2D tests" begin
    @test Mesh2D.renumber_nodes(3, 3) ==
          [1, 2, 5, 11, 3, 4, 6, 10, 8, 7, 9, 12, 15, 14, 13, 16]

    let mesh = Mesh2D.Mesh{Int,Float64}(1, 1, 1.0, 1.0, 1.0, 1.0)
        I, J = Mesh2D.sparsity(mesh)
        refJ = [j for i = 1:32 for j = 1:i]
        @test J == refJ
        refI = [1, 2, 4, 7, 11, 16, 22, 29, 37, 46, 56, 67, 79, 92, 106, 121, 137,
            154, 172, 191, 211, 232, 254, 277, 301, 326, 352, 379, 407, 436,
            466, 497, 529]
        @test I == refI
    end
end