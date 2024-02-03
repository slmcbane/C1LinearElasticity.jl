using Test

using C1LinearElasticity
using C1LinearElasticity.Mesh2D

@test Mesh2D.renumber_nodes(3, 3) ==
    [1, 2, 5, 11, 3, 4, 6, 10, 8, 7, 9, 12, 15, 14, 13, 16]