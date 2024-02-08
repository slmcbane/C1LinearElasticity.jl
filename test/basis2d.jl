using .Basis2D
using .Quadrature

@testset "2D basis tests" begin
    let λ = zeros(Float64, 16), points = five_point_rule_2d(Float64).points
        λ[1:4:end] .= 1
        map = Basis2D.map_to_points(points)
        mapped = map * λ
        println(mapped)
        @test all(≈(1.0), mapped)
    end
end