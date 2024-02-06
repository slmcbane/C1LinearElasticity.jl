using HCubature: hcubature, hquadrature
using ForwardDiff

using .Quadrature
using .Basis2D

function my_1d_quadrature(f)
    ws = Quadrature.five_point_weights(Float64)
    xs = Quadrature.five_point_points(Float64)
    sum(w * f(x) for (w, x) in zip(ws, xs))
end

function my_2d_quadrature(f)
    rule = five_point_rule_2d(Float64)
    sum(w * f(x) for (w, x) in zip(rule.weights, rule.points))
end

@testset "Quadrature tests" begin
    let f = x -> x^9 + 4 * x^5 - x^4 + 3
        @test hquadrature(f, -1, 1)[1] ≈ my_1d_quadrature(f)
    end

    let f = x -> Basis2D.terms[1](x[1], x[2])
        @test hcubature(f, [-1, -1], [1, 1])[1] ≈ my_2d_quadrature(f)
    end

    ∂x(f) = (x, y) -> ForwardDiff.derivative(z -> f(z, y), x)
    let f = x -> ∂x(Basis2D.terms[1])(x[1], x[2]) * ∂x(Basis2D.terms[2])(x[1], x[2])
        @test hcubature(f, [-1, -1], [1, 1])[1] ≈ my_2d_quadrature(f)
    end
end
