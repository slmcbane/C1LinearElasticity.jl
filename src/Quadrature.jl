module Quadrature

export five_point_rule_2d

struct Rule{T}
    points::Vector{Tuple{T, T}}
    weights::Vector{T}
end

"""
The points for a five point Gauss-Legendre rule (integrates over [-1, 1]) in one dimension.
This rule integrates polynomials of degree 9 or less exactly which is enough to evaluate
the entries of the stiffness (resp. mass) matrix exactly.
"""
function five_point_points(T::DataType)
    pt1 = ((5 - 2√(10 * one(T) / 7)) |> sqrt) / 3
    pt2 = ((5 + 2√(10 * one(T) / 7)) |> sqrt) / 3
    ( zero(T), pt1, -pt1, pt2, -pt2 )
end

"""
Weights for five point Gauss-Legendre rule in 1 dimension.
"""
function five_point_weights(T::DataType)
    w1 = (322 + 13√(70 * one(T))) / 900
    w2 = (322 - 13√(70 * one(T))) / 900
    ((one(T) * 128) / 225, w1, w1, w2, w2)
end

"""
Generate a 25-point rule to integrate in two dimensions given the type
(probably should be Float32 or Float64)
"""
function five_point_rule_2d(T::DataType)
    w = five_point_weights(T)
    x = five_point_points(T)
    Rule{T}(
        [(x[i], x[j]) for i = 1:5 for j = 1:5],
        [w[i] * w[j] for i = 1:5 for j = 1:5])
end

end # module Quadrature