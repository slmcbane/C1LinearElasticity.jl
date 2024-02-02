module Stiffness2D

using ..Quadrature: five_point_rule_2d
using ..Basis2D: map_to_points, evaluate_partials

"""
See notes for definition of K1. Returns the flattened, lower diagonal part of K1 / lambda
at point (x, y). Since K1 is symmetric it would be wasteful to compute the full thing
(we only need 528 entries vs. 1024).
"""
function flat_K1_at_point(x, y)
    T = promote_type(typeof(x), typeof(y))
    dϕ = evaluate_partials(x, y)
    out = Vector{T}(undef, 528)
    out_index = 1

    for col = 1:32
        for row = col:32
            i = (row - 1) ÷ 2 + 1
            j = (col - 1) ÷ 2 + 1
            p_i = (row + 1) % 2 + 1
            p_j = (col + 1) % 2 + 1
            out[out_index] = dϕ[i, p_i] * dϕ[j, p_j]
            out_index += 1
        end
    end
    out
end

function K1_map(T::DataType)
    rule = five_point_rule_2d(BigFloat)
    Q = Matrix{BigFloat}(undef, 528, 25)
    for j in 1:25
        Q[:, j] .= rule.weights[j] .* flat_K1_at_point(rule.points[j]...)
    end
    L = map_to_points(rule.points)
    T.(Q * L)
end

end # Stiffness2D
