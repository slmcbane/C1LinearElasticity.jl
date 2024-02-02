module Basis2D

export map_to_points, evaluate_partials

using LinearAlgebra: I, Matrix

"""
The individual polynomial terms in each of 16 basis functions.
"""
const terms = [
    (x, y) -> x^3 * y^3,
    (x, y) -> x^3 * y^2,
    (x, y) -> y^3 * x^2,
    (x, y) -> x^2 * y^2,
    (x, y) -> x^3 * y,
    (x, y) -> x * y^3,
    (x, y) -> x^2 * y,
    (x, y) -> x * y^2,
    (x, y) -> x^3,
    (x, y) -> y^3,
    (x, y) -> x * y,
    (x, y) -> x^2,
    (x, y) -> y^2,
    (x, y) -> x,
    (x, y) -> y,
    (x, y) -> 1
]

zero_fn(x, y) = zero(x)
one_fn(x, y) = one(x)

"""
Partial derivative with respect to x of each term in `terms`.
"""
const dx = [
    (x, y) -> 3x^2 * y^3,
    (x, y) -> 3x^2 * y^2,
    (x, y) -> 2x * y^3,
    (x, y) -> 2x * y^2,
    (x, y) -> 3x^2 * y,
    (x, y) -> y^3,
    (x, y) -> 2x * y,
    (x, y) -> y^2,
    (x, y) -> 3x^2,
    zero_fn,
    (x, y) -> y,
    (x, y) -> 2x,
    zero_fn,
    one_fn,
    zero_fn,
    zero_fn
]

"""
Partial derivative with respect to y of each term in `terms`.
"""
const dy = [
    (x, y) -> x^3 * 3y^2,
    (x, y) -> x^3 * 2y,
    (x, y) -> x^2 * 3y^2,
    (x, y) -> x^2 * 2y,
    (x, y) -> x^3,
    (x, y) -> x * 3y^2,
    (x, y) -> x^2,
    (x, y) -> x * 2y,
    zero_fn,
    (x, y) -> 3y^2,
    (x, y) -> x,
    zero_fn,
    (x, y) -> 2y,
    zero_fn,
    one_fn,
    zero_fn
]

"""
Mixed partial derivative d/dx(d/dy(f)) for each term in `terms`.
"""
const dxdy = [
    (x, y) -> 9x^2 * y^2,
    (x, y) -> 6x^2 * y,
    (x, y) -> 6x * y^2, (x, y) -> 4x * y, (x, y) -> 3x^2,
    (x, y) -> 3y^2,
    (x, y) -> 2x,
    (x, y) -> 2y,
    zero_fn,
    zero_fn,
    one_fn,
    zero_fn,
    zero_fn,
    zero_fn,
    zero_fn,
    zero_fn
]

"""
constraint_values(i, T)

For the i-th term in `terms`, evaluate its value, the value of its x partial derivative,
the value of its y partial derivative, and the value of its mixed partial derivative at each
vertex of the reference square [-1,-1] x [1, 1]
"""
constraint_values(i, T::DataType) =
    let verts = ((-one(T), -one(T)), (-one(T), one(T)), (one(T), one(T)), (one(T), -one(T)))
        [ (terms[i](v...) for v in verts)...,
          (dx[i](v...) for v in verts)...,
          (dy[i](v...) for v in verts)...,
          (dxdy[i](v...) for v in verts)... ]
    end

constraint_matrix(T::DataType) = vcat((constraint_values(i, T)' for i in 1:16)...)

"""
Solve for the coefficients of the basis functions given the degrees of freedom documented
in `constraint_values`. Result is a 16 x 16 matrix with the coefficients of the i-th basis
function in its i-th row.
"""
coefficients(T::DataType) = constraint_matrix(T) \ I

"""
Given a vector of points (which should be 2-tuple shaped) returns a matrix whose action
on an input vector of coefficients of a function in the basis defined in this module
maps it to the values of the function at the given points. Useful for quadrature.
"""
function map_to_points(points)
    T = points |> eltype |> eltype
    N = length(points)
    map = Matrix{T}(undef, N, 16)

    for j in 1:16
        for i in 1:N
            map[i, j] = terms[j](points[i]...)
        end
    end

    map
end

function evaluate_partials(x, y)
    T = promote_type(typeof(x), typeof(y))
    coeffs = coefficients(T)
    dxs = T[ dx[i](x, y) for i in 1:16 ]
    dys = T[ dy[i](x, y) for i in 1:16 ]
    coeffs * hcat(dxs, dys)
end

end # module Basis