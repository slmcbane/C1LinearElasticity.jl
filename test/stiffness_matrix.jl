using .Basis2D: terms, coefficients
using HCubature: hcubature
using LinearAlgebra: Symmetric
using .Stiffness2D: K1_map
using ForwardDiff

@testset "2D stiffness tests" begin

    basis_coefficients = coefficients(Float64)

    @test size(basis_coefficients) == (16, 16)

    naive_basis = [(x, y) -> sum(basis_coefficients[i, j] * terms[j](x, y) for j in 1:16) for i in 1:16]

    let x = Quadrature.five_point_rule_2d(Float64).points[7]
        @test naive_basis[1](x...) + naive_basis[2](x...) + naive_basis[3](x...) + naive_basis[4](x...) ≈ 1.0
    end

    ∂x(f) = (x, y) -> ForwardDiff.derivative(z -> f(z, y), x)
    ∂y(f) = (x, y) -> ForwardDiff.derivative(z -> f(x, z), y)

    function K1_element(i, j, lambda)
        idir = i % 2 == 0 ? :y : :x
        jdir = j % 2 == 0 ? :y : :x
        i = (i - 1) ÷ 2 + 1
        j = (j - 1) ÷ 2 + 1
        ϕᵢ = naive_basis[i]
        ϕⱼ = naive_basis[j]

        f1 = idir == :x ? ∂x(ϕᵢ) : ∂y(ϕᵢ)
        f2 = jdir == :x ? ∂x(ϕⱼ) : ∂y(ϕⱼ)

        integrand = x -> begin
            (ξ, η) = (x[1], x[2])
            lambda(ξ, η) * f1(ξ, η) * f2(ξ, η)
        end

        I, _ = hcubature(integrand, [-1, -1], [1, 1])
        I
    end

    function naive_K1(lambda)
        K = Matrix{Float64}(undef, 32, 32)
        for j = 1:32
            for i = j:32
                K[i, j] = K1_element(i, j, lambda)
            end
        end
        Symmetric(K, :L)
    end

    function reshape_to_symmetric(K::Vector)
        @assert length(K) == 528
        dest = Matrix{Float64}(undef, 32, 32)
        ind = 1
        for j = 1:32
            for i = j:32
                dest[i, j] = K[ind]
                ind += 1
            end
        end
        Symmetric(dst, :L)
    end

    let M = K1_map(Float64)
        lambda = (x, y) -> 1.0
        L = zeros(16)
        L[1:4] .= 1
        #Kref = naive_K1(lambda)
        K = M * L
        @test K[1] ≈ K1_element(1, 1, lambda)
        @test K[2] ≈ K1_element(2, 1, lambda)
    end
end