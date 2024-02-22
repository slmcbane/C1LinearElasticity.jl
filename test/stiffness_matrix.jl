using .Basis2D: terms, coefficients
using HCubature: hcubature
using LinearAlgebra: Symmetric
using .Stiffness2D: K1_map, K2_map
using ForwardDiff
using ProgressMeter: @showprogress

@testset "Test the K1 component of the stiffness matrix" begin
    basis_coefficients = coefficients(Float64)

    @test size(basis_coefficients) == (16, 16)

    naive_basis = [(x, y) -> sum(basis_coefficients[i, j] * terms[j](x, y) for j in 1:16) for i in 1:16]

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

        I, _ = hcubature(integrand, [-1, -1], [1, 1], atol=1e-10)
        I
    end

    let M = K1_map(Float64)
        lambda = (x, y) -> 1.0
        L = zeros(16)
        L[1:4:end] .= 1
        i = 1
        j = 1
        K = M * L
        @showprogress desc = "Checking elements of K1..." for k = 1:528
            @test isapprox(K[k], K1_element(i, j, lambda), atol=1e-10)
            i += 1
            if i == 33
                j += 1
                i = j
            end
        end
    end

    function K2_element(i, j, mu)
        idir = i % 2 == 0 ? :y : :x
        jdir = j % 2 == 0 ? :y : :x
        i = (i - 1) ÷ 2 + 1
        j = (j - 1) ÷ 2 + 1
        ϕᵢ = naive_basis[i]
        ϕⱼ = naive_basis[j]

        if idir == jdir
            f1 = idir == :x ? ∂x(ϕᵢ) : ∂y(ϕᵢ)
            f2 = jdir == :x ? ∂x(ϕⱼ) : ∂y(ϕⱼ)
            g1 = idir == :x ? ∂y(ϕᵢ) : ∂x(ϕᵢ)
            g2 = idir == :x ? ∂y(ϕⱼ) : ∂x(ϕⱼ)
            integrand = x -> begin
                (ξ, η) = (x[1], x[2])
                mu(ξ, η) * (2 * f1(ξ, η) * f2(ξ, η) + g1(ξ, η) * g2(ξ, η))
            end
        else
            f1 = jdir == :x ? ∂x(ϕᵢ) : ∂y(ϕᵢ)
            f2 = idir == :x ? ∂x(ϕⱼ) : ∂y(ϕⱼ)
            integrand = x -> begin
                (ξ, η) = (x[1], x[2])
                mu(ξ, η) * f1(ξ, η) * f2(ξ, η)
            end
        end
        I, _ = hcubature(integrand, [-1, -1], [1, 1], atol=1e-10)
        I
    end

    let M = K2_map(Float64)
        mu = (x, y) -> 1.0
        L = zeros(16)
        L[1:4:end] .= 1
        i = 1
        j = 1
        K = M * L
        @showprogress desc = "Checking elements of K2..." for k = 1:528
            @test isapprox(K[k], K2_element(i, j, mu), atol=1e-10)
            i += 1
            if i == 33
                j += 1
                i = j
            end
        end
    end
end