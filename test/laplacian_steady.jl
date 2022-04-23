# Validation test
#
# Solves equation
#
#       -2∇²u(x) = 3     x ∈ [0, 1]
#           u(0) = 1
#           u(1) = 1
#
# Which is known to be equal to
#
#     u(x) = 1 + 3/4*x*(1 - x)
#

include("../src/integrators/time_integrator_steady.jl")

using Test

function solve_simple_laplacian(nelems::Integer, polynomial_order::Integer, n_gauss::Integer)::Vector{Float64}
    # Problem-specific settings
    length = 1.0                        # Size of the domain
    elements = ("Laplacian", polynomial_order, nelems)
    left_bc  = ("Dirichlet", 1.0)       # Left boundary condition
    right_bc = ("Dirichlet", 1.0)       # Right boundary condition
    s(t, x) = 3                         # Source term
    μ(t, x) = 2                         # Diffusivity constant

    mesh = generate_mesh(length, elements, left_bc, right_bc)

    space_integrator = SpaceIntegrator(mesh, n_gauss)
    time_integrator = TimeIntegratorSteady(space_integrator)

    return integrate(time_integrator; s=s, μ=μ)
end

analytical(x::Float64)::Float64 = x*(1 - x)*3/4 + 1

@testset "Steady state" begin
    # First order polynomials, 2 elements
    U = solve_simple_laplacian(2, 1, 1)
    @test size(U, 1) == 3
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=1e-8
    end

    # First order polynomials, 3 elements
    U = solve_simple_laplacian(3, 1, 1)
    @test size(U, 1) == 4
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=1e-8
    end

    # Second order polynomials,1 element
    U = solve_simple_laplacian(1, 2, 2)
    @test size(U, 1) == 3
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=1e-8
    end

    # Second order polynomials,2 elements
    U = solve_simple_laplacian(2, 2, 2)
    @test size(U, 1) == 5
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=1e-8
    end

    # Third order polynomials, 1 element
    U = solve_simple_laplacian(1, 3, 3)
    @test size(U, 1) == 4
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=1e-8
    end

    # Third order polynomials, 2 elements
    U = solve_simple_laplacian(2, 3, 3)
    @test size(U, 1) == 7
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=1e-8
    end
end