# Validation test
#
# Solves equation
#
#       ∂u/∂t - 2∇²u(x) = 3     x ∈ [0, 1]
#           u(0) = 1
#           u(1) = 1
#
# Which is known to converge to
#
#     u(x) = 1 + 3/4*x*(1 - x)
#
# By t=10 the solution has converged to this value

include("../src/integrators/time_integrator_theta_method.jl")

using Test

function solve_simple_laplacian_theta(polynomial_order::Integer, n_steps::Integer, θ::Float64)::Vector{Float64}
    # Problem-specific settings
    length = 1.0                        # Size of the domain
    t_end  = 20.0
    nelems = Int(12 / polynomial_order)
    elements = ("Laplacian", polynomial_order, nelems)
    left_bc  = ("Dirichlet", 1.0)       # Left boundary condition
    right_bc = ("Dirichlet", 1.0)       # Right boundary condition
    s(t, x) = 3                         # Source term
    μ(t, x) = 2                         # Diffusivity constant
    u0(x) = 1                           # Initial condition

    mesh = generate_mesh(length, elements, left_bc, right_bc)

    space_integrator = SpaceIntegrator(mesh, polynomial_order+1)
    time_integrator = TimeIntegratorThetaMethod(space_integrator, t_end, n_steps, θ)

    return integrate(time_integrator; s=s, μ=μ, u0=u0)
end

analytical(x::Float64)::Float64 = x*(1 - x)*3/4 + 1

@testset "Theta method" begin
    # First order polynomials,
    U = solve_simple_laplacian_theta(1, 20, 0.5)
    @test size(U, 1) == 13
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=0.01
    end

    U = solve_simple_laplacian_theta(1, 20, 1.0)
    @test size(U, 1) == 13
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=0.001
    end

    # Second order polynomials
    U = solve_simple_laplacian_theta(2, 20, 0.5)
    @test size(U, 1) == 13
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=0.01
    end

    U = solve_simple_laplacian_theta(2, 20, 1.0)
    @test size(U, 1) == 13
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=0.01
    end

    # Third order polynomials,
    U = solve_simple_laplacian_theta(3, 20, 0.5)
    @test size(U, 1) == 13
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=0.01
    end

    U = solve_simple_laplacian_theta(3, 20, 1.0)
    @test size(U, 1) == 13
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=0.01
    end
end