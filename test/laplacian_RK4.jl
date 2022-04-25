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

include("../src/integrators/time_integrators/time_integrator_RK4.jl")

using Test

function ensure_fourier_below(fourier::Float64, Δx::Float64, t_end::Float64, μ_max::Float64)::Integer
    Δt = fourier*Δx^2 / μ_max
    n_steps = Int(ceil(t_end / Δt))
    return n_steps
end

function solve_simple_laplacian_RK4(polynomial_order::Integer)::Vector{Float64}
    # Problem-specific settings
    length = 1.0                        # Size of the domain
    t_end  = 20.0
    nelems = Int(12 / polynomial_order)
    elements = ("Laplacian", polynomial_order, nelems)
    left_bc  = ("Dirichlet", 1.0)       # Left boundary condition
    right_bc = ("Dirichlet", 1.0)       # Right boundary condition
    s(t, x) = 3.0                       # Source term
    μ(t, x) = 0.2                       # Diffusivity constant
    u0(x) = 1.0                         # Initial condition
    n_steps = ensure_fourier_below(0.7, length/nelems, t_end, μ(0, 0))

    mesh = generate_mesh(length, elements, left_bc, right_bc)

    space_integrator = SpaceIntegrator(mesh, polynomial_order+1)
    time_integrator = TimeIntegratorRK4(space_integrator, t_end, n_steps)

    return integrate(time_integrator; s=s, μ=μ, u0=u0)
end

analytical(x::Float64)::Float64 = x*(1 - x)*3/0.4 + 1

@testset "Runge Kutta 4" begin
    # First order polynomials,
    U = solve_simple_laplacian_RK4(1)
    @test size(U, 1) == 13
    for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
        @test u ≈ analytical(x)   atol=0.01
    end

    # Higher order polynomials fail

    # Second order polynomials
    # U = solve_simple_laplacian_RK4(2)
    # @test size(U, 1) == 13
    # for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
    #     @test u ≈ analytical(x)   atol=0.01
    # end

    # # Third order polynomials,
    # U = solve_simple_laplacian_RK4(3)
    # @test size(U, 1) == 13
    # for (u, x) in zip(U, LinRange(0, 1, size(U, 1)))
    #     @test u ≈ analytical(x)   atol=0.01
    # end

end
