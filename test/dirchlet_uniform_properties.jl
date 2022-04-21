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

include("../src/math/gauss.jl")
include("../src/math/shape_fun.jl")
include("../src/mesh/mesh.jl")
include("../src/math/system_of_equations.jl")

using Test

function solve_simple_laplacian(nelems::Integer, polynomial_order::Integer, n_gauss::Integer)::Vector{Float64}
    # Problem-specific settings
    length = 1.0                        # Size of the domain
    left_bc = DIRICHLET, 1.0            # Left boundary condition
    right_bc = DIRICHLET, 1.0           # Right boundary condition
    source(x) = 3                       # Source term f
    diffusivity(x) = 2                  # Diffusivity constant k

    mesh = generate_mesh(nelems, polynomial_order, length, left_bc, right_bc)
    gauss_data = get_gauss_quadrature(n_gauss)
    shape_functions = compute_shape_functions(polynomial_order, gauss_data)
    system = build(mesh, shape_functions, gauss_data, diffusivity, source)
    return solve(system, mesh)
end

analytical(x::Float64)::Float64 = x*(1 - x)*3/4 + 1

@testset "Validation: dirchlet_uniform_properties" begin
    # First order polynomials, 2 elements
    u = solve_simple_laplacian(2, 1, 1)
    @test size(u,1) == 3
    @test u[1] ≈ analytical(0.0)   atol=1e-8
    @test u[2] ≈ analytical(0.5)   atol=1e-8
    @test u[3] ≈ analytical(1.0)   atol=1e-8

    # First order polynomials, 3 elements
    u = solve_simple_laplacian(3, 1, 1)
    @test size(u,1) == 4
    @test u[1] ≈ analytical(0.0)   atol=1e-8
    @test u[2] ≈ analytical(1/3)   atol=1e-8
    @test u[3] ≈ analytical(2/3)   atol=1e-8
    @test u[4] ≈ analytical(1.0)   atol=1e-8

    # Second order polynomials,1 element
    u = solve_simple_laplacian(1, 2, 2)
    @test size(u,1) == 3
    @test u[1] ≈ analytical(0.0)   atol=1e-8
    @test u[2] ≈ analytical(0.5)   atol=1e-8
    @test u[3] ≈ analytical(1.0)   atol=1e-8

    # Second order polynomials,2 elements
    u = solve_simple_laplacian(2, 2, 2)
    @test size(u,1) == 5
    @test u[1] ≈ analytical(0.0)   atol=1e-8
    @test u[2] ≈ analytical(1/4)   atol=1e-8
    @test u[3] ≈ analytical(1/2)   atol=1e-8
    @test u[4] ≈ analytical(3/4)   atol=1e-8
    @test u[5] ≈ analytical(1.0)   atol=1e-8

    # Third order polynomials, 1 element
    u = solve_simple_laplacian(1, 3, 3)
    @test size(u,1) == 4
    @test u[1] ≈ analytical(0.0)   atol=1e-8
    @test u[2] ≈ analytical(1/3)   atol=1e-8
    @test u[3] ≈ analytical(2/3)   atol=1e-8
    @test u[4] ≈ analytical(1.0)   atol=1e-8

    # Third order polynomials, 2 elements
    u = solve_simple_laplacian(2, 3, 3)
    @test size(u,1) == 7
    @test u[1] ≈ analytical(0.0)   atol=1e-8
    @test u[2] ≈ analytical(1/6)   atol=1e-8
    @test u[3] ≈ analytical(1/3)   atol=1e-8
    @test u[4] ≈ analytical(1/2)   atol=1e-8
    @test u[5] ≈ analytical(2/3)   atol=1e-8
    @test u[6] ≈ analytical(5/6)   atol=1e-8
    @test u[7] ≈ analytical(1.0)   atol=1e-8
end