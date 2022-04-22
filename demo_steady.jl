#!/bin/bash
#=
exec julia -i --compile=min "${BASH_SOURCE[0]}" "$@"
=#
# Finite element program
#-----------------------
# Solves 1D equations of type
#
#       -∇·k∇u = f
#
# k is the diffusivity
# f is the source term
# u is the unknown
# They're all functions of x

include("src/all.jl")

function main()
    println("Starting Finite Element program")
    ## Settings
    # Numerical settings
    nelems = 10                         # Number of elements
    polynomial_order = 3                # Polynomial order of the elements
    n_gauss_numerical = 3               # Gauss interpolation for integration
    n_gauss_plotting = 7                # Gauss interpolation for plotting solution

    # Domain settings
    length = 1.0                        # Size of the domain
    left_bc = "Neumann", -1.0             # Left boundary condition
    right_bc = "Dirichlet", -2.0          # Right boundary condition

    # Physical settings
    source(x) = - 100 * cos.(3*pi*x)    # Source term f
    diffusivity(x) = 1                  # Diffusivity constant k

    ## Meshing
    mesh = generate_mesh(length, ("Laplacian", polynomial_order, nelems), left_bc, right_bc)
    println("Meshing completed")

    # Precomputing data
    gauss_data = get_gauss_quadrature(n_gauss_numerical)
    shape_functions = compute_shape_functions(polynomial_order, gauss_data)
    bns = BuilderAndSolver(mesh)
    println("Preliminaries completed")


    # Assembly
    build(bns, shape_functions, gauss_data, diffusivity, source)
    println("Assembly completed")

    # Solution
    u = solve(bns)
    println("Solving completed")

    # Output
    plotting_gauss = get_gauss_quadrature(n_gauss_plotting)
    plotting_shape_fun = compute_shape_functions(polynomial_order, plotting_gauss)
    p = plot_solution(mesh, u, plotting_shape_fun, plotting_gauss; title = "'Solution'")

    return (mesh, u, p)
end

(mesh, solution, p) = main()
display(p)