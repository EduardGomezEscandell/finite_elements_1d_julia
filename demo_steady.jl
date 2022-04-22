# Finite element program
#-----------------------
# Solves 1D equations of type
#
#       -∇·μ∇u = s
#
# μ is the diffusivity
# s is the source term
# u is the unknown
# They're all functions of x

include("src/all.jl")

function demo_steady_main()
    println("Starting Finite Element program")
    ## Settings
    # Numerical settings
    nelems = 10                         # Number of elements
    polynomial_order = 3                # Polynomial order of the elements
    n_gauss_numerical = 3               # Gauss interpolation for integration
    n_gauss_plotting = 7                # Gauss interpolation for plotting solution

    # Domain settings
    length = 1.0                        # Size of the domain
    left_bc = "Neumann", -1.0           # Left boundary condition
    right_bc = "Dirichlet", -2.0        # Right boundary condition

    # Physical settings
    s(x) = - 100 * cos.(3*pi*x)         # Source term
    μ(x) = 1                            # Diffusivity term

    ## Meshing
    mesh = generate_mesh(length, ("Laplacian", polynomial_order, nelems), left_bc, right_bc)
    println("Meshing completed")

    # Precomputing data
    gauss_data = get_gauss_quadrature(n_gauss_numerical)
    shape_functions = compute_shape_functions(polynomial_order, gauss_data)
    bns = BuilderAndSolver(mesh)
    println("Preliminaries completed")

    # Assembly
    build(bns, shape_functions, gauss_data, μ=μ, s=s)
    println("Assembly completed")

    # Solution
    u = solve(bns)
    println("Solving completed")

    # Output
    plotting_gauss = get_gauss_quadrature(n_gauss_plotting)
    plotting_shape_fun = compute_shape_functions(polynomial_order, plotting_gauss)
    p = plot_solution(mesh, u, plotting_shape_fun, plotting_gauss; title = "'Solution'")

    display(p)
end

demo_steady_main()