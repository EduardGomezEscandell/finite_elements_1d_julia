# Finite element program
#-----------------------
# Solves 1D equations of type
#
#     ∂u/∂t - ∇·μ∇u = s,    u ∈ [0, 1],  t ∈ [0, 10]
#            u(0,t) = 0
#          μ∇u(1,t) = 0
#            u(x,0) = sin(19*π/2*x)
#
# μ is the diffusivity
# s is the source term
# u is the unknown
# They're all functions of x
#
#
# Note: if you change the boundary conditions or or initial conditions.
# ensure the initial condition is consistent with the boundary conditions,
# otherwise the solution will oscilate.


include("../src/all.jl")

function demo_unsteady_main()
    println("Starting Finite Element program")

    ## Settings
    # Time settings
    t_start = 0.0                       # Start time
    t_end   = 10.0                      # Start time
    n_steps = 100                       # Number of time steps

    # Space settings
    length = 1.0                        # Size of the domain
    nelems = 10                         # Number of elements
    polynomial_order = 3                # Polynomial order of the elements
    n_gauss_numerical = 4               # Gauss interpolation for integration

    # Plotting settings
    n_gauss_plotting = 7                # Gauss interpolation for plotting solution
    wallclock_wait_time = 1/12          # Time between frames

    # Domain settings
    left_bc  = "Dirichlet", 0.0         # Left boundary condition
    right_bc = "Neumann", 0.0           # Right boundary condition

    # Physical settings
    μ(t, x) = 1e-3                      # Diffusivity constant
    s(t, x) = 0.1/(1 + t)               # Source term
    u0(x) = sin.(9.5*pi*x)              # Initial condition

    ## Meshing
    mesh = generate_mesh(length, ("UnsteadyLaplacian", polynomial_order, nelems), left_bc, right_bc)
    println("Meshing completed")

    ## Precomputing data
    gauss_data = get_gauss_quadrature(n_gauss_numerical)
    shape_functions = compute_shape_functions(polynomial_order, gauss_data)
    builder_and_solver = BuilderAndSolver(mesh)
    Δt = (t_end - t_start) / n_steps
    u = u0([node.x for node in mesh.nodes]) .* ones(size(mesh.nodes))
    println("Preliminaries completed")

    for (step, time) in enumerate(LinRange(t_start, t_end, n_steps))
        wallclock_start = time_ns()/1e9

        formatted_time = string(floor(1000 * time) / 1000)
        println("STEP $(step)  t = $(formatted_time)s  Δt = $(Δt)s")

        # Assembly
        build(builder_and_solver, shape_functions, gauss_data; μ=μ, s=s, t=time, Δt=Δt, u_old=u)

        # Solution
        u = solve(builder_and_solver)

        # Output
        plotting_gauss = get_gauss_quadrature(n_gauss_plotting)
        plotting_shape_fun = compute_shape_functions(polynomial_order, plotting_gauss)
        p = plot_solution(mesh, u, plotting_shape_fun, plotting_gauss; title = "'Solution at t=$(formatted_time)s'", yrange=(-1.1, 1.1))

        wallclock_time_elapsed = time_ns()/1e9 - wallclock_start
        sleep(max(wallclock_wait_time - wallclock_time_elapsed, 0))
        display(p)
    end

    return nothing
end


demo_unsteady_main()