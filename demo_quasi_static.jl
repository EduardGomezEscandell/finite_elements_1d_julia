#!/bin/bash
#=
exec julia -i --compile=min "${BASH_SOURCE[0]}" "$@"
=#
# Finite element program
#-----------------------
# Solves 1D equations of type
#
#       -∇·μ∇u = f
#
# μ is the diffusivity
# f is the source term
# u is the unknown
# They're all functions of x

include("src/math/gauss.jl")
include("src/math/shape_fun.jl")
include("src/mesh/mesh.jl")
include("src/math/system_of_equations.jl")
include("src/post_process.jl")

function main()
    println("Starting Finite Element program")

    ## Settings
    # Time settings
    t_start = 0.0                       # Start time
    t_end = 10.0                        # Start time
    n_steps = 100                       # Number of time steps

    # Space settings
    length = 1.0                        # Size of the domain
    nelems = 10                         # Number of elements
    polynomial_order = 3                # Polynomial order of the elements
    n_gauss_numerical = 3               # Gauss interpolation for integration

    # Plotting settings
    n_gauss_plotting = 7                # Gauss interpolation for plotting solution
    wallclock_wait_time = 0.1           # Time between frames

    # Domain settings
    left_bc = DIRICHLET, 0.0            # Left boundary condition
    right_bc = DIRICHLET, 0.0           # Right boundary condition

    # Physical settings
    source(t, x) = - 100 * cos.(3*pi*x) * sin(t)    # Source term f
    diffusivity(t, x) = 1                           # Diffusivity constant μ

    ## Meshing
    mesh = generate_mesh(nelems, polynomial_order, length, left_bc, right_bc)
    println("Meshing completed")

    ## Precomputing data
    gauss_data = get_gauss_quadrature(n_gauss_numerical)
    shape_functions = compute_shape_functions(polynomial_order, gauss_data)
    println("Preliminaries completed")

    # Time-dependent stuff
    curry(f, x) = (xs...) -> f(x, xs...)

    # Assembly
    for (step, time) in enumerate(LinRange(t_start, t_end, n_steps))
        wallclock_start = time_ns()/1e9

        formatted_time = string(floor(1000 * time) / 1000)
        println("STEP $(step)  t = $(formatted_time)s")

        μ = curry(diffusivity, time)
        s = curry(source, time)

        system = build(mesh, shape_functions, gauss_data, μ, s)

        # Solution
        solve(system, mesh)

        # Output
        plotting_gauss = get_gauss_quadrature(n_gauss_plotting)
        plotting_shape_fun = compute_shape_functions(polynomial_order, plotting_gauss)
        p = plot_solution(mesh, system.sol, plotting_shape_fun, plotting_gauss; title = "'Solution at t=$(formatted_time)s'",yrange=(-2, 2))

        wallclock_time_elapsed = time_ns()/1e9 - wallclock_start
        sleep(max(wallclock_wait_time - wallclock_time_elapsed, 0))
        display(p)
    end

    return nothing
end


main()