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

include("../src/integrators/time_integrator_quasistatic.jl")
include("../src/post_process/post_process.jl")

function demo_quasi_static_main()::Nothing
    @info "Starting Finite Element program"

    ## Settings
    # Time settings
    t_end = 10.0                        # End time
    n_steps = 100                       # Number of time steps

    # Space settings
    length = 1.0                        # Size of the domain
    nelems = 10                         # Number of elements
    polynomial_order = 3                # Polynomial order of the elements
    n_gauss_numerical = 3               # Gauss interpolation for integration

    # Plotting settings
    n_gauss_plotting = 7                # Gauss interpolation for plotting solution
    wallclock_wait_time = 1/30          # Time between frames

    # Domain settings
    left_bc  = "Dirichlet", 0.0         # Left boundary condition
    right_bc = "Dirichlet", 0.0         # Right boundary condition

    # Physical settings
    s(t, x) = - 100 * cos.(3*pi*x) * sin(t)    # Source term
    μ(t, x) = 1                                # Diffusivity constant

    ## Meshing
    mesh = generate_mesh(length, ("Laplacian", polynomial_order, nelems), left_bc, right_bc)
    @info "Meshing completed"

    # Chosing tools
    space_integrator = SpaceIntegrator(mesh, n_gauss_numerical)
    time_integrator = TimeIntegratorQuasiStatic(space_integrator, t_end, n_steps)
    plotter = Plotter(mesh, n_gauss_plotting)

    end_of_step_hook = (u; kwargs...) -> begin
        step = kwargs[:step]
        time = kwargs[:time]
        @info "TimeIntegratorQuasiStatic: Solved step $(step) t=$(time)"
        display(plot_step(plotter, u; title =  "'Solution at t=$(time)s'",yrange=(-2, 2)))
        sleep(wallclock_wait_time)
    end

    # Solving
    integrate(time_integrator, end_of_step_hook; s=s, μ=μ)

    @info "Solved"

    return nothing
end


demo_quasi_static_main()