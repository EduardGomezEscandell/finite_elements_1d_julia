# Finite element program
#-----------------------
# Solves 1D equations of type
#
#  ∂u/∂t - ∇·μ∇u + u·∇u = s,    u ∈ [0, 5],  t ∈ [0, 10]
#
# μ is the diffusivity
# c is the convection speed
# u is the unknown
# They're all functions of t, x
#
# Note: if you change the boundary conditions or or initial conditions.
# ensure the initial condition is consistent with the boundary conditions,
# otherwise the solution will oscilate.

include("../src/integrators/time_integrator_theta_method.jl")
include("../src/integrators/time_integrator_RK4.jl")
include("../src/post_process/post_process.jl")

function demo_unsteady_main()
    @info "Starting Finite Element program"

    ## Settings
    # Time settings
    θ = 0.5                             # Time integration point: θ ∈ [0, 1]
    t_end   = 10.0                      # Start time
    n_steps = 100                       # Number of time steps

    # Space settings
    length = 5.0                        # Size of the domain
    nelems = 50                         # Number of elements
    polynomial_order = 1                # Polynomial order of the elements
    n_gauss_numerical = 4               # Gauss interpolation for integration

    # Plotting settings
    n_gauss_plotting = 7                # Gauss interpolation for plotting solution
    wallclock_wait_time = 1/12          # Time between frames

    # Domain settings
    left_bc  = "Neumann", 0.0           # Left boundary condition
    right_bc = "Neumann", 0.0           # Right boundary condition

    # Physical settings
    μ(t, x) = 0.04                       # Diffusivity constant
    s(t, x) = 0.0                        # Source term
    u0(x) = 0.2*(1.0 .- x) .+ 0.5        # Initial condition

    ## Meshing
    mesh = generate_mesh(length, ("Burgers", polynomial_order, nelems), left_bc, right_bc)
    @info "Meshing completed"

    # Chosing tools
    space_integrator = SpaceIntegrator(mesh, n_gauss_numerical)
    time_integrator = TimeIntegratorRK4(space_integrator, t_end, n_steps)
    plotter = Plotter(mesh, n_gauss_plotting)

    end_of_step_hook = (u; kwargs...) -> begin
        time = kwargs[:time]
        display(plot_step(plotter, u; title = "'Solution at t=$(floor(1000*time)/1000)s'", yrange=(-0.5, 1.1)))
        sleep(wallclock_wait_time)
    end

    # Solving
    integrate(time_integrator, end_of_step_hook; s=s, μ=μ, u0=u0)
    @info "Solved"

    return nothing
end


demo_unsteady_main()