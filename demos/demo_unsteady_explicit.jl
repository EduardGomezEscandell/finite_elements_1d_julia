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

include("../src/integrators/time_integrator_RK4.jl")
include("../src/integrators/time_integrator_forward_euler.jl")
include("../src/post_process/post_process.jl")

@enum Scheme FORWARD_EULER RUNGE_KUTTA_4

function demo_unsteady_main()
    @info "Starting Finite Element program"

    ## Settings
    # Time settings
    scheme::Scheme = RUNGE_KUTTA_4
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
    mesh = generate_mesh(length, ("Laplacian", polynomial_order, nelems), left_bc, right_bc)
    @info "Meshing completed"

    # Validating
    Δx = mesh.nodes[2].x -mesh.nodes[1].x
    Δt = t_end / n_steps
    μ_max = maximum([μ(0, node.x) for node in mesh.nodes])
    fourier = μ_max * Δt / Δx^2
    max_fourier = scheme == FORWARD_EULER ? 0.5 : 0.7  # Experimental data
    if fourier < max_fourier
        @info "Fourier number: $(fourier)"
    else
        @warn "Fourier number beyond stability threshold!: $(fourier) > $(max_fourier)"
    end

    # Chosing tools
    space_integrator = SpaceIntegrator(mesh, n_gauss_numerical)

    if scheme == FORWARD_EULER
        time_integrator = TimeIntegratorForwardEuler(space_integrator, t_end, n_steps)
    elseif scheme == RUNGE_KUTTA_4
        time_integrator = TimeIntegratorRK4(space_integrator, t_end, n_steps)
    end

    plotter = Plotter(mesh, n_gauss_plotting)

    end_of_step_hook = (u; kwargs...) -> begin
        time = kwargs[:time]
        display(plot_step(plotter, u; title = "'Solution at t=$(floor(1000*time)/1000)s'", yrange=(-1.1, 1.1)))
        sleep(wallclock_wait_time)
    end

    # Solving
    integrate(time_integrator, end_of_step_hook; s=s, μ=μ, u0=u0)
    @info "Solved"

    return nothing
end


demo_unsteady_main()