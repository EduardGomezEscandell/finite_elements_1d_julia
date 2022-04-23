# Finite element program
#-----------------------
# Solves 1D equations of type
#
#     ∂u/∂t - ∇·μ∇u + c·∇u = s,    u ∈ [0, 1],  t ∈ [0, 10]
#
# μ is the diffusivity
# c is the convection speed
# s is the source term
# u is the unknown
# They're all functions of t, x
#
# Note: if you change the boundary conditions or or initial conditions.
# ensure the initial condition is consistent with the boundary conditions,
# otherwise the solution will oscilate.

include("../src/integrators/time_integrator_theta_method.jl")
include("../src/post_process/post_process.jl")

function demo_unsteady_main()
    @info "Starting Finite Element program"

    ## Settings
    # Time settings
    θ = 0.5                             # Time integration point: θ ∈ [0, 1]
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
    store_to_disk = false

    # Domain settings
    left_bc  = "Neumann", 0.0           # Left boundary condition
    right_bc = "Neumann", 0.0           # Right boundary condition

    # Physical settings
    μ(t, x) = 1e-3                      # Diffusivity constant
    c(t, x) = 0.1                       # Convection velocity
    s(t, x) = 0.0                       # Source term
    u0(x) = (abs.(x.-0.2) .< 0.05)*1.0  # Initial condition

    ## Meshing
    mesh = generate_mesh(length, ("ConvectionDiffusion", polynomial_order, nelems), left_bc, right_bc)
    @info "Meshing completed"

    # Validating
    Δx = mesh.nodes[2].x -mesh.nodes[1].x
    μ_max = maximum([μ(0, node.x) for node in mesh.nodes])
    c_max = maximum([c(0, node.x) for node in mesh.nodes])
    peclet = c_max * Δx / μ_max
    @info "Peclet number: $(peclet)"

    # Chosing tools
    space_integrator = SpaceIntegrator(mesh, n_gauss_numerical)
    time_integrator = TimeIntegratorThetaMethod(space_integrator, t_end, n_steps, θ)
    plotter = Plotter(mesh, n_gauss_plotting)

    end_of_step_hook = (u; kwargs...) -> begin
        time = kwargs[:time]
        step = kwargs[:step]
        p = plot_step(plotter, u; title = "'Convection-diffusion equation. t=$(floor(1000*time)/1000)s'", yrange=(-0.5, 1.1))
        display(p)
        if store_to_disk
            save(term = "png", saveopts = "size 600,400", output="results/convdiff_$(step).gif")
        end
        sleep(wallclock_wait_time)
    end

    # Solving
    integrate(time_integrator, end_of_step_hook; s=s, μ=μ, u0=u0, c=c)
    @info "Solved"

    return nothing
end


demo_unsteady_main()