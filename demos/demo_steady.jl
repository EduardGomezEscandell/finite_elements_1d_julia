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

include("../src/integrators/time_integrator_factory.jl")
include("../src/post_process/post_process.jl")

function demo_steady_main()::Nothing
    @info "Starting Finite Element program"
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
    s(t, x) = - 100 * cos.(3*pi*x)      # Source term
    μ(t, x) = 1                         # Diffusivity term

    ## Meshing
    mesh = generate_mesh(length, ("Laplacian", polynomial_order, nelems), left_bc, right_bc)
    @info "Meshing completed"

    # Chosing tools
    space_integrator = SpaceIntegrator(mesh, n_gauss_numerical)
    time_integrator = time_integrator_factory("steady", space_integrator)
    plotter = Plotter(mesh, n_gauss_plotting)

    end_of_step_hook = (u; kwargs...) -> display(plot_step(plotter, u; title = "'Solution'"))

    # Solving
    integrate(time_integrator, end_of_step_hook; s=s, μ=μ)
    @info "Solved"

    return nothing

end

demo_steady_main()