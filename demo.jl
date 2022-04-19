# Finite element program
# Solves 1D equations of type ∇·k∇u = f
#
# k is the diffusivity
# f is the source term
# u is the unknown
# They're all functions of x

include("src/gauss.jl")
include("src/shape_fun.jl")
include("src/mesh.jl")
include("src/build_and_solve.jl")

## Settings
# Numerical settings
n_gauss = 4
polynomial_order = 3
nelems = 2

# Domain settings
length = 2.0
left_bc = DIRICHLET, 0.0
right_bc = DIRICHLET, 0.0

# Physical settings
source(x) = 1.0 .* ones(size(x))
diffusivity = 1.0

## Meshing
mesh = generate_mesh(nelems, polynomial_order, length)
println("Meshing completed")

# Precomputing data
gauss_data = get_gauss_quadrature(n_gauss)
shape_functions = compute_shape_functions(polynomial_order, gauss_data)
println("Preliminaries completed")

# Assembly
system = build(mesh, shape_functions, gauss_data, diffusivity, source)
println("Assembly completed")

# Solution
solve(system)
println("Solving completed")

# Output
println("\nSolution:")
display(system.sol)
println()