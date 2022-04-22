include("../math/system_of_equations.jl")
include("../elements/element.jl")
include("../conditions/condition.jl")
include("../mesh/mesh.jl")

mutable struct BuilderAndSolver
    mesh::Mesh
    system::SystemOfEquations
end

function BuilderAndSolver(mesh::Mesh)
    return BuilderAndSolver(
        mesh, SystemOfEquations()
    )
end

# Takes local system of equations (from element and condition) and assembles it into the global system
function local_assembly(self::BuilderAndSolver, nodes::Vector{Node}, K::Matrix{Float64}, F::Vector{Float64}, U::Vector{Float64})::Nothing
    for (i::Int64, dof_i::Dof) in enumerate(map(m -> m.dof, nodes))

        if !dof_i.free
            if size(U, 1) > 0
                self.system.u_l[dof_i.id] += U[i]
            end
            continue
        end

        self.system.b_f[dof_i.id] += F[i]

        for (j::Int64, dof_j::Dof) in enumerate(map(m -> m.dof, nodes))

            if dof_j.free
                push!(self.system.A_ff, dof_i.id, dof_j.id, K[i, j])
            else
                push!(self.system.A_fl, dof_i.id, dof_j.id, K[i, j])
            end

        end
    end
end


function build(self::BuilderAndSolver, shape_functions::ShapeFunctions, gauss_data::GaussData, diffusivity::Function, source::Function)::Nothing
    (free_counter, lock_counter) = setup_dofs(self.mesh)
    self.system = SystemOfEquations(free_counter, lock_counter)

    for e in self.mesh.elems
        (K, F, U) = local_system(e, self.system, shape_functions, gauss_data, diffusivity, source)
        local_assembly(self, e.nodes, K, F, U)
    end

    for c in self.mesh.conds
        (K, F, U) = local_system(c)
        local_assembly(self, [c.node], K, F, U)
    end
    return nothing
end

function solve(self::BuilderAndSolver)::Vector{Float64}
    solve(self.system)

    free_dofs   = map(n -> n.id, filter(n ->  n.dof.free, self.mesh.nodes))
    locked_dofs = map(n -> n.id, filter(n -> !n.dof.free, self.mesh.nodes))

    u = zeros(Float64, size(self.mesh.nodes))
    u[free_dofs] = self.system.u_f
    u[locked_dofs] = self.system.u_l
    return u
end