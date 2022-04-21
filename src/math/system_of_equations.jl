include("../mesh/element.jl")
include("../mesh/condition.jl")
include("sparse.jl")

mutable struct SystemOfEquations
    A_ff::DokMatrix
    A_fl::DokMatrix
    b_f::Vector{Float64}
    sol::Vector{Float64}
end

function SystemOfEquations(free_dof::Integer, locked_dof::Integer)::SystemOfEquations
    return SystemOfEquations(
        DokMatrix(free_dof, free_dof),
        DokMatrix(free_dof, locked_dof),
        zeros(Float64, free_dof),
        zeros(Float64, free_dof + locked_dof))
end

function setup_dofs(self::Mesh)::Tuple{Int64, Int64}
    for n in self.nodes
        n.dof.free = true
    end

    for c in self.conds
        lock_dofs(c)
    end

    free_counter = 0
    lock_counter = 0

    for n in self.nodes
        if n.dof.free
            free_counter += 1
            n.dof.id = free_counter
        else
            lock_counter += 1
            n.dof.id = lock_counter
        end
    end

    return (free_counter, lock_counter)
end

# Takes local system of equations (from element and condition) and assembles it into the global system
function local_assembly(self::SystemOfEquations, nodes::Vector{Node}, K::Matrix{Float64}, F::Vector{Float64}, U::Vector{Float64})::Nothing
    for (i::Int64, node::Node) in enumerate(nodes)

        if !node.dof.free
            if size(U, 1) > 0
                self.sol[node.id] += U[i]
            end
            continue
        end

        self.b_f[node.dof.id] += F[i]

        for (j::Int64, dof::Dof) in enumerate(map(m -> m.dof, nodes))

            if dof.free
                push!(self.A_ff, node.dof.id, dof.id, K[i, j])
            else
                push!(self.A_fl, node.dof.id, dof.id, K[i, j])
            end

        end
    end
end


function build(mesh::Mesh, shape_functions::ShapeFunctions, gauss_data::GaussData, diffusivity::Function, source::Function)::SystemOfEquations
    (free_counter, lock_counter) = setup_dofs(mesh)
    system = SystemOfEquations(free_counter, lock_counter)

    for e in mesh.elems
        (K, F, U) = local_system(e, shape_functions, gauss_data, diffusivity, source)
        local_assembly(system, e.nodes, K, F, U)
    end

    for c in mesh.conds
        (K, F, U) = local_system(c)
        local_assembly(system, [c.node], K, F, U)
    end

    return system
end

function solve(self::SystemOfEquations, mesh::Mesh)::Vector{Float64}
    locked_nodes = map(n -> n.id, filter(n -> !n.dof.free, mesh.nodes))
    u_l = view(self.sol, locked_nodes)

    lhs = to_csc(self.A_ff)
    rhs = self.b_f - to_csc(self.A_fl) * u_l

    free_nodes = map(n -> n.id, filter(n ->  n.dof.free, mesh.nodes))
    self.sol[free_nodes] = lhs \ rhs

    return self.sol
end