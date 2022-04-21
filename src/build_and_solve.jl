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

function local_system(self::Element, shape_fun::ShapeFunctions, gauss_data::GaussData, k::Function, f::Function)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    nnodes = size(self.nodes, 1)
    L = self.nodes[nnodes].x - self.nodes[1].x
    jacobian = L / 2

    x = self.nodes[1].x .+ L * (gauss_data.points .+ 1) / 2

    B = shape_fun.B ./ jacobian
    f_gauss = f(x)
    k_gauss = k(x)

    if isa(f_gauss, Number)
        f_gauss = f_gauss .* ones(size(x))
    end
    if isa(k_gauss, Number)
        k_gauss = k_gauss .* ones(size(x))
    end

    Klocal = reduce(+, B[:, i] .* k_gauss[i] .* transpose(B[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)
    Flocal = reduce(+, shape_fun.N[:, i] .* f_gauss[i] .* gauss_data.weights[i] for i=1:gauss_data.size)

    Klocal .*= jacobian
    Flocal .*= jacobian

    return (Klocal, Flocal, zeros(Float64,0))
end

function local_system(self::Condition)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    if self.type == DIRICHLET
        return (zeros(Float64, 0, 0), zeros(Float64,0), [self.value])
    elseif self.type == NEUMANN
        return (zeros(Float64, 1, 1), [self.normal * self.value], zeros(Float64,0))
    end
end

function lock_dofs(self::Condition)::Nothing
    if self.type == DIRICHLET
        self.node.dof.free = false
    end
    return nothing
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