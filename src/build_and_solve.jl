include("sparse.jl")

mutable struct SystemOfEquations
    mat::DokMatrix
    vec::Vector{Float64}
    sol::Vector{Float64}
    locked_dofs::Vector{Integer}
    size::Integer
end

function SystemOfEquations(size::Integer)::SystemOfEquations
    return SystemOfEquations(DokMatrix(size, size), zeros(Float64, size), zeros(Float64, size), [], size)
end

function assemble(self::Element, system::SystemOfEquations, shape_fun::ShapeFunctions, gauss_data::GaussData, k::Function, f::Function)::Nothing
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

    i = 1
    for n in self.nodes
        j = 1
        for m in self.nodes
            add(system.mat, n.id, m.id, Klocal[i, j])
            j += 1
        end
        system.vec[n.id] += Flocal[i]
        i += 1
    end
    return nothing
end

function assemble(self::Condition, system::SystemOfEquations)::Nothing
    if self.type == DIRICHLET
        push!(system.locked_dofs, self.node.id)
        system.sol[self.node.id] = self.value
    elseif self.type == NEUMANN
        Flocal = self.normal * self.value
        system.vec[self.node.id] += Flocal
    end
    return nothing
end

function build(mesh::Mesh, shape_functions::ShapeFunctions, gauss_data::GaussData, diffusivity::Function, source::Function)::SystemOfEquations
    nnodes = size(mesh.nodes)[1]
    system = SystemOfEquations(nnodes)
    for e in mesh.elems
        assemble(e, system, shape_functions, gauss_data, diffusivity, source)
    end
    for c in mesh.conds
        assemble(c, system)
    end
    return system
end

function solve(self::SystemOfEquations)::Vector{Float64}
    free_dofs = setdiff(1:self.size, self.locked_dofs)
    A = to_csc(self.mat, self.size, self.size)
    Aff = view(A, free_dofs, free_dofs)
    Afl = view(A, free_dofs, self.locked_dofs)
    ul = view(self.sol, self.locked_dofs)
    bf = view(self.vec, free_dofs)

    rhs = bf - Afl * ul
    lhs = SparseMatrixCSC(Aff) # Should not be a copy, but Julia can't do lu!(view, vector)

    self.sol[free_dofs] = lhs \ rhs

    return self.sol
end