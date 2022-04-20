mutable struct SystemOfEquations
    mat::Matrix{Float64}
    vec::Vector{Float64}
    sol::Vector{Float64}
    locked_dofs::Vector{Integer}
end

function SystemOfEquations(nnodes::Int)
    return SystemOfEquations(zeros(Float64, nnodes, nnodes), zeros(Float64, nnodes), zeros(Float64, nnodes), [])
end

function assemble(self::Element, system::SystemOfEquations, shape_fun::ShapeFunctions, gauss_data::GaussData, k::Float64, f::Function)
    nnodes = size(self.nodes, 1)
    L = self.nodes[nnodes].x - self.nodes[1].x
    jacobian = L / 2

    x = self.nodes[1].x .+ L * (gauss_data.points .+ 1) / 2

    B = shape_fun.B ./ jacobian
    f_gauss = f(x)

    Klocal = k .* reduce(+, B[:, i] .* transpose(B[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)
    Flocal = reduce(+, f_gauss[i] .* shape_fun.N[:, i]  .* gauss_data.weights[i] for i=1:gauss_data.size)

    Klocal .*= jacobian
    Flocal .*= jacobian

    i = 1
    for n in self.nodes
        j = 1
        for m in self.nodes
            system.mat[n.id, m.id] += Klocal[i, j]
            j += 1
        end
        system.vec[n.id] += Flocal[i]
        i += 1
    end
end

function assemble(self::Condition, system::SystemOfEquations, shape_fun::ShapeFunctions, gauss_data::GaussData, k::Float64)

    if self.type == DIRICHLET
        push!(system.locked_dofs, self.node.id)
        system.sol[self.node.id] = self.value
    elseif self.type == NEUMANN
        Flocal = self.normal * self.value
        system.vec[self.node.id] += Flocal
    end
end


function build(mesh::Mesh, shape_functions::ShapeFunctions, gauss_data::GaussData, diffusivity::Float64, source::Function)
    nnodes = size(mesh.nodes)[1]
    system = SystemOfEquations(nnodes)
    for e in mesh.elems
        assemble(e, system, shape_functions, gauss_data, diffusivity, source)
    end
    for c in mesh.conds
        assemble(c, system, shape_functions, gauss_data, diffusivity)
    end
    return system
end


function solve(self::SystemOfEquations)
    free_dofs = setdiff(1:size(self.mat)[1], self.locked_dofs)
    Aff = view(self.mat, free_dofs, free_dofs)
    Afl = view(self.mat, free_dofs, self.locked_dofs)
    ul = view(self.sol, self.locked_dofs)
    bf = view(self.vec, free_dofs)

    lhs = Aff
    rhs = bf - Afl * ul

    self.sol[free_dofs] = lhs \ rhs
end