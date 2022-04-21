mutable struct Dof
    id::Int64
    free::Bool
end
mutable struct Node
    id::Int64
    x::Float64
    dof::Dof
end

mutable struct Element
    id::Int64
    nodes::Vector{Node}
end

@enum BoundaryTypes DIRICHLET NEUMANN

mutable struct Condition
    id::Int64
    node::Node
    type::BoundaryTypes
    value::Float64
    normal::Float64
end

struct Mesh
    nodes::Vector{Node}
    elems::Vector{Element}
    conds::Vector{Condition}
end

function generate_mesh(
        nelems::Integer,
        polynomial_order::Integer,
        length::Float64,
        left_bc::Tuple{BoundaryTypes, Float64},
        right_bc::Tuple{BoundaryTypes, Float64}
    )::Mesh
    nnodes = nelems*polynomial_order + 1
    nodes = [Node(i, (i-1)*length/(nnodes-1), Dof(0, true)) for i=1:nnodes]
    elems = [Element(el, nodes[(el-1)*polynomial_order+1:el*polynomial_order+1]) for el=1:nelems]
    conds = [Condition(1, nodes[1],      left_bc[1],  left_bc[2], -1),
            Condition(2, nodes[nnodes], right_bc[1], right_bc[2], 1)]
    return Mesh(nodes, elems, conds)
end