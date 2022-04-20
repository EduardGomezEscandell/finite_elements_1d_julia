struct DofNode
    id::Integer
    x::Float64
end

struct Element
    id::Integer
    nodes::Vector{DofNode}
end

@enum BoundaryTypes DIRICHLET NEUMANN

struct Condition
    id::Integer
    node::DofNode
    type::BoundaryTypes
    value::Float64
end

struct Mesh
    nodes::Vector{DofNode}
    elems::Vector{Element}
    conds::Vector{Condition}
end

function generate_mesh(
        nelems::Integer,
        polynomial_order::Integer,
        length::Float64,
        left_bc::Tuple{BoundaryTypes, Float64},
        right_bc::Tuple{BoundaryTypes, Float64}
    )
    nnodes = nelems*polynomial_order + 1
    nodes = [DofNode(i, (i-1)*length/(nnodes-1)) for i=1:nnodes]
    elems = [Element(el, nodes[(el-1)*polynomial_order+1:el*polynomial_order+1]) for el=1:nelems]
    conds = [Condition(1, nodes[1],      left_bc[1],  left_bc[2]),
            Condition(2, nodes[nnodes], right_bc[1], right_bc[2])]
    return Mesh(nodes, elems, conds)
end