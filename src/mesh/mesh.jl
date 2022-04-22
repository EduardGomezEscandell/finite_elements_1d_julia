include("node.jl")
include("../elements/factory.jl")
include("../conditions/factory.jl")

mutable struct Mesh
    nodes::Vector{Node}
    elems::Vector{Element}
    conds::Vector{Condition}
    shape_functions::ShapeFunctions
end

function generate_mesh(
        length::Float64,
        element::Tuple{String, Integer, Integer},
        left_bc::Tuple{String, Float64},
        right_bc::Tuple{String, Float64},
    )::Mesh
    element_name, polynomial_order, nelems = element
    left_bc_name, left_bc_value = left_bc
    right_bc_name, right_bc_value = right_bc

    nnodes = nelems*polynomial_order + 1
    nodes = [Node(i, (i-1)*length/(nnodes-1), Dof(0, true)) for i=1:nnodes]
    elems = [ElementFactory(element_name, el, nodes[(el-1)*polynomial_order+1:el*polynomial_order+1]) for el=1:nelems]
    conds = [ConditionFactory(left_bc_name,  1, nodes[1],       left_bc_value, -1.0),
             ConditionFactory(right_bc_name, 2, nodes[nnodes], right_bc_value,  1.0)]
    return Mesh(nodes, elems, conds, ShapeFunctions(polynomial_order))
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