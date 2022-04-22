include("condition.jl")
include("dirichlet_condition.jl")
include("neumann_condition.jl")


function ConditionFactory(
    name::String,
    id::Int64,
    node::Node,
    value::Float64,
    normal::Float64)::Condition

    if name == "Dirichlet"
        return DirichletCondition(id, node, value, normal)
    elseif name == "Neumann"
        return NeumannCondition(id, node, value, normal)
    end

    error("Unknown condition type: $(name)")
end