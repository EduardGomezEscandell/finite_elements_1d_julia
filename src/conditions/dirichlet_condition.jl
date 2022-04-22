include("condition.jl")

mutable struct DirichletCondition <: Condition
    id::Int64
    node::Node
    value::Float64
    normal::Float64
end

function local_system(self::DirichletCondition)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    return (zeros(Float64, 0, 0), zeros(Float64,0), [self.value])
end

function lock_dofs(self::DirichletCondition)::Nothing
    self.node.dof.free = false
    return nothing
end
