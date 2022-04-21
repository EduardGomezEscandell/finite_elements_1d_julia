include("node.jl")

@enum ConditionType DIRICHLET NEUMANN

mutable struct Condition
    id::Int64
    node::Node
    type::ConditionType
    value::Float64
    normal::Float64
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
