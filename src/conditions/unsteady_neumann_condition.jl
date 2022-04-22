include("condition.jl")

mutable struct NeumannCondition <: Condition
    id::Int64
    node::Node
    value::Float64
    normal::Float64
end

function local_system(self::NeumannCondition; kwargs...)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    return (zeros(Float64, 1, 1), [self.normal * self.value], zeros(Float64,0))
end

function lock_dofs(self::NeumannCondition)::Nothing
    return nothing
end
