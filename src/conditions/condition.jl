include("../mesh/node.jl")

abstract type Condition end

#=

Fields
---------

id::Int64
node::Node
value::Float64
normal::Float64

Interface
---------
=#
function Condition(id::Int64, node::Node, value::Float64, normal::Float64)::Condition
    error("Calling abstract base class Condition constructor")
end

function lock_dofs(self::Condition)::Nothing
    error("Calling abstract base class lock_dofs")
end

function local_system(self::Condition; kwargs...)::Tuple{Matrix{Float64}, Vector{Float64}}
    error("Calling abstract base class local_system")
end