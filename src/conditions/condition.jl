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
function Condition(id::Int64, node::Node, value::Float64, normal::Float64)::Condition
function lock_dofs(self::Condition)::Nothing
function local_system(self::Condition)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
=#