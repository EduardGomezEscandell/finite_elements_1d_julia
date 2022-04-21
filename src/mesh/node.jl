mutable struct Dof
    id::Int64
    free::Bool
end

mutable struct Node
    id::Int64
    x::Float64
    dof::Dof
end