include("../mesh/node.jl")
include("../math/system_of_equations.jl")
include("../math/shape_functions.jl")
include("../math/gauss_data.jl")

abstract type Element end

#=
Fields
---------
id:Int64
nodes::Vector{Node}


Interface
---------
=#
function Element(id::Int64, node::Vector{Node})
    error("Calling base class constructor")
end

function local_system(self::Element, system::SystemOfEquations, shape_fun::ShapeFunctions, gauss_data::GaussData, k::Function, f::Function)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    error("Calling base class local_system")
end