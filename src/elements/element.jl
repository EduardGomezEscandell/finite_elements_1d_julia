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
    error("Calling abstract base class constructor")
end

function local_system(self::Element, system::SystemOfEquations, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    error("Calling abstract base class local_system")
end