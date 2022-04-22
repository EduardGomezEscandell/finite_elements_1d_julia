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

#=
For an element

    ∂u/∂t + L(u)*u = f

discretized with shape functions w as:

    (w, ∂u/∂t) - (w, L(u)*u) = (w, f)

- calculate_M: returns (w, u)
- calculate_L: returns (w, L(u))
- calculate_F: returns (w, F)
=#

function calculate_L(self::Element, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Matrix{Float64}
    error("Calling abstract base class local_system")
end

function calculate_M(self::Element, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Matrix{Float64}
    error("Calling abstract base class local_system")
end

function calculate_F(self::Element, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Vector{Float64}
    error("Calling abstract base class local_system")
end