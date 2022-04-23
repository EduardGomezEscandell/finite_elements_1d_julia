include("gauss_data.jl")

mutable struct ShapeFunctions
    order::Integer
    w::Function
    ∂w::Function
    N::Matrix{Float64}
    B::Matrix{Float64}
end

function ShapeFunctions(other::ShapeFunctions)::ShapeFunctions
    return ShapeFunctions(other.order, other.w, other.∂w, zeros(Float64,0,0), zeros(Float64,0,0))
end

function ShapeFunctions(interp_order::Integer)::ShapeFunctions
    if interp_order == 1
        return ShapeFunctions(1,
            ξ -> [(1 - ξ) / 2;    (1 + ξ) / 2],
            ξ -> [-1/2;                  1/2],
            zeros(Float64,0,0), zeros(Float64,0,0))
    elseif interp_order == 2
        return ShapeFunctions(2,
            ξ -> [-0.5*ξ*(1-ξ);   1 - ξ^2;   0.5*ξ*(1+ξ)],
            ξ -> [-0.5+ξ;           -2*ξ;         0.5+ξ],
            zeros(Float64,0,0), zeros(Float64,0,0))
    elseif interp_order == 3
        return ShapeFunctions(3,
            ξ -> [- 9/16*(ξ-1)*(ξ^2 - 1/9);    +27/16*(ξ^2 - 1)*(ξ - 1/3);    -27/16*(ξ^2 - 1)*(ξ + 1/3);    + 9/16*(ξ+1)*(ξ^2 - 1/9)],
            ξ -> [(1 + 18*ξ - 27*ξ^2)/16;        9/16*(-3 - 2*ξ + 9*ξ^2);    - 9/16*(-3 + 2*ξ + 9*ξ^2);      (-1 + 18*ξ + 27*ξ^2)/16],
            zeros(Float64,0,0), zeros(Float64,0,0))
    end
    error("Not implemented for interp_order = ", interp_order)
end

# Converts analytical shape functions to pre-calculated ones
function cache(self::ShapeFunctions, gauss_data::GaussData)::Nothing
    self.N = mapreduce(permutedims, hcat, transpose([self.w(ξ) for ξ=gauss_data.points]))
    self.B = mapreduce(permutedims, hcat, transpose([self.∂w(ξ) for ξ=gauss_data.points]))
    return nothing
end
