include("gauss_data.jl")

struct ShapeFunctions
    N::Matrix{Float64}
    B::Matrix{Float64}
end

# Converts analytical shape functions to pre-calculated ones
function ShapeFunctions(w::Function, dw::Function, gauss_data::GaussData)::ShapeFunctions
    return ShapeFunctions(
        mapreduce(permutedims, hcat, transpose([w(ξ) for ξ=gauss_data.points])),
        mapreduce(permutedims, hcat, transpose([dw(ξ) for ξ=gauss_data.points])),
    )
end

function compute_shape_functions(interp_order::Integer, gauss_data::GaussData)::ShapeFunctions
    if interp_order == 1
        w1(ξ) = [(1 - ξ) / 2;
                 (1 + ξ) / 2]
        dw1(ξ) = [-1/2;
                   1/2]
        return ShapeFunctions(w1, dw1, gauss_data)
    elseif interp_order == 2
        w2(ξ) = [-0.5*ξ*(1-ξ);
                 1 - ξ^2;
                 0.5*ξ*(1+ξ)]
        dw2(ξ) = [-0.5+ξ;
                 -2*ξ;
                 0.5+ξ]
        return ShapeFunctions(w2, dw2, gauss_data)
    elseif interp_order == 3
        w3(ξ) = [- 9/16*(ξ-1)*(ξ^2 - 1/9);
                 +27/16*(ξ^2 - 1)*(ξ - 1/3);
                 -27/16*(ξ^2 - 1)*(ξ + 1/3);
                 + 9/16*(ξ+1)*(ξ^2 - 1/9)]
        dw3(ξ) = [ (1 + 18*ξ - 27*ξ^2)/16;
                   9/16*(-3 - 2*ξ + 9*ξ^2);
                 - 9/16*(-3 + 2*ξ + 9*ξ^2);
                 (-1 + 18*ξ + 27*ξ^2)/16]
        return ShapeFunctions(w3, dw3, gauss_data)
    end
    error("Not implemented for interp_order = ", interp_order)
end