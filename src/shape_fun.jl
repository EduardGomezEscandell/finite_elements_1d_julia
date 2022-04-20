struct ShapeFunctions
    N::Matrix{Float64}
    B::Matrix{Float64}
end

function compute_shape_functions(interp_order::Integer, gauss_data::GaussData)
    if interp_order == 1
        w1(ξ) = [(1 - ξ) / 2;
                 (1 + ξ) / 2]
        dw1(ξ) = [-1/2;
                   1/2]
        return generate_N(w1, dw1, gauss_data)
    elseif interp_order == 2
        w2(ξ) = [-0.5*ξ*(1-ξ);
                 1 - ξ^2;
                 0.5*ξ*(1+ξ)]
        dw2(ξ) = [-0.5+ξ;
                 -2*ξ;
                 0.5+ξ]
        return generate_N(w2, dw2, gauss_data)
    elseif interp_order == 3
        w3(ξ) = [- 9/16*(ξ-1)*(ξ^2 - 1/9);
                 +27/16*(ξ^2 - 1)*(ξ - 1/3);
                 -27/16*(ξ^2 - 1)*(ξ + 1/3);
                 + 9/16*(ξ+1)*(ξ^2 - 1/9)]
        dw3(ξ) = [ (1 + 18*ξ - 27*ξ^2)/16;
                   9/16*(-3 - 2*ξ + 9*ξ^2);
                 - 9/16*(-3 + 2*ξ + 9*ξ^2);
                 (-1 + 18*ξ + 27*ξ^2)/16]
        return generate_N(w3, dw3, gauss_data)
    end
    error("Not implemented for interp_order = ", interp_order)
end

function generate_N(w::Function, dw::Function, gauss_data::GaussData)
    return ShapeFunctions(
        mapreduce(permutedims, hcat, transpose([w(ξ) for ξ=gauss_data.points])),
        mapreduce(permutedims, hcat, transpose([dw(ξ) for ξ=gauss_data.points])),
    )
end