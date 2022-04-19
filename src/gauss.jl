struct GaussData
    size :: Int
    weights::Vector{Float64}
    points::Vector{Float64}
end

function get_gauss_quadrature(npoints::Integer)
    if npoints == 2
        return GaussData(2,
            [ 1.0000000000000000,  1.0000000000000000],
            [-0.5773502691896257,  0.5773502691896257]
        )
    elseif npoints == 3
        return GaussData(3,
            [ 0.8888888888888888,  0.5555555555555556,  0.5555555555555556],
            [ 0.0000000000000000, -0.7745966692414834,  0.7745966692414834]
        )
    elseif npoints == 4
        return GaussData(4,
            [ 0.6521451548625461,  0.6521451548625461,  0.3478548451374538,  0.3478548451374538],
            [-0.3399810435848563,  0.3399810435848563, -0.8611363115940526,  0.8611363115940526]
        )
    elseif npoints == 5
        return GaussData(5,
            [ 0.5688888888888889,  0.4786286704993665,  0.4786286704993665,  0.2369268850561891,  0.2369268850561891],
            [ 0.0000000000000000, -0.5384693101056831,  0.5384693101056831, -0.9061798459386640,  0.9061798459386640]
        )
    end
    error("Gauss quadrature not implemented for npoints = ", npoints)
end
