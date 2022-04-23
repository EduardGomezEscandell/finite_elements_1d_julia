struct GaussData
    size :: Int
    weights::Vector{Float64}
    points::Vector{Float64}
end

function get_gauss_quadrature(npoints::Integer)::GaussData
    if npoints == 1
        return GaussData(1,
            [ 2.0000000000000000 ],
            [ 0.0000000000000000 ]
        )
    elseif npoints == 2
        return GaussData(2,
            [ 1.0000000000000000,  1.0000000000000000],
            [-0.5773502691896257,  0.5773502691896257]
        )
    elseif npoints == 3
        return GaussData(3,
            [ 0.5555555555555556,  0.8888888888888888,  0.5555555555555556],
            [-0.7745966692414834,  0.0000000000000000,  0.7745966692414834]
        )
    elseif npoints == 4
        return GaussData(4,
            [ 0.3478548451374538, 0.6521451548625461,  0.6521451548625461,  0.3478548451374538],
            [-0.8611363115940526, -0.3399810435848563,  0.3399810435848563, 0.8611363115940526]
        )
    elseif npoints == 5
        return GaussData(5,
            [ 0.2369268850561891,  0.4786286704993665,  0.5688888888888889,  0.4786286704993665,  0.2369268850561891],
            [-0.9061798459386640, -0.5384693101056831,  0.0000000000000000,  0.5384693101056831,  0.9061798459386640]
        )
    elseif npoints == 6
        return GaussData(6,
            [ 0.1713244923791704,  0.3607615730481386,  0.4679139345726910, 0.3607615730481386, 0.4679139345726910, 0.1713244923791704],
            [-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.6612093864662645, 0.2386191860831969, 0.9324695142031521]
        )
    elseif npoints == 7
        return GaussData(7,
            [ 0.1294849661688697,  0.2797053914892766,  0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697],
            [-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0000000000000000, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585]
        )
    end
    error("Gauss quadrature not implemented for npoints = ", npoints)
end

function get_lobatto_quadrature(npoints::Integer)::GaussData
    midpoint_mass = 2 / (npoints - 1)
    endpoint_mass = 1 / (npoints - 1)

    weights = zeros(Float64, npoints)
    weights[1] = endpoint_mass
    weights[npoints] = endpoint_mass
    weights[2:(npoints - 1)] .= midpoint_mass

    return GaussData(npoints,
        weights,
        LinRange(-1.0, 1.0, npoints),
    )
end