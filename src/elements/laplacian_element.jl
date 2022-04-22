include("element.jl")

mutable struct LaplacianElement <: Element
    id::Int64
    nodes::Vector{Node}
end

function calculate_L(self::LaplacianElement, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Matrix{Float64}
    μ = kwargs[:μ]
    t = kwargs[:t]

    nnodes = size(self.nodes, 1)
    L = self.nodes[nnodes].x - self.nodes[1].x
    jacobian = L / 2

    x = self.nodes[1].x .+ L * (gauss_data.points .+ 1) / 2

    B = shape_fun.B ./ jacobian
    μ_gauss = μ(t, x)

    if isa(μ_gauss, Number)
        μ_gauss = μ_gauss .* ones(size(x))
    end

    Klocal = jacobian * reduce(+, B[:, i] .* μ_gauss[i] .* transpose(B[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)

    return Klocal
end

function calculate_F(self::LaplacianElement, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Vector{Float64}
    s = kwargs[:s]
    t = kwargs[:t]

    nnodes = size(self.nodes, 1)
    L = self.nodes[nnodes].x - self.nodes[1].x
    jacobian = L / 2

    x = self.nodes[1].x .+ L * (gauss_data.points .+ 1) / 2

    s_gauss = s(t, x)

    if isa(s_gauss, Number)
        s_gauss = s_gauss .* ones(size(x))
    end

    Flocal = jacobian * reduce(+, shape_fun.N[:, i] .* s_gauss[i] .* gauss_data.weights[i] for i=1:gauss_data.size)

    return Flocal
end

function calculate_M(self::LaplacianElement, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Matrix{Float64}
    nnodes = size(self.nodes, 1)
    L = self.nodes[nnodes].x - self.nodes[1].x
    jacobian = L / 2

    N = shape_fun.N

    return jacobian * reduce(+, N[:, i] .* transpose(N[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)
end