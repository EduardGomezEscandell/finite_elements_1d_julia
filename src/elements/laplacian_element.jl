include("element.jl")

mutable struct LaplacianElement <: Element
    id::Int64
    nodes::Vector{Node}
end

function local_system(self::LaplacianElement, system::SystemOfEquations, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    s = kwargs[:s]
    μ = kwargs[:μ]

    nnodes = size(self.nodes, 1)
    L = self.nodes[nnodes].x - self.nodes[1].x
    jacobian = L / 2

    x = self.nodes[1].x .+ L * (gauss_data.points .+ 1) / 2

    B = shape_fun.B ./ jacobian
    s_gauss = s(x)
    μ_gauss = μ(x)

    if isa(s_gauss, Number)
        s_gauss = s_gauss .* ones(size(x))
    end
    if isa(μ_gauss, Number)
        μ_gauss = μ_gauss .* ones(size(x))
    end

    Klocal = jacobian * reduce(+, B[:, i] .* μ_gauss[i] .* transpose(B[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)
    Flocal = jacobian * reduce(+, shape_fun.N[:, i] .* s_gauss[i] .* gauss_data.weights[i] for i=1:gauss_data.size)

    return (Klocal, Flocal, zeros(Float64,0))
end