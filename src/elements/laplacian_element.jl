include("element.jl")

mutable struct LaplacianElement <: Element
    id::Int64
    nodes::Vector{Node}
end

function local_system(self::LaplacianElement, system::SystemOfEquations, shape_fun::ShapeFunctions, gauss_data::GaussData, k::Function, f::Function)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    nnodes = size(self.nodes, 1)
    L = self.nodes[nnodes].x - self.nodes[1].x
    jacobian = L / 2

    x = self.nodes[1].x .+ L * (gauss_data.points .+ 1) / 2

    B = shape_fun.B ./ jacobian
    f_gauss = f(x)
    k_gauss = k(x)

    if isa(f_gauss, Number)
        f_gauss = f_gauss .* ones(size(x))
    end
    if isa(k_gauss, Number)
        k_gauss = k_gauss .* ones(size(x))
    end

    Klocal = jacobian * reduce(+, B[:, i] .* k_gauss[i] .* transpose(B[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)
    Flocal = jacobian * reduce(+, shape_fun.N[:, i] .* f_gauss[i] .* gauss_data.weights[i] for i=1:gauss_data.size)

    return (Klocal, Flocal, zeros(Float64,0))
end