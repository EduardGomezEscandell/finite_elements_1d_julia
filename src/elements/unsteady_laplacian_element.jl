include("element.jl")

# Solves for:
#
#      ∂u/∂t - ∇·μ∇u = s
#
# with a Crank-Nicolson scheme.

mutable struct UnsteadyLaplacianElement <: Element
    id::Int64
    nodes::Vector{Node}
end

function local_system(self::UnsteadyLaplacianElement, system::SystemOfEquations, shape_fun::ShapeFunctions, gauss_data::GaussData; kwargs...)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    s = kwargs[:s]
    μ = kwargs[:μ]
    t = kwargs[:t]
    Δt = kwargs[:Δt]
    u_old = [kwargs[:u_old][node.id] for node in self.nodes]

    nnodes = size(self.nodes, 1)
    L = self.nodes[nnodes].x - self.nodes[1].x
    jacobian = L / 2

    x = self.nodes[1].x .+ L * (gauss_data.points .+ 1) / 2

    B = shape_fun.B ./ jacobian
    N = shape_fun.N

    s_gauss = (s(t, x) + s(t + Δt, x)) / 2
    μ_gauss_old = μ(t, x)
    μ_gauss_new = μ(t + Δt, x)

    if isa(s_gauss, Number)
        s_gauss = s_gauss .* ones(size(x))
    end
    if isa(μ_gauss_old, Number)
        μ_gauss_old = μ_gauss_old .* ones(size(x))
    end
    if isa(μ_gauss_new, Number)
        μ_gauss_new = μ_gauss_new .* ones(size(x))
    end

    M     = jacobian * reduce(+, N[:, i] .*                   transpose(N[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)
    K_old = jacobian * reduce(+, B[:, i] .* μ_gauss_old[i] .* transpose(B[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)
    K_new = jacobian * reduce(+, B[:, i] .* μ_gauss_new[i] .* transpose(B[:, i]) .* gauss_data.weights[i] for i=1:gauss_data.size)

    F     = jacobian * reduce(+, N[:, i] .* s_gauss[i] .* gauss_data.weights[i] for i=1:gauss_data.size)

    A =  M + Δt/2 * K_new
    b = (M - Δt/2 * K_old) * u_old + Δt*F

    return (A, b, zeros(Float64,0))
end
