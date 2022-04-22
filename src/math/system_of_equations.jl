include("sparse.jl")

mutable struct SystemOfEquations
    LHS_ff::SparseMatrixCSC{Float64}
    LHS_fl::SparseMatrixCSC{Float64}
    RHS_f::Vector{Float64}
    u_f::Vector{Float64}
    u_l::Vector{Float64}
end

function solve(self::SystemOfEquations)::Nothing
    lhs = self.LHS_ff
    rhs = self.RHS_f - self.LHS_fl * self.u_l
    self.u_f = lhs \ rhs
    return nothing
end