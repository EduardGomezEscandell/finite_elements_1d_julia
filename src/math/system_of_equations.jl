include("sparse.jl")

mutable struct SystemOfEquations
    A_ff::DokMatrix
    A_fl::DokMatrix
    b_f::Vector{Float64}
    u_f::Vector{Float64}
    u_l::Vector{Float64}
end

function SystemOfEquations()::SystemOfEquations
    return SystemOfEquations(
        DokMatrix(0, 0),
        DokMatrix(0, 0),
        zeros(Float64, 0),
        zeros(Float64, 0),
        zeros(Float64, 0)
        )
end

function SystemOfEquations(free_dof::Integer, locked_dof::Integer)::SystemOfEquations
    return SystemOfEquations(
        DokMatrix(free_dof, free_dof),
        DokMatrix(free_dof, locked_dof),
        zeros(Float64, free_dof),
        zeros(Float64, free_dof),
        zeros(Float64, locked_dof)
        )
end

function solve(self::SystemOfEquations)::Nothing
    lhs = to_csc(self.A_ff)
    rhs = self.b_f - to_csc(self.A_fl) * self.u_l
    self.u_f = lhs \ rhs
    return nothing
end