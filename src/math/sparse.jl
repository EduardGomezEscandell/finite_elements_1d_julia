using SparseArrays

mutable struct DokMatrix
    I::Array{Int64,1}
    J::Array{Int64,1}
    X::Array{Float64,1}
    nrows::Int64
    ncols::Int64
    ε::Float64              # Numbers smaller in absolute value than ε are considered zero
end

function DokMatrix(nrows::Integer, ncols::Integer, ε::Float64=1e-15)::DokMatrix
    return DokMatrix(
        Array{Int64,1}(),
        Array{Int64,1}(),
        Array{Float64,1}(),
        nrows,
        ncols,
        ε
    )
end

function Base.push!(self::DokMatrix, i::Int64, j::Int64, x::Float64)::Nothing
    if(i < 1 || i > self.nrows || j < 1 || j > self.ncols)
        error("Index out of bounds")
    end

    if abs(x) < self.ε
        return nothing
    end

    push!(self.I, i)
    push!(self.J, j)
    push!(self.X, x)
    return nothing
end

function to_csc(self::DokMatrix)::SparseMatrixCSC{Float64}
    return sparse(self.I, self.J, self.X, self.nrows, self.ncols, +)
end