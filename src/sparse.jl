using SparseArrays

mutable struct DokMatrix
    I::Array{Int64,1}
    J::Array{Int64,1}
    X::Array{Float64,1}
    nrows::Int64
    ncols::Int64
end

function DokMatrix(nrows::Integer, ncols::Integer)::DokMatrix
    return DokMatrix(
        Array{Int64,1}(),
        Array{Int64,1}(),
        Array{Float64,1}(),
        nrows,
        ncols
    )
end

function Base.push!(self::DokMatrix, i::Int64, j::Int64, x::Float64)::Nothing
    if(i < 1 || i > self.nrows || j < 1 || j > self.ncols)
        error("Index out of bounds")
    end
    push!(self.I, i)
    push!(self.J, j)
    push!(self.X, x)
    return nothing
end

function to_csc(self::DokMatrix, nrows::Integer, ncols::Integer)::SparseMatrixCSC{Float64}
    return sparse(self.I, self.J, self.X, nrows, ncols, +)
end