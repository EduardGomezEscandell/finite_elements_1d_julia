include("../math/system_of_equations.jl")
include("../elements/element.jl")
include("../conditions/condition.jl")
include("../mesh/mesh.jl")

mutable struct SpaceIntegrator
    mesh::Mesh
    gauss_data::GaussData
    free_counter::Integer
    lock_counter::Integer
end

function SpaceIntegrator(mesh::Mesh, n_gauss::Integer)
    return SpaceIntegrator(
        mesh, get_gauss_quadrature(n_gauss), 0, 0
    )
end

# Takes local matrix  (from element and condition) and assembles it into the global matrices
function local_assembly(nodes::Vector{Node}, mat::Matrix{Float64}, output_ff::DokMatrix, output_fl::DokMatrix)::Nothing
    for (i::Int64, dof_i::Dof) in enumerate(map(m -> m.dof, nodes))

        if !dof_i.free
            continue
        end

        for (j::Int64, dof_j::Dof) in enumerate(map(m -> m.dof, nodes))

            if dof_j.free
                push!(output_ff, dof_i.id, dof_j.id, mat[i, j])
            else
                push!(output_fl, dof_i.id, dof_j.id, mat[i, j])
            end
        end
    end
end

# Takes local vector (from element and condition) and assembles it into the global vectors
function local_assembly(nodes::Vector{Node}, vec::Vector{Float64}, output_f::Vector{Float64}, output_l::Vector{Float64})::Nothing
    for (i::Int64, dof_i::Dof) in enumerate(map(m -> m.dof, nodes))
        if dof_i.free
            output_f[dof_i.id] += vec[i]
        else
            output_l[dof_i.id] += vec[i]
        end
    end
end

function initialize(self::SpaceIntegrator)::Nothing
    cache(self.mesh.shape_functions, self.gauss_data)
    (self.free_counter, self.lock_counter) = setup_dofs(self.mesh)
    return nothing
end

function build_diferential_and_source(
        self::SpaceIntegrator; kwargs...
    )::Tuple{SparseMatrixCSC{Float64}, SparseMatrixCSC{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}

    L_ff = DokMatrix(self.free_counter, self.free_counter)
    L_fl = DokMatrix(self.free_counter, self.lock_counter)
    F_f = zeros(Float64, self.free_counter)
    F_l = zeros(Float64, self.lock_counter)
    U_f = zeros(Float64, self.free_counter)
    U_l = zeros(Float64, self.lock_counter)

    for e in self.mesh.elems
        L = calculate_L(e, self.mesh.shape_functions, self.gauss_data; kwargs...)
        F = calculate_F(e, self.mesh.shape_functions, self.gauss_data; kwargs...)
        local_assembly(e.nodes, L, L_ff, L_fl)
        local_assembly(e.nodes, F, F_f,  F_l)
    end

    for c in self.mesh.conds
        (L, F_or_U) = local_system(c; kwargs...)
        local_assembly([c.node], L,      L_ff, L_fl)
        local_assembly([c.node], F_or_U, F_f,  U_l)
    end
    return (to_csc(L_ff), to_csc(L_fl), F_f, U_f, U_l)
end

function build_mass_matrix(self::SpaceIntegrator; kwargs...)::Tuple{SparseMatrixCSC{Float64}, SparseMatrixCSC{Float64}}
    M_ff = DokMatrix(self.free_counter, self.free_counter)
    M_fl = DokMatrix(self.free_counter, self.lock_counter)

    for e in self.mesh.elems
        M = calculate_M(e, self.mesh.shape_functions, self.gauss_data; kwargs...)
        local_assembly(e.nodes, M, M_ff, M_fl)
    end
    return to_csc(M_ff), to_csc(M_fl)
end

function build_lumped_mass_matrix(self::SpaceIntegrator; kwargs...)::Tuple{SparseMatrixCSC{Float64}, SparseMatrixCSC{Float64}}
    M_ff = DokMatrix(self.free_counter, self.free_counter)
    M_fl = DokMatrix(self.free_counter, self.lock_counter)

    shape_fun = ShapeFunctions(self.mesh.shape_functions)
    gauss_data = get_lobatto_quadrature(shape_fun.order+1)
    cache(shape_fun, gauss_data)


    for e in self.mesh.elems
        M = calculate_M(e, shape_fun, gauss_data; kwargs...)
        local_assembly(e.nodes, M, M_ff, M_fl)
    end
    return to_csc(M_ff), to_csc(M_fl)
end

function reconstruct_solution(self::SpaceIntegrator, system::SystemOfEquations)::Vector{Float64}
    reconstruct_solution(self, system.u_f, system.u_l)
end

function reconstruct_solution(self::SpaceIntegrator, u_f::Vector{Float64}, u_l::Vector{Float64})::Vector{Float64}

    free_dofs   = map(n -> n.id, filter(n ->  n.dof.free, self.mesh.nodes))
    locked_dofs = map(n -> n.id, filter(n -> !n.dof.free, self.mesh.nodes))

    u = zeros(Float64, size(self.mesh.nodes))
    u[free_dofs] = u_f
    u[locked_dofs] = u_l
    return u
end
